#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(writexl)
  library(RColorBrewer)
  library(BiocParallel)

set.seed(123)
})

pdf.options(family = "Helvetica")
options(stringsAsFactors = FALSE)

n_cores <- min(parallel::detectCores() - 1, 8)
register(MulticoreParam(workers = n_cores, progressbar = TRUE))

cat("================================================================================\n")
cat("SCENIC Publication Enhancement Analysis (PARALLEL VERSION)\n")
cat(paste("Using", n_cores, "CPU cores for parallel processing\n"))
cat("Strategy: Unbiased Data-Driven Discovery\n")
cat("================================================================================\n\n")

cat("PART 1: Loading existing SCENIC results...\n")
cat("----------------------------------------\n")

scenicOptions <- readRDS("int/scenicOptions.rds")
regulonAUC <- readRDS("int/3.4_regulonAUC.Rds")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

cellInfo <- readRDS("int/cellInfo.Rds")

regulon_targets <- read.table("output/Step2_regulonTargetsInfo.tsv",
                              header = TRUE, sep = "\t", quote = "")

cat(paste("OK Loaded", nrow(regulonAUC), "regulons\n"))
cat(paste("OK Loaded", ncol(regulonAUC), "cells\n"))
cat(paste("OK Cell types:", paste(unique(cellInfo$celltype), collapse = ", "), "\n\n"))

cat("PART 2: Statistical significance analysis...\n")
cat("----------------------------------------\n")

auc_matrix <- getAUC(regulonAUC)
cell_types <- unique(cellInfo$celltype)

cat(paste("Performing statistical tests for", nrow(auc_matrix), "regulons across",
          length(cell_types), "cell types...\n"))

stat_results <- data.frame()

for (regulon in rownames(auc_matrix)) {

  auc_by_type <- split(auc_matrix[regulon, ], cellInfo$celltype)

  kw_test <- kruskal.test(auc_matrix[regulon, ] ~ cellInfo$celltype)

  pairwise_results <- list()
  for (i in 1:(length(cell_types)-1)) {
    for (j in (i+1):length(cell_types)) {
      type1 <- cell_types[i]
      type2 <- cell_types[j]

      wilcox_result <- wilcox.test(auc_by_type[[type1]], auc_by_type[[type2]])

      n1 <- length(auc_by_type[[type1]])
      n2 <- length(auc_by_type[[type2]])
      U <- wilcox_result$statistic
      r <- 1 - (2*U)/(n1*n2)

      pairwise_results[[paste0(type1, "_vs_", type2)]] <- list(
        p_value = wilcox_result$p.value,
        effect_size = r
      )
    }
  }

  mean_auc <- sapply(auc_by_type, mean)

  top_celltype <- names(which.max(mean_auc))

  result_row <- data.frame(
    Regulon = regulon,
    KW_pvalue = kw_test$p.value,
    Top_CellType = top_celltype,
    Top_Mean_AUC = max(mean_auc),
    stringsAsFactors = FALSE
  )

  for (ct in cell_types) {
    result_row[[paste0("Mean_AUC_", ct)]] <- mean_auc[ct]
  }

  for (comparison in names(pairwise_results)) {
    result_row[[paste0("Pvalue_", comparison)]] <- pairwise_results[[comparison]]$p_value
    result_row[[paste0("EffectSize_", comparison)]] <- pairwise_results[[comparison]]$effect_size
  }

  stat_results <- rbind(stat_results, result_row)
}

stat_results$KW_FDR <- p.adjust(stat_results$KW_pvalue, method = "BH")

stat_results$Significant <- stat_results$KW_FDR < 0.1
sig_count <- sum(stat_results$Significant)

cat(paste("OK Found", sig_count, "regulons with significant activity differences (FDR < 0.1)\n"))

write.table(stat_results,
            "publication_results/statistical_tests/regulon_statistical_tests.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("OK Saved: publication_results/statistical_tests/regulon_statistical_tests.tsv\n\n")

cat("PART 3: PARALLEL genome-wide enrichment analysis...\n")
cat("----------------------------------------\n")
cat(paste("Strategy: Multi-threaded enrichment using", n_cores, "CPU cores\n"))
cat("          Unbiased data-driven discovery (GO/KEGG)\n\n")

unique_regulons <- unique(regulon_targets$TF)

cat(paste("Performing enrichment for", length(unique_regulons), "regulons...\n"))
cat(paste("Estimated time: ~", round(length(unique_regulons) / (n_cores * 2.5)), "minutes\n\n"))

perform_enrichment <- function(regulon, regulon_targets) {
  results <- list()

  target_genes <- regulon_targets$gene[regulon_targets$TF == regulon]

  if (length(target_genes) < 5) return(NULL)

  tryCatch({
    gene_entrez <- bitr(target_genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    if (nrow(gene_entrez) < 5) return(NULL)

    ego_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                      universe = NULL,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

    if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
      ego_bp_df <- as.data.frame(ego_bp)
      ego_bp_df$Regulon <- regulon
      ego_bp_df$Ontology <- "BP"
      results[[paste0(regulon, "_BP")]] <- ego_bp_df
    }

    ego_mf <- enrichGO(gene = gene_entrez$ENTREZID,
                      universe = NULL,
                      OrgDb = org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

    if (!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
      ego_mf_df <- as.data.frame(ego_mf)
      ego_mf_df$Regulon <- regulon
      ego_mf_df$Ontology <- "MF"
      results[[paste0(regulon, "_MF")]] <- ego_mf_df
    }

    ego_cc <- enrichGO(gene = gene_entrez$ENTREZID,
                      universe = NULL,
                      OrgDb = org.Hs.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

    if (!is.null(ego_cc) && nrow(ego_cc@result) > 0) {
      ego_cc_df <- as.data.frame(ego_cc)
      ego_cc_df$Regulon <- regulon
      ego_cc_df$Ontology <- "CC"
      results[[paste0(regulon, "_CC")]] <- ego_cc_df
    }

    ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                       organism = 'hsa',
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      ekegg_df <- as.data.frame(ekegg)
      ekegg_df$Regulon <- regulon
      results[[paste0(regulon, "_KEGG")]] <- ekegg_df
    }

    return(results)

  }, error = function(e) {
    return(NULL)
  })
}

enrichment_results <- bplapply(unique_regulons, perform_enrichment,
                               regulon_targets = regulon_targets,
                               BPPARAM = MulticoreParam(workers = n_cores, progressbar = TRUE))

enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]

all_results <- unlist(enrichment_results, recursive = FALSE)

go_results <- all_results[grepl("_(BP|MF|CC)$", names(all_results))]
kegg_results <- all_results[grepl("_KEGG$", names(all_results))]

cat("\n")

if (length(go_results) > 0) {
  go_combined <- do.call(rbind, go_results)
  dir.create("publication_results/enrichment_analysis", recursive = TRUE, showWarnings = FALSE)
  write.table(go_combined,
              "publication_results/enrichment_analysis/GO_enrichment_complete.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("OK GO enrichment:", nrow(go_combined), "significant terms across all ontologies\n"))

  go_summary <- go_combined %>%
    group_by(Ontology) %>%
    summarise(N_Terms = n(), N_Regulons = n_distinct(Regulon))
  cat("\n  GO Enrichment Summary:\n")
  for (i in 1:nrow(go_summary)) {
    cat(paste("    ", go_summary$Ontology[i], ":",
              go_summary$N_Terms[i], "terms from",
              go_summary$N_Regulons[i], "regulons\n"))
  }
}

if (length(kegg_results) > 0) {
  kegg_combined <- do.call(rbind, kegg_results)
  write.table(kegg_combined,
              "publication_results/enrichment_analysis/KEGG_enrichment_complete.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("\nOK KEGG enrichment:", nrow(kegg_combined), "significant pathways\n"))
  cat(paste("  From", length(unique(kegg_combined$Regulon)), "regulons\n"))
}

cat("\n")

cat("PART 4: Generating publication-quality figures...\n")
cat("----------------------------------------\n")

cat("Generating Figure 1: Regulon activity heatmap...\n")

top_sig_regulons <- stat_results %>%
  filter(Significant == TRUE) %>%
  arrange(KW_pvalue) %>%
  head(30) %>%
  pull(Regulon)

if (length(top_sig_regulons) > 0) {

  regulon_activity_by_type <- sapply(cell_types, function(ct) {
    cells <- rownames(cellInfo)[cellInfo$celltype == ct]
    rowMeans(auc_matrix[top_sig_regulons, cells, drop = FALSE])
  })

  regulon_activity_scaled <- t(scale(t(regulon_activity_by_type)))

  regulon_activity_scaled <- regulon_activity_scaled[!apply(is.na(regulon_activity_scaled), 1, all), , drop = FALSE]

  top_sig_regulons <- rownames(regulon_activity_scaled)

  sig_annotation <- data.frame(
    FDR = cut(stat_results$KW_FDR[match(top_sig_regulons, stat_results$Regulon)],
              breaks = c(0, 0.01, 0.05, 0.1, 1),
              labels = c("***", "**", "*", "ns")),
    row.names = top_sig_regulons
  )

  ann_colors <- list(
    FDR = c("***" = "#D73027", "**" = "#FC8D59", "*" = "#FEE090", "ns" = "#E0E0E0")
  )

  data_range <- range(regulon_activity_scaled, na.rm = TRUE)
  break_min <- max(data_range[1], -3)
  break_max <- min(data_range[2], 3)
  color_breaks <- seq(break_min, break_max, length.out = 101)

  pdf("publication_results/figures/Fig1_Regulon_Activity_Heatmap.pdf",
      width = 8, height = 12)
  pheatmap(regulon_activity_scaled,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = color_breaks,
           annotation_row = sig_annotation,
           annotation_colors = ann_colors,
           fontsize_row = 7,
           fontsize_col = 10,
           cluster_cols = FALSE,
           border_color = NA,
           main = "Top 30 Significantly Differential Regulons\n(Scaled Activity by Cell Type)")
  dev.off()

  cat("OK Saved: publication_results/figures/Fig1_Regulon_Activity_Heatmap.pdf\n")
}

cat("Generating Figure 2: Regulon specificity scores...\n")

if (file.exists("output/regulon_specificity_scores.rds")) {
  rss <- readRDS("output/regulon_specificity_scores.rds")

  top_specific <- lapply(cell_types, function(ct) {
    ct_rss <- sort(rss[, ct], decreasing = TRUE)
    names(ct_rss)[1:min(5, length(ct_rss))]
  })
  top_specific_regulons <- unique(unlist(top_specific))

  rss_scaled <- t(scale(t(rss), center = FALSE, scale = TRUE))

  valid_rows <- !apply(is.na(rss_scaled), 1, all)
  rss_scaled <- rss_scaled[valid_rows, , drop = FALSE]

  top_specific_regulons <- intersect(top_specific_regulons, rownames(rss_scaled))

  if (length(top_specific_regulons) > 0) {

    data_range <- range(rss_scaled[top_specific_regulons, ], na.rm = TRUE)
    color_breaks <- seq(data_range[1], data_range[2], length.out = 101)

    pdf("publication_results/figures/Fig2_Regulon_Specificity_Scores.pdf",
        width = 8, height = 10)
    pheatmap(rss_scaled[top_specific_regulons, ],
             color = colorRampPalette(c("white", "orange", "red"))(100),
             breaks = color_breaks,
             fontsize_row = 7,
             border_color = NA,
             main = "Cell Type-Specific Regulons\n(Regulon Specificity Scores)")
    dev.off()

    cat("OK Saved: publication_results/figures/Fig2_Regulon_Specificity_Scores.pdf\n")
  }
}

if (exists("go_combined") && nrow(go_combined) > 0) {
  cat("Generating Figure 3: Top GO enrichment terms...\n")

  top_regulons_for_viz <- stat_results %>%
    filter(Significant == TRUE) %>%
    arrange(KW_pvalue) %>%
    head(20) %>%
    pull(Regulon)

  plot_data <- go_combined %>%
    filter(Regulon %in% top_regulons_for_viz, Ontology == "BP") %>%
    group_by(Regulon) %>%
    slice_min(order_by = p.adjust, n = 3) %>%
    ungroup()

  if (nrow(plot_data) > 0) {

    plot_data$Description_short <- substr(plot_data$Description, 1, 50)

    p <- ggplot(plot_data, aes(x = Regulon, y = Description_short,
                                size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      theme_minimal(base_size = 9, base_family = "Helvetica") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.text.y = element_text(size = 7)) +
      labs(title = "Top GO Terms per Regulon (Biological Process)",
           subtitle = "Data-driven discovery - Top 20 significant regulons",
           x = "Regulon",
           y = "GO Term",
           size = "Gene Count",
           color = "Adjusted P-value")

    ggsave("publication_results/figures/Fig3_GO_Enrichment_Dotplot.pdf",
           plot = p, width = 14, height = 10)

    cat("OK Saved: publication_results/figures/Fig3_GO_Enrichment_Dotplot.pdf\n")
  }
}

if (exists("kegg_combined") && nrow(kegg_combined) > 0) {
  cat("Generating Figure 4: KEGG pathway enrichment overview...\n")

  top_kegg <- kegg_combined %>%
    group_by(Description) %>%
    summarise(
      N_Regulons = n_distinct(Regulon),
      Min_Pvalue = min(p.adjust),
      Total_Count = sum(Count)
    ) %>%
    arrange(desc(N_Regulons), Min_Pvalue) %>%
    head(30)

  if (nrow(top_kegg) > 0) {
    p <- ggplot(top_kegg, aes(x = reorder(Description, N_Regulons),
                               y = N_Regulons, fill = -log10(Min_Pvalue))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient(low = "lightblue", high = "darkred") +
      theme_minimal(base_size = 10, base_family = "Helvetica") +
      labs(title = "Top KEGG Pathways Enriched Across Regulons",
           subtitle = "Data-driven discovery",
           x = "KEGG Pathway",
           y = "Number of Regulons",
           fill = "-log10(min FDR)")

    ggsave("publication_results/figures/Fig4_KEGG_Pathway_Overview.pdf",
           plot = p, width = 10, height = 12)

    cat("OK Saved: publication_results/figures/Fig4_KEGG_Pathway_Overview.pdf\n")
  }
}

cat("\n")

cat("PART 5: Generating supplementary tables...\n")
cat("----------------------------------------\n")

supp_table_s1 <- stat_results %>%
  select(Regulon, KW_pvalue, KW_FDR, Significant, Top_CellType, Top_Mean_AUC,
         matches("Mean_AUC_")) %>%
  arrange(KW_FDR)

target_counts <- regulon_targets %>%
  group_by(TF) %>%
  summarise(N_Targets = n())

supp_table_s1 <- supp_table_s1 %>%
  left_join(target_counts, by = c("Regulon" = "TF"))

supp_table_s2 <- stat_results

if (exists("go_combined")) {
  supp_table_s3 <- go_combined %>%
    select(Regulon, Ontology, ID, Description, GeneRatio, BgRatio,
           pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Ontology, p.adjust)
} else {
  supp_table_s3 <- data.frame(Note = "No significant GO enrichment found")
}

if (exists("kegg_combined")) {
  supp_table_s4 <- kegg_combined %>%
    select(Regulon, ID, Description, GeneRatio, BgRatio,
           pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(p.adjust)
} else {
  supp_table_s4 <- data.frame(Note = "No significant KEGG enrichment found")
}

supp_tables <- list(
  "S1_Regulon_Statistics" = supp_table_s1,
  "S2_Pairwise_Tests" = supp_table_s2,
  "S3_GO_Enrichment" = supp_table_s3,
  "S4_KEGG_Enrichment" = supp_table_s4
)

write_xlsx(supp_tables,
           "publication_results/supplementary_tables/SCENIC_Supplementary_Tables.xlsx")

cat("OK Saved: publication_results/supplementary_tables/SCENIC_Supplementary_Tables.xlsx\n")
cat("  - Table S1: Regulon statistics and cell type specificity\n")
cat("  - Table S2: Pairwise statistical tests (Wilcoxon)\n")
cat("  - Table S3: Complete GO enrichment (BP, MF, CC)\n")
cat("  - Table S4: Complete KEGG pathway enrichment\n\n")

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Output files generated:\n")
cat("\n1. Statistical Tests:\n")
cat("   - publication_results/statistical_tests/regulon_statistical_tests.tsv\n")

cat("\n2. Enrichment Analysis (PARALLEL MULTI-THREADED):\n")
if (exists("go_combined")) {
  cat("   - publication_results/enrichment_analysis/GO_enrichment_complete.tsv\n")
  cat(paste("     (", nrow(go_combined), "GO terms across BP/MF/CC ontologies)\n"))
}
if (exists("kegg_combined")) {
  cat("   - publication_results/enrichment_analysis/KEGG_enrichment_complete.tsv\n")
  cat(paste("     (", nrow(kegg_combined), "KEGG pathways)\n"))
}

cat("\n3. Publication Figures:\n")
cat("   - publication_results/figures/Fig1_Regulon_Activity_Heatmap.pdf\n")
if (file.exists("publication_results/figures/Fig2_Regulon_Specificity_Scores.pdf")) {
  cat("   - publication_results/figures/Fig2_Regulon_Specificity_Scores.pdf\n")
}
if (exists("go_combined") && nrow(go_combined) > 0) {
  cat("   - publication_results/figures/Fig3_GO_Enrichment_Dotplot.pdf\n")
}
if (exists("kegg_combined") && nrow(kegg_combined) > 0) {
  cat("   - publication_results/figures/Fig4_KEGG_Pathway_Overview.pdf\n")
}

cat("\n4. Supplementary Tables:\n")
cat("   - publication_results/supplementary_tables/SCENIC_Supplementary_Tables.xlsx\n")

cat("\n================================================================================\n")
cat(paste("PARALLEL PROCESSING:", n_cores, "cores used\n"))
cat("Analysis Strategy: UNBIASED DATA-DRIVEN DISCOVERY\n")
cat("No pre-defined pathway lists - all enrichments discovered from data\n")
cat(paste("Completed:", format(Sys.time()), "\n"))
cat("================================================================================\n")

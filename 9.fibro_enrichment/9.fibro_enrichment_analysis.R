rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(GOplot)
  library(ggthemes)
  library(msigdbr)
  library(GSVA)
  library(tidyverse)
  library(pheatmap)
  library(grid)
  library(gridExtra)
  library(reshape2)
})

dir.create("enrichment_result/GO", recursive = TRUE, showWarnings = FALSE)
dir.create("enrichment_result/KEGG", recursive = TRUE, showWarnings = FALSE)
dir.create("enrichment_result/GSVA", recursive = TRUE, showWarnings = FALSE)
dir.create("enrichment_result/combined", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

write_log <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp, msg))
  cat(paste(timestamp, msg, "\n"),
      file = "logs/enrichment_analysis.log",
      append = TRUE)
}

write_log("========================================")
write_log("Fibroblast Enrichment Analysis Started")
write_log("========================================")

write_log("\n1. Loading DEG data...")

subtypes <- c("ProInflamF", "SecRetF", "SecPapF")

deg_data_all <- list()
for (subtype in subtypes) {
  file_path <- paste0("../8.fibro_DEG/data/DEG_", subtype, "_CF_vs_Normal.csv")
  if (file.exists(file_path)) {
    deg_data_all[[subtype]] <- read.csv(file_path)
    write_log(sprintf("  OK Loaded %s: %d genes", subtype, nrow(deg_data_all[[subtype]])))
  } else {
    write_log(sprintf("  ERROR File not found: %s", file_path))
  }
}

write_log("\n2. Filtering DEGs...")

p_adj_cutoff <- 0.01
log2fc_cutoff <- 0.25

deg_gene_lists <- list()

for (subtype in subtypes) {
  deg <- deg_data_all[[subtype]]

  up_genes <- deg %>%
    filter(p_val_adj < p_adj_cutoff, avg_log2FC > log2fc_cutoff) %>%
    pull(gene)

  down_genes <- deg %>%
    filter(p_val_adj < p_adj_cutoff, avg_log2FC < -log2fc_cutoff) %>%
    pull(gene)

  deg_gene_lists[[paste0(subtype, "_Up")]] <- up_genes
  deg_gene_lists[[paste0(subtype, "_Down")]] <- down_genes

  write_log(sprintf("  %s: Up=%d, Down=%d",
                    subtype, length(up_genes), length(down_genes)))
}

write_log("\n3. Converting gene IDs...")

all_genes <- unique(unlist(deg_gene_lists))
write_log(sprintf("  Total unique genes: %d", length(all_genes)))

entrezIDs <- bitr(all_genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db,
                  drop = TRUE)

write_log(sprintf("  Successfully converted: %d genes", nrow(entrezIDs)))
write_log(sprintf("  Failed to convert: %d genes", length(all_genes) - nrow(entrezIDs)))

deg_entrez_lists <- lapply(deg_gene_lists, function(genes) {
  matched_ids <- entrezIDs[entrezIDs$SYMBOL %in% genes, "ENTREZID"]
  return(matched_ids)
})

deg_entrez_lists <- deg_entrez_lists[sapply(deg_entrez_lists, length) > 0]

write_log(sprintf("  Final gene lists: %d groups", length(deg_entrez_lists)))
for (name in names(deg_entrez_lists)) {
  write_log(sprintf("    - %s: %d genes", name, length(deg_entrez_lists[[name]])))
}

write_log("\n========================================")
write_log("4. GO Enrichment Analysis (BP)")
write_log("========================================")

go_results <- compareCluster(
  deg_entrez_lists,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

write_log(sprintf("  GO enrichment completed: %d terms found", nrow(go_results@compareClusterResult)))

write.csv(go_results@compareClusterResult,
          "enrichment_result/GO/GO_enrichment_results.csv",
          row.names = FALSE)
write_log("  OK Saved: enrichment_result/GO/GO_enrichment_results.csv")

saveRDS(go_results, "enrichment_result/GO/GO_results.rds")

write_log("\n  Creating GO visualizations...")
pdf("enrichment_result/GO/GO_dotplot.pdf", width = 14, height = 10)
p_go_dot <- dotplot(go_results, showCategory = 10) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "GO Biological Process Enrichment",
       subtitle = "CF vs Normal in Fibroblast Subtypes")
print(p_go_dot)
dev.off()
write_log("  OK Saved: enrichment_result/GO/GO_dotplot.pdf")

pdf("enrichment_result/GO/GO_barplot.pdf", width = 12, height = 8)

go_summary <- go_results@compareClusterResult %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  summarise(Count = n())

p_go_bar <- ggplot(go_summary, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(title = "Number of Enriched GO Terms",
       subtitle = "p.adjust < 0.05",
       x = "Group", y = "Count of GO Terms")
print(p_go_bar)
dev.off()
write_log("  OK Saved: enrichment_result/GO/GO_barplot.pdf")

write_log("\n========================================")
write_log("5. KEGG Pathway Enrichment")
write_log("========================================")

kegg_results <- compareCluster(
  deg_entrez_lists,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

write_log(sprintf("  KEGG enrichment completed: %d pathways found",
                  nrow(kegg_results@compareClusterResult)))

write.csv(kegg_results@compareClusterResult,
          "enrichment_result/KEGG/KEGG_enrichment_results.csv",
          row.names = FALSE)
write_log("  OK Saved: enrichment_result/KEGG/KEGG_enrichment_results.csv")

saveRDS(kegg_results, "enrichment_result/KEGG/KEGG_results.rds")

write_log("\n  Creating KEGG visualizations...")
pdf("enrichment_result/KEGG/KEGG_dotplot.pdf", width = 14, height = 10)
p_kegg_dot <- dotplot(kegg_results, showCategory = 10) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "KEGG Pathway Enrichment",
       subtitle = "CF vs Normal in Fibroblast Subtypes")
print(p_kegg_dot)
dev.off()
write_log("  OK Saved: enrichment_result/KEGG/KEGG_dotplot.pdf")

write_log("\n========================================")
write_log("6. GSVA Analysis")
write_log("========================================")

write_log("  Loading Seurat object...")
sce <- readRDS("../7.0.fibro_extrction/07_fibroblast/fibroblast_annotated.rds")
write_log(sprintf("  OK Loaded %d cells", ncol(sce)))

write_log("  Loading MSigDB C2 gene sets...")
genesets <- msigdbr(species = "Homo sapiens", category = "C2")
genesets <- subset(genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
write_log(sprintf("  OK Loaded %d gene sets", length(genesets)))

write_log("  Calculating average expression...")

sce$group <- paste0(sce$celltype, "_", sce$condition)
Idents(sce) <- sce$group

expr <- AverageExpression(sce, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0, ]
write_log(sprintf("  OK Expression matrix: %d genes x %d groups",
                  nrow(expr), ncol(expr)))

write_log("  Running GSVA (this may take several minutes)...")

gsva_param <- gsvaParam(as.matrix(expr), genesets)

gsva_res <- gsva(gsva_param, verbose = FALSE)
write_log(sprintf("  OK GSVA completed: %d gene sets", nrow(gsva_res)))

saveRDS(gsva_res, "enrichment_result/GSVA/GSVA_results.rds")
gsva_df <- data.frame(Genesets = rownames(gsva_res), gsva_res, check.names = FALSE)
write.csv(gsva_df, "enrichment_result/GSVA/GSVA_scores.csv", row.names = FALSE)
write_log("  OK Saved: enrichment_result/GSVA/GSVA_scores.csv")

write_log("\n  Renaming columns for visualization...")
colname_mapping <- c(
  "pro-inflammatory-CF" = "ProInflamF_CF",
  "pro-inflammatory-Normal" = "ProInflamF_Normal",
  "secretory-reticular-CF" = "SecRetF_CF",
  "secretory-reticular-Normal" = "SecRetF_Normal",
  "secretory-papillary-CF" = "SecPapF_CF",
  "secretory-papillary-Normal" = "SecPapF_Normal",
  "mesenchymal-CF" = "MesF_CF",
  "mesenchymal-Normal" = "MesF_Normal"
)

for (old_name in names(colname_mapping)) {
  if (old_name %in% colnames(gsva_res)) {
    colnames(gsva_res)[colnames(gsva_res) == old_name] <- colname_mapping[old_name]
  }
}
write_log(sprintf("  OK Renamed columns: %s", paste(colname_mapping, collapse = ", ")))

write_log("\n  Creating GSVA visualizations...")
gsva_var <- apply(gsva_res, 1, var)
top_gs <- names(sort(gsva_var, decreasing = TRUE)[1:30])
gsva_top <- gsva_res[top_gs, ]

pdf("enrichment_result/GSVA/GSVA_heatmap.pdf", width = 12, height = 10)
pheatmap(gsva_top,
         show_colnames = TRUE,
         show_rownames = TRUE,
         scale = "row",
         angle_col = 45,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "GSVA Enrichment Score - Top 30 Variable Gene Sets")
dev.off()
write_log("  OK Saved: enrichment_result/GSVA/GSVA_heatmap.pdf")

gsva_long <- melt(gsva_top)
colnames(gsva_long) <- c("Geneset", "Group", "Score")

pdf("enrichment_result/GSVA/GSVA_bubble.pdf", width = 14, height = 10)
p_gsva_bubble <- ggplot(gsva_long, aes(x = Group, y = Geneset, size = abs(Score), color = Score)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 6)) +
  scale_color_gradient2(low = "#427183", mid = "white", high = "#D2D470", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8)) +
  labs(title = "GSVA Enrichment Score - Bubble Plot",
       subtitle = "Top 30 Variable Gene Sets",
       x = "Group", y = "Gene Set", size = "Abs(GSVA Score)", color = "GSVA Score")
print(p_gsva_bubble)
dev.off()
write_log("  OK Saved: enrichment_result/GSVA/GSVA_bubble.pdf")

write_log("\n  Calculating differential GSVA scores...")
diff_gsva_list <- list()

celltype_mapping <- list(
  "pro-inflammatory" = "ProInflamF",
  "secretory-reticular" = "SecRetF",
  "secretory-papillary" = "SecPapF",
  "mesenchymal" = "MesF"
)

for (celltype_full in names(celltype_mapping)) {
  subtype <- celltype_mapping[[celltype_full]]
  cf_col <- paste0(celltype_full, "-CF")
  normal_col <- paste0(celltype_full, "-Normal")

  if (cf_col %in% colnames(gsva_res) && normal_col %in% colnames(gsva_res)) {
    diff_gsva <- gsva_res[, cf_col] - gsva_res[, normal_col]
    diff_gsva_list[[subtype]] <- diff_gsva
    write_log(sprintf("    OK %s: %s vs %s", subtype, cf_col, normal_col))
  }
}

diff_gsva_df <- as.data.frame(diff_gsva_list)
diff_gsva_df$Geneset <- rownames(diff_gsva_df)
write.csv(diff_gsva_df,
          "enrichment_result/GSVA/GSVA_differential.csv",
          row.names = FALSE)
write_log("  OK Saved: enrichment_result/GSVA/GSVA_differential.csv")

write_log("\n========================================")
write_log("7. Creating Summary Report")
write_log("========================================")

summary_stats <- data.frame(
  Analysis = c("GO Enrichment", "KEGG Pathways", "GSVA Gene Sets"),
  Total_Terms = c(
    nrow(go_results@compareClusterResult),
    nrow(kegg_results@compareClusterResult),
    nrow(gsva_res)
  ),
  Significant = c(
    sum(go_results@compareClusterResult$p.adjust < 0.05),
    sum(kegg_results@compareClusterResult$p.adjust < 0.05),
    NA
  )
)

write.csv(summary_stats,
          "enrichment_result/combined/enrichment_summary.csv",
          row.names = FALSE)
write_log("  OK Saved: enrichment_result/combined/enrichment_summary.csv")

write_log("\nEnrichment Analysis Summary:")
write_log(sprintf("  GO Terms: %d total, %d significant (p.adj < 0.05)",
                  summary_stats$Total_Terms[1], summary_stats$Significant[1]))
write_log(sprintf("  KEGG Pathways: %d total, %d significant (p.adj < 0.05)",
                  summary_stats$Total_Terms[2], summary_stats$Significant[2]))
write_log(sprintf("  GSVA Gene Sets: %d scored", summary_stats$Total_Terms[3]))

write_log("\n========================================")
write_log(" Enrichment Analysis Completed Successfully!")
write_log("========================================")
write_log("\nOutput files generated:")
write_log("  - enrichment_result/GO/")
write_log("  - enrichment_result/KEGG/")
write_log("  - enrichment_result/GSVA/")
write_log("  - enrichment_result/combined/")
write_log("  - logs/enrichment_analysis.log")
write_log("========================================")

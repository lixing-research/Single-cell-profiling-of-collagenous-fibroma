#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

log_file <- "logs/04_leading_edge.log"
cat("", file = log_file)
log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s\n", timestamp, msg)
  cat(log_line)
  cat(log_line, file = log_file, append = TRUE)
}

log_msg("================================================================================")
log_msg("04:Leading Edge")
log_msg("================================================================================")

log_msg("\n1. GSEA...")

gsea_hallmark <- readRDS("input.data.gsea_results/gsea_hallmark_results.rds")
log_msg(sprintf("   Hallmark: %d", nrow(gsea_hallmark)))

gsea_reactome <- readRDS("input.data.gsea_results/gsea_reactome_results.rds")
log_msg(sprintf("   Reactome: %d", nrow(gsea_reactome)))

pseudobulk <- readRDS("input.data.gsea_results/pseudobulk_overall.rds")
log_msg("   Pseudo-bulk")

log_msg("\n2. (Leading Edge)...")

key_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_ANGIOGENESIS"
)

log_msg(sprintf("   %dHallmark", length(key_pathways)))

log_msg("\n3. Leading Edge...")

leading_edge_list <- list()
all_le_genes <- character(0)

for (pathway in key_pathways) {
  pw_data <- gsea_hallmark[gsea_hallmark$pathway == pathway, ]

  if (nrow(pw_data) == 0) {
    log_msg(sprintf("   :  %s", pathway))
    next
  }

  le_genes <- unlist(pw_data$leadingEdge)
  leading_edge_list[[pathway]] <- le_genes

  pathway_short <- gsub("HALLMARK_", "", pathway)
  log_msg(sprintf("   %s:", pathway_short))
  log_msg(sprintf("     NES=%.3f, FDR=%.4f", pw_data$NES, pw_data$padj))
  log_msg(sprintf("     Leading Edge genes: %d", length(le_genes)))

  all_le_genes <- c(all_le_genes, le_genes)
}

unique_le_genes <- unique(all_le_genes)
log_msg(sprintf("\n   Leading Edge(): %d", length(unique_le_genes)))

log_msg("\n4. Leading EdgeCF vs Normal...")

CF_expr <- pseudobulk$CF
Normal_expr <- pseudobulk$Normal

le_genes_in_data <- intersect(unique_le_genes, names(CF_expr))
log_msg(sprintf("   %d/%dLE", length(le_genes_in_data), length(unique_le_genes)))

le_fc <- log2((CF_expr[le_genes_in_data] + 0.01) / (Normal_expr[le_genes_in_data] + 0.01))
le_fc_sorted <- sort(le_fc, decreasing = TRUE)

log_msg(sprintf("   LElog2FC: %.3f to %.3f", min(le_fc), max(le_fc)))
log_msg(sprintf("   (log2FC>0.5): %d", sum(le_fc > 0.5)))
log_msg(sprintf("   (log2FC<-0.5): %d", sum(le_fc < -0.5)))

log_msg("\n5. Leading Edge...")

oxphos_genes <- leading_edge_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]
if (!is.null(oxphos_genes)) {
  oxphos_df <- data.frame(
    gene = oxphos_genes,
    log2FC = le_fc[oxphos_genes],
    CF_expr = CF_expr[oxphos_genes],
    Normal_expr = Normal_expr[oxphos_genes]
  )
  oxphos_df <- oxphos_df[order(-oxphos_df$log2FC), ]
  write.csv(oxphos_df, "results/leading_edge_oxphos.csv", row.names = FALSE)
  log_msg(sprintf("   OXPHOS leading edge: %d genes", nrow(oxphos_df)))
}

upr_genes <- leading_edge_list[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]]
if (!is.null(upr_genes)) {
  upr_df <- data.frame(
    gene = upr_genes,
    log2FC = le_fc[upr_genes],
    CF_expr = CF_expr[upr_genes],
    Normal_expr = Normal_expr[upr_genes]
  )
  upr_df <- upr_df[order(-upr_df$log2FC), ]
  write.csv(upr_df, "results/leading_edge_upr.csv", row.names = FALSE)
  log_msg(sprintf("   UPR leading edge: %d genes", nrow(upr_df)))
}

emt_genes <- leading_edge_list[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]
if (!is.null(emt_genes)) {
  emt_df <- data.frame(
    gene = emt_genes,
    log2FC = le_fc[emt_genes],
    CF_expr = CF_expr[emt_genes],
    Normal_expr = Normal_expr[emt_genes]
  )
  emt_df <- emt_df[order(-emt_df$log2FC), ]
  write.csv(emt_df, "results/leading_edge_emt.csv", row.names = FALSE)
  log_msg(sprintf("   EMT leading edge: %d genes", nrow(emt_df)))
}

glycolysis_genes <- leading_edge_list[["HALLMARK_GLYCOLYSIS"]]
if (!is.null(glycolysis_genes)) {
  glyc_df <- data.frame(
    gene = glycolysis_genes,
    log2FC = le_fc[glycolysis_genes],
    CF_expr = CF_expr[glycolysis_genes],
    Normal_expr = Normal_expr[glycolysis_genes]
  )
  glyc_df <- glyc_df[order(-glyc_df$log2FC), ]
  write.csv(glyc_df, "results/leading_edge_glycolysis.csv", row.names = FALSE)
  log_msg(sprintf("   GLYCOLYSIS leading edge: %d genes", nrow(glyc_df)))
}

log_msg("\n6. Leading Edge...")

le_summary <- data.frame()

for (pathway_name in names(leading_edge_list)) {
  genes <- leading_edge_list[[pathway_name]]
  pathway_short <- gsub("HALLMARK_", "", pathway_name)

  for (gene in genes) {
    if (gene %in% names(le_fc)) {
      le_summary <- rbind(le_summary, data.frame(
        pathway = pathway_short,
        gene = gene,
        log2FC = le_fc[gene],
        CF_expr = CF_expr[gene],
        Normal_expr = Normal_expr[gene]
      ))
    }
  }
}

le_summary <- le_summary[order(-abs(le_summary$log2FC)), ]
write.csv(le_summary, "results/leading_edge_genes_summary.csv", row.names = FALSE)
log_msg(sprintf("   : %d(x)", nrow(le_summary)))
log_msg("   : results/leading_edge_genes_summary.csv")

log_msg("\n7. qPCR...")

qpcr_candidates <- data.frame()

if (!is.null(oxphos_genes)) {
  oxphos_top <- head(oxphos_df, 10)
  oxphos_top$pathway <- "OXPHOS"
  oxphos_top$priority <- "High"
  qpcr_candidates <- rbind(qpcr_candidates, oxphos_top[, c("pathway", "gene", "log2FC", "priority")])
}

if (!is.null(upr_genes)) {
  upr_top <- head(upr_df, 8)
  upr_top$pathway <- "UPR"
  upr_top$priority <- "High"
  qpcr_candidates <- rbind(qpcr_candidates, upr_top[, c("pathway", "gene", "log2FC", "priority")])
}

if (!is.null(emt_genes)) {
  emt_top <- head(emt_df, 8)
  emt_top$pathway <- "EMT"
  emt_top$priority <- "High"
  qpcr_candidates <- rbind(qpcr_candidates, emt_top[, c("pathway", "gene", "log2FC", "priority")])
}

if (!is.null(glycolysis_genes)) {
  glyc_top <- head(glyc_df, 6)
  glyc_top$pathway <- "GLYCOLYSIS"
  glyc_top$priority <- "Medium"
  qpcr_candidates <- rbind(qpcr_candidates, glyc_top[, c("pathway", "gene", "log2FC", "priority")])
}

qpcr_candidates <- qpcr_candidates[!duplicated(qpcr_candidates$gene), ]
qpcr_candidates <- qpcr_candidates[order(-qpcr_candidates$log2FC), ]

write.csv(qpcr_candidates, "results/qpcr_candidate_genes.csv", row.names = FALSE)
log_msg(sprintf("   qPCR: %d()", nrow(qpcr_candidates)))
log_msg("   : results/qpcr_candidate_genes.csv")

log_msg("\n8. Leading Edge...")

heatmap_genes <- c()
if (!is.null(oxphos_genes)) heatmap_genes <- c(heatmap_genes, head(oxphos_df$gene, 10))
if (!is.null(upr_genes)) heatmap_genes <- c(heatmap_genes, head(upr_df$gene, 10))
if (!is.null(emt_genes)) heatmap_genes <- c(heatmap_genes, head(emt_df$gene, 10))
if (!is.null(glycolysis_genes)) heatmap_genes <- c(heatmap_genes, head(glyc_df$gene, 10))

heatmap_genes <- unique(heatmap_genes)
log_msg(sprintf("   : %d", length(heatmap_genes)))

expr_matrix <- cbind(
  CF = CF_expr[heatmap_genes],
  Normal = Normal_expr[heatmap_genes]
)

expr_matrix_log <- log2(expr_matrix + 0.01)
expr_matrix_scaled <- t(scale(t(expr_matrix_log)))

gene_annotation <- data.frame(
  gene = heatmap_genes,
  pathway = sapply(heatmap_genes, function(g) {
    if (g %in% oxphos_genes) return("OXPHOS")
    if (g %in% upr_genes) return("UPR")
    if (g %in% emt_genes) return("EMT")
    if (g %in% glycolysis_genes) return("GLYCOLYSIS")
    return("Other")
  })
)
rownames(gene_annotation) <- gene_annotation$gene

pdf("results/leading_edge_heatmap.pdf", width = 6, height = 12)
pheatmap(expr_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),
         annotation_row = gene_annotation["pathway"],
         main = "Leading Edge Genes Expression\n(CF vs Normal)",
         border_color = NA)
dev.off()

log_msg("   : results/leading_edge_heatmap.pdf")

png("results/leading_edge_heatmap.png", width = 1800, height = 3600, res = 300)
pheatmap(expr_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),
         annotation_row = gene_annotation["pathway"],
         main = "Leading Edge Genes Expression (CF vs Normal)",
         border_color = NA)
dev.off()

log_msg("   PNG: results/leading_edge_heatmap.png")

log_msg("\n9. Leading Edge...")

pathway_stats <- data.frame()
for (pathway_name in names(leading_edge_list)) {
  genes <- leading_edge_list[[pathway_name]]
  pathway_short <- gsub("HALLMARK_", "", pathway_name)

  genes_in_data <- intersect(genes, names(le_fc))
  mean_fc <- mean(le_fc[genes_in_data])
  median_fc <- median(le_fc[genes_in_data])
  upregulated <- sum(le_fc[genes_in_data] > 0.5)

  pathway_stats <- rbind(pathway_stats, data.frame(
    pathway = pathway_short,
    n_genes = length(genes),
    n_in_data = length(genes_in_data),
    mean_log2FC = mean_fc,
    median_log2FC = median_fc,
    n_upregulated = upregulated
  ))
}

pathway_stats <- pathway_stats[order(-pathway_stats$mean_log2FC), ]
write.csv(pathway_stats, "results/leading_edge_pathway_stats.csv", row.names = FALSE)
log_msg("   : results/leading_edge_pathway_stats.csv")

log_msg("\n10. Leading Edge...")

log_msg("\n   === OXPHOS Leading Edge (Top 10) ===")
if (!is.null(oxphos_genes)) {
  for (i in 1:min(10, nrow(oxphos_df))) {
    log_msg(sprintf("   %2d. %s (log2FC=%.3f)", i, oxphos_df$gene[i], oxphos_df$log2FC[i]))
  }
}

log_msg("\n   === UPR Leading Edge (Top 8) ===")
if (!is.null(upr_genes)) {
  for (i in 1:min(8, nrow(upr_df))) {
    log_msg(sprintf("   %2d. %s (log2FC=%.3f)", i, upr_df$gene[i], upr_df$log2FC[i]))
  }
}

log_msg("\n   === EMT Leading Edge (Top 8) ===")
if (!is.null(emt_genes)) {
  for (i in 1:min(8, nrow(emt_df))) {
    log_msg(sprintf("   %2d. %s (log2FC=%.3f)", i, emt_df$gene[i], emt_df$log2FC[i]))
  }
}

log_msg("\n11. session...")
session_file <- "logs/session_info.txt"
writeLines(capture.output(sessionInfo()), session_file)
log_msg(sprintf("   Session info: %s", session_file))

log_msg("\n================================================================================")
log_msg("04")
log_msg("================================================================================")
log_msg(sprintf("OK %d", length(leading_edge_list)))
log_msg(sprintf("OK Leading Edge: %d()", length(unique_le_genes)))
log_msg(sprintf("OK qPCR: %d", nrow(qpcr_candidates)))
log_msg(sprintf("OK : %d", length(heatmap_genes)))
log_msg("")
log_msg(":")
log_msg("  - results/leading_edge_genes_summary.csv")
log_msg("  - results/leading_edge_oxphos.csv")
log_msg("  - results/leading_edge_upr.csv")
log_msg("  - results/leading_edge_emt.csv")
log_msg("  - results/leading_edge_glycolysis.csv")
log_msg("  - results/leading_edge_heatmap.pdf/.png")
log_msg("  - results/qpcr_candidate_genes.csv")
log_msg("  - results/leading_edge_pathway_stats.csv")
log_msg("")
log_msg(": 05 - ")
log_msg("================================================================================")

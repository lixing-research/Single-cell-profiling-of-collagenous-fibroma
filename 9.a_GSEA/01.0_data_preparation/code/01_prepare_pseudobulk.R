#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

log_file <- "logs/01_data_prep.log"
cat("", file = log_file)
log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s\n", timestamp, msg)
  cat(log_line)
  cat(log_line, file = log_file, append = TRUE)
}

log_msg("================================================================================")
log_msg("01:Pseudo-bulk")
log_msg("================================================================================")

log_msg("\n1. ...")

input_file <- "input.data.fibroblasts/fibroblast_annotated.rds"
if (!file.exists(input_file)) {
}

sce <- readRDS(input_file)
log_msg(sprintf("   : %s", input_file))
log_msg(sprintf("   : %d", ncol(sce)))
log_msg(sprintf("   : %d", nrow(sce)))

log_msg("\n2. ...")

meta_cols <- colnames(sce@meta.data)
log_msg(sprintf("   : %s", paste(meta_cols, collapse=", ")))

if ("subtype" %in% meta_cols) {
  log_msg("   'subtype'")
} else if ("celltype_abbr" %in% meta_cols) {
  log_msg("   'celltype_abbr'")
  sce$subtype <- sce$celltype_abbr
} else if ("celltype" %in% meta_cols) {
  log_msg("   'celltype'")
  sce$subtype <- sce$celltype
} else if ("seurat_clusters" %in% meta_cols) {
  log_msg("   'seurat_clusters'")
  sce$subtype <- sce$seurat_clusters
} else {
}

log_msg("\n   :")
if ("sample" %in% meta_cols) {
  sample_table <- table(sce$sample)
  for (i in 1:length(sample_table)) {
    log_msg(sprintf("     %s: %d cells", names(sample_table)[i], sample_table[i]))
  }
}

log_msg("\n   :")
subtype_table <- table(sce$subtype)
for (i in 1:length(subtype_table)) {
  log_msg(sprintf("     %s: %d cells", names(subtype_table)[i], subtype_table[i]))
}

if ("group" %in% meta_cols) {
  group_col <- "group"
  log_msg(sprintf("   '%s'CF/Normal", group_col))
  sce$analysis_group <- sce$group
} else if ("condition" %in% meta_cols) {
  group_col <- "condition"
  log_msg(sprintf("   '%s'CF/Normal", group_col))
  sce$analysis_group <- sce$condition
} else if ("sample" %in% meta_cols) {
  group_col <- "sample"
  log_msg(sprintf("   '%s'CF/Normal", group_col))
  sce$analysis_group <- ifelse(grepl("CF|Tumor", sce$sample, ignore.case = TRUE),
                               "CF", "Normal")
} else {
}

group_counts <- table(sce$analysis_group)
log_msg(sprintf("   CF cells: %d", group_counts["CF"]))
log_msg(sprintf("   Normal cells: %d", group_counts["Normal"]))

log_msg("\n3. Pseudo-bulk...")

if ("SCT" %in% names(sce@assays)) {
  expr_data <- GetAssayData(sce, assay = "SCT", slot = "data")
  log_msg("   SCT")
} else if ("RNA" %in% names(sce@assays)) {
  expr_data <- GetAssayData(sce, assay = "RNA", slot = "data")
  log_msg("   RNA")
} else {
  expr_data <- GetAssayData(sce, slot = "data")
  log_msg("   ")
}

log_msg(sprintf("   : %d genes x %d cells", nrow(expr_data), ncol(expr_data)))

CF_cells <- which(sce$analysis_group == "CF")
Normal_cells <- which(sce$analysis_group == "Normal")

log_msg(sprintf("   CF: %d", length(CF_cells)))
log_msg(sprintf("   Normal: %d", length(Normal_cells)))

CF_bulk <- Matrix::rowMeans(expr_data[, CF_cells])
Normal_bulk <- Matrix::rowMeans(expr_data[, Normal_cells])

log_msg(sprintf("   CF pseudo-bulk: %d genes", length(CF_bulk)))
log_msg(sprintf("   Normal pseudo-bulk: %d genes", length(Normal_bulk)))

min_expr <- 0.01
genes_keep <- (CF_bulk > min_expr) | (Normal_bulk > min_expr)
log_msg(sprintf("   (min_expr=%.3f): %d -> %d genes",
                min_expr, length(genes_keep), sum(genes_keep)))

CF_bulk <- CF_bulk[genes_keep]
Normal_bulk <- Normal_bulk[genes_keep]

pseudobulk_overall <- list(
  CF = CF_bulk,
  Normal = Normal_bulk,
  genes = names(CF_bulk),
  n_cells = c(CF = length(CF_cells), Normal = length(Normal_cells))
)

saveRDS(pseudobulk_overall, "results/pseudobulk_overall.rds")
log_msg("   : results/pseudobulk_overall.rds")

log_msg("\n4. Log2 Fold Change...")

pseudocount <- 0.01
log2FC <- log2((CF_bulk + pseudocount) / (Normal_bulk + pseudocount))

ranked_genes <- sort(log2FC, decreasing = TRUE)

log_msg(sprintf("   %dlog2FC", length(log2FC)))
log_msg(sprintf("   log2FC: %.3f to %.3f", min(log2FC), max(log2FC)))
log_msg(sprintf("   (log2FC>0.5): %d", sum(log2FC > 0.5)))
log_msg(sprintf("   (log2FC<-0.5): %d", sum(log2FC < -0.5)))

saveRDS(ranked_genes, "results/ranked_genes_overall.rds")
log_msg("   : results/ranked_genes_overall.rds")

ranked_df <- data.frame(
  gene = names(ranked_genes),
  log2FC = ranked_genes,
  CF_expr = CF_bulk[names(ranked_genes)],
  Normal_expr = Normal_bulk[names(ranked_genes)]
)
write.csv(ranked_df, "results/ranked_genes_overall.csv", row.names = FALSE)
log_msg("   : results/ranked_genes_overall.csv")

log_msg("\n5. Pseudo-bulk...")

subtypes <- unique(sce$subtype)
log_msg(sprintf("   %d: %s", length(subtypes), paste(subtypes, collapse=", ")))

pseudobulk_by_subtype <- list()

for (st in subtypes) {
  log_msg(sprintf("\n   : %s", st))

  st_cells <- which(sce$subtype == st)
  log_msg(sprintf("     : %d", length(st_cells)))

  st_CF_cells <- intersect(st_cells, CF_cells)
  st_Normal_cells <- intersect(st_cells, Normal_cells)

  log_msg(sprintf("     CF cells: %d", length(st_CF_cells)))
  log_msg(sprintf("     Normal cells: %d", length(st_Normal_cells)))

  if (length(st_CF_cells) < 50 || length(st_Normal_cells) < 50) {
    log_msg(sprintf("     :,"))
    next
  }

  st_CF_bulk <- Matrix::rowMeans(expr_data[, st_CF_cells])
  st_Normal_bulk <- Matrix::rowMeans(expr_data[, st_Normal_cells])

  st_genes_keep <- (st_CF_bulk > min_expr) | (st_Normal_bulk > min_expr)
  st_CF_bulk <- st_CF_bulk[st_genes_keep]
  st_Normal_bulk <- st_Normal_bulk[st_genes_keep]

  log_msg(sprintf("     : %d", length(st_CF_bulk)))

  pseudobulk_by_subtype[[st]] <- list(
    CF = st_CF_bulk,
    Normal = st_Normal_bulk,
    genes = names(st_CF_bulk),
    n_cells = c(CF = length(st_CF_cells), Normal = length(st_Normal_cells))
  )
}

saveRDS(pseudobulk_by_subtype, "results/pseudobulk_by_subtype.rds")
log_msg(sprintf("\n   %dpseudo-bulk", length(pseudobulk_by_subtype)))
log_msg("   : results/pseudobulk_by_subtype.rds")

log_msg("\n6. ...")

qc_report <- "results/data_qc_report.txt"
cat("", file = qc_report)

write_qc <- function(text) {
  cat(text, "\n", file = qc_report, append = TRUE)
}

write_qc("================================================================================")
write_qc("GSEA - ")
write_qc("================================================================================")
write_qc(sprintf(": %s", Sys.time()))
write_qc("")

write_qc("1. ")
write_qc(sprintf("   : %s", input_file))
write_qc(sprintf("   : %d", ncol(sce)))
write_qc(sprintf("   : %d", nrow(sce)))
write_qc("")

write_qc("2. ")
write_qc(sprintf("   CF: %d (%.1f%%)", group_counts["CF"],
                group_counts["CF"]/sum(group_counts)*100))
write_qc(sprintf("   Normal: %d (%.1f%%)", group_counts["Normal"],
                group_counts["Normal"]/sum(group_counts)*100))
write_qc("")

write_qc("3. ")
for (st in names(subtype_table)) {
  write_qc(sprintf("   %s: %d cells (%.1f%%)", st, subtype_table[st],
                  subtype_table[st]/sum(subtype_table)*100))
}
write_qc("")

write_qc("4. Pseudo-bulk")
write_qc(sprintf("   : %d", length(CF_bulk)))
write_qc(sprintf("   log2FC: %.3f to %.3f", min(log2FC), max(log2FC)))
write_qc(sprintf("   (log2FC>0.5): %d", sum(log2FC > 0.5)))
write_qc(sprintf("   (log2FC<-0.5): %d", sum(log2FC < -0.5)))
write_qc("")

write_qc("5. pseudo-bulk")
for (st in names(pseudobulk_by_subtype)) {
  pb <- pseudobulk_by_subtype[[st]]
  write_qc(sprintf("   %s:", st))
  write_qc(sprintf("     CF cells: %d", pb$n_cells["CF"]))
  write_qc(sprintf("     Normal cells: %d", pb$n_cells["Normal"]))
  write_qc(sprintf("     : %d", length(pb$genes)))
}
write_qc("")

write_qc("6. (Top 10)")
top_up <- head(ranked_genes, 10)
for (i in 1:length(top_up)) {
  write_qc(sprintf("   %d. %s (log2FC=%.3f)", i, names(top_up)[i], top_up[i]))
}
write_qc("")

write_qc("7. (Top 10)")
top_down <- tail(ranked_genes, 10)
for (i in length(top_down):1) {
  write_qc(sprintf("   %d. %s (log2FC=%.3f)", 11-i, names(top_down)[i], top_down[i]))
}
write_qc("")

write_qc("================================================================================")
write_qc("")
write_qc("================================================================================")

log_msg(sprintf("   : %s", qc_report))

log_msg("\n7. QC...")

pdf("results/qc_plots.pdf", width = 12, height = 8)

par(mfrow = c(2, 2))

hist(log2FC, breaks = 100, main = "Log2 Fold Change Distribution",
     xlab = "log2(CF/Normal)", col = "lightblue", border = "white")
abline(v = c(-0.5, 0.5), col = "red", lty = 2, lwd = 2)
text(-2, par("usr")[4]*0.9, sprintf("Down: %d", sum(log2FC < -0.5)), pos = 4)
text(1, par("usr")[4]*0.9, sprintf("Up: %d", sum(log2FC > 0.5)), pos = 4)

plot(log10(Normal_bulk + 0.01), log10(CF_bulk + 0.01),
     pch = 16, cex = 0.3, col = rgb(0,0,0,0.3),
     xlab = "log10(Normal expression)", ylab = "log10(CF expression)",
     main = "Pseudo-bulk Expression: CF vs Normal")
abline(0, 1, col = "red", lty = 2, lwd = 2)

barplot(group_counts, main = "Cell Counts by Group",
        ylab = "Number of Cells", col = c("CF" = "#E41A1C", "Normal" = "#377EB8"))

barplot(subtype_table, main = "Cell Counts by Subtype",
        ylab = "Number of Cells", las = 2, col = rainbow(length(subtype_table)))

dev.off()

log_msg("   QC: results/qc_plots.pdf")

log_msg("\n8. session...")

session_file <- "logs/session_info.txt"
writeLines(capture.output(sessionInfo()), session_file)
log_msg(sprintf("   Session info: %s", session_file))

log_msg("\n================================================================================")
log_msg("01")
log_msg("================================================================================")
log_msg(sprintf("OK : %d cells, %d genes", ncol(sce), nrow(sce)))
log_msg(sprintf("OK CF vs Normal: %d vs %d cells", group_counts["CF"], group_counts["Normal"]))
log_msg(sprintf("OK Pseudo-bulk: %d", length(CF_bulk)))
log_msg(sprintf("OK : %d", length(pseudobulk_by_subtype)))
log_msg("")
log_msg(":")
log_msg("  - results/pseudobulk_overall.rds")
log_msg("  - results/pseudobulk_by_subtype.rds")
log_msg("  - results/ranked_genes_overall.rds")
log_msg("  - results/ranked_genes_overall.csv")
log_msg("  - results/data_qc_report.txt")
log_msg("  - results/qc_plots.pdf")
log_msg("")
log_msg(": 02 - GSEA()")
log_msg("================================================================================")

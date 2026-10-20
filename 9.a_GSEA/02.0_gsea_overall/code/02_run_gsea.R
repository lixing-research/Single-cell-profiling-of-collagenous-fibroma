#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

log_file <- "logs/02_gsea.log"
cat("", file = log_file)
log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s\n", timestamp, msg)
  cat(log_line)
  cat(log_line, file = log_file, append = TRUE)
}

log_msg("================================================================================")
log_msg("02:GSEA - ")
log_msg("================================================================================")

log_msg("\n1. ranked genes...")

ranked_genes <- readRDS("input.data.pseudobulk/ranked_genes_overall.rds")
log_msg(sprintf("   : %d", length(ranked_genes)))
log_msg(sprintf("   log2FC: %.3f to %.3f", min(ranked_genes), max(ranked_genes)))

log_msg("\n2. ...")

log_msg("   Hallmark...")
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
log_msg(sprintf("   Hallmark: %d", length(unique(hallmark_gene_sets$gs_name))))

hallmark_list <- split(hallmark_gene_sets$gene_symbol, hallmark_gene_sets$gs_name)
log_msg(sprintf("   list: %d", length(hallmark_list)))

log_msg("\n   Reactome...")
reactome_gene_sets <- msigdbr(species = "Homo sapiens",
                              category = "C2",
                              subcategory = "CP:REACTOME")
log_msg(sprintf("   Reactome: %d", length(unique(reactome_gene_sets$gs_name))))

reactome_list <- split(reactome_gene_sets$gene_symbol, reactome_gene_sets$gs_name)
log_msg(sprintf("   list: %d", length(reactome_list)))

log_msg("\n3. fgsea - Hallmark...")
log_msg("   : nperm=10000, minSize=15, maxSize=500")

set.seed(42)

start_time <- Sys.time()
gsea_hallmark <- fgsea(pathways = hallmark_list,
                       stats = ranked_genes,
                       nperm = 10000,
                       minSize = 15,
                       maxSize = 500)
end_time <- Sys.time()

log_msg(sprintf("   ,: %.1f", as.numeric(difftime(end_time, start_time, units = "secs"))))
log_msg(sprintf("   %dHallmark", nrow(gsea_hallmark)))

gsea_hallmark <- gsea_hallmark[order(NES, decreasing = TRUE)]
saveRDS(gsea_hallmark, "results/gsea_hallmark_results.rds")

gsea_hallmark_csv <- gsea_hallmark %>%
  select(-leadingEdge) %>%
  as.data.frame()
write.csv(gsea_hallmark_csv, "results/gsea_hallmark_results.csv", row.names = FALSE)
log_msg("   : results/gsea_hallmark_results.rds/.csv")

sig_hallmark <- gsea_hallmark[padj < 0.05]
log_msg(sprintf("   (FDR<0.05): %d/%d", nrow(sig_hallmark), nrow(gsea_hallmark)))

sig_up <- sig_hallmark[NES > 0]
sig_down <- sig_hallmark[NES < 0]
log_msg(sprintf("     : %d", nrow(sig_up)))
log_msg(sprintf("     : %d", nrow(sig_down)))

log_msg("\n4. fgsea - Reactome...")
log_msg(sprintf("   %dReactome()...", length(reactome_list)))

start_time <- Sys.time()
gsea_reactome <- fgsea(pathways = reactome_list,
                       stats = ranked_genes,
                       nperm = 10000,
                       minSize = 15,
                       maxSize = 500)
end_time <- Sys.time()

log_msg(sprintf("   ,: %.1f", as.numeric(difftime(end_time, start_time, units = "secs"))))

gsea_reactome <- gsea_reactome[order(NES, decreasing = TRUE)]
saveRDS(gsea_reactome, "results/gsea_reactome_results.rds")

gsea_reactome_csv <- gsea_reactome %>%
  select(-leadingEdge) %>%
  as.data.frame()
write.csv(gsea_reactome_csv, "results/gsea_reactome_results.csv", row.names = FALSE)
log_msg("   : results/gsea_reactome_results.rds/.csv")

sig_reactome <- gsea_reactome[padj < 0.05]
log_msg(sprintf("   (FDR<0.05): %d/%d", nrow(sig_reactome), nrow(gsea_reactome)))

sig_reactome_up <- sig_reactome[NES > 0]
sig_reactome_down <- sig_reactome[NES < 0]
log_msg(sprintf("     : %d", nrow(sig_reactome_up)))
log_msg(sprintf("     : %d", nrow(sig_reactome_down)))

log_msg("\n5. ...")

all_significant <- rbind(
  sig_hallmark %>% mutate(database = "Hallmark"),
  sig_reactome %>% mutate(database = "Reactome")
)

all_significant <- all_significant[order(-abs(NES))]

all_sig_csv <- all_significant %>%
  select(-leadingEdge) %>%
  as.data.frame()
write.csv(all_sig_csv, "results/gsea_significant_pathways.csv", row.names = FALSE)

saveRDS(all_significant, "results/gsea_significant_pathways.rds")

log_msg(sprintf("   : %d", nrow(all_significant)))
log_msg("   : results/gsea_significant_pathways.csv/.rds")

log_msg("\n6. ...")

log_msg("\n   === Hallmark(Top 10) ===")
top_hallmark_up <- head(sig_up, 10)
for (i in 1:min(10, nrow(top_hallmark_up))) {
  pw <- top_hallmark_up[i, ]
  log_msg(sprintf("   %d. %s", i, gsub("HALLMARK_", "", pw$pathway)))
  log_msg(sprintf("      NES=%.3f, p=%.4f, FDR=%.4f, size=%d",
                  pw$NES, pw$pval, pw$padj, pw$size))
}

if (nrow(sig_down) > 0) {
  log_msg("\n   === Hallmark(Top 10) ===")
  top_hallmark_down <- head(sig_down[order(NES)], 10)
  for (i in 1:min(10, nrow(top_hallmark_down))) {
    pw <- top_hallmark_down[i, ]
    log_msg(sprintf("   %d. %s", i, gsub("HALLMARK_", "", pw$pathway)))
    log_msg(sprintf("      NES=%.3f, p=%.4f, FDR=%.4f, size=%d",
                    pw$NES, pw$pval, pw$padj, pw$size))
  }
}

log_msg("\n   === Reactome(Top 10) ===")
top_reactome_up <- head(sig_reactome_up[order(-NES)], 10)
for (i in 1:min(10, nrow(top_reactome_up))) {
  pw <- top_reactome_up[i, ]
  pathway_short <- gsub("REACTOME_", "", pw$pathway)
  pathway_short <- substr(pathway_short, 1, 60)
  log_msg(sprintf("   %d. %s", i, pathway_short))
  log_msg(sprintf("      NES=%.3f, p=%.4f, FDR=%.4f", pw$NES, pw$pval, pw$padj))
}

log_msg("\n7. GSEA...")

report_file <- "results/gsea_summary_report.txt"
cat("", file = report_file)

write_report <- function(text) {
  cat(text, "\n", file = report_file, append = TRUE)
}

write_report("================================================================================")
write_report("GSEA -  (CF vs Normal)")
write_report("================================================================================")
write_report(sprintf(": %s", Sys.time()))
write_report("")

write_report("1. ")
write_report(sprintf("   Ranked genes: %d", length(ranked_genes)))
write_report(sprintf("   log2FC: %.3f to %.3f", min(ranked_genes), max(ranked_genes)))
write_report("")

write_report("2. ")
write_report(sprintf("   Hallmark: %d", length(hallmark_list)))
write_report(sprintf("   Reactome: %d", length(reactome_list)))
write_report("")

write_report("3. GSEA")
write_report("   Permutations: 10,000")
write_report("   Gene set size: 15-500")
write_report("   Significance: FDR < 0.05")
write_report("")

write_report("4. Hallmark")
write_report(sprintf("   : %d", nrow(gsea_hallmark)))
write_report(sprintf("   (FDR<0.05): %d (%.1f%%)", nrow(sig_hallmark),
                    nrow(sig_hallmark)/nrow(gsea_hallmark)*100))
write_report(sprintf("     : %d", nrow(sig_up)))
write_report(sprintf("     : %d", nrow(sig_down)))
write_report("")

write_report("5. Reactome")
write_report(sprintf("   : %d", nrow(gsea_reactome)))
write_report(sprintf("   (FDR<0.05): %d (%.1f%%)", nrow(sig_reactome),
                    nrow(sig_reactome)/nrow(gsea_reactome)*100))
write_report(sprintf("     : %d", nrow(sig_reactome_up)))
write_report(sprintf("     : %d", nrow(sig_reactome_down)))
write_report("")

write_report("6. Hallmark(Top 15, )")
if (nrow(sig_up) > 0) {
  for (i in 1:min(15, nrow(sig_up))) {
    pw <- sig_up[i, ]
    write_report(sprintf("   %2d. %s", i, gsub("HALLMARK_", "", pw$pathway)))
    write_report(sprintf("       NES=%.3f, FDR=%.4f, size=%d", pw$NES, pw$padj, pw$size))
  }
} else {
  write_report("   ")
}
write_report("")

write_report("7. Hallmark(Top 10, )")
if (nrow(sig_down) > 0) {
  for (i in 1:min(10, nrow(sig_down[order(NES)]))) {
    pw <- sig_down[order(NES)][i, ]
    write_report(sprintf("   %2d. %s", i, gsub("HALLMARK_", "", pw$pathway)))
    write_report(sprintf("       NES=%.3f, FDR=%.4f, size=%d", pw$NES, pw$padj, pw$size))
  }
} else {
  write_report("   ")
}
write_report("")

write_report("8. Reactome(Top 20, )")
if (nrow(sig_reactome_up) > 0) {
  for (i in 1:min(20, nrow(sig_reactome_up))) {
    pw <- sig_reactome_up[i, ]
    pathway_name <- gsub("REACTOME_", "", pw$pathway)
    pathway_name <- substr(pathway_name, 1, 70)
    write_report(sprintf("   %2d. %s", i, pathway_name))
    write_report(sprintf("       NES=%.3f, FDR=%.4f", pw$NES, pw$padj))
  }
} else {
  write_report("   ")
}
write_report("")

write_report("================================================================================")
write_report("GSEA")
write_report("================================================================================")

log_msg(sprintf("   : %s", report_file))

log_msg("\n8. session...")
session_file <- "logs/session_info.txt"
writeLines(capture.output(sessionInfo()), session_file)
log_msg(sprintf("   Session info: %s", session_file))

log_msg("\n================================================================================")
log_msg("02")
log_msg("================================================================================")
log_msg(sprintf("OK Hallmark GSEA: %d,%d", nrow(gsea_hallmark), nrow(sig_hallmark)))
log_msg(sprintf("OK Reactome GSEA: %d,%d", nrow(gsea_reactome), nrow(sig_reactome)))
log_msg(sprintf("OK : %d", nrow(all_significant)))
log_msg("")
log_msg(":")
log_msg("  - results/gsea_hallmark_results.rds/.csv")
log_msg("  - results/gsea_reactome_results.rds/.csv")
log_msg("  - results/gsea_significant_pathways.csv")
log_msg("  - results/gsea_summary_report.txt")
log_msg("")
log_msg(": 03 - GSEA")
log_msg("================================================================================")

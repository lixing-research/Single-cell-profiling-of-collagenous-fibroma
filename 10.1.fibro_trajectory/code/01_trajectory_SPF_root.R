#!/usr/bin/env Rscript

rm(list=ls())
set.seed(23456)

cat("========================================\n")
cat("Fibroblast Trajectory Analysis - SPF as Root\n")
cat("========================================\n\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/data", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

DEBUG_MODE <- FALSE

CHECKPOINT_DIR <- "results/data"
CDS_CHECKPOINT <- file.path(CHECKPOINT_DIR, "cds_ordered.rds")
GSVA_CHECKPOINT <- file.path(CHECKPOINT_DIR, "gsva_matrix.rds")

log_file <- "logs/trajectory_analysis.log"
write_log <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp, msg))
  cat(paste(timestamp, msg, "\n"), file = log_file, append = TRUE)
}

write_log("Analysis started - SPF as root state")

write_log("Loading libraries...")
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(monocle)
  library(ggplot2)
  library(cowplot)
  library(GSVA)
  library(msigdbr)
  library(pheatmap)
  library(viridis)
  library(RColorBrewer)
  library(parallel)
  library(doParallel)
  library(BiocParallel)
})

n_cores <- min(80, parallel::detectCores() - 10)
registerDoParallel(cores = n_cores)
register(MulticoreParam(workers = n_cores))
write_log(paste("Using", n_cores, "cores"))

pdf.options(family = "Helvetica")
theme_set(theme_minimal(base_family = "Helvetica", base_size = 11))
FIGURE_DPI <- 300

write_log("Step 1: Loading fibroblast data...")

data_path <- 'input.data/fibro_sce.rds'
sce.all <- readRDS(data_path)
write_log(paste("Loaded", ncol(sce.all), "cells"))

if("celltype_abbr" %in% colnames(sce.all@meta.data)) {
  Idents(sce.all) <- sce.all$celltype_abbr
  celltype_field <- "celltype_abbr"
} else {
  Idents(sce.all) <- sce.all$celltype
  celltype_field <- "celltype"
}

cat("Cell types:\n")
print(table(Idents(sce.all)))

allCells <- names(Idents(sce.all))
n_cells <- length(allCells)
if(n_cells > 10000) {
  set.seed(123)
  sample_cells <- sample(allCells, 10000)
  cg_sce <- sce.all[, sample_cells]
} else {
  cg_sce <- sce.all
}

write_log(paste("Using", ncol(cg_sce), "cells for analysis"))

write_log("Step 2: Creating Monocle CDS object...")

Mono_matrix <- as(as.matrix(GetAssayData(cg_sce, slot = "counts")), 'sparseMatrix')

feature_ann <- data.frame(
  gene_id = rownames(Mono_matrix),
  gene_short_name = rownames(Mono_matrix)
)
rownames(feature_ann) <- rownames(Mono_matrix)
Mono_fd <- new("AnnotatedDataFrame", data = feature_ann)

sample_ann <- cg_sce@meta.data
Mono_pd <- new("AnnotatedDataFrame", data = sample_ann)

Mono.cds <- newCellDataSet(
  Mono_matrix,
  phenoData = Mono_pd,
  featureData = Mono_fd,
  expressionFamily = negbinomial.size()
)

write_log(paste("CDS created:", ncol(Mono.cds), "cells,", nrow(Mono.cds), "genes"))

write_log("Step 3: Preprocessing...")

Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)

disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
write_log(paste("Selected", nrow(unsup_clustering_genes), "genes"))

Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)

if (DEBUG_MODE && file.exists(CDS_CHECKPOINT)) {
  write_log("Loading CDS checkpoint (skipping DDRTree computation)...")
  cds <- readRDS(CDS_CHECKPOINT)

  root_state <- attr(cds, "root_state")
  if (is.null(root_state)) root_state <- 1
  write_log(paste("Loaded CDS from checkpoint. Root state:", root_state))
} else {
  write_log("Step 4: Running DDRTree dimension reduction...")

  Mono.cds <- reduceDimension(
    Mono.cds,
    max_components = 2,
    method = 'DDRTree',
    verbose = TRUE
  )

  Mono.cds <- orderCells(Mono.cds)

  write_log("Step 5: Setting SPF as root state...")

  state_df <- pData(Mono.cds) %>%
    group_by(State) %>%
    summarize(n_cells = n(), .groups = 'drop')

  for(ct in unique(pData(Mono.cds)[[celltype_field]])) {
    state_df[[paste0("pct_", ct)]] <- sapply(state_df$State, function(s) {
      state_cells <- pData(Mono.cds)$State == s
      ct_cells <- pData(Mono.cds)[[celltype_field]] == ct
      sum(state_cells & ct_cells) / sum(state_cells) * 100
    })
  }

  cat("\nState composition:\n")
  print(state_df)

  spf_col <- grep("SPF|SecPapF|papillary", colnames(state_df), value = TRUE, ignore.case = TRUE)
  if(length(spf_col) > 0) {
    root_state <- state_df$State[which.max(state_df[[spf_col[1]]])]
    write_log(paste("Root state set to", root_state, "(highest SPF proportion:",
                    round(max(state_df[[spf_col[1]]]), 1), "%)"))
  } else {

    spf_counts <- sapply(state_df$State, function(s) {
      state_cells <- pData(Mono.cds)$State == s
      ct_cells <- grepl("SPF|SecPapF|papillary", pData(Mono.cds)[[celltype_field]], ignore.case = TRUE)
      sum(state_cells & ct_cells)
    })
    root_state <- state_df$State[which.max(spf_counts)]
    write_log(paste("Root state set to", root_state, "(most SPF cells:", max(spf_counts), ")"))
  }

  cds <- orderCells(Mono.cds, root_state = root_state)
  cds <- detectGenes(cds, min_expr = 0.1)

  write_log("Trajectory ordering complete with SPF as root")

  attr(cds, "root_state") <- root_state
  saveRDS(cds, CDS_CHECKPOINT)
  write_log(paste("CDS checkpoint saved:", CDS_CHECKPOINT))
}

write_log("Step 6: Loading gene sets...")

hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()

fibroblast_hallmarks <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_P53_PATHWAY"
)

hallmark_filtered <- hallmark_sets %>% filter(gs_name %in% fibroblast_hallmarks)

kegg_sets <- tryCatch({
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    as.data.frame()
}, error = function(e) {
  c2_sets <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    as.data.frame()
  c2_sets[grepl("KEGG_", c2_sets$gs_name), ]
})

kegg_relevant <- c(
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_FOCAL_ADHESION",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_WNT_SIGNALING_PATHWAY",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_CELL_CYCLE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
)

kegg_filtered <- kegg_sets %>% filter(grepl(paste(kegg_relevant, collapse = "|"), gs_name))

reactome_sets <- tryCatch({
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    as.data.frame()
}, error = function(e) {
  c2_sets <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    as.data.frame()
  c2_sets[grepl("REACTOME_", c2_sets$gs_name), ]
})

reactome_relevant <- c("REACTOME_COLLAGEN", "REACTOME_EXTRACELLULAR_MATRIX")
reactome_filtered <- reactome_sets %>% filter(grepl(paste(reactome_relevant, collapse = "|"), gs_name))

all_genesets <- rbind(hallmark_filtered, kegg_filtered, reactome_filtered)
geneset_list <- split(all_genesets$gene_symbol, all_genesets$gs_name)

write_log(paste("Loaded", length(geneset_list), "gene sets"))

if (DEBUG_MODE && file.exists(GSVA_CHECKPOINT)) {
  write_log("Loading GSVA checkpoint (skipping GSVA computation)...")
  gsva_matrix <- readRDS(GSVA_CHECKPOINT)
  write_log(paste("Loaded GSVA matrix from checkpoint:", nrow(gsva_matrix), "gene sets x", ncol(gsva_matrix), "cells"))
} else {
  write_log("Step 7: Running GSVA analysis...")

  expr_matrix <- as.matrix(exprs(cds))

  geneset_list_filtered <- lapply(geneset_list, function(genes) {
    intersect(genes, rownames(expr_matrix))
  })
  geneset_list_filtered <- geneset_list_filtered[sapply(geneset_list_filtered, length) >= 5]

  write_log(paste("Running GSVA on", length(geneset_list_filtered), "gene sets"))

  gsva_version <- packageVersion("GSVA")
  if (gsva_version >= "1.50.0") {
    if (exists("gsvaParam", mode = "function")) {
      gsva_param <- gsvaParam(
        exprData = expr_matrix,
        geneSets = geneset_list_filtered,
        kcdf = "Poisson",
        minSize = 5,
        maxSize = 500
      )
      gsva_results <- gsva(gsva_param, verbose = TRUE, BPPARAM = MulticoreParam(workers = n_cores))

      if (class(gsva_results)[1] == "matrix") {
        gsva_matrix <- gsva_results
      } else {
        gsva_matrix <- SummarizedExperiment::assay(gsva_results)
      }
    } else {
      gsva_matrix <- GSVA::gsva(
        expr = expr_matrix,
        gset.idx.list = geneset_list_filtered,
        kcdf = "Poisson",
        min.sz = 5,
        max.sz = 500,
        parallel.sz = n_cores,
        verbose = TRUE
      )
    }
  } else {
    gsva_matrix <- gsva(
      expr_matrix,
      geneset_list_filtered,
      method = "gsva",
      kcdf = "Poisson",
      min.sz = 5,
      max.sz = 500,
      parallel.sz = n_cores,
      verbose = TRUE
    )
  }

  write_log("GSVA complete")

  saveRDS(gsva_matrix, GSVA_CHECKPOINT)
  write_log(paste("GSVA checkpoint saved:", GSVA_CHECKPOINT))
}

write_log("Step 8: Calculating pseudotime correlations...")

pseudotime_gsva <- data.frame(
  Pseudotime = pData(cds)$Pseudotime,
  State = pData(cds)$State,
  Celltype = pData(cds)[[celltype_field]],
  t(gsva_matrix)
)

gsva_correlations <- cor(pseudotime_gsva$Pseudotime,
                         pseudotime_gsva[,4:ncol(pseudotime_gsva)],
                         use = "complete.obs")

gsva_cor_df <- data.frame(
  GeneSet = colnames(gsva_correlations),
  Correlation = as.numeric(gsva_correlations),
  AbsCorrelation = abs(as.numeric(gsva_correlations))
) %>%
  arrange(desc(AbsCorrelation))

write.csv(gsva_cor_df, "results/data/geneset_pseudotime_correlations_SPF_root.csv", row.names = FALSE)
write_log("Saved: geneset_pseudotime_correlations_SPF_root.csv")

cat("\n=== Top 10 Gene Sets Correlated with Pseudotime (SPF root) ===\n")
print(head(gsva_cor_df[,1:2], 10))

write_log("Step 9: Calculating pseudotime statistics by celltype...")

pseudotime_by_celltype <- pData(cds) %>%
  group_by(!!sym(celltype_field)) %>%
  summarize(
    mean_pseudotime = mean(Pseudotime, na.rm = TRUE),
    median_pseudotime = median(Pseudotime, na.rm = TRUE),
    sd_pseudotime = sd(Pseudotime, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  arrange(mean_pseudotime)

write.csv(pseudotime_by_celltype, "results/data/pseudotime_by_celltype_SPF_root.csv", row.names = FALSE)

cat("\n=== Pseudotime by Cell Type (SPF root) ===\n")
print(pseudotime_by_celltype)

write_log("Step 10: Generating visualizations...")

trajectory_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)
)

p1 <- plot_cell_trajectory(cds, color_by = celltype_field, size = 1, show_backbone = TRUE) +
  ggtitle("Trajectory by Subtype (SPF root)") +
  trajectory_theme +
  scale_color_brewer(palette = "Set2", name = "Cell Type")

p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE) +
  ggtitle("Pseudotime (SPF root)") +
  trajectory_theme +
  scale_color_viridis(name = "Pseudotime", option = "C")

p3 <- plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  ggtitle("Trajectory States") +
  trajectory_theme +
  scale_color_brewer(palette = "Set1", name = "State")

p4 <- plot_cell_trajectory(cds, color_by = "orig.ident", size = 1, show_backbone = TRUE) +
  ggtitle("By Sample Origin") +
  trajectory_theme +
  scale_color_manual(values = c("CF" = "#E41A1C", "Normal" = "#377EB8"), name = "Sample")

combined_basic <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv")
ggsave("results/figures/trajectory_combined_SPF_root.pdf", combined_basic,
       width = 16, height = 16, dpi = FIGURE_DPI)
write_log("Saved: trajectory_combined_SPF_root.pdf")

ggsave("results/figures/trajectory_by_celltype.pdf", p1, width = 8, height = 8, dpi = FIGURE_DPI)
ggsave("results/figures/trajectory_by_pseudotime.pdf", p2, width = 8, height = 8, dpi = FIGURE_DPI)

write_log("Step 11: Saving results...")

saveRDS(cds, "results/data/fibroblast_trajectory_cds_SPF_root.rds")

trajectory_data <- pData(cds) %>%
  dplyr::select(all_of(c("orig.ident", celltype_field, "State", "Pseudotime"))) %>%
  mutate(
    trajectory_stage = case_when(
      Pseudotime < quantile(Pseudotime, 0.33, na.rm = TRUE) ~ "Early",
      Pseudotime < quantile(Pseudotime, 0.67, na.rm = TRUE) ~ "Mid",
      TRUE ~ "Late"
    )
  )

write.csv(trajectory_data, "results/data/trajectory_metadata_SPF_root.csv", row.names = TRUE)
write_log("Saved: trajectory_metadata_SPF_root.csv")

write_log("========================================")
write_log("Analysis Complete!")
write_log("========================================")

cat("\n\n========== SUMMARY ==========\n")
cat("Root state:", root_state, "(SPF-enriched)\n")
cat("\nPseudotime by cell type (SPF root):\n")
print(pseudotime_by_celltype)
cat("\nTop 5 positively correlated gene sets:\n")
print(head(gsva_cor_df[gsva_cor_df$Correlation > 0, 1:2], 5))
cat("\nTop 5 negatively correlated gene sets:\n")
print(head(gsva_cor_df[gsva_cor_df$Correlation < 0, 1:2], 5))

write_log("All outputs saved to results/")

rm(list=ls())
options(stringsAsFactors = F)
source('../scRNA_scripts/lib.R')
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(harmony)

cat("========================================\n")
cat("Step 02: Batch Integration")
cat("========================================\n\n")

sce.all <- readRDS("../1.QC/01_QC/sce_qc.rds")
cat("Loaded QC data from 01_QC/sce_qc.rds\n")
cat(sprintf("Total cells: %d\n", ncol(sce.all)))
cat(sprintf("Total genes: %d\n", nrow(sce.all)))

sce.all.filt = sce.all
cat(sprintf("\nUsing QC data: %d cells, %d genes\n", ncol(sce.all.filt), nrow(sce.all.filt)))

cat("\n========================================\n")
cat("Step 2: Normalization and Harmony Integration\n")
cat("========================================\n\n")

set.seed(10086)
table(sce.all.filt$orig.ident)

output_dir <- "02_harmony"
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

source('../../scRNA_scripts/harmony.R')
sce.all.int = run_harmony(sce.all.filt)

cat("\nHarmony integration completed\n")
setwd('..')

cat("\n========================================\n")
cat("Step 3: Saving Integrated Data\n")
cat("========================================\n\n")

table(sce.all.int$RNA_snn_res.0.2)

sce.all.int <- SetIdent(sce.all.int, value = "RNA_snn_res.0.2")
sce.all.int$seurat_clusters <- sce.all.int$RNA_snn_res.0.2

cat(sprintf("Number of clusters (res=0.2): %d\n", length(unique(sce.all.int$seurat_clusters))))

saveRDS(sce.all.int, file = "sce_integrated.rds")
cat("\nIntegrated data saved to: 02_harmony/sce_integrated.rds\n")

integration_summary <- list(
  date = Sys.Date(),
  n_cells = ncol(sce.all.int),
  n_genes = nrow(sce.all.int),
  n_clusters = length(unique(sce.all.int$seurat_clusters)),
  samples = unique(sce.all.int$orig.ident),
  resolution = 0.2
)

saveRDS(integration_summary, "integration_summary.rds")
cat("Integration summary saved to: 02_harmony/integration_summary.rds\n")

setwd('..')

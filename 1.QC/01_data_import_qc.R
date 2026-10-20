rm(list=ls())
options(stringsAsFactors = F)

cat("========================================\n")
cat("Step 1.1: Data Import")
cat("========================================\n\n")

cat("Loading required libraries...\n")
source('../scRNA_scripts/lib.R')

cat("\n========================================\n")
cat("Step 1: Importing Pre-filtered Data")
cat("========================================\n")

upstream_dir <- '../data'
samples <- c('CF_sample', 'Normal_sample')
data_paths <- file.path(upstream_dir, samples, 'output', 'filter_matrix')
names(data_paths) <- samples

dir.create("01_QC", showWarnings = FALSE)

cat("Reading DNBelab C4 pre-filtered data:\n")
for(i in seq_along(data_paths)) {
  cat(sprintf("  %s: %s\n", samples[i], data_paths[i]))
}

sceList <- lapply(samples, function(pro) {
  cat(paste("Loading sample:", pro, "...\n"))

  tmp <- Read10X(data.dir = data_paths[pro], gene.column = 1)

  if(length(tmp) == 2) {
    ct <- tmp[[1]]
    cat("  Multi-modal data detected, using Gene Expression matrix\n")
  } else {
    ct <- tmp
  }

  sce <- CreateSeuratObject(
    counts = ct,
    project = pro,
    min.cells = 3,
    min.features = 200
  )

  cat(paste("  Cells:", ncol(sce), "| Genes:", nrow(sce), "\n"))
  return(sce)
})

cat("\n========================================\n")
cat("Sample Dimensions Summary:\n")
cat("========================================\n")
sample_dims <- do.call(rbind, lapply(sceList, dim))
rownames(sample_dims) <- samples
colnames(sample_dims) <- c("Genes", "Cells")
print(sample_dims)

setwd("01_QC")

cat("\n========================================\n")
cat("Step 2: Merging Samples")
cat("========================================\n")

sce.all <- merge(
  x = sceList[[1]],
  y = sceList[-1],
  add.cell.ids = samples
)

cat("\nMerging completed:\n")
cat(paste("Total cells:", ncol(sce.all), "\n"))
cat(paste("Total genes:", nrow(sce.all), "\n"))

if("RNA" %in% names(sce.all@assays)) {
  cat("\nSeurat v5 layer information:\n")
  cat("Layers found:", paste(names(sce.all@assays$RNA@layers), collapse = ", "), "\n")

  sce.all <- JoinLayers(sce.all)
  cat("Layers joined successfully.\n")
}

cat("\n========================================\n")
cat("Step 3: Adding Metadata")
cat("========================================\n")

phe <- sce.all@meta.data

cat("\nOriginal sample identities:\n")
print(table(phe$orig.ident))

phe$patient <- str_split(phe$orig.ident, '[_]', simplify = T)[,1]
cat("\nPatient/Condition groups:\n")
print(table(phe$patient))

phe$sample <- phe$patient
cat("\nSample groups:\n")
print(table(phe$sample))

sce.all@meta.data <- phe

cat("\n========================================\n")
cat("Step 4: Calculating QC Metrics")
cat("========================================\n")

sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
cat("Mitochondrial gene percentage calculated.\n")

sce.all[["percent.ribo"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
cat("Ribosomal gene percentage calculated.\n")

sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all, pattern = "^HB[^(P)]")
cat("Hemoglobin gene percentage calculated.\n")

cat("\nQC Metrics Summary:\n")
cat("------------------\n")
cat(sprintf("nFeature_RNA: %d - %d (median: %d)\n",
            min(sce.all$nFeature_RNA),
            max(sce.all$nFeature_RNA),
            median(sce.all$nFeature_RNA)))
cat(sprintf("nCount_RNA: %d - %d (median: %d)\n",
            min(sce.all$nCount_RNA),
            max(sce.all$nCount_RNA),
            median(sce.all$nCount_RNA)))
cat(sprintf("percent.mt: %.2f%% - %.2f%% (median: %.2f%%)\n",
            min(sce.all$percent.mt),
            max(sce.all$percent.mt),
            median(sce.all$percent.mt)))
cat(sprintf("percent.ribo: %.2f%% - %.2f%% (median: %.2f%%)\n",
            min(sce.all$percent.ribo),
            max(sce.all$percent.ribo),
            median(sce.all$percent.ribo)))

cat("\n========================================\n")
cat("Step 5: Saving Data")
cat("========================================\n")

saveRDS(sce.all, file = "sce_qc.rds")
cat("QC data saved to: sce_qc.rds (NOT normalized)\n")

qc_thresholds <- data.frame(
  metric = c("nFeature_RNA_min", "nFeature_RNA_max",
             "nCount_RNA_min", "nCount_RNA_max",
             "percent.mt_max", "percent.ribo_min"),
  recommended_value = c(
    max(200, quantile(sce.all$nFeature_RNA, 0.01)),
    quantile(sce.all$nFeature_RNA, 0.99),
    max(500, quantile(sce.all$nCount_RNA, 0.01)),
    quantile(sce.all$nCount_RNA, 0.99),
    median(sce.all$percent.mt) + 3*mad(sce.all$percent.mt),
    5
  ),
  description = c(
    "Minimum genes per cell (remove low-quality cells)",
    "Maximum genes per cell (remove potential doublets)",
    "Minimum UMIs per cell (remove low-quality cells)",
    "Maximum UMIs per cell (remove potential doublets)",
    "Maximum mitochondrial percentage (remove dying cells)",
    "Minimum ribosomal percentage (remove ambient RNA)"
  )
)

saveRDS(qc_thresholds, file = "qc_thresholds.rds")
cat("\nRecommended QC thresholds:\n")
print(qc_thresholds)

data_summary <- list(
  date = Sys.Date(),
  samples = samples,
  n_samples = length(samples),
  total_cells = ncol(sce.all),
  total_genes = nrow(sce.all),
  cells_per_sample = table(sce.all$orig.ident),
  qc_metrics = data.frame(
    metric = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    min = c(min(sce.all$nFeature_RNA), min(sce.all$nCount_RNA),
            min(sce.all$percent.mt), min(sce.all$percent.ribo)),
    median = c(median(sce.all$nFeature_RNA), median(sce.all$nCount_RNA),
              median(sce.all$percent.mt), median(sce.all$percent.ribo)),
    max = c(max(sce.all$nFeature_RNA), max(sce.all$nCount_RNA),
            max(sce.all$percent.mt), max(sce.all$percent.ribo))
  )
)

saveRDS(data_summary, file = "data_import_summary.rds")
cat("Data summary saved to: data_import_summary.rds\n")

sink("data_import_report.txt")
cat("========================================\n")
cat("Data Import Report\n")
cat("========================================\n\n")
cat(paste("Date:", Sys.Date(), "\n"))
cat(paste("Samples imported:", length(samples), "\n"))
cat(paste("Sample names:", paste(samples, collapse = ", "), "\n\n"))

cat("Dimensions:\n")
cat(paste("Total cells:", ncol(sce.all), "\n"))
cat(paste("Total genes:", nrow(sce.all), "\n\n"))

cat("Cells per sample:\n")
print(table(sce.all$orig.ident))
cat("\n")

cat("QC Metrics Summary:\n")
print(data_summary$qc_metrics)
cat("\n")

cat("Recommended QC Thresholds:\n")
print(qc_thresholds)
cat("\n")

cat("Next steps:\n")
cat("1. Review QC plots (violin plots, scatter plots) to validate thresholds\n")
cat("2. Run QC filtering based on recommended thresholds or adjust manually\n")
cat("3. Run Harmony integration for batch effect correction\n")
cat("4. Proceed with downstream analysis (clustering, marker identification)\n")
sink()

cat("\n========================================\n")
cat("Step 6: Generating QC Plots")
cat("========================================\n")

cat("Creating violin plots...\n")
p1 <- VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
              ncol = 2, pt.size = 0, group.by = "sample") +
  plot_annotation(title = "QC Metrics by Sample")

ggsave("qc_violin_plots.pdf", p1, width = 12, height = 10)
cat("  Saved: qc_violin_plots.pdf\n")

cat("Creating scatter plots...\n")
p2 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "sample") +
  geom_hline(yintercept = qc_thresholds$recommended_value[1], linetype = "dashed", color = "red") +
  geom_hline(yintercept = qc_thresholds$recommended_value[2], linetype = "dashed", color = "red")

p3 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "percent.mt",
                     group.by = "sample") +
  geom_hline(yintercept = qc_thresholds$recommended_value[5], linetype = "dashed", color = "red")

p_scatter <- p2 | p3
ggsave("qc_scatter_plots.pdf", p_scatter, width = 14, height = 6)
cat("  Saved: qc_scatter_plots.pdf\n")

cat("\nSample-specific QC summary:\n")
qc_by_sample <- sce.all@meta.data %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_mt = median(percent.mt),
    median_ribo = median(percent.ribo)
  )
print(qc_by_sample)

saveRDS(qc_by_sample, file = "qc_by_sample_summary.rds")
cat("\nSaved: qc_by_sample_summary.rds\n")

cat("\n========================================\n")
cat("OK Data import completed successfully!\n")
cat("OK Results saved in 01_QC/\n")
cat("========================================\n")

setwd("..")

gc()

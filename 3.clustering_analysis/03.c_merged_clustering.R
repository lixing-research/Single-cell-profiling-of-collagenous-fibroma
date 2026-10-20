rm(list = ls())
set.seed(123)
gc()

cat("========================================\n")
cat("Clustering Analysis - Part C: Merging\n")
cat("========================================\n\n")

cat("Loading required libraries...\n")

current_dir <- basename(getwd())
if(current_dir != "3.clustering_analysis") {
  stop("Error: Please run this script from the 3.clustering_analysis directory!\n",
       "Current directory: ", getwd(), "\n",
       "Expected: .../3.clustering_analysis")
}

source('../scRNA_scripts/load_packages.R')
source('../scRNA_scripts/lib.R')

required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "ComplexHeatmap")
load_packages(required_packages)

output_dir <- "03c_result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(paste("Results will be saved to:", file.path(getwd(), output_dir), "\n"))

cat("\n========================================\n")
cat("Step 1: Loading Data")
cat("========================================\n")

input_files <- c(
  "03b_result/sce_clustered_final.rds",
  "03a_analysis/sce_clustered_annotated.rds"
)

input_file <- NULL
for(f in input_files) {
  if(file.exists(f)) {
    input_file <- f
    break
  }
}

if(is.null(input_file)) {
  stop("Cannot find clustered data. Please run 03.a or 03.b first.\n",
       "Looking for files in: 03b_result/ or 03a_analysis/\n")
}

cat(paste("Loading data from:", input_file, "\n"))
sce.all <- readRDS(input_file)

cat(paste("  Cells:", ncol(sce.all), "\n"))
cat(paste("  Original clusters:", length(unique(sce.all$seurat_clusters)), "\n"))

if("predicted_celltype" %in% colnames(sce.all@meta.data)) {
  cat("OK Found cell type annotations from Part B\n")
  cat(paste("  Cell types:", length(unique(sce.all$predicted_celltype)), "\n"))
} else {
  cat("WARNING No cell type annotations found. Run Part B first for better results.\n")
}

cat("\n========================================\n")
cat("Step 2: Defining Merging Strategy")
cat("========================================\n")

sce.all@meta.data$original_clusters <- sce.all$seurat_clusters

cat("\nCurrent cluster distribution (17 clusters from Resolution=0.2):\n")
cluster_counts <- table(sce.all$seurat_clusters)
for(i in names(cluster_counts)) {
  cat(sprintf("  Cluster %2s: %5d cells\n", i, cluster_counts[i]))
}

cat("\nWARNING IMPORTANT: Merging strategy must be defined based on marker gene analysis!\n")
cat("Please review the marker genes from Part A and Part B before defining merges.\n\n")

cat("Example merging strategies (adjust based on actual data):\n")
cat("  Option 1: Keep all 17 clusters if they represent distinct cell types\n")
cat("  Option 2: Merge similar epithelial subtypes (e.g., basal, ciliated, secretory)\n")
cat("  Option 3: Merge fibroblast subtypes if they show similar signatures\n")
cat("  Option 4: Merge immune cell subtypes (T cells, B cells, macrophages)\n\n")

cat("Current setting: NO MERGING - Keeping all 17 original clusters\n")
cat("To merge clusters, edit the 'cluster_merge_map' below.\n\n")

cluster_merge_map <- setNames(
  as.character(0:16),
  as.character(0:16)
)

merged_celltype_mapping <- c(
  "0" = "Fib",
  "1" = "Fib",
  "2" = "Endo",
  "3" = "PeriSMC",
  "4" = "Epi",
  "5" = "Epi",
  "6" = "Epi",
  "7" = "Macro",
  "8" = "Epi",
  "9" = "Tcell",
  "10" = "Schwann",
  "11" = "GlandEpi",
  "12" = "GlandEpi",
  "13" = "Endo",
  "14" = "Neuro",
  "15" = "Mast"
)

celltype_order <- c("Fib", "Endo", "PeriSMC", "Epi", "Macro", "Tcell",
                    "Schwann", "GlandEpi", "Neuro", "Mast", "MyoepiSMC")

if("16" %in% unique(sce.all$seurat_clusters)) {
  merged_celltype_mapping["16"] <- "MyoepiSMC"
}

cat("\nOK Keeping all clusters for analysis\n")

sce.all@meta.data$merged_celltype <- merged_celltype_mapping[as.character(sce.all$seurat_clusters)]

available_celltypes <- celltype_order[celltype_order %in% unique(sce.all@meta.data$merged_celltype)]
sce.all@meta.data$merged_celltype <- factor(
  sce.all@meta.data$merged_celltype,
  levels = available_celltypes
)

cat("\n========================================\n")
cat("Step 3: Merging Statistics")
cat("========================================\n")

original_stats <- table(sce.all$original_clusters)
cat("\nOriginal clusters (17 total):\n")
for(i in names(sort(original_stats))) {
  cat(sprintf("  Cluster %s: %d cells (%.1f%%)\n",
              i, original_stats[i],
              original_stats[i]/ncol(sce.all)*100))
}

merged_stats <- table(sce.all$merged_celltype)
cat("\nMerged cell types (11 total):\n")
for(ct in available_celltypes) {
  if(ct %in% names(merged_stats)) {
    cat(sprintf("  %s: %d cells (%.1f%%)\n",
                ct, merged_stats[ct],
                merged_stats[ct]/ncol(sce.all)*100))
  }
}

merge_summary <- data.frame(
  Original_Cluster = names(merged_celltype_mapping),
  Merged_CellType = unname(merged_celltype_mapping),
  stringsAsFactors = FALSE
)

merge_summary$NCells <- sapply(merge_summary$Original_Cluster, function(x) {
  sum(sce.all$original_clusters == x)
})

merge_summary$Percentage <- round(merge_summary$NCells / ncol(sce.all) * 100, 2)

cat("\nMerge mapping summary:\n")
print(merge_summary)

write.csv(merge_summary, file.path(output_dir, "cluster_merge_mapping.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("Step 4: Validating Merge Rationale")
cat("========================================\n")

merge_validation_markers <- list(
  Fibroblast = list(
    cluster0 = c("CFD", "APOD", "IGFBP6"),
    cluster1 = c("MMP11", "POSTN", "COL5A2"),
    shared = c("COL1A1", "DCN", "LUM")
  ),
  Endothelial = list(
    cluster2 = c("CCL14", "PECAM1", "CLEC14A"),
    cluster13 = c("TFF3", "CCL21", "PROX1"),
    shared = c("VWF", "CD31", "PECAM1")
  ),
  Epithelial = list(
    cluster4 = c("DEFB1", "KRT23", "S100A2"),
    cluster6 = c("CDH1", "CDH3", "TIAM1"),
    cluster8 = c("KRT15", "KRT17", "CALML3"),
    shared = c("EPCAM", "KRT8", "KRT18")
  )
)

cat("\n========================================\n")
cat("Step 5: Creating Comparison Visualizations")
cat("========================================\n")

pdf(file.path(output_dir, "cluster_merging_comparison.pdf"), width = 20, height = 10)

p1 <- DimPlot(sce.all,
              group.by = "original_clusters",
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
  ggtitle("Original 17 Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- DimPlot(sce.all,
              group.by = "merged_celltype",
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
  ggtitle("Merged 11 Cell Types (Abbreviated)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

combined_plot <- p1 | p2
print(combined_plot)

library(reshape2)

comparison_df <- data.frame(
  Category = c(rep("Original", length(original_stats)),
               rep("Merged", length(merged_stats))),
  Type = c(names(original_stats), names(merged_stats)),
  Count = c(as.numeric(original_stats), as.numeric(merged_stats))
)

p_bar <- ggplot(comparison_df, aes(x = Type, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Category, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Distribution: Original vs Merged",
       x = "Cluster/Cell Type",
       y = "Number of Cells") +
  scale_fill_manual(values = c("Original" = "#377EB8", "Merged" = "#E41A1C"))

print(p_bar)

merge_flow <- merge_summary %>%
  group_by(Merged_CellType) %>%
  summarise(
    Original_Clusters = paste(Original_Cluster, collapse = "+"),
    Total_Cells = sum(NCells),
    Total_Percentage = sum(Percentage)
  )

merge_flow$Merged_CellType <- factor(merge_flow$Merged_CellType, levels = available_celltypes)
merge_flow <- merge_flow %>% arrange(Merged_CellType)

cat("\nMerge flow summary:\n")
print(merge_flow)

dev.off()

cat("\n========================================\n")
cat("Step 6: Marker Expression Validation")
cat("========================================\n")

pdf(file.path(output_dir, "merged_celltype_markers.pdf"), width = 16, height = 12)

merged_markers <- list(
  Fib = c("COL1A1", "DCN", "LUM", "CFD", "MMP11"),
  Endo = c("PECAM1", "VWF", "CCL14", "CCL21", "LYVE1"),
  PeriSMC = c("RGS5", "ACTA2", "TAGLN", "MYH11"),
  Epi = c("EPCAM", "KRT8", "KRT18", "CDH1", "KRT17"),
  Macro = c("CD68", "CD163", "C1QA", "C1QB"),
  Tcell = c("CD3D", "CD3E", "IL32", "TRBC2"),
  Schwann = c("PTPRG", "FTX", "DCLK1", "NAV3"),
  GlandEpi = c("AZGP1", "PIP", "AQP5", "KRT19", "SCGB2A2", "CRISP3"),
  Neuro = c("CNTN5", "DPP10", "RBFOX1"),
  Mast = c("CPA3", "TPSAB1", "TPSB2", "CTSG"),
  MyoepiSMC = c("PCP4", "SYNM", "KCNMB1")
)

merged_markers <- merged_markers[available_celltypes]

all_markers <- unique(unlist(merged_markers))
available_markers <- all_markers[all_markers %in% rownames(sce.all)]

cat(paste("Testing", length(available_markers), "marker genes\n"))

if(length(available_markers) > 0) {
  p_dot <- DotPlot(sce.all,
                   features = available_markers,
                   group.by = "merged_celltype",
                   scale = TRUE,
                   cols = c("lightgrey", "red"),
                   dot.scale = 6,
                   scale.min = 0,
                   scale.max = 100) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("Merged Cell Type Marker Expression (Abbreviated Names)")
  print(p_dot)
}

dev.off()

cat("\n========================================\n")
cat("Step 7: Saving Merged Results")
cat("========================================\n")

saveRDS(sce.all, file.path(output_dir, "sce_merged_clusters.rds"))
cat("Saved merged Seurat object to sce_merged_clusters.rds\n")

celltype_summary <- data.frame(
  CellType = names(merged_stats),
  NCells = as.numeric(merged_stats),
  Percentage = round(as.numeric(merged_stats)/ncol(sce.all)*100, 2),
  Original_Clusters = sapply(names(merged_stats), function(ct) {
    clusters <- names(merged_celltype_mapping)[merged_celltype_mapping == ct]
    paste(clusters, collapse = ",")
  })
)

celltype_summary$CellType <- factor(celltype_summary$CellType, levels = available_celltypes)
celltype_summary <- celltype_summary[order(celltype_summary$CellType),]
write.csv(celltype_summary, file.path(output_dir, "merged_celltype_summary.csv"), row.names = FALSE)

cat("\nFinal cell type summary:\n")
print(celltype_summary)

cat("\n========================================\n")
cat("Step 8: Preparing for Downstream Analysis")
cat("========================================\n")

cat("\nFinding markers for merged cell types...\n")
Idents(sce.all) <- "merged_celltype"
merged_markers_de <- FindAllMarkers(sce.all,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.5,
                                   test.use = "wilcox",
                                   verbose = FALSE)

write.csv(merged_markers_de, file.path(output_dir, "merged_celltype_markers.csv"), row.names = FALSE)

top10_merged <- merged_markers_de %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top10_merged, file.path(output_dir, "top10_merged_markers.csv"), row.names = FALSE)

cat("\nTop 5 markers for each merged cell type:\n")
for(ct in available_celltypes) {
  if(ct %in% unique(top10_merged$cluster)) {
    cat(paste0("\n", ct, ": "))
    top5_genes <- top10_merged %>%
      filter(cluster == ct) %>%
      head(5) %>%
      pull(gene)
    cat(paste(top5_genes, collapse = ", "))
  }
}

cat("\n\n========================================\n")
cat("Part C: Cluster merging completed successfully!\n")
cat("========================================\n")

cat("\nSummary:\n")
cat(paste("  - Original clusters: 17 (0-16)\n"))
cat(paste("  - Merged cell types: 11 (abbreviated names)\n"))
cat(paste("  - Total cells:", ncol(sce.all), "\n"))
cat("\nAbbreviated cell types created:\n")
cat("  - Fib, Endo, PeriSMC, Epi, Macro, Tcell, Schwann, GlandEpi, Neuro, Mast, MyoepiSMC\n")
cat("\nResults saved:\n")
cat("  - sce_merged_clusters.rds (Seurat object with abbreviated annotations)\n")
cat("  - cluster_merge_mapping.csv (Merge mapping with abbreviated names)\n")
cat("  - merged_celltype_summary.csv (Cell type statistics with abbreviated names)\n")
cat("  - merged_celltype_markers.csv (Marker genes for abbreviated types)\n")
cat("  - Visualization PDFs (with abbreviated labels)\n")
cat("\nReady for downstream analysis with simplified cell type names!\n")
cat("========================================\n")

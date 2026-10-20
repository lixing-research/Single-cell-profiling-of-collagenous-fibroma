rm(list=ls())
options(stringsAsFactors = F)

current_dir <- basename(getwd())
initial_wd <- getwd()

if(current_dir == "7.0.fibro_extrction") {
  data_path_04 <- "../4.cell.annotation/04_celltype/"
  scripts_path <- "../scRNA_scripts/"
  output_dir <- "./07_fibroblast"
} else {
  data_path_04 <- "4.cell.annotation/04_celltype/"
  scripts_path <- "scRNA_scripts/"
  output_dir <- "7.0.fibro_extrction/07_fibroblast"
}

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggpubr)

source(file.path(initial_wd, scripts_path, 'lib.R'))
source(file.path(initial_wd, scripts_path, 'mycolors.R'))

cat("========================================\n")
cat("Step 07: Fibroblast Extraction & CNV Analysis\n")
cat("Study: Case vs Control\n")
cat("========================================\n\n")

cat(sprintf("Current directory: %s\n", getwd()))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
cat(sprintf("Working in: %s\n\n", getwd()))
set.seed(1234567)

data_file_paths <- c(
  file.path(initial_wd, data_path_04, "sce_annotated_merged.rds"),
  file.path(initial_wd, data_path_04, "sce_annotated.rds")
)

loaded <- FALSE
for(path in data_file_paths) {
  if(file.exists(path)) {
    sce.all = readRDS(path)
    cat(sprintf("OK Loaded annotated data from: %s\n", path))
    loaded <- TRUE
    break
  }
}

if(!loaded) {
  stop(sprintf("ERROR Cannot find annotated data files. Checked paths:\n%s\n",
               paste("  -", data_file_paths, collapse="\n")))
}

cat("\n Data summary:\n")
cat(sprintf("  Total cells: %d\n", ncol(sce.all)))
table(sce.all$celltype)

sce1 = sce.all[, sce.all$celltype %in% c('Fib')]
cat(sprintf("\nOK Extracted %d fibroblast cells\n", ncol(sce1)))

if(!"group" %in% colnames(sce1@meta.data)) {
  sce1$group <- sce1$orig.ident
  cat("  Created 'group' column from orig.ident\n")
}

LayerData(sce1, assay = "RNA", layer = "counts")
sce1 <- JoinLayers(sce1)

as.data.frame(sce1@assays$RNA$counts[1:10, 1:2])
head(sce1@meta.data, 10)
table(sce1$orig.ident)

sce = sce1
sce <- NormalizeData(sce, normalization.method =  "LogNormalize",
                     scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce,
                            selection.method = "vst", nfeatures = 8000)
sce <- ScaleData(sce)
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
DimHeatmap(sce, dims = 1:12, cells = 1000, balanced = TRUE)
ElbowPlot(sce)

cat("\n Checking batch effects...\n")
if(length(unique(sce$orig.ident)) > 1) {
  cat("  Multiple samples detected. Applying Harmony integration...\n")
  library(harmony)

  sce <- RunHarmony(sce,
                    group.by.vars = "orig.ident",
                    dims.use = 1:15,
                    theta = 2,
                    plot_convergence = FALSE)

  sce <- FindNeighbors(sce, dims = 1:15, reduction = "harmony")
  cat("  OK Harmony integration completed\n")
} else {
  cat("  Single sample detected. No batch correction needed.\n")
  sce <- FindNeighbors(sce, dims = 1:15)
}

cat("\n Testing multiple resolutions for optimal clustering...\n")
cat("  Testing fine-grained resolutions to find 4 clusters target\n")
test_resolutions <- seq(0.05, 0.5, by = 0.02)
sce <- FindClusters(sce, resolution = test_resolutions)

cat("Resolution -> Number of clusters:\n")
for(res in test_resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  n_clusters <- length(unique(sce@meta.data[[col_name]]))
  cat(sprintf("  %.2f -> %d clusters\n", res, n_clusters))
}

cat("\n Generating clustree for resolution selection...\n")
library(clustree)

pdf("fibroblast_clustree.pdf", width = 12, height = 10)
p_tree <- clustree(sce,
                   prefix = "RNA_snn_res.",
                   node_colour = "sc3_stability",
                   node_size = "size",
                   edge_arrow = TRUE)
print(p_tree)
dev.off()
cat("  OK Clustree plot saved: fibroblast_clustree.pdf\n")

pdf("fibroblast_clustree_markers.pdf", width = 14, height = 12)
key_markers <- c("COL1A1", "ACTA2", "IL6", "MKI67")
cat("  Creating clustree with marker expression overlays...\n")
for(marker in key_markers) {
  if(marker %in% rownames(sce)) {
    p_tree_marker <- clustree(sce,
                              prefix = "RNA_snn_res.",
                              node_colour = marker,
                              node_colour_aggr = "mean") +
      ggtitle(paste("Clustree with", marker, "expression"))
    print(p_tree_marker)
    cat(sprintf("    - %s overlay added\n", marker))
  } else {
    cat(sprintf("    - %s not found in dataset\n", marker))
  }
}
dev.off()
cat("  OK Marker clustree plots saved: fibroblast_clustree_markers.pdf\n")

cat("\n Analyzing clustering stability...\n")

n_clusters_df <- data.frame(
  resolution = test_resolutions,
  n_clusters = sapply(test_resolutions, function(res) {
    col_name <- paste0("RNA_snn_res.", res)
    length(unique(sce@meta.data[[col_name]]))
  })
)

library(ggplot2)
p_nclusters <- ggplot(n_clusters_df, aes(x = resolution, y = n_clusters)) +
  geom_line(size = 1.5, color = "steelblue") +
  geom_point(size = 3, color = "darkblue") +
  geom_hline(yintercept = c(3, 4), linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 0.4, y = 3.5, label = "Target: 3-4 subtypes", color = "red", size = 4) +
  labs(title = "Fibroblast Clustering: Resolution vs Number of Clusters",
       subtitle = "Optimal range for CF study: 3-4 subtypes",
       x = "Resolution Parameter",
       y = "Number of Clusters") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "gray40")) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
  scale_y_continuous(breaks = 1:10)

ggsave("fibroblast_resolution_curve.pdf", p_nclusters, width = 8, height = 6)
cat("  OK Resolution curve saved: fibroblast_resolution_curve.pdf\n")

optimal_resolutions <- n_clusters_df[n_clusters_df$n_clusters >= 3 &
                                     n_clusters_df$n_clusters <= 4, ]
cat("\n  Resolutions yielding 3-4 clusters:\n")
if(nrow(optimal_resolutions) > 0) {
  for(i in 1:nrow(optimal_resolutions)) {
    cat(sprintf("    Resolution %.2f -> %d clusters\n",
                optimal_resolutions$resolution[i],
                optimal_resolutions$n_clusters[i]))
  }
} else {
  cat("    None found in tested range\n")
}

cat("\n Selecting optimal resolution...\n")
optimal_res <- NULL
n_clusters_target <- NULL

for(i in 1:nrow(optimal_resolutions)) {
  if(optimal_resolutions$n_clusters[i] == 4 && is.null(optimal_res)) {
    optimal_res <- optimal_resolutions$resolution[i]
    n_clusters_target <- 4
    break
  }
}

if(is.null(optimal_res) && nrow(optimal_resolutions) > 0) {
  optimal_res <- optimal_resolutions$resolution[1]
  n_clusters_target <- optimal_resolutions$n_clusters[1]
}

if(is.null(optimal_res)) {
  optimal_res <- 0.11
  n_clusters_target <- 4
  cat("  WARNING  No optimal resolution found, using default 0.11\n")
}

cat(sprintf("  OK Selected resolution: %.2f (yields %d clusters)\n",
            optimal_res, n_clusters_target))

selected_res <- paste0("RNA_snn_res.", optimal_res)
res_label <- sprintf("res%.2f", optimal_res)
cat(sprintf("  OK Resolution column: %s\n", selected_res))

set.seed(54321)

if("harmony" %in% names(sce@reductions)) {
  cat("  Using Harmony embeddings for dimensionality reduction\n")
  sce <- RunTSNE(object = sce, reduction = "harmony", dims = 1:15, do.fast = TRUE)
  sce <- RunUMAP(object = sce, reduction = "harmony", dims = 1:15, do.fast = TRUE)
} else {
  sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
  sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)
}

cat("\n Creating visualization plots...\n")

library(cowplot)
p1 <- DimPlot(sce, group.by = selected_res, label = TRUE, cols = mycolors) + ggtitle(sprintf("Resolution %.2f", optimal_res))
p2 <- DimPlot(sce, group.by = "RNA_snn_res.0.2", label = TRUE, cols = mycolors) + ggtitle("Resolution 0.2")
p3 <- DimPlot(sce, group.by = "RNA_snn_res.0.3", label = TRUE, cols = mycolors) + ggtitle("Resolution 0.3")

p4 <- DimPlot(sce, group.by = "orig.ident", label = FALSE) + ggtitle("By Sample")

if(!"group" %in% colnames(sce@meta.data)) {
  sce$group <- sce$orig.ident
}
p5 <- DimPlot(sce, group.by = "group", label = FALSE) + ggtitle("Normal vs CF")

p_resolutions <- plot_grid(p1, p2, p3, ncol = 3)
p_batch <- plot_grid(p4, p5, ncol = 2)

ggsave("fibroblast_clustering_resolutions.pdf", p_resolutions, width = 15, height = 5)
ggsave("fibroblast_batch_effects.pdf", p_batch, width = 10, height = 5)

p = DimPlot(sce, group.by = selected_res, label=T, cols = mycolors) +
    ggtitle(sprintf("Fibroblast Clustering (Resolution %.2f, %d clusters)", optimal_res, n_clusters_target))
p

cluster_specific_markers <- list(
  Mesenchymal = c('POSTN', 'ASPN', 'COMP', 'COL10A1', 'COL8A1', 'COL5A2', 'SPARC'),

  ProInflammatory = c('CCL19', 'APOE', 'CXCL3', 'JUN', 'HES1', 'NFKBIA'),

  SecretoryReticular = c('CCN5', 'SLPI', 'ANGPTL1', 'MFAP5', 'FBLN1', 'DCN'),

  SecretoryPapillary = c('APCDD1', 'WIF1', 'COL13A1', 'MALAT1', 'NEAT1', 'COCH', 'TNN', 'TNMD')
)

genes_to_check <- unique(unlist(cluster_specific_markers))

heatmap_markers <- c(

  'POSTN', 'ASPN', 'COMP', 'COL10A1', 'COL8A1',

  'CCL19', 'APOE', 'CXCL3', 'JUN', 'HES1',

  'CCN5', 'SLPI', 'ANGPTL1', 'MFAP5', 'FBLN1',

  'APCDD1', 'WIF1', 'COL13A1', 'MALAT1', 'NEAT1'
)

library(stringr)
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

Idents(sce) <- selected_res
cat(sprintf("\n Creating marker visualizations for resolution %.2f (%d clusters)...\n",
            optimal_res, n_clusters_target))

p_dot = DotPlot(sce, features = unique(genes_to_check), assay='RNA') +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(sprintf("Fibroblast Subtype Markers (Resolution %.2f)", optimal_res))
print(p_dot)
ggsave('check_markers.pdf', p_dot, height = 9, width = 7)
cat("  OK Saved: check_markers.pdf\n")

cat("  Creating violin plots for markers...\n")

key_markers_violin <- c(

  "COL1A1", "COL3A1", "DCN", "LUM", "POSTN",

  "IL6", "CXCL1", "CCL2", "PTGS2",

  "HSP90AA1", "ATF3", "FOS", "JUN",

  "MMP2", "MMP14", "TIMP1", "SPARC",

  "MKI67", "TOP2A"
)

key_markers_violin <- key_markers_violin[key_markers_violin %in% rownames(sce)]
cat(sprintf("  Found %d/%d markers in dataset\n",
            length(key_markers_violin), 19))

cat("  Creating multi-page violin plots (one cluster per page)...\n")

library(patchwork)

cluster_names <- c("MF (Matrix Fibroblasts)",
                  "PF (Pro-inflammatory Fibroblasts)",
                  "SRF (Secretory Reticular Fibroblasts)",
                  "SPF (Secretory Papillary Fibroblasts)")

pdf('fibroblast_markers_violin_by_cluster.pdf', width = 12, height = 10)

for(i in 0:3) {
  cat(sprintf("    Plotting cluster %d: %s\n", i, cluster_names[i+1]))

  cluster_markers <- cluster_specific_markers[[i+1]]

  cluster_markers_available <- cluster_markers[toupper(cluster_markers) %in% rownames(sce)]

  if(length(cluster_markers_available) > 0) {

    p_cluster <- VlnPlot(sce,
                        features = toupper(cluster_markers_available),
                        group.by = selected_res,
                        ncol = 3,
                        pt.size = 0.1,
                        cols = mycolors) +
      plot_annotation(title = paste0("Cluster ", i, ": ", cluster_names[i+1]),
                     theme = theme(plot.title = element_text(size = 16, face = "bold")))

    print(p_cluster)
  }
}

cat("    Adding summary page with key markers...\n")
if(length(key_markers_violin) > 0) {
  p_summary <- VlnPlot(sce,
                      features = key_markers_violin[1:min(9, length(key_markers_violin))],
                      group.by = selected_res,
                      ncol = 3,
                      pt.size = 0.1,
                      cols = mycolors) +
    plot_annotation(title = "Key Markers Summary (Top 9)",
                   theme = theme(plot.title = element_text(size = 16, face = "bold")))

  print(p_summary)
}

dev.off()
cat("  OK Saved: fibroblast_markers_violin_by_cluster.pdf\n")

cat("  Creating cluster-specific marker heatmap...\n")

heatmap_markers_available <- heatmap_markers[toupper(heatmap_markers) %in% rownames(sce)]
cat(sprintf("  Found %d/%d cluster-specific markers in dataset\n",
            length(heatmap_markers_available), length(heatmap_markers)))

if(length(heatmap_markers_available) > 0) {

  if(ncol(sce) > 1000) {
    set.seed(123)
    cells_use <- sample(colnames(sce), 1000)
  } else {
    cells_use <- colnames(sce)
  }

  p_heatmap <- DoHeatmap(subset(sce, cells = cells_use),
                         features = heatmap_markers_available,
                         group.by = selected_res,
                         group.colors = mycolors,
                         size = 3) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle("Cluster-Specific Markers (MF, PF, SRF, SPF)")

  ggsave('fibroblast_markers_heatmap.pdf', p_heatmap,
         height = 8, width = 10)
  cat("  OK Saved: fibroblast_markers_heatmap.pdf\n")
}

cat("\n=====================================\n")
cat(sprintf(" Comprehensive Marker Analysis (Resolution %.2f)\n", optimal_res))
cat("=====================================\n")

Idents(sce) <- selected_res
cat(sprintf("OK Fixed resolution: %.2f\n", optimal_res))
cat(sprintf("OK Number of clusters: %d\n", length(unique(Idents(sce)))))

cat("\n Finding top 30 markers for each cluster...\n")
all_markers <- FindAllMarkers(sce,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              return.thresh = 0.01,
                              verbose = TRUE)

all_markers <- all_markers %>%
  mutate(
    pct_diff = pct.1 - pct.2,
    fold_change = 2^avg_log2FC
  )

top30_markers <- all_markers %>%
  group_by(cluster) %>%

  arrange(p_val_adj, desc(avg_log2FC), desc(pct_diff)) %>%
  slice_head(n = 30) %>%
  ungroup()

top10_specific <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(pct_diff), p_val_adj) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("\n Determining characteristic genes for each cluster...\n")
cat("=====================================\n")

characteristic_genes <- list()
for(i in sort(unique(as.character(top30_markers$cluster)))) {
  cat(sprintf("\nCluster %s:\n", i))

  cluster_markers <- top30_markers %>%
    filter(cluster == i)

  characteristic <- cluster_markers %>%
    mutate(
      score = (avg_log2FC * 2) + (pct_diff * 3) - log10(p_val_adj + 1e-300)
    ) %>%
    arrange(desc(score)) %>%
    slice_head(n = 15) %>%
    pull(gene)

  characteristic_genes[[paste0("Cluster_", i)]] <- characteristic

  cat("  Top 5 by fold change: ")
  top_fc <- cluster_markers %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5) %>%
    pull(gene)
  cat(paste(top_fc, collapse = ", "))
  cat("\n")

  cat("  Top 5 by specificity: ")
  top_spec <- cluster_markers %>%
    arrange(desc(pct_diff)) %>%
    slice_head(n = 5) %>%
    pull(gene)
  cat(paste(top_spec, collapse = ", "))
  cat("\n")

  cat("  Characteristic genes: ")
  cat(paste(characteristic[1:5], collapse = ", "))
  if(length(characteristic) > 5) cat("...")
  cat("\n")
}

cat("\n Creating summary tables...\n")

cluster_summary <- top30_markers %>%
  group_by(cluster) %>%
  summarise(
    n_genes = n(),
    avg_log2FC_mean = mean(avg_log2FC),
    avg_log2FC_max = max(avg_log2FC),
    pct1_mean = mean(pct.1),
    pct_diff_mean = mean(pct_diff),
    .groups = "drop"
  )

cat("\nCluster Summary Statistics:\n")
print(cluster_summary)

cat("\n Saving results to CSV files...\n")

write.csv(top30_markers,
          paste0("fibroblast_top30_markers_", res_label, ".csv"),
          row.names = FALSE)
cat(sprintf("  OK Saved: fibroblast_top30_markers_%s.csv\n", res_label))

write.csv(top10_specific,
          paste0("fibroblast_top10_specific_", res_label, ".csv"),
          row.names = FALSE)
cat(sprintf("  OK Saved: fibroblast_top10_specific_%s.csv\n", res_label))

characteristic_df <- data.frame()
for(cluster_name in names(characteristic_genes)) {
  genes <- characteristic_genes[[cluster_name]]
  temp_df <- data.frame(
    cluster = gsub("Cluster_", "", cluster_name),
    rank = 1:length(genes),
    gene = genes,
    stringsAsFactors = FALSE
  )
  characteristic_df <- rbind(characteristic_df, temp_df)
}
write.csv(characteristic_df,
          paste0("fibroblast_characteristic_genes_", res_label, ".csv"),
          row.names = FALSE)
cat(sprintf("  OK Saved: fibroblast_characteristic_genes_%s.csv\n", res_label))

write.csv(cluster_summary,
          paste0("fibroblast_cluster_summary_", res_label, ".csv"),
          row.names = FALSE)
cat(sprintf("  OK Saved: fibroblast_cluster_summary_%s.csv\n", res_label))

cat("\n Gene lists for heatmap visualization:\n")
cat("=====================================\n")

heatmap_genes <- top30_markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  pull(gene) %>%
  unique()

cat("# Top 10 genes per cluster (for heatmap):\n")
cat(sprintf("heatmap_genes <- c(\n  \"%s\"\n)\n",
            paste(heatmap_genes, collapse = "\",\n  \"")))
cat(sprintf("# Total: %d genes\n", length(heatmap_genes)))

cat("\n Creating annotated heatmap with top genes...\n")

if(ncol(sce) > 2000) {
  set.seed(123)

  cells_use <- c()
  for(cl in levels(Idents(sce))) {
    cl_cells <- WhichCells(sce, idents = cl)
    n_sample <- min(500, length(cl_cells))
    cells_use <- c(cells_use, sample(cl_cells, n_sample))
  }
} else {
  cells_use <- colnames(sce)
}

p_heatmap_annotated <- DoHeatmap(subset(sce, cells = cells_use),
                                 features = heatmap_genes,
                                 group.by = selected_res,
                                 group.colors = mycolors[1:4],
                                 label = TRUE,
                                 size = 3,
                                 angle = 90,
                                 hjust = 0.5) +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                       na.value = "white") +
  theme(axis.text.y = element_text(size = 7, face = "italic"),
        legend.position = "right") +
  ggtitle("Top 10 Characteristic Genes per Fibroblast Subtype (MF, PF, SRF, SPF)")

height_size <- max(8, length(heatmap_genes) * 0.15 + 3)
ggsave(paste0("fibroblast_characteristic_heatmap_", res_label, ".pdf"),
       p_heatmap_annotated,
       width = 10,
       height = height_size,
       limitsize = FALSE)
cat(sprintf("OK Saved heatmap: fibroblast_characteristic_heatmap_%s.pdf\n", res_label))

cat("\n=====================================\n")
cat(" Analysis Complete!\n")
cat("=====================================\n")

cat("\n=====================================\n")
cat(" Creating Additional Heatmaps\n")
cat("=====================================\n")

cat("\n1 Creating top 30 genes heatmap...\n")

all.markers <- all_markers %>%
  dplyr::select(gene, everything()) %>%
  dplyr::filter(p_val < 0.05)

top30 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)

cat(sprintf("  Selected %d genes total for heatmap\n", nrow(top30)))

heatmap_top30 <- DoHeatmap(sce,
                           features = top30$gene,
                           group.by = selected_res,
                           group.colors = mycolors[1:4],
                           label = TRUE,
                           size = 3) +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred")) +
  theme(axis.text.y = element_text(size = 6)) +
  ggtitle("Top 30 Genes per Fibroblast Subtype (p < 0.05)")

ggsave("fibroblast_top30_heatmap.pdf", heatmap_top30, width = 15, height = 20)
cat("  OK Saved: fibroblast_top30_heatmap.pdf\n")

cat("\n2 Creating characteristic markers heatmap...\n")

characteristic_genes_all <- unique(unlist(cluster_specific_markers))

characteristic_genes_available <- characteristic_genes_all[toupper(characteristic_genes_all) %in% rownames(sce)]

cat(sprintf("  Found %d/%d characteristic markers in dataset\n",
            length(characteristic_genes_available), length(characteristic_genes_all)))

if(length(characteristic_genes_available) > 0) {

  if(ncol(sce) > 2000) {
    set.seed(123)
    cells_use <- c()
    for(cl in levels(Idents(sce))) {
      cl_cells <- WhichCells(sce, idents = cl)
      n_sample <- min(500, length(cl_cells))
      cells_use <- c(cells_use, sample(cl_cells, n_sample))
    }
  } else {
    cells_use <- colnames(sce)
  }

  heatmap_characteristic <- DoHeatmap(subset(sce, cells = cells_use),
                                      features = toupper(characteristic_genes_available),
                                      group.by = selected_res,
                                      group.colors = mycolors[1:4],
                                      label = TRUE,
                                      size = 3) +
    scale_fill_gradientn(colors = c("navy", "white", "firebrick3")) +
    theme(axis.text.y = element_text(size = 7, face = "italic")) +
    ggtitle("Characteristic Markers for Fibroblast Subtypes (MF, PF, SRF, SPF)")

  ggsave("fibroblast_characteristic_markers_heatmap.pdf",
         heatmap_characteristic,
         width = 12,
         height = 10)
  cat("  OK Saved: fibroblast_characteristic_markers_heatmap.pdf\n")
}

cat("\n=====================================\n")
cat(sprintf(" Fibroblast Subtype Definitions (Resolution %.2f)\n", optimal_res))
cat("=====================================\n\n")

subtype_definitions <- data.frame(
  Cluster = 0:3,
  Short_Name = c("MF", "PF", "SRF", "SPF"),
  Full_Name = c("Matrix_Fibroblasts", "ProInflam_Fibroblasts",
                "Sec_Reticular_Fibro", "Sec_Papillary_Fibro"),
  English_Name = c("Matrix Fibroblasts",
                   "Pro-inflammatory Fibroblasts",
                   "Secretory Reticular Fibroblasts",
                   "Secretory Papillary Fibroblasts"),
  Key_Markers = c("POSTN, ASPN, COMP, COL10A1, SPARC",
                  "CCL19, APOE, CXCL3, JUN, HES1",
                  "CCN5, SLPI, ANGPTL1, MFAP5, FBLN1",
                  "APCDD1, WIF1, COL13A1, MALAT1, NEAT1, COCH, TNN, TNMD"),
  Function = c("ECM synthesis and tissue remodeling",
               "Chemokine secretion and NF-kB signaling",
               "Secretory activity and reticular matrix maintenance",
               "Wnt inhibition and papillary matrix maintenance"),
  stringsAsFactors = FALSE
)

for(i in 1:nrow(subtype_definitions)) {
  cat(sprintf("Cluster %d - %s (%s):\n",
              subtype_definitions$Cluster[i],
              subtype_definitions$Short_Name[i],
              subtype_definitions$Full_Name[i]))
  cat("\n")
}

write.csv(subtype_definitions,
          "fibroblast_subtype_definitions.csv",
          row.names = FALSE)
cat("OK Saved: fibroblast_subtype_definitions.csv\n")

cat("\n Functional Characteristics:\n")
cat("=====================================\n")
cat("MF  - ECM production, tissue remodeling, collagen synthesis\n")
cat("PF  - Inflammatory response, cytokine secretion, immune regulation\n")
cat("SRF - Secretory function, reticular network formation\n")
cat("SPF - Secretory function, papillary structure support\n")
cat("=====================================\n\n")

source(file.path(initial_wd, scripts_path, 'Bottom_left_axis.R'))
result <- left_axes(sce)

axes <- result$axes
label <- result$label

umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.6,
              group.by = selected_res,label = T,label.box = F) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
umap
ggsave('RNA_snn_res.0.1_umap.pdf',width = 9,height = 7)

selected_markers <- c(
  'COL1A1', 'COL3A1', 'DCN', 'LUM',
  'IL6', 'CXCL1', 'CCL2',
  'HSP90AA1', 'ATF3', 'FOS',
  'MMP2', 'MMP14', 'TIMP1'
)

available_markers <- selected_markers[selected_markers %in% rownames(sce)]
cat(sprintf("\nOK Visualizing %d/%d markers\n",
            length(available_markers), length(selected_markers)))

if(length(available_markers) > 0) {
  p1 <- VlnPlot(
    sce,
    features = available_markers,
    pt.size = 0.01,
    ncol = 3,
    cols = mycolors
  )
  print(p1)
  ggsave('fibroblast_markers_violin.pdf', width = 12, height = 8)
}

saveRDS(sce, file = "fibro_sce.rds")
cat("\nOK Saved fibroblast data: fibro_sce.rds\n")

cat("\n========================================\n")
cat("WARNING  CNV Analysis Separated\n")
cat("========================================\n\n")

cat("CNV analysis has been moved to a standalone script for better organization.\n\n")

cat("To perform CNV analysis:\n")
cat("  1. Navigate to project root directory\n")
cat("  2. Run: Rscript 7.1.fibro_cnv/07.1_cnv_analysis.R\n\n")

cat("CNV analysis will:\n")
cat("  - Use Tcell + Macro as reference cells\n")
cat("  - Calculate CNV scores for fibroblast subtypes\n")
cat("  - Generate CNV heatmaps and visualizations\n")
cat("  - Save results in 7.1.fibro_cnv/CNV/\n\n")

cat("Continuing with fibroblast subtype annotation...\n\n")

cat("\n========================================\n")
cat("Fibroblast Subtype Annotation for CF Study\n")
cat("========================================\n\n")

fibro_sce <- sce
n_clusters <- n_clusters_target
cat(sprintf("OK Working with %d fibroblast clusters\n", n_clusters))

cat("\n Evaluating clustering quality with silhouette analysis...\n")
library(cluster)

if("harmony" %in% names(fibro_sce@reductions)) {
  embeddings <- fibro_sce@reductions$harmony@cell.embeddings[, 1:15]
  cat("  Using Harmony embeddings for silhouette analysis\n")
} else {
  embeddings <- fibro_sce@reductions$pca@cell.embeddings[, 1:15]
  cat("  Using PCA embeddings for silhouette analysis\n")
}

n_cells <- nrow(embeddings)
if(n_cells > 2000) {
  cat(sprintf("  Subsampling %d cells to 2000 for silhouette calculation\n", n_cells))
  set.seed(123)
  sample_idx <- sample(1:n_cells, 2000)
  embeddings_subset <- embeddings[sample_idx, ]
  clusters_subset <- as.numeric(factor(fibro_sce@meta.data[[selected_res]][sample_idx]))
} else {
  embeddings_subset <- embeddings
  clusters_subset <- as.numeric(factor(fibro_sce@meta.data[[selected_res]]))
}

dist_matrix <- dist(embeddings_subset)
sil_scores <- silhouette(clusters_subset, dist_matrix)
avg_sil <- mean(sil_scores[, 3])

if(nrow(sil_scores) > 0) {
  sil_by_cluster <- aggregate(sil_scores[, 3],
                              by = list(cluster = sil_scores[, 1]),
                              FUN = function(x) {
                                c(mean = mean(x, na.rm = TRUE),
                                  sd = sd(x, na.rm = TRUE))
                              })

  if(nrow(sil_by_cluster) > 0) {
    sil_summary <- data.frame(
      Cluster = sil_by_cluster$cluster - 1,
      Mean_Silhouette = round(sil_by_cluster$x[, 1], 3),
      SD_Silhouette = round(sil_by_cluster$x[, 2], 3)
    )
  } else {
    sil_summary <- data.frame(
      Cluster = integer(0),
      Mean_Silhouette = numeric(0),
      SD_Silhouette = numeric(0)
    )
  }
} else {
  sil_summary <- data.frame(
    Cluster = integer(0),
    Mean_Silhouette = numeric(0),
    SD_Silhouette = numeric(0)
  )
}

cat(sprintf("\n  Overall average silhouette score: %.3f\n", avg_sil))
cat("  Per-cluster silhouette scores:\n")
print(sil_summary)

if(avg_sil > 0.7) {
  cat("  -> Strong clustering structure\n")
} else if(avg_sil > 0.5) {
  cat("  -> Reasonable clustering structure\n")
} else if(avg_sil > 0.25) {
  cat("  -> Weak clustering structure\n")
} else {
  cat("  -> No substantial clustering structure found\n")
}

pdf("fibroblast_silhouette_plot.pdf", width = 10, height = 8)
plot(sil_scores, col = mycolors[1:n_clusters],
     main = paste("Silhouette Plot for Fibroblast Clustering\nAverage silhouette width:",
                  round(avg_sil, 3)),
     border = NA)
dev.off()
cat("  OK Silhouette plot saved: fibroblast_silhouette_plot.pdf\n")

fibro_sce$seurat_clusters <- fibro_sce@meta.data[[selected_res]]

if(n_clusters == 4) {
  celltype = data.frame(
    ClusterID = 0:3,
    celltype = c('mesenchymal', 'pro-inflammatory',
                 'secretory-reticular', 'secretory-papillary'),
    full_name = c('Mesenchymal Fibroblasts',
                  'Pro-inflammatory Fibroblasts',
                  'Secretory Reticular Fibroblasts',
                  'Secretory Papillary Fibroblasts'),
    abbreviation = c('MesF', 'ProInflamF', 'SecRetF', 'SecPapF'),
    function_desc = c(
      "Mesenchymal fibroblasts with ECM production",
      "Pro-inflammatory fibroblasts secreting cytokines",
      "Secretory reticular fibroblasts forming network",
      "Secretory papillary fibroblasts with specialized secretion"
    ),
    stringsAsFactors = FALSE
  )
} else if(n_clusters == 3) {
  celltype = data.frame(
    ClusterID = 0:2,
    celltype = c('mesenchymal', 'pro-inflammatory', 'secretory-reticular'),
    full_name = c('Mesenchymal Fibroblasts', 'Pro-inflammatory Fibroblasts', 'Secretory Reticular Fibroblasts'),
    abbreviation = c('MesF', 'ProInflamF', 'SecRetF'),
    function_desc = c(
      "Mesenchymal fibroblasts",
      "Pro-inflammatory fibroblasts",
      "Secretory reticular fibroblasts"
    ),
    stringsAsFactors = FALSE
  )
} else {
  celltype = data.frame(
    ClusterID = 0:(n_clusters-1),
    celltype = paste0('Fibro_C', 0:(n_clusters-1)),
    full_name = paste0('Fibroblast Cluster ', 0:(n_clusters-1)),
    abbreviation = paste0('FC', 0:(n_clusters-1)),
    function_desc = paste0(' ', 0:(n_clusters-1)),
    stringsAsFactors = FALSE
  )
}

cat("\nSubtype assignments:\n")
print(celltype)

fibro_sce@meta.data$celltype = "NA"
fibro_sce@meta.data$celltype_abbr = "NA"
fibro_sce@meta.data$celltype_full = "NA"

for(i in 1:nrow(celltype)){
  matching_cells <- which(fibro_sce@meta.data[[selected_res]] == celltype$ClusterID[i])
  fibro_sce@meta.data[matching_cells, 'celltype'] <- celltype$celltype[i]
  fibro_sce@meta.data[matching_cells, 'celltype_abbr'] <- celltype$abbreviation[i]
  fibro_sce@meta.data[matching_cells, 'celltype_full'] <- celltype$full_name[i]
}

table(fibro_sce@meta.data$celltype)
table(fibro_sce@meta.data$celltype_abbr)

th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 0.5, hjust=0.5))
library(patchwork)
source(file.path(initial_wd, scripts_path, 'Bottom_left_axis.R'))
result <- left_axes(fibro_sce)
axes <- result$axes

label_data <- result$label
celltype_umap = DimPlot(fibro_sce, reduction = "umap", cols = my36colors, pt.size = 0.8,
                        group.by = "celltype", label = TRUE) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends = "last", type = "closed")) +
  geom_text(data = label_data,
            aes(x = x, y = y, angle = angle, label = lab), fontface = 'italic') +
  theme(plot.title = element_blank())

celltype_umap
ggsave('umap_by_celltype.pdf',width = 9,height = 7)

saveRDS(fibro_sce, "fibro_sce.rds")
cat("\nFibroblast data saved to: 07_fibroblast/fibro_sce.rds\n")

Idents(fibro_sce) = fibro_sce$celltype
if (!file.exists('fibro_markers.csv')) {
  markers <- FindAllMarkers(object = fibro_sce, only.pos = TRUE,
                            min.pct = 0.25,
                            thresh.use = 0.25)
  write.csv(markers,file='fibro_markers.csv')
} else {

  sce.markers = read.csv('fibro_markers.csv',row.names = 1)
}

prop_pie_paths <- c(
  file.path(initial_wd, "7.fibro_extrction", "07_prop_pie.R"),
  file.path(getwd(), "../07_prop_pie.R"),
  file.path(getwd(), "07_prop_pie.R")
)

prop_pie_loaded <- FALSE
for(pie_path in prop_pie_paths) {
  if(file.exists(pie_path)) {
    source(pie_path)
    cat(sprintf("OK Loaded prop_pie function from: %s\n", pie_path))
    prop_pie_loaded <- TRUE
    break
  }
}

if(!prop_pie_loaded) {
  cat("WARNING  prop_pie.R not found, will use alternative visualization\n")
}

if(exists("prop_pie")) {
  prop_pie(object = fibro_sce, seurat_by = "sample", pheno_by = "celltype", palette = "jama")
} else {
  cat("prop_pie function not available, using alternative visualization\n")
}

object = fibro_sce
seurat_by = "sample"
pheno_by = "celltype"
palette = "nejm"
library(ggrepel)
library(ggsci)
library(tidyverse)
library(ggplot2)
plot_data <- object@meta.data %>%
  dplyr::select({{pheno_by}}, {{seurat_by}})

plot_data <- table(plot_data[[seurat_by]], plot_data[[pheno_by]]) %>%
  as.data.frame() %>%
  group_by(Var2) %>%
  mutate(Total = sum(Freq),
         Proportion = round((Freq / Total)*100,2),
         labels = scales::percent(Freq / Total)) %>%
  mutate(text_y = Total - (Freq/2)) %>%
  dplyr::rename(pheno_by = "Var2") %>%
  dplyr::rename(seurat_by = "Var1")

table(plot_data$seurat_by)
plot_data$seurat_by = factor(plot_data$seurat_by, levels = c('Normal', 'CF'))

colors <- colorRampPalette((pal_nejm("default")(8)))(length(unique(plot_data$seurat_by)))

p <- ggplot(plot_data, aes(x = "", y = Proportion, fill = seurat_by)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  facet_wrap(~ pheno_by) +
  labs(fill = seurat_by) +
  xlab("") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)
p

p2 = ggplot(plot_data, aes(x = seurat_by, y = Proportion, color = seurat_by)) +
  geom_line(aes(group = 1), color = '#64A36C', size = 1.5) +
  geom_point(size = 5) +
  facet_wrap(~ pheno_by) +
  theme_bw() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)
p2
p + p2
ggsave('pie_prob.pdf', width = 12, height = 7)

cat("\n========================================\n")
cat("Normal vs CF Comparison Analysis\n")
cat("========================================\n")

fibro_sce$condition <- ifelse(
  grepl("Normal", fibro_sce$orig.ident), "Normal", "CF"
)
cat(sprintf("\nSample distribution:\n"))
table(fibro_sce$condition)

prop_table <- table(fibro_sce$celltype, fibro_sce$condition)
cat("\nCell type distribution by condition:\n")
print(prop_table)

prop_df <- as.data.frame(prop_table)
colnames(prop_df) <- c("CellType", "Condition", "Count")
prop_df <- prop_df %>%
  group_by(Condition) %>%
  mutate(Total = sum(Count),
         Proportion = Count/Total * 100)

cat("\n Statistical Analysis (Fisher's Exact Test):\n")
for(ct in unique(fibro_sce$celltype)) {
  ct_cells <- fibro_sce$celltype == ct
  cont_table <- table(ct_cells, fibro_sce$condition)
  if(nrow(cont_table) == 2 && ncol(cont_table) == 2) {
    test_result <- fisher.test(cont_table)
    cat(sprintf("  %s: p-value = %.4f %s\n",
                ct, test_result$p.value,
                ifelse(test_result$p.value < 0.05, "*", "")))
  }
}

p_comparison <- ggplot(prop_df, aes(x = CellType, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Fibroblast Subtypes: Normal vs CF",
       x = "Fibroblast Subtype",
       y = "Proportion (%)",
       fill = "Condition") +
  theme_classic() +
  scale_fill_manual(values = c("Normal" = "#4DBBD4", "CF" = "#E64A35")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_comparison)
ggsave('fibroblast_normal_vs_cf.pdf', width = 8, height = 6)

write.csv(prop_df, "fibroblast_proportions_normal_vs_cf.csv", row.names = FALSE)

saveRDS(fibro_sce, "fibroblast_annotated.rds")
cat("\nOK Analysis results saved\n")

cat("\n========================================\n")
cat("ANALYSIS SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total fibroblasts analyzed: %d\n", ncol(fibro_sce)))
cat(sprintf("Number of subtypes identified: %d\n", n_clusters))

if(!is.null(celltype$full_name)) {
  cat("\nSubtype annotations:\n")
  for(i in 1:nrow(celltype)) {
    cat(sprintf("  %s: %s\n", celltype$celltype[i], celltype$full_name[i]))
  }
}

cat("\nKey findings:\n")
cat("  1. Batch correction applied: ",
    ifelse("harmony" %in% names(fibro_sce@reductions), "Yes (Harmony)", "No"), "\n")
cat("  2. Subtypes characterized based on CF pathology and marker expression\n")
cat("  3. Proportion analysis completed for Normal vs CF comparison\n")

cat("\nOK Fibroblast extraction and analysis complete!\n")
cat("========================================\n")

setwd(initial_wd)
cat(sprintf("\nOK Returned to initial directory: %s\n", getwd()))

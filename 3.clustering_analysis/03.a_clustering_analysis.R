rm(list = ls())
set.seed(123)
gc()

cat("========================================\n")
cat("Clustering Analysis - Part A\n")
cat("========================================\n\n")

cat("Loading required libraries...\n")

source('../scRNA_scripts/load_packages.R')
source('../scRNA_scripts/lib.R')

required_packages <- c("Seurat", "ggplot2", "dplyr", "clustree", "patchwork",
                      "ComplexHeatmap", "circlize")
load_packages(required_packages)

current_dir <- basename(getwd())
if(current_dir != "3.clustering_analysis") {
  stop("Error: Please run this script from the 3.clustering_analysis directory!\n",
       "Current directory: ", getwd(), "\n",
       "Expected: .../3.clustering_analysis")
}

cat(paste("Output directory:", getwd(), "\n"))

cat("\n========================================\n")
cat("Step 1: Loading Data")
cat("========================================\n")

data_paths <- c(
  "../2.harmony/sce_integrated.rds",
  "../2.harmony/02_harmony/sce.all_int.rds"
)

sce.all <- NULL
for(data_path in data_paths) {
  if(file.exists(data_path)) {
    cat(paste("Loading integrated data from:", data_path, "\n"))
    sce.all <- readRDS(data_path)
    cat(paste("  OK Successfully loaded!\n"))
    break
  }
}

if(is.null(sce.all)) {
  stop("Cannot find integrated data. Please run 02_batch_integration.R first.\n")
}

cat(paste("Total cells:", ncol(sce.all), "\n"))
cat(paste("Total genes:", nrow(sce.all), "\n"))

if("seurat_clusters" %in% colnames(sce.all@meta.data)) {
  cat(paste("Current clusters:", length(unique(sce.all$seurat_clusters)), "\n"))
  cat("Cluster distribution:\n")
  print(table(sce.all$seurat_clusters))
} else {
  cat("No existing clustering found.\n")
}

cat("\n========================================\n")
cat("Step 2: Checking Current Clustering")
cat("========================================\n")

if("seurat_clusters" %in% colnames(sce.all@meta.data)) {
  current_clusters <- length(unique(sce.all$seurat_clusters))
  cat(paste("OK Current number of clusters:", current_clusters, "\n"))

  cluster_table <- table(sce.all$seurat_clusters)
  cat("\nCluster size distribution:\n")
  for(i in names(cluster_table)) {
    cat(sprintf("  Cluster %2s: %5d cells (%.1f%%)\n",
                i, cluster_table[i],
                cluster_table[i]/ncol(sce.all)*100))
  }

  if(current_clusters == 17) {
    cat("\nOK Optimal clustering detected (17 clusters from Resolution=0.2)!\n")
    cat("  This matches the Harmony integration output.\n")
    cat("  Proceeding with marker gene identification and validation.\n")
  } else if(current_clusters > 20) {
    cat("\nWARNING Warning: Too many clusters detected (>20)!\n")
    cat("This may lead to over-clustering and difficulty in cell type annotation.\n")
  } else {
    cat("\nOK Clustering appears reasonable (<=20 clusters).\n")
  }

  has_existing_clustering <- TRUE
} else {
  cat("No existing clustering found. Will perform clustering analysis.\n")
  has_existing_clustering <- FALSE
}

cat("\n========================================\n")
cat("Step 3: Testing Different Resolution Parameters")
cat("========================================\n")

if(!"RNA_snn" %in% names(sce.all@graphs)) {
  cat("Computing neighbor graph...\n")
  if("harmony" %in% names(sce.all@reductions)) {
    cat("OK Using harmony reduction for neighbor graph (dims 1:15)\n")
    sce.all <- FindNeighbors(sce.all, reduction = "harmony", dims = 1:15)
  } else if("pca" %in% names(sce.all@reductions)) {
    cat("Using PCA reduction for neighbor graph (dims 1:30)\n")
    sce.all <- FindNeighbors(sce.all, reduction = "pca", dims = 1:30)
  } else {
    stop("No suitable reduction found for neighbor graph computation!")
  }
} else {
  cat("OK Neighbor graph already exists.\n")
}

skip_reclustering <- FALSE
if(has_existing_clustering && current_clusters == 17) {
  cat("\nOK Optimal clustering already present (17 clusters from Resolution=0.2)\n")
  cat("  Testing limited resolutions for comparison only.\n")
  skip_reclustering <- TRUE
  resolutions <- c(0.1, 0.15, 0.2, 0.25, 0.3)
} else {
  cat("Testing multiple resolutions to find optimal clustering.\n")
  resolutions <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8)
}

cat("Testing resolutions:", paste(resolutions, collapse = ", "), "\n\n")

for(res in resolutions) {
  cat(sprintf("  Resolution %.2f ... ", res))
  sce.all <- FindClusters(sce.all, resolution = res, verbose = FALSE)
  col_name <- paste0("RNA_snn_res.", res)
  n_clusters <- length(unique(sce.all@meta.data[[col_name]]))
  cat(sprintf("%2d clusters", n_clusters))
  if(res == 0.2) cat(" <- Recommended")
  cat("\n")
}

cat("\n========================================\n")
cat("Step 4: Clustering Tree Analysis")
cat("========================================\n")

suppressPackageStartupMessages(library(clustree))

pdf("clustree_analysis.pdf", width = 12, height = 10)
p_tree <- clustree(sce.all, prefix = "RNA_snn_res.")
print(p_tree)
dev.off()

cat("Clustree plot saved. Check for stable clustering patterns.\n")

cat("\n========================================\n")
cat("Step 5: Evaluating Clustering Quality")
cat("========================================\n")

clustering_stats <- data.frame(
  resolution = resolutions,
  n_clusters = numeric(length(resolutions)),
  min_cluster_size = numeric(length(resolutions)),
  median_cluster_size = numeric(length(resolutions)),
  silhouette_score = numeric(length(resolutions))
)

if(!"umap" %in% names(sce.all@reductions)) {
  cat("Computing UMAP...\n")
  sce.all <- RunUMAP(sce.all, dims = 1:30, verbose = FALSE)
}

library(cluster)
for(i in 1:length(resolutions)) {
  res <- resolutions[i]
  col_name <- paste0("RNA_snn_res.", res)
  clusters <- sce.all@meta.data[[col_name]]

  clustering_stats$n_clusters[i] <- length(unique(clusters))
  cluster_sizes <- table(clusters)
  clustering_stats$min_cluster_size[i] <- min(cluster_sizes)
  clustering_stats$median_cluster_size[i] <- median(cluster_sizes)

  if(ncol(sce.all) > 5000) {
    set.seed(123)
    sample_idx <- sample(1:ncol(sce.all), 5000)
  } else {
    sample_idx <- 1:ncol(sce.all)
  }

  umap_coords <- Embeddings(sce.all, "umap")[sample_idx,]
  clusters_sample <- as.numeric(clusters[sample_idx])

  if(length(unique(clusters_sample)) > 1) {
    sil <- silhouette(clusters_sample, dist(umap_coords))
    clustering_stats$silhouette_score[i] <- mean(sil[,3])
  }
}

cat("\nClustering Statistics:\n")
print(clustering_stats)

write.csv(clustering_stats, "clustering_statistics.csv", row.names = FALSE)

cat("\n========================================\n")
cat("Step 6: Recommending Optimal Resolution")
cat("========================================\n")

optimal_res <- 0.2
cat("\nOK Using user-specified resolution: 0.2 for detailed initial clustering\n")
cat("  This will enable fine-grained cell type identification\n")

res_02_stats <- clustering_stats[clustering_stats$resolution == 0.2,]
if(nrow(res_02_stats) > 0) {
  cat(paste("  - Number of clusters:", res_02_stats$n_clusters, "\n"))
  cat(paste("  - Minimum cluster size:", res_02_stats$min_cluster_size, "\n"))
  cat(paste("  - Median cluster size:", res_02_stats$median_cluster_size, "\n"))
  cat(paste("  - Silhouette score:", round(res_02_stats$silhouette_score, 3), "\n"))
}

cat("\n========================================\n")
cat("Step 7: Applying Optimal Clustering")
cat("========================================\n")

sce.all <- FindClusters(sce.all, resolution = optimal_res, verbose = FALSE)
cat(paste("Applied resolution", optimal_res, "clustering\n"))
cat(paste("Final number of clusters:", length(unique(sce.all$seurat_clusters)), "\n"))

cat("\n========================================\n")
cat("Step 8: Identifying Marker Genes")
cat("========================================\n")

cat("Finding marker genes for each cluster...\n")
cluster_markers <- FindAllMarkers(sce.all,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.5,
                                 test.use = "wilcox")

write.csv(cluster_markers, "all_cluster_markers.csv", row.names = FALSE)
cat(paste("Total marker genes identified:", nrow(cluster_markers), "\n"))

top20_markers_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top20_markers_per_cluster, "top20_markers_per_cluster.csv", row.names = FALSE)

cat("\nSelecting top 30 markers for each cluster...\n")
top30_markers_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top30_markers_per_cluster, "top30_markers_per_cluster.csv", row.names = FALSE)
cat(paste("Top 30 markers saved for", length(unique(top30_markers_per_cluster$cluster)), "clusters\n"))

cat("\nTop 10 markers for each cluster:\n")
for(i in sort(unique(cluster_markers$cluster))) {
  cat(paste0("\nCluster ", i, ":\n"))
  cluster_genes <- cluster_markers %>%
    filter(cluster == i) %>%
    arrange(desc(avg_log2FC)) %>%
    head(10) %>%
    pull(gene)
  cat(paste(cluster_genes, collapse = ", "), "\n")
}

cat("\n========================================\n")
cat("Step 9: Creating Heatmap of Top Genes")
cat("========================================\n")

heatmap_genes <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

cat(paste("Preparing heatmap with", length(heatmap_genes), "genes\n"))

if(!"scale.data" %in% names(sce.all@assays$RNA)) {
  cat("Scaling data for heatmap...\n")
  sce.all <- ScaleData(sce.all, features = heatmap_genes, verbose = FALSE)
}

expr_mat <- GetAssayData(sce.all, slot = "scale.data")
heatmap_genes_available <- heatmap_genes[heatmap_genes %in% rownames(expr_mat)]
expr_mat <- expr_mat[heatmap_genes_available, ]

cell_order <- order(sce.all$seurat_clusters)
expr_mat <- expr_mat[, cell_order]

cluster_annotation <- sce.all$seurat_clusters[cell_order]
sample_annotation <- sce.all$orig.ident[cell_order]

library(ComplexHeatmap)
library(circlize)

col_ann <- HeatmapAnnotation(
  Cluster = cluster_annotation,
  Sample = sample_annotation,
  col = list(
    Cluster = setNames(rainbow(length(unique(cluster_annotation))),
                       sort(unique(cluster_annotation))),
    Sample = c("CF" = "#E41A1C", "Normal" = "#377EB8")
  ),
  annotation_name_side = "left"
)

gene_cluster_assignment <- cluster_markers %>%
  filter(gene %in% heatmap_genes_available) %>%
  group_by(gene) %>%
  top_n(1, wt = avg_log2FC) %>%
  select(gene, cluster) %>%
  distinct()

gene_order <- gene_cluster_assignment %>%
  arrange(cluster, desc(gene)) %>%
  pull(gene)

expr_mat <- expr_mat[gene_order, ]

pdf("cluster_top_genes_heatmap.pdf", width = 20, height = 16)
ht <- Heatmap(expr_mat,
        name = "Expression",
        top_annotation = col_ann,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_split = cluster_annotation,
        row_split = gene_cluster_assignment$cluster[match(gene_order, gene_cluster_assignment$gene)],
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        column_title = "Cells grouped by cluster",
        row_title = "Top 10 genes per cluster",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"))
draw(ht)
dev.off()

cat("Heatmap saved to cluster_top_genes_heatmap.pdf\n")

cat("\n========================================\n")
cat("Step 10: Cluster Statistics Summary")
cat("========================================\n")

cluster_table <- table(sce.all$seurat_clusters)
cluster_stats <- data.frame(
  Cluster = names(cluster_table),
  NCells = as.numeric(cluster_table),
  Percentage = round(as.numeric(cluster_table)/ncol(sce.all)*100, 2)
)
cluster_stats <- cluster_stats[order(cluster_stats$NCells, decreasing = TRUE),]

cat("\nCluster distribution:\n")
print(cluster_stats)

cat("\n\nTop 5 marker genes per cluster (for reference):\n")
cat("================================================\n")
for(i in sort(unique(cluster_markers$cluster))) {
  cat(paste0("Cluster ", i, ": "))
  top5_genes <- cluster_markers %>%
    filter(cluster == i) %>%
    arrange(desc(avg_log2FC)) %>%
    head(5) %>%
    pull(gene)
  cat(paste(top5_genes, collapse = ", "), "\n")
}

cat("\n========================================\n")
cat("Saving intermediate results")
cat("========================================\n")

saveRDS(sce.all, "sce_clustered_annotated.rds")
saveRDS(cluster_markers, "cluster_markers_full.rds")
saveRDS(top30_markers_per_cluster, "top30_markers.rds")
saveRDS(cluster_stats, "cluster_stats.rds")

cat("\nPart A analysis completed!\n")
cat("Results saved. Ready for Part B visualization.\n")
cat("========================================\n")

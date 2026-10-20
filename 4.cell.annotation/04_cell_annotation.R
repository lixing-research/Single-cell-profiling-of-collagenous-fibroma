rm(list=ls())

current_dir <- basename(getwd())
initial_wd <- getwd()
if(current_dir == "4.cell.annotation") {
  source('../scRNA_scripts/lib.R')
  source('../scRNA_scripts/mycolors.R')
  base_path <- "../3.clustering_analysis/"
  scripts_path <- "../scRNA_scripts/"
} else {
  source('scRNA_scripts/lib.R')
  source('scRNA_scripts/mycolors.R')
  base_path <- "3.clustering_analysis/"
  scripts_path <- "scRNA_scripts/"
}

cat("========================================\n")
cat("Step 04: Cell Type Annotation")
cat("========================================\n\n")

data_file_paths <- c(
  paste0(base_path, "03c_result/sce_merged_clusters.rds"),
  paste0(base_path, "03b_result/sce_clustered_final.rds")
)

loaded <- FALSE
for(path in data_file_paths) {
  if(file.exists(path)) {
    sce.all.int = readRDS(path)
    cat(sprintf("OK Loaded data from: %s\n", path))
    loaded <- TRUE
    break
  }
}

if(!loaded) {
  stop(sprintf("ERROR Cannot find clustering results. Checked paths:\n%s\n",
               paste("  -", data_file_paths, collapse="\n")))
}

cat(sprintf("\n Data summary:\n"))
cat(sprintf("  Total cells: %d\n", ncol(sce.all.int)))
if("merged_celltype" %in% colnames(sce.all.int@meta.data)) {
  cat(sprintf("  Merged cell types: %d\n", length(unique(sce.all.int$merged_celltype))))
  cat(sprintf("  Cell types: %s\n", paste(sort(unique(sce.all.int$merged_celltype)), collapse=", ")))
}
if("original_clusters" %in% colnames(sce.all.int@meta.data)) {
  cat(sprintf("  Original clusters: %d\n", length(unique(sce.all.int$original_clusters))))
}

if(current_dir == "4.cell.annotation") {
  output_dir <- "./04_celltype"
} else {
  output_dir <- "4.cell.annotation/04_celltype"
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
scRNA=sce.all.int

genes_to_check = c(
  'COL1A1', 'DCN', 'LUM', 'CFD', 'PTGDS', 'OGN', 'HAPLN1', 'STEAP1',

  'PECAM1', 'VWF', 'CCL14', 'ECSCR', 'CCL21', 'CLDN5', 'SOX18', 'CLEC14A',

  'RGS5', 'ACTA2', 'TAGLN', 'MUSTN1', 'HIGD1B', 'COX4I2', 'FAM162B',

  'EPCAM', 'KRT8', 'KRT18', 'DEFB1', 'CCL27', 'LGALS7', 'LGALS7B', 'KRT16P6',

  'CD68', 'C1QA', 'AIF1', 'MS4A7', 'FOLR2', 'LYZ', 'FCGR3A',

  'CD3D', 'CD3E', 'CD52', 'TRBC2', 'LTB', 'CD48', 'S100B',

  'PDZRN4', 'SGCZ', 'NAV3', 'DCLK1', 'ROBO2', 'PTPRG',

  'PIP', 'SCGB2A2', 'DCD', 'ROPN1B', 'SNORC', 'AQP5',

  'DPP10', 'CNTN5', 'GRID2', 'KSR2', 'RBFOX1',

  'TPSB2', 'MS4A2', 'TPSAB1', 'CPA3', 'CTSG',

  'PCP4', 'SYNM', 'ACTG2', 'TCEAL2', 'KCNMB1'
)

library(stringr)
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

if("merged_celltype" %in% colnames(scRNA@meta.data)) {
  cat("\nUsing merged cell types for visualization\n")
  group_var <- "merged_celltype"
} else {
  cat("\nUsing original clusters for visualization\n")
  group_var <- "seurat_clusters"
}

p = DotPlot(scRNA, features = unique(genes_to_check),
            group.by = group_var,
            assay='RNA'  )  + coord_flip()
p
ggsave('check_last_markers.pdf',height = 12,width = 8)

source(file.path(initial_wd, scripts_path, 'Bottom_left_axis.R'))
result <- left_axes(scRNA)
axes <- result$axes
label <- result$label

cluster_var <- if("original_clusters" %in% colnames(scRNA@meta.data)) "original_clusters" else "seurat_clusters"

umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
              group.by = cluster_var,label = T,label.box = F) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank()) +
  ggtitle(paste("UMAP by", cluster_var))
umap
ggsave('original_clusters_umap.pdf',width = 9,height = 7)

sample_umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                     group.by = "sample") +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
sample_umap
ggsave('sample_umap.pdf',width = 9,height = 7)

umap+sample_umap

celltype_order <- c("Fib", "Endo", "PeriSMC", "Epi", "Macro", "Tcell",
                   "Schwann", "GlandEpi", "Neuro", "Mast", "MyoepiSMC")

if("merged_celltype" %in% colnames(scRNA@meta.data)) {
  cat("\nOK Using existing merged cell types from 03c analysis\n")
  scRNA@meta.data$celltype <- as.character(scRNA@meta.data$merged_celltype)

  scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype,
                                     levels = celltype_order[celltype_order %in% unique(scRNA@meta.data$celltype)])

  cat("\nCell type distribution:\n")
  celltype_table <- table(scRNA@meta.data$celltype)
  for(ct in levels(scRNA@meta.data$celltype)) {
    cat(sprintf("  %s: %d cells (%.1f%%)\n",
                ct, celltype_table[ct],
                celltype_table[ct]/ncol(scRNA)*100))
  }
} else {
  cat("\nWARNING merged_celltype not found, using seurat_clusters\n")
  scRNA@meta.data$celltype <- as.character(scRNA@meta.data$seurat_clusters)
}

th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 0.5, hjust=0.5))
library(patchwork)
celltype_colors <- setNames(my36colors[1:length(levels(scRNA@meta.data$celltype))],
                           levels(scRNA@meta.data$celltype))

celltype_umap =DimPlot(scRNA, reduction = "umap",
                       cols = celltype_colors,
                       pt.size = 0.8,
                       group.by = "celltype",label = T) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank()) +
  ggtitle("UMAP by Cell Types (Merged)")
celltype_umap
ggsave('umap_by_merged_celltype.pdf',width = 9,height = 7)

umap + sample_umap+celltype_umap
ggsave('combine_umap_merged.pdf',width = 21,height = 7)

cat("\n========================================\n")
cat("Saving annotated data...\n")
saveRDS(scRNA, "sce_annotated_merged.rds")
cat("OK Saved annotated Seurat object to sce_annotated_merged.rds\n")

celltype_info <- data.frame(
  Cell = colnames(scRNA),
  CellType = scRNA@meta.data$celltype,
  Sample = scRNA@meta.data$orig.ident
)
write.csv(celltype_info, "celltype_annotations.csv", row.names = FALSE)
cat("OK Exported cell type annotations to celltype_annotations.csv\n")

cat("\n========================================\n")
cat(" Cell type annotation completed successfully!\n")
cat("========================================\n")

setwd('../')

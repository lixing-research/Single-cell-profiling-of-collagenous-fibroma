#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
})

set.seed(123)

seurat_path     <- "4.cell.annotation/04_celltype/sce_annotated_merged.rds"
gene_order_file <- "gene_order_file.txt"
out_dir         <- "infercnv_output"

if (!file.exists(seurat_path))     stop("Seurat object not found: ", seurat_path)
if (!file.exists(gene_order_file)) stop("Gene order file not found: ", gene_order_file)

full <- readRDS(seurat_path)

pick <- function(obj, label, n) {
  cells <- colnames(obj)[obj$celltype == label]
  if (length(cells) < n) stop("Fewer than ", n, " cells for ", label)
  sample(cells, n)
}

cf_cells     <- pick(subset(full, orig.ident == "CF"),     "Fib",  400)
normal_cells <- pick(subset(full, orig.ident == "Normal"), "Fib",  400)
t_cells      <- pick(full, "Tcell", 300)
macro_cells  <- pick(full, "Macro", 300)

cells  <- c(cf_cells, normal_cells, t_cells, macro_cells)
counts <- GetAssayData(full, layer = "counts")[, cells, drop = FALSE]

annotations <- data.frame(
  cell_type = c(rep("CF_Fib", length(cf_cells)),
                rep("Normal_Fib", length(normal_cells)),
                rep("Tcell", length(t_cells)),
                rep("Macro", length(macro_cells))),
  row.names = cells
)
write.table(annotations, "infercnv_annotations.txt",
            sep = "\t", quote = FALSE, col.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = "infercnv_annotations.txt",
  gene_order_file   = gene_order_file,
  ref_group_names   = c("Tcell", "Macro")
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = out_dir,
  denoise = TRUE,
  HMM = TRUE,
  HMM_type = "i6",
  analysis_mode = "subclusters",
  cluster_by_groups = TRUE,
  plot_steps = FALSE
)

saveRDS(infercnv_obj, file.path(out_dir, "infercnv_obj.rds"))

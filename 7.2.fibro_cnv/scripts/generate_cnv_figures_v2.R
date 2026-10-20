#!/usr/bin/env Rscript

options(bitmapType = "cairo")

library(RColorBrewer)
library(grDevices)

expr_obj_path <- "CNV_results/infercnv_output/14_invert_log_transformHMMi6.infercnv_obj"
hmm_obj_path <- "CNV_results/infercnv_output/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.infercnv_obj"
output_dir <- "figures"

plot_title <- "CNV Analysis: CF vs Normal Fibroblasts"
pdf_width <- 16
png_res <- 300

cat("========================================\n")
cat("InferCNV Heatmap Fix - Version 7\n")
cat("(Fixed color key Y-axis numbers)\n")
cat("========================================\n\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

get_chr_order <- function(chr_names) {
  chr_names <- as.character(chr_names)
  chr_ids <- gsub("^chr", "", chr_names)
  is_numeric <- grepl("^[0-9]+$", chr_ids)
  numeric_chrs <- chr_names[is_numeric]
  other_chrs <- chr_names[!is_numeric]

  if (length(numeric_chrs) > 0) {
    numeric_nums <- as.numeric(gsub("^chr", "", numeric_chrs))
    numeric_order <- numeric_chrs[order(numeric_nums)]
  } else {
    numeric_order <- character(0)
  }

  if (length(other_chrs) > 0) {
    other_ids <- gsub("^chr", "", other_chrs)
    special_order <- c("X", "Y", "M", "MT")
    other_ranks <- match(other_ids, special_order)
    other_ranks[is.na(other_ranks)] <- 100 + rank(other_ids[is.na(other_ranks)])
    other_order <- other_chrs[order(other_ranks)]
  } else {
    other_order <- character(0)
  }

  c(numeric_order, other_order)
}

calculate_legend_dims <- function(group_names, plot_width = 16, cex = 0.85) {
  if (length(group_names) == 0) return(list(ncol = 1, nrow = 0, height = 0.5))
  max_chars <- max(nchar(group_names))
  n_groups <- length(group_names)
  char_width <- 0.08 * cex
  item_width <- max_chars * char_width + 0.4
  avail_width <- plot_width - 3
  ncol <- max(1, min(6, floor(avail_width / item_width)))
  nrow <- ceiling(n_groups / ncol)
  row_height <- 0.35 * cex
  height <- nrow * row_height + 0.6
  list(ncol = ncol, nrow = nrow, height = height, cex = cex)
}

generate_cnv_colors <- function(n = 100) {
  colorRampPalette(c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
  ))(n)
}

merge_dendrograms <- function(dend_list) {
  if (length(dend_list) == 0) return(NULL)
  if (length(dend_list) == 1) return(dend_list[[1]])
  merged <- dend_list[[1]]
  for (i in 2:length(dend_list)) {
    merged <- merge(merged, dend_list[[i]], height = max(attr(merged, "height"),
                                                          attr(dend_list[[i]], "height")) + 0.5)
  }
  return(merged)
}

cat("Loading expression data object...\n")
expr_obj <- readRDS(expr_obj_path)
cat("Loading HMM object (for dendrograms)...\n")
hmm_obj <- readRDS(hmm_obj_path)

expr_data <- expr_obj@expr.data
cat(sprintf("Expression matrix: %d genes x %d cells\n", nrow(expr_data), ncol(expr_data)))

gene_order <- expr_obj@gene_order

cat("\nFixing chromosome order...\n")

current_chrs <- unique(as.character(gene_order$chr))
correct_chr_order <- get_chr_order(current_chrs)
cat(sprintf("Chromosome order: %s\n", paste(correct_chr_order, collapse = ", ")))

gene_order$chr <- factor(gene_order$chr, levels = correct_chr_order)
gene_order <- gene_order[order(gene_order$chr, gene_order$start), ]
expr_data <- expr_data[rownames(gene_order), ]

n_genes <- nrow(expr_data)
n_cells <- ncol(expr_data)
chromosomes <- correct_chr_order
chr_labels <- gsub("^chr", "", chromosomes)
cat(sprintf("Simplified labels: %s\n", paste(chr_labels, collapse = ", ")))

ref_groups <- expr_obj@reference_grouped_cell_indices
obs_groups <- expr_obj@observation_grouped_cell_indices
ref_names <- names(ref_groups)
obs_names <- names(obs_groups)
cat(sprintf("Reference groups: %d (%s)\n", length(ref_names), paste(ref_names, collapse = ", ")))
cat(sprintf("Observation groups: %d (%s)\n", length(obs_names), paste(obs_names, collapse = ", ")))

cat("\nExtracting dendrograms...\n")

if (!is.null(hmm_obj@tumor_subclusters) && !is.null(hmm_obj@tumor_subclusters$hc)) {
  obs_dendrograms <- list()
  obs_cell_order_list <- list()

  for (grp in obs_names) {
    if (!is.null(hmm_obj@tumor_subclusters$hc[[grp]])) {
      hc <- hmm_obj@tumor_subclusters$hc[[grp]]
      obs_dendrograms[[grp]] <- as.dendrogram(hc)
      obs_cell_order_list[[grp]] <- hc$labels[hc$order]
      cat(sprintf("  Found dendrogram for %s (%d cells)\n", grp, length(hc$order)))
    } else {
      grp_cells <- names(obs_groups[[grp]])
      if (is.null(grp_cells)) {
        grp_cells <- colnames(expr_data)[obs_groups[[grp]]]
      }
      obs_cell_order_list[[grp]] <- grp_cells
      cat(sprintf("  No dendrogram for %s, using original order\n", grp))
    }
  }

  obs_cell_order_final <- unlist(obs_cell_order_list)
  has_dendrogram <- length(obs_dendrograms) > 0

  if (has_dendrogram) {
    obs_merged_dend <- merge_dendrograms(obs_dendrograms)
    cat(sprintf("Merged %d observation dendrograms\n", length(obs_dendrograms)))
  }
} else {
  cat("No tumor_subclusters$hc found in HMM object\n")
  obs_cell_order_final <- unlist(lapply(obs_groups, function(x) colnames(expr_data)[x]))
  has_dendrogram <- FALSE
}

ref_cell_order <- unlist(lapply(ref_groups, function(x) colnames(expr_data)[x]))

cat(sprintf("\nTotal cells: ref=%d, obs=%d\n", length(ref_cell_order), length(obs_cell_order_final)))

cat("\nPreparing expression data...\n")

ref_expr <- expr_data[, ref_cell_order, drop = FALSE]
obs_expr <- expr_data[, obs_cell_order_final, drop = FALSE]

all_expr <- c(as.vector(ref_expr), as.vector(obs_expr))
x_center <- 1.0
quantiles <- quantile(all_expr, c(0.01, 0.99), na.rm = TRUE)
max_dev <- max(abs(quantiles - x_center), 0.15)
x_range <- c(x_center - max_dev, x_center + max_dev)
cat(sprintf("Expression range for coloring: [%.3f, %.3f]\n", x_range[1], x_range[2]))

cnv_colors <- generate_cnv_colors(100)

chr_factor <- factor(gene_order$chr, levels = chromosomes)
chr_boundaries <- which(diff(as.numeric(chr_factor)) != 0)
chr_boundaries <- c(0, chr_boundaries, n_genes)
chr_mids <- (chr_boundaries[-length(chr_boundaries)] + chr_boundaries[-1]) / 2 / n_genes

cat("\nCreating cell color bars (CF vs Normal only)...\n")

cf_normal_palette <- c("CF" = "#E41A1C", "Normal" = "#377EB8")

obs_group_type <- sapply(obs_names, function(x) {
  if (grepl("^CF", x)) return("CF")
  if (grepl("^Normal", x)) return("Normal")
  return("Unknown")
})
names(obs_group_type) <- obs_names

obs_cell_colors <- character(length(obs_cell_order_final))
obs_cell_types <- character(length(obs_cell_order_final))

cell_idx <- 1
for (grp in obs_names) {
  grp_cells <- obs_cell_order_list[[grp]]
  n_cells_in_grp <- length(grp_cells)
  cell_type <- obs_group_type[grp]
  obs_cell_colors[cell_idx:(cell_idx + n_cells_in_grp - 1)] <- cf_normal_palette[cell_type]
  obs_cell_types[cell_idx:(cell_idx + n_cells_in_grp - 1)] <- cell_type
  cell_idx <- cell_idx + n_cells_in_grp
}

ref_palette <- c("ref-1" = "#FFA500", "ref-2" = "#87CEEB")
ref_cell_colors <- character(length(ref_cell_order))
idx <- 1
for (grp in ref_names) {
  n <- length(ref_groups[[grp]])
  ref_cell_colors[idx:(idx + n - 1)] <- ref_palette[grp]
  idx <- idx + n
}

n_cf <- sum(obs_cell_types == "CF")
n_normal <- sum(obs_cell_types == "Normal")
cat(sprintf("CF cells: %d, Normal cells: %d\n", n_cf, n_normal))

cat("\nCalculating layout dimensions...\n")

obs_legend_names <- c("CF", "Normal")
obs_legend_dims <- calculate_legend_dims(obs_legend_names, pdf_width, cex = 1.0)
ref_legend_dims <- calculate_legend_dims(ref_names, pdf_width, cex = 1.0)

n_ref_cells <- length(ref_cell_order)
n_obs_cells <- length(obs_cell_order_final)
ref_heatmap_height <- max(2, 2.5 * n_ref_cells / 600)
obs_heatmap_height <- max(5, 6 * n_obs_cells / 800)

total_height <- 1.8 + ref_heatmap_height + 0.2 + obs_heatmap_height + 0.5 +
                ref_legend_dims$height + 0.2 + obs_legend_dims$height + 0.5

cat(sprintf("Total figure height: %.2f inches\n", total_height))

plot_infercnv_heatmap_v7 <- function() {

  layout_mat <- matrix(c(
    1,  1,  1,
    2,  3,  4,
    5,  6,  7,
    0,  8,  8,
    9,  9,  9,
    10, 10, 10
  ), nrow = 6, ncol = 3, byrow = TRUE)

  layout_heights <- c(
    1.5,
    ref_heatmap_height,
    obs_heatmap_height,
    1.0,
    ref_legend_dims$height,
    obs_legend_dims$height + 0.3
  )

  layout_widths <- c(1.5, 0.4, 12)

  layout(layout_mat, heights = layout_heights, widths = layout_widths)

  par(mar = c(2, 4, 1, 5))
  z <- seq(x_range[1], x_range[2], length.out = 100)
  image(z = matrix(z, ncol = 1), col = cnv_colors, axes = FALSE)
  axis(1, at = seq(0, 1, length.out = 5),
       labels = round(seq(x_range[1], x_range[2], length.out = 5), 2), cex.axis = 0.8)
  mtext("Modified Expression", side = 1, line = 2, cex = 0.8)
  title(main = "Distribution of Expression", cex.main = 0.9, font.main = 1)

  h <- hist(all_expr, breaks = 50, plot = FALSE)
  max_count <- max(h$counts)
  hx <- (h$mids - x_range[1]) / diff(x_range)
  hy <- h$counts / max_count * 0.8
  par(new = TRUE)
  plot(hx, hy, type = "s", col = "blue", lwd = 1,
       axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))

  pretty_counts <- pretty(c(0, max_count), n = 5)
  pretty_y_pos <- pretty_counts / max_count * 0.8
  pretty_labels <- formatC(pretty_counts, format = "e", digits = 1)

  axis(4, at = pretty_y_pos,
       labels = pretty_labels,
       cex.axis = 0.8, las = 1)
  mtext("Count", side = 4, line = 3.2, cex = 0.8)

  par(mar = c(0, 0, 0.5, 0))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))

  par(mar = c(0, 0.1, 0.5, 0.1))
  n_ref <- length(ref_cell_colors)
  image(z = matrix(1:n_ref, nrow = 1), col = ref_cell_colors,
        axes = FALSE, xlab = "", ylab = "")

  par(mar = c(0, 0, 0.5, 0.5))
  ref_mat <- as.matrix(ref_expr)
  ref_mat[ref_mat < x_range[1]] <- x_range[1]
  ref_mat[ref_mat > x_range[2]] <- x_range[2]
  image(z = ref_mat, col = cnv_colors, zlim = x_range,
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  for (b in chr_boundaries[-c(1, length(chr_boundaries))]) {
    abline(v = b / n_genes, col = "black", lwd = 0.5)
  }
  mtext("References (Cells)", side = 4, line = 0.2, cex = 0.6, las = 0)

  par(mar = c(0, 0, 0, 0))
  if (has_dendrogram && !is.null(obs_merged_dend)) {
    plot(obs_merged_dend, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  } else {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1))
  }

  par(mar = c(0, 0.1, 0, 0.1))
  n_obs <- length(obs_cell_colors)
  image(z = matrix(1:n_obs, nrow = 1), col = obs_cell_colors,
        axes = FALSE, xlab = "", ylab = "")

  par(mar = c(0, 0, 0, 0.5))
  obs_mat <- as.matrix(obs_expr)
  obs_mat[obs_mat < x_range[1]] <- x_range[1]
  obs_mat[obs_mat > x_range[2]] <- x_range[2]
  image(z = obs_mat, col = cnv_colors, zlim = x_range,
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  for (b in chr_boundaries[-c(1, length(chr_boundaries))]) {
    abline(v = b / n_genes, col = "black", lwd = 0.5)
  }
  mtext("Observations (Cells)", side = 4, line = 0.2, cex = 0.6, las = 0)

  par(mar = c(1.5, 0.1, 0.2, 0.5))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  for (i in seq_along(chr_labels)) {
    text(chr_mids[i], 0.5, labels = chr_labels[i],
         cex = 0.7, srt = 0, adj = c(0.5, 0.5))
  }
  mtext("Genomic Region", side = 1, line = 0.5, cex = 0.7)

  par(mar = c(0.3, 1, 0.3, 1))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend("center", legend = ref_names, fill = ref_palette[ref_names],
         ncol = ref_legend_dims$ncol, bty = "n", cex = ref_legend_dims$cex,
         title = "Reference Groups", title.adj = 0.5, xpd = TRUE)

  par(mar = c(0.5, 1, 0.3, 1))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend("center", legend = c("CF", "Normal"),
         fill = cf_normal_palette[c("CF", "Normal")],
         ncol = 2, bty = "n", cex = 1.2,
         title = "Observation Groups", title.adj = 0.5, xpd = TRUE)

  mtext(plot_title, side = 3, outer = TRUE, line = -1.5, cex = 1.2, font = 2)
}

cat("\n========================================\n")
cat("Generating outputs...\n")
cat("========================================\n\n")

cat("Generating PDF...\n")
pdf_file <- file.path(output_dir, "infercnv_heatmap.pdf")
cairo_pdf(pdf_file, width = pdf_width, height = total_height)
plot_infercnv_heatmap_v7()
dev.off()
cat(sprintf("PDF saved: %s (%.1f MB)\n", basename(pdf_file), file.info(pdf_file)$size / 1024 / 1024))

cat("\nGenerating PNG (300 dpi)...\n")
png_file <- file.path(output_dir, "infercnv_heatmap.png")
png(png_file, width = pdf_width, height = total_height,
    units = "in", res = png_res, type = "cairo")
plot_infercnv_heatmap_v7()
dev.off()
cat(sprintf("PNG saved: %s (%.1f MB)\n", basename(png_file), file.info(png_file)$size / 1024 / 1024))

cat("\n========================================\n")
cat("Output Summary\n")
cat("========================================\n\n")

output_files <- list.files(output_dir, full.names = FALSE)
for (f in output_files) {
  file_size <- file.info(file.path(output_dir, f))$size
  cat(sprintf("  - %s (%.1f MB)\n", f, file_size / 1024 / 1024))
}

cat("\nDone!\n")

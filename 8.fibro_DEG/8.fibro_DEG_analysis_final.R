rm(list = ls())
options(stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

log_file <- "logs/DEG_analysis_log.txt"
write_log <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- paste0("[", timestamp, "] ", message)
  cat(log_message, "\n")
  cat(log_message, "\n", file = log_file, append = TRUE)
}

write_log("========================================")
write_log("Fibroblast DEG Analysis: CF vs Normal")
write_log("========================================")

write_log("\nLoading required packages...")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

write_log("\nLoading fibroblast data...")
fibro_path <- "../7.0.fibro_extrction/07_fibroblast/fibroblast_annotated.rds"
if(!file.exists(fibro_path)) {
  stop("Fibroblast data not found at: ", fibro_path)
}

fibro_sce <- readRDS(fibro_path)
write_log(sprintf("OK Loaded %d fibroblast cells", ncol(fibro_sce)))

write_log("\nCell type distribution:")
celltype_table <- table(fibro_sce$celltype_abbr, fibro_sce$condition)
write_log(capture.output(print(celltype_table)))

subtypes <- c("ProInflamF", "SecRetF", "SecPapF")
write_log(sprintf("\nAnalyzing %d fibroblast subtypes: %s",
                  length(subtypes), paste(subtypes, collapse=", ")))

all_deg_results <- list()
deg_summary <- data.frame()

write_log("\n========================================")
write_log("Differential Expression Analysis")
write_log("========================================")

for(subtype in subtypes) {
  write_log(sprintf("\n--- Analyzing %s ---", subtype))

  cells_subset <- subset(fibro_sce, subset = celltype_abbr == subtype)
  n_cells <- ncol(cells_subset)

  n_cf <- sum(cells_subset$condition == "CF")
  n_normal <- sum(cells_subset$condition == "Normal")

  write_log(sprintf("  Total cells: %d (CF: %d, Normal: %d)", n_cells, n_cf, n_normal))

  if(n_cf < 3 || n_normal < 3) {
    write_log(sprintf("  WARNING  WARNING: Insufficient cells for %s (need >=3 per group)", subtype))
    write_log(sprintf("  Skipping %s", subtype))
    next
  }

  Idents(cells_subset) <- "condition"

  write_log(sprintf("  Running FindMarkers (CF vs Normal)..."))
  tryCatch({
    deg_markers <- FindMarkers(
      cells_subset,
      ident.1 = "CF",
      ident.2 = "Normal",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.25,
      verbose = FALSE
    )

    deg_markers$gene <- rownames(deg_markers)
    deg_markers$neg_log10_padj <- -log10(deg_markers$p_val_adj + 1e-300)
    deg_markers$significance <- ifelse(deg_markers$p_val_adj < 0.01,
                                       "Significant (p<0.01)",
                                       "Not Significant")
    deg_markers$subtype <- subtype

    all_deg_results[[subtype]] <- deg_markers

    n_total <- nrow(deg_markers)
    n_sig <- sum(deg_markers$p_val_adj < 0.01)
    n_up <- sum(deg_markers$p_val_adj < 0.01 & deg_markers$avg_log2FC > 0)
    n_down <- sum(deg_markers$p_val_adj < 0.01 & deg_markers$avg_log2FC < 0)

    write_log(sprintf("  OK Found %d DEGs (p_adj < 0.01): %d up, %d down", n_sig, n_up, n_down))

    output_file <- sprintf("data/DEG_%s_CF_vs_Normal.csv", subtype)
    deg_export <- deg_markers
    deg_export$p_val <- ifelse(deg_export$p_val == 0, "<1e-300",
                               formatC(deg_export$p_val, format = "e", digits = 3))
    deg_export$p_val_adj <- ifelse(deg_export$p_val_adj == 0, "<1e-300",
                                   formatC(deg_export$p_val_adj, format = "e", digits = 3))
    write.csv(deg_export, output_file, row.names = FALSE)
    write_log(sprintf("  OK Saved: %s", output_file))

    deg_summary <- rbind(deg_summary, data.frame(
      Subtype = subtype,
      N_cells_CF = n_cf,
      N_cells_Normal = n_normal,
      Total_genes_tested = n_total,
      Significant_DEGs = n_sig,
      Upregulated_in_CF = n_up,
      Downregulated_in_CF = n_down
    ))

  }, error = function(e) {
    write_log(sprintf("  ERROR ERROR in %s: %s", subtype, e$message))
  })
}

write.csv(deg_summary, "data/DEG_summary_all_subtypes.csv", row.names = FALSE)
write_log("\nOK Saved summary: data/DEG_summary_all_subtypes.csv")

write_log("\nDEG Summary:")
write_log(capture.output(print(deg_summary)))

write_log("\n========================================")
write_log("Creating Volcano Plots (Figure 4A style)")
write_log("========================================")

all_deg_combined <- data.frame()

for(subtype in names(all_deg_results)) {
  deg_data <- all_deg_results[[subtype]]

  deg_data <- deg_data %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))

  top10_genes <- head(deg_data, 10)$gene

  deg_data$label <- ifelse(deg_data$gene %in% top10_genes, deg_data$gene, "")

  all_deg_combined <- rbind(all_deg_combined, deg_data)

  write_log(sprintf("  %s: Top 10 genes - %s", subtype,
                    paste(head(top10_genes, 5), collapse=", ")))
}

all_deg_combined$subtype <- factor(all_deg_combined$subtype,
                                   levels = c("ProInflamF", "SecRetF", "SecPapF"))

all_deg_combined <- all_deg_combined %>%
  mutate(
    regulation = case_when(
      p_val_adj < 0.01 & avg_log2FC > 0 ~ "Upregulated (CF)",
      p_val_adj < 0.01 & avg_log2FC < 0 ~ "Downregulated (CF)",
      TRUE ~ "Not Significant"
    )
  )

all_deg_combined$regulation <- factor(all_deg_combined$regulation,
                                      levels = c("Upregulated (CF)",
                                                "Downregulated (CF)",
                                                "Not Significant"))

write_log("\nRegulation summary:")
write_log(capture.output(print(table(all_deg_combined$regulation))))

write_log("\nCreating faceted volcano plot with differential colors...")

p_volcano <- ggplot(all_deg_combined, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = regulation),
             alpha = 0.6,
             size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated (CF)" = "#E64B35",
      "Downregulated (CF)" = "#4DBBD5",
      "Not Significant" = "gray70"
    ),
    name = "Regulation"
  ) +
  geom_text_repel(
    aes(label = label),
    size = 2.5,
    max.overlaps = 12,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2,
    min.segment.length = 0
  ) +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed",
             color = "gray40",
             linewidth = 0.5) +
  geom_vline(xintercept = 0,
             linetype = "solid",
             color = "gray40",
             linewidth = 0.5) +
  facet_wrap(~ subtype, nrow = 1, scales = "free_x") +
  labs(
    title = "Fibroblast Subtypes: CF vs Normal",
    subtitle = "Differential Gene Expression Analysis (p_adj < 0.01)",
    x = "log2 Fold Change (CF vs Normal)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_classic(base_size = 11) +
  theme(

    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 11),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),

    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),

    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0.8, "lines"),

    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9)
  )

ggsave("figures/DEG_volcano_CF_vs_Normal_final.pdf", p_volcano, width = 14, height = 6)
write_log("OK Saved: figures/DEG_volcano_CF_vs_Normal_final.pdf")

write_log("\nCreating summary barplot...")

if(nrow(deg_summary) > 0) {
  p_summary <- ggplot(deg_summary, aes(x = Subtype)) +
    geom_bar(aes(y = Upregulated_in_CF, fill = "Upregulated"),
             stat = "identity", position = position_dodge(width = 0.8)) +
    geom_bar(aes(y = -Downregulated_in_CF, fill = "Downregulated"),
             stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5")) +
    geom_hline(yintercept = 0, color = "black") +
    labs(title = "Number of DEGs per Fibroblast Subtype",
         subtitle = "CF vs Normal (p_adj < 0.01)",
         x = "Fibroblast Subtype",
         y = "Number of DEGs",
         fill = "Regulation") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("figures/DEG_summary_barplot.pdf", p_summary, width = 8, height = 6)
  write_log("OK Saved: figures/DEG_summary_barplot.pdf")
}

write_log("\n========================================")
write_log("Creating Jitter Volcano Plot (scatter style)")
write_log("========================================")

all_deg_jitter <- data.frame()

for(subtype in names(all_deg_results)) {
  deg_data <- all_deg_results[[subtype]]
  deg_data$gene <- rownames(deg_data)
  deg_data$subtype <- subtype
  all_deg_jitter <- rbind(all_deg_jitter, deg_data)
}

all_deg_filtered <- all_deg_jitter
write_log(sprintf("  Total genes for jitter plot: %d genes", nrow(all_deg_filtered)))

all_deg_filtered$sig_label <- ifelse(all_deg_filtered$p_val_adj < 0.01,
                                     "adjust P-val<0.01",
                                     "adjust P-val>=0.01")

top_genes <- all_deg_filtered %>%
  group_by(subtype) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 25) %>%
  ungroup()

write_log(sprintf("  Selected Top 25 genes per subtype: %d genes total", nrow(top_genes)))
write_log(sprintf("  Unique genes among Top 25: %d", length(unique(top_genes$gene))))

gene_counts <- table(top_genes$gene)
duplicated_genes <- gene_counts[gene_counts > 1]
if(length(duplicated_genes) > 0) {
  write_log(sprintf("  %d genes appear in multiple subtypes:", length(duplicated_genes)))
  for(gene in names(duplicated_genes)) {
    subtypes <- top_genes$subtype[top_genes$gene == gene]
    write_log(sprintf("    - %s: %s", gene, paste(unique(subtypes), collapse = ", ")))
  }
}

all_deg_filtered$subtype <- factor(all_deg_filtered$subtype,
                                   levels = c("ProInflamF", "SecRetF", "SecPapF"))
top_genes$subtype <- factor(top_genes$subtype,
                            levels = c("ProInflamF", "SecRetF", "SecPapF"))

all_deg_filtered$is_top <- ifelse(all_deg_filtered$gene %in% top_genes$gene, 2, 1)
dt_non_top <- filter(all_deg_filtered, is_top == 1)

set.seed(123)
jitter_width <- 0.4
top_genes$x_jitter <- as.numeric(top_genes$subtype) + runif(nrow(top_genes), -jitter_width, jitter_width)
dt_non_top$x_jitter <- as.numeric(dt_non_top$subtype) + runif(nrow(dt_non_top), -jitter_width, jitter_width)

write_log(sprintf("  Added jitter coordinates (seed=123, width=%.1f)", jitter_width))

max_y <- ceiling(max(abs(all_deg_filtered$avg_log2FC)) + 1)
min_y <- -max_y

subtypes <- c("ProInflamF", "SecRetF", "SecPapF")

dfbar <- data.frame(x = factor(subtypes, levels = subtypes),
                    y = c(max_y, max_y, max_y))
dfbar1 <- data.frame(x = factor(subtypes, levels = subtypes),
                     y = c(min_y, min_y, min_y))

dfcol <- data.frame(x = factor(subtypes, levels = subtypes),
                    y = 0,
                    label = subtypes)
mycol <- c("#E64B35", "#4DBBD5", "#F39B7F")

write_log("\nCreating jitter scatter volcano plot...")

p_jitter <- ggplot() +

  geom_col(data = dfbar, aes(x = x, y = y),
           fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = dfbar1, aes(x = x, y = y),
           fill = "#dcdcdc", alpha = 0.6) +

  geom_point(data = dt_non_top,
             aes(x = x_jitter, y = avg_log2FC, color = sig_label),
             size = 0.85, alpha = 0.6) +

  geom_point(data = top_genes,
             aes(x = x_jitter, y = avg_log2FC, color = sig_label),
             size = 1.8, alpha = 0.9) +

  geom_tile(data = dfcol, aes(x = x, y = y),
            height = 0.8, color = "black",
            fill = mycol, alpha = 0.7,
            show.legend = FALSE) +

  geom_text_repel(data = top_genes,
                  aes(x = x_jitter, y = avg_log2FC, label = gene),
                  size = 2.3,
                  force = 3,
                  max.overlaps = 100,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.size = 0.25,
                  segment.color = "gray30",
                  min.segment.length = 0,
                  arrow = arrow(length = unit(0.01, "npc"),
                               type = "open", ends = "last")) +

  scale_color_manual(name = NULL,
                     values = c("adjust P-val<0.01" = "red",
                               "adjust P-val>=0.01" = "black")) +

  scale_x_discrete() +

  labs(title = "Fibroblast Subtypes: CF vs Normal",
       subtitle = "Jitter Volcano Plot - Top 25 genes per subtype (p_adj < 0.01)",
       x = "Fibroblast Subtype",
       y = "average log2FC (CF vs Normal)") +

  geom_text(data = dfcol, aes(x = x, y = y, label = label),
            size = 4, color = "white", fontface = "bold") +
  theme_minimal() +
  theme(

    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),

    axis.title = element_text(size = 13, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", linewidth = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),

    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 12)
  )

ggsave("figures/DEG_jitter_volcano_CF_vs_Normal.pdf", p_jitter, width = 14, height = 10)
write_log("OK Saved: figures/DEG_jitter_volcano_CF_vs_Normal.pdf")

write_log(sprintf("\nJitter plot statistics:"))
write_log(sprintf("  - Total genes plotted: %d", nrow(all_deg_filtered)))
write_log(sprintf("  - Significant (p_adj < 0.01): %d", sum(all_deg_filtered$sig_label == "adjust P-val<0.01")))
write_log(sprintf("  - Top 25 genes per subtype: %d total", nrow(top_genes)))
write_log(sprintf("  - Labels point to exact jittered positions (seed=123)"))

write_log("\n========================================")
write_log("Analysis Complete!")
write_log("========================================")

write_log("\nOutput files generated:")
write_log("  - data/DEG_summary_all_subtypes.csv")
write_log("  - figures/DEG_volcano_CF_vs_Normal_final.pdf (traditional faceted)")
write_log("  - figures/DEG_jitter_volcano_CF_vs_Normal.pdf (scatter style)")
write_log("  - figures/DEG_summary_barplot.pdf")
write_log(sprintf("  - %d individual DEG CSV files", length(all_deg_results)))

write_log("\nKey findings:")
for(i in 1:nrow(deg_summary)) {
  write_log(sprintf("  %s: %d significant DEGs (%d up, %d down)",
                    deg_summary$Subtype[i],
                    deg_summary$Significant_DEGs[i],
                    deg_summary$Upregulated_in_CF[i],
                    deg_summary$Downregulated_in_CF[i]))
}

write_log("\nColor scheme:")
write_log("   Red   - Upregulated in CF (p_adj < 0.01)")
write_log("   Blue  - Downregulated in CF (p_adj < 0.01)")
write_log("   Gray  - Not significant")

write_log("\n Differential expression analysis completed successfully!")
write_log("========================================\n")

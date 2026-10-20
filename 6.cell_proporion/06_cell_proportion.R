rm(list=ls())
options(stringsAsFactors = F)

current_dir <- basename(getwd())
initial_wd <- getwd()

if(current_dir == "6.cell_proporion") {
  data_path_04 <- "../4.cell.annotation/04_celltype/"
  data_path_03c <- "../3.clustering_analysis/03c_result/"
  scripts_path <- "../scRNA_scripts/"
  output_dir <- "./06_proportion"
} else {
  data_path_04 <- "4.cell.annotation/04_celltype/"
  data_path_03c <- "3.clustering_analysis/03c_result/"
  scripts_path <- "scRNA_scripts/"
  output_dir <- "6.cell_proporion/06_proportion"
}

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(celldex)
library(singleseqgset)
library(devtools)
library(reshape2)
library(tidyr)
library(ggthemes)
library(gplots)
library(patchwork)

cat("========================================\n")
cat("Step 06: Cell Proportion Analysis")
cat("========================================\n\n")

celltype_order <- c("Fib", "Endo", "PeriSMC", "Epi", "Macro", "Tcell",
                   "Schwann", "GlandEpi", "Neuro", "Mast", "MyoepiSMC")

cat(sprintf("Current directory: %s\n", getwd()))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

source(file.path(initial_wd, scripts_path, 'mycolors.R'))

data_file_paths <- c(
  file.path(initial_wd, data_path_04, "sce_annotated_merged.rds"),
  file.path(initial_wd, data_path_03c, "sce_merged_clusters.rds")
)

loaded <- FALSE
for(path in data_file_paths) {
  if(file.exists(path)) {
    sce.all = readRDS(path)
    cat(sprintf("OK Loaded data from: %s\n", path))
    loaded <- TRUE
    break
  }
}

if(!loaded) {
  stop(sprintf("ERROR Cannot find data files. Checked paths:\n%s\n",
               paste("  -", data_file_paths, collapse="\n")))
}

cat("\n Data summary:\n")
cat(sprintf("  Total cells: %d\n", ncol(sce.all)))
cat(sprintf("  Samples: %s\n", paste(unique(sce.all$orig.ident), collapse=", ")))

if(!"celltype" %in% colnames(sce.all@meta.data)) {
  if("merged_celltype" %in% colnames(sce.all@meta.data)) {
    sce.all$celltype <- sce.all$merged_celltype
    cat("  Using merged_celltype as celltype\n")
  } else {
    stop("ERROR No celltype information found in data\n")
  }
}

cat(sprintf("  Cell types: %d\n", length(unique(sce.all$celltype))))
print(table(sce.all$celltype))

cat("\n========================================\n")
cat("Preparing data for analysis...\n")
cat("========================================\n")

if(!"sample" %in% colnames(sce.all@meta.data)) {
  sce.all$sample <- sce.all$orig.ident
}

unique_samples <- unique(sce.all$sample)
cat(sprintf("  Detected samples: %s\n", paste(unique_samples, collapse=", ")))

if(length(unique_samples) != 2) {
  warning(sprintf("WARNING  Expected 2 samples, found %d: %s\n",
                  length(unique_samples),
                  paste(unique_samples, collapse=", ")))
}

if(all(c('Normal', 'CF') %in% unique_samples)) {
  sce.all$sample = factor(sce.all$sample, levels = c('Normal','CF'))
  cat("  OK Sample order set: Normal -> CF\n")
} else {
  sce.all$sample = factor(sce.all$sample, levels = sort(unique_samples))
  cat(sprintf("   Sample order set: %s\n", paste(levels(sce.all$sample), collapse=" -> ")))
}

available_celltypes <- intersect(celltype_order, unique(sce.all$celltype))
sce.all$celltype <- factor(sce.all$celltype, levels = available_celltypes)

celltype_colors <- setNames(my36colors[1:length(available_celltypes)], available_celltypes)

tb <- table(sce.all$sample, sce.all$celltype)
cat("\nOK Cell count table:\n")
print(tb)

bar_data <- as.data.frame(tb)
colnames(bar_data) <- c("Sample", "CellType", "Count")

bar_per <- bar_data %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count)) %>%
  mutate(Percent = Count / Total * 100)

cat("\nOK Cell proportions calculated\n")
print(bar_per)

write.csv(bar_per, file = "celltype_proportions.csv", row.names = FALSE)
cat("\n========================================\n")
cat("Statistical Analysis")
cat("========================================\n")

statistical_results <- data.frame()

for(ct in available_celltypes) {
  normal_count <- tb["Normal", ct]
  cf_count <- tb["CF", ct]
  normal_total <- sum(tb["Normal", ])
  cf_total <- sum(tb["CF", ])

  contingency_table <- matrix(c(normal_count, cf_count,
                                normal_total - normal_count,
                                cf_total - cf_count),
                              nrow = 2, byrow = TRUE)
  fisher_result <- fisher.test(contingency_table)

  normal_prop <- normal_count / normal_total
  cf_prop <- cf_count / cf_total

  statistical_results <- rbind(statistical_results,
    data.frame(
      CellType = ct,
      Normal_Count = normal_count,
      Normal_Percent = round(normal_prop * 100, 2),
      CF_Count = cf_count,
      CF_Percent = round(cf_prop * 100, 2),
      Fold_Change = round(cf_prop / normal_prop, 2),
      P_Value = fisher_result$p.value,
      Significance = ifelse(fisher_result$p.value < 0.001, "***",
                           ifelse(fisher_result$p.value < 0.01, "**",
                                 ifelse(fisher_result$p.value < 0.05, "*", "ns")))
    )
  )
}

cat("\nOK Statistical test results:\n")
print(statistical_results)

write.csv(statistical_results, file = "celltype_statistics.csv", row.names = FALSE)

cat("\n========================================\n")
cat("Visualization")
cat("========================================\n")

p1 <- ggplot(bar_per, aes(x = Percent, y = Sample)) +
  geom_bar(aes(fill = CellType), stat = "identity") +
  coord_flip() +
  labs(y = "Sample", x = "Relative proportion (%)", fill = "Cell Type") +
  scale_fill_manual(values = celltype_colors) +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  ggtitle("Cell Type Proportions (Stacked)")

p2 <- ggplot(bar_per, aes(x = CellType, y = Percent, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cell Type", y = "Proportion (%)", fill = "Sample") +
  scale_fill_manual(values = c("Normal" = "#4DBBD4", "CF" = "#E64A35")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  ggtitle("Cell Type Proportions Comparison")

p1 + p2
ggsave(filename = "proportion_comparison.pdf", width = 14, height = 6)
cat("OK Saved: proportion_comparison.pdf\n")

heatmap_data <- reshape2::dcast(bar_per, CellType ~ Sample, value.var = "Percent")
rownames(heatmap_data) <- heatmap_data$CellType
heatmap_data <- heatmap_data[, -1]

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(0, max(heatmap_data, na.rm = TRUE)/2, max(heatmap_data, na.rm = TRUE)),
                      c("white", "#FDB863", "#B2182B"))

ht <- Heatmap(as.matrix(heatmap_data),
              name = "Proportion (%)",
              col = col_fun,
              row_names_side = "left",
              column_names_side = "bottom",
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.1f", heatmap_data[i, j]),
                         x, y, gp = gpar(fontsize = 8))
              },
              column_title = "Cell Type Proportions Heatmap",
              heatmap_legend_param = list(title = "Proportion (%)")
)

pdf("proportion_heatmap.pdf", width = 6, height = 8)
draw(ht)
dev.off()
cat("OK Saved: proportion_heatmap.pdf\n")

p3 <- ggplot(statistical_results, aes(x = CellType)) +
  geom_point(aes(y = Normal_Percent), color = "#4DBBD4", size = 3) +
  geom_point(aes(y = CF_Percent), color = "#E64A35", size = 3) +
  geom_segment(aes(xend = CellType, y = Normal_Percent, yend = CF_Percent),
               color = "gray50", linetype = "dashed") +
  geom_text(aes(y = pmax(Normal_Percent, CF_Percent) + 2,
                label = Significance),
            size = 4) +
  labs(x = "Cell Type", y = "Proportion (%)",
       title = "Cell Type Proportion Changes",
       subtitle = "Blue: Normal, Red: CF") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

ggsave("proportion_dotplot.pdf", plot = p3, width = 10, height = 6)
cat("OK Saved: proportion_dotplot.pdf\n")

library(ggrepel)

pie_data_list <- split(bar_per, bar_per$Sample)

pie_plots <- list()
for(sample_name in names(pie_data_list)) {
  pie_data <- pie_data_list[[sample_name]]

  pie_plots[[sample_name]] <- ggplot(pie_data, aes(x = "", y = Percent, fill = CellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = celltype_colors) +
    theme_void() +
    theme(legend.position = "right") +
    ggtitle(paste("Cell Composition -", sample_name)) +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
                    position = position_stack(vjust = 0.5),
                    size = 3)
}

pie_combined <- pie_plots[[1]] + pie_plots[[2]]
ggsave("proportion_pie.pdf", plot = pie_combined, width = 14, height = 6)
cat("OK Saved: proportion_pie.pdf\n")

p4 <- ggplot(bar_per, aes(x = CellType, y = Percent, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Normal" = "#4DBBD4", "CF" = "#E64A35")) +
  labs(x = "Cell Type", y = "Proportion (%)",
       title = "Cell Type Distribution") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  facet_wrap(~ Sample, ncol = 1, scales = "free_x")

ggsave("proportion_facet.pdf", plot = p4, width = 10, height = 8)
cat("OK Saved: proportion_facet.pdf\n")

cat("\n========================================\n")
cat("Summary Report")
cat("========================================\n")

summary_text <- paste0(
  "Cell Proportion Analysis Summary\n",
  "================================\n\n",
  "Date: ", Sys.Date(), "\n",
  "Total cells analyzed: ", ncol(sce.all), "\n",
  "  - Normal: ", sum(tb["Normal", ]), " cells\n",
  "  - CF: ", sum(tb["CF", ]), " cells\n\n",
  "Cell types identified: ", length(available_celltypes), "\n\n",
  "Major findings:\n"
)

sig_changes <- statistical_results[statistical_results$Significance != "ns", ]
if(nrow(sig_changes) > 0) {
  summary_text <- paste0(summary_text,
    "\nSignificantly changed cell types (p < 0.05):\n")
  for(i in 1:nrow(sig_changes)) {
    summary_text <- paste0(summary_text,
      sprintf("  - %s: %.1f%% (Normal) vs %.1f%% (CF), FC=%.2f, p=%.3e %s\n",
              sig_changes$CellType[i],
              sig_changes$Normal_Percent[i],
              sig_changes$CF_Percent[i],
              sig_changes$Fold_Change[i],
              sig_changes$P_Value[i],
              sig_changes$Significance[i]))
  }
} else {
  summary_text <- paste0(summary_text,
    "\nNo significantly changed cell types found.\n")
}

top_normal <- statistical_results[order(statistical_results$Normal_Percent, decreasing = TRUE), ][1:3, ]
top_cf <- statistical_results[order(statistical_results$CF_Percent, decreasing = TRUE), ][1:3, ]

summary_text <- paste0(summary_text,
  "\nTop 3 cell types in Normal:\n")
for(i in 1:3) {
  summary_text <- paste0(summary_text,
    sprintf("  %d. %s (%.1f%%)\n", i, top_normal$CellType[i], top_normal$Normal_Percent[i]))
}

summary_text <- paste0(summary_text,
  "\nTop 3 cell types in CF:\n")
for(i in 1:3) {
  summary_text <- paste0(summary_text,
    sprintf("  %d. %s (%.1f%%)\n", i, top_cf$CellType[i], top_cf$CF_Percent[i]))
}

cat(summary_text)
writeLines(summary_text, "proportion_summary.txt")
cat("\nOK Summary report saved to proportion_summary.txt\n")

cat("\n========================================\n")
cat(" Cell proportion analysis completed successfully!\n")
cat("========================================\n\n")

cat("Output files generated:\n")
cat("  - celltype_proportions.csv: Detailed proportion data\n")
cat("  - celltype_statistics.csv: Statistical test results\n")
cat("  - proportion_comparison.pdf: Stacked and grouped bar plots\n")
cat("  - proportion_heatmap.pdf: Proportion heatmap\n")
cat("  - proportion_dotplot.pdf: Dot plot with significance\n")
cat("  - proportion_pie.pdf: Pie charts for each sample\n")
cat("  - proportion_facet.pdf: Faceted bar plots\n")
cat("  - proportion_summary.txt: Text summary report\n")

setwd('../')

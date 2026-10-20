#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)

set.seed(123)
})

cat("============================================================================\n")
cat("FIBROBLAST SUBTYPE COMMUNICATION ANALYSIS\n")
cat("============================================================================\n\n")

source("utils/resource_monitor.R")
source("utils/parallel_config.R")

dir.create("08.fibro_subtype_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("08.fibro_subtype_analysis/CF", showWarnings = FALSE, recursive = TRUE)
dir.create("08.fibro_subtype_analysis/Normal", showWarnings = FALSE, recursive = TRUE)
dir.create("08.fibro_subtype_analysis/comparison", showWarnings = FALSE, recursive = TRUE)
dir.create("08.fibro_subtype_analysis/inter_subtype", showWarnings = FALSE, recursive = TRUE)

setup_parallel_computing(n_cores = 80)

perf_log <- "logs/fibro_subtype_performance.log"
cat("", file = perf_log)

log_step <- function(step, start_time = NULL) {
  if (is.null(start_time)) {
    return(Sys.time())
  } else {
    elapsed <- difftime(Sys.time(), start_time, units = "secs")
    ram <- get_ram_usage()
    msg <- sprintf("%s completed in %.2f seconds (RAM: %.2f GB)\n",
                   step, elapsed, ram)
    cat(msg)
    cat(msg, file = perf_log, append = TRUE)
  }
}

cat("\nSTEP 1: LOADING FIBROBLAST DATA\n")
cat("----------------------------------------\n")
t1 <- log_step()

cat("Loading data from 01.input_data/input.Rdata...\n")
load("01.input_data/input.Rdata")
cat(sprintf("  OK Loaded full dataset: %d cells\n", ncol(sce.all)))

fibro <- subset(sce.all, subset = celltype == "Fib")
cat(sprintf("  OK Extracted %d fibroblast cells\n", ncol(fibro)))

fibro_subtypes <- fibro.annotated
cat(sprintf("  OK Loaded subtype annotations for %d cells\n", ncol(fibro_subtypes)))

subtype_col <- "celltype"
if (!subtype_col %in% colnames(fibro_subtypes@meta.data)) {
  subtype_col <- colnames(fibro_subtypes@meta.data)[grepl("celltype", colnames(fibro_subtypes@meta.data))][1]
}

common_cells <- intersect(colnames(fibro), colnames(fibro_subtypes))
cat(sprintf("\nCommon cells: %d\n", length(common_cells)))

fibro <- fibro[, common_cells]
fibro_subtypes <- fibro_subtypes[, common_cells]

subtype_map <- fibro_subtypes@meta.data[[subtype_col]]
names(subtype_map) <- colnames(fibro_subtypes)
fibro$subtype <- subtype_map

if ("group" %in% colnames(fibro@meta.data)) {

  cat("Using existing 'group' column\n")
} else if ("condition" %in% colnames(fibro@meta.data)) {
  fibro$group <- fibro$condition
  cat("Using 'condition' as group\n")
} else if ("orig.ident" %in% colnames(fibro@meta.data)) {
  fibro$group <- fibro$orig.ident
  cat("Using 'orig.ident' as group\n")
} else {
  stop("No group/condition/orig.ident column found in metadata!")
}

cat("\nSubtype distribution:\n")
print(table(fibro$subtype))

cat("\nSubtype by group:\n")
print(table(fibro$subtype, fibro$group))

subtype_names <- c(
  "mesenchymal" = "MesF",
  "pro-inflammatory" = "ProInflamF",
  "secretory-reticular" = "SecRetF",
  "secretory-papillary" = "SecPapF"
)

fibro@meta.data$subtype_short <- subtype_names[as.character(fibro$subtype)]

log_step("Data Loading", t1)

cat("\n\nSTEP 2: SUBTYPE-SPECIFIC ANALYSIS\n")
cat("========================================\n")

analyze_subtype <- function(subtype_name, subtype_short, fibro_data) {
  cat(sprintf("\n--- Analyzing %s (%s) ---\n", subtype_name, subtype_short))
  t_sub <- Sys.time()

  sub_cells <- subset(fibro_data, subset = subtype == subtype_name)

  cf_n <- sum(sub_cells$group == "CF")
  normal_n <- sum(sub_cells$group == "Normal")
  cat(sprintf("  Sample sizes: CF=%d, Normal=%d\n", cf_n, normal_n))

  if (cf_n < 50 || normal_n < 50) {
    cat(sprintf("  WARNING WARNING: Small sample size. Results should be interpreted cautiously.\n"))
  }

  if (cf_n < 10 || normal_n < 10) {
    cat(sprintf("  ERROR SKIPPED: Sample size too small for reliable CellChat analysis.\n"))
    return(NULL)
  }

  sub_cf <- subset(sub_cells, subset = group == "CF")
  sub_normal <- subset(sub_cells, subset = group == "Normal")

  prepare_cellchat_data <- function(seurat_obj, group_name) {

    if (length(seurat_obj@assays$RNA@layers) > 1) {
      seurat_obj <- JoinLayers(seurat_obj)
    }

    data.input <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
    meta <- seurat_obj@meta.data

    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subtype")

    CellChatDB <- CellChatDB.human
    cellchat@DB <- CellChatDB

    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)

    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    return(cellchat)
  }

  cat(sprintf("  Running CellChat for %s-CF...\n", subtype_short))
  cellchat_cf <- tryCatch({
    prepare_cellchat_data(sub_cf, "CF")
  }, error = function(e) {
    cat(sprintf("  Error in CF analysis: %s\n", e$message))
    return(NULL)
  })

  cat(sprintf("  Running CellChat for %s-Normal...\n", subtype_short))
  cellchat_normal <- tryCatch({
    prepare_cellchat_data(sub_normal, "Normal")
  }, error = function(e) {
    cat(sprintf("  Error in Normal analysis: %s\n", e$message))
    return(NULL)
  })

  if (is.null(cellchat_cf) || is.null(cellchat_normal)) {
    cat("  ERROR Analysis failed for this subtype\n")
    return(NULL)
  }

  saveRDS(cellchat_cf, sprintf("08.fibro_subtype_analysis/CF/%s_cellchat.rds", subtype_short))
  saveRDS(cellchat_normal, sprintf("08.fibro_subtype_analysis/Normal/%s_cellchat.rds", subtype_short))

  cat("  Comparing CF vs Normal...\n")

  cellchat_list <- list(CF = cellchat_cf, Normal = cellchat_normal)
  cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

  saveRDS(cellchat_merged, sprintf("08.fibro_subtype_analysis/comparison/%s_merged.rds", subtype_short))

  lr_cf <- subsetCommunication(cellchat_cf)
  lr_normal <- subsetCommunication(cellchat_normal)

  if (nrow(lr_cf) > 0) {
    lr_cf$group <- "CF"
    write.csv(lr_cf, sprintf("08.fibro_subtype_analysis/CF/%s_LR_pairs.csv", subtype_short), row.names = FALSE)
  }

  if (nrow(lr_normal) > 0) {
    lr_normal$group <- "Normal"
    write.csv(lr_normal, sprintf("08.fibro_subtype_analysis/Normal/%s_LR_pairs.csv", subtype_short), row.names = FALSE)
  }

  pdf(sprintf("08.fibro_subtype_analysis/comparison/%s_comparison.pdf", subtype_short), width = 12, height = 8)

  tryCatch({
    gg <- compareInteractions(cellchat_merged, show.legend = TRUE, group = c(1,2))
    print(gg)
  }, error = function(e) cat("  Warning: Could not generate interaction comparison\n"))

  tryCatch({
    gg <- compareInteractions(cellchat_merged, show.legend = TRUE, group = c(1,2), measure = "weight")
    print(gg)
  }, error = function(e) cat("  Warning: Could not generate weight comparison\n"))

  tryCatch({
    netVisual_heatmap(cellchat_merged, comparison = c(1, 2))
  }, error = function(e) cat("  Warning: Could not generate heatmap\n"))

  dev.off()

  elapsed <- difftime(Sys.time(), t_sub, units = "secs")
  cat(sprintf("  OK %s analysis completed in %.2f seconds\n", subtype_short, elapsed))

  return(list(
    cf = cellchat_cf,
    normal = cellchat_normal,
    merged = cellchat_merged,
    lr_cf = lr_cf,
    lr_normal = lr_normal
  ))
}

subtype_results <- list()

for (st in names(subtype_names)) {
  st_short <- subtype_names[st]
  result <- analyze_subtype(st, st_short, fibro)

  if (!is.null(result)) {
    subtype_results[[st_short]] <- result
  }
}

cat("\n\nSTEP 3: INTER-SUBTYPE COMMUNICATION\n")
cat("========================================\n")
cat("Analyzing communication between fibroblast subtypes within CF...\n\n")
t3 <- log_step()

fibro_cf <- subset(fibro, subset = group == "CF")

if (length(fibro_cf@assays$RNA@layers) > 1) {
  fibro_cf <- JoinLayers(fibro_cf)
}

data_cf <- GetAssayData(fibro_cf, slot = "data", assay = "RNA")
meta_cf <- fibro_cf@meta.data

cellchat_cf_subtypes <- createCellChat(object = data_cf, meta = meta_cf, group.by = "subtype_short")

CellChatDB <- CellChatDB.human
cellchat_cf_subtypes@DB <- CellChatDB

cat("  Preprocessing...\n")
cellchat_cf_subtypes <- subsetData(cellchat_cf_subtypes)
cellchat_cf_subtypes <- identifyOverExpressedGenes(cellchat_cf_subtypes)
cellchat_cf_subtypes <- identifyOverExpressedInteractions(cellchat_cf_subtypes)

cat("  Computing communication probabilities...\n")
cellchat_cf_subtypes <- computeCommunProb(cellchat_cf_subtypes, type = "triMean")
cellchat_cf_subtypes <- filterCommunication(cellchat_cf_subtypes, min.cells = 10)

cat("  Computing pathway activities...\n")
cellchat_cf_subtypes <- computeCommunProbPathway(cellchat_cf_subtypes)
cellchat_cf_subtypes <- aggregateNet(cellchat_cf_subtypes)

cat("  Computing network centrality...\n")
cellchat_cf_subtypes <- netAnalysis_computeCentrality(cellchat_cf_subtypes)

saveRDS(cellchat_cf_subtypes, "08.fibro_subtype_analysis/inter_subtype/CF_all_subtypes_cellchat.rds")

lr_inter <- subsetCommunication(cellchat_cf_subtypes)
write.csv(lr_inter, "08.fibro_subtype_analysis/inter_subtype/inter_subtype_LR_pairs.csv", row.names = FALSE)

cat(sprintf("  OK Found %d inter-subtype L-R interactions\n", nrow(lr_inter)))

cat("\nKey inter-subtype interactions:\n")
if (nrow(lr_inter) > 0) {
  key_pairs <- lr_inter %>%
    group_by(source, target) %>%
    summarise(
      n_interactions = n(),
      avg_prob = mean(prob),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_prob)) %>%
    head(10)

  print(key_pairs)
  write.csv(key_pairs, "08.fibro_subtype_analysis/inter_subtype/key_subtype_pairs.csv", row.names = FALSE)
}

log_step("Inter-subtype Analysis", t3)

cat("\n\nSTEP 4: GENERATING INTER-SUBTYPE VISUALIZATIONS\n")
cat("========================================\n")
t4 <- log_step()

pdf("08.fibro_subtype_analysis/inter_subtype/network_visualizations.pdf", width = 12, height = 10)

cat("  Generating circle plot...\n")
tryCatch({
  groupSize <- as.numeric(table(cellchat_cf_subtypes@idents))
  netVisual_circle(cellchat_cf_subtypes@net$count,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = "Number of interactions between subtypes")
}, error = function(e) cat("    Warning: Circle plot failed\n"))

tryCatch({
  netVisual_circle(cellchat_cf_subtypes@net$weight,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = "Interaction strength between subtypes")
}, error = function(e) cat("    Warning: Strength circle plot failed\n"))

cat("  Generating heatmap...\n")
tryCatch({
  netVisual_heatmap(cellchat_cf_subtypes, measure = "count")
}, error = function(e) cat("    Warning: Count heatmap failed\n"))

tryCatch({
  netVisual_heatmap(cellchat_cf_subtypes, measure = "weight")
}, error = function(e) cat("    Warning: Weight heatmap failed\n"))

cat("  Analyzing signaling roles...\n")
tryCatch({
  cellchat_cf_subtypes <- netAnalysis_signalingRole_network(cellchat_cf_subtypes,
                                                             signaling = NULL,
                                                             width = 12, height = 6)
}, error = function(e) cat("    Warning: Signaling role analysis failed\n"))

tryCatch({
  gg <- netAnalysis_signalingRole_scatter(cellchat_cf_subtypes)
  print(gg)
}, error = function(e) cat("    Warning: Signaling role scatter failed\n"))

dev.off()

log_step("Visualization", t4)

cat("\n\nSTEP 5: KEY PATHWAY ANALYSIS\n")
cat("========================================\n")
t5 <- log_step()

all_pathways <- cellchat_cf_subtypes@netP$pathways

pathway_strength <- sapply(all_pathways, function(pw) {
  sum(cellchat_cf_subtypes@net$weight[,,pw], na.rm = TRUE)
})

pathway_strength_sorted <- sort(pathway_strength, decreasing = TRUE)
top_n <- min(20, length(pathway_strength_sorted))
top_pathways <- names(pathway_strength_sorted)[1:top_n]

cat(sprintf("\nDetected %d active pathways\n", length(all_pathways)))
cat(sprintf("Analyzing top %d pathways by communication strength:\n\n", top_n))
print(head(pathway_strength_sorted, 20))

pathway_stats <- data.frame(
  Pathway = names(pathway_strength_sorted),
  Strength = as.numeric(pathway_strength_sorted),
  Rank = 1:length(pathway_strength_sorted)
)
write.csv(pathway_stats,
          "08.fibro_subtype_analysis/inter_subtype/pathway_strength_stats.csv",
          row.names = FALSE)
cat("\nOK Saved pathway statistics to pathway_strength_stats.csv\n\n")

pdf("08.fibro_subtype_analysis/inter_subtype/pathway_analysis.pdf", width = 14, height = 10)

visualized_count <- 0
for (pathway in top_pathways) {
  cat(sprintf("  Visualizing %s (strength: %.3f)...\n",
              pathway, pathway_strength[pathway]))

  tryCatch({

    pathway_sum <- sum(cellchat_cf_subtypes@net$weight[,,pathway], na.rm = TRUE)

    if (pathway_sum > 0.01) {

      netVisual_aggregate(cellchat_cf_subtypes, signaling = pathway,
                        layout = "hierarchy")

      netVisual_aggregate(cellchat_cf_subtypes, signaling = pathway,
                        layout = "circle")

      netVisual_heatmap(cellchat_cf_subtypes, signaling = pathway)

      visualized_count <- visualized_count + 1
      cat(sprintf("    OK Successfully visualized %s\n", pathway))
    } else {
      cat(sprintf("    WARNING Skipped %s (strength too low: %.4f)\n", pathway, pathway_sum))
    }
  }, error = function(e) {
    cat(sprintf("    ERROR Error visualizing %s: %s\n", pathway, e$message))
  })
}

dev.off()

cat(sprintf("\nOK Successfully visualized %d/%d pathways\n", visualized_count, top_n))

log_step("Pathway Analysis", t5)

cat("\n\nSTEP 6: GENERATING SUMMARY STATISTICS\n")
cat("========================================\n")

summary_list <- list()

for (st_short in names(subtype_results)) {
  res <- subtype_results[[st_short]]

  summary_list[[st_short]] <- data.frame(
    Subtype = st_short,
    CF_Cells = ncol(res$cf@data.signaling),
    Normal_Cells = ncol(res$normal@data.signaling),
    CF_Interactions = nrow(res$lr_cf),
    Normal_Interactions = nrow(res$lr_normal),
    CF_Pathways = length(res$cf@netP$pathways),
    Normal_Pathways = length(res$normal@netP$pathways),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, "08.fibro_subtype_analysis/subtype_summary_statistics.csv", row.names = FALSE)

cat("\nSubtype Summary:\n")
print(summary_df)

inter_summary <- data.frame(
  Analysis = "Inter-subtype (CF)",
  Total_Cells = ncol(cellchat_cf_subtypes@data.signaling),
  Subtypes = length(unique(cellchat_cf_subtypes@idents)),
  Total_Interactions = nrow(lr_inter),
  Active_Pathways = length(cellchat_cf_subtypes@netP$pathways)
)

write.csv(inter_summary, "08.fibro_subtype_analysis/inter_subtype/summary_statistics.csv", row.names = FALSE)

cat("\nInter-subtype Summary:\n")
print(inter_summary)

cat("\n\n============================================================================\n")
cat("FIBROBLAST SUBTYPE ANALYSIS COMPLETE\n")
cat("============================================================================\n\n")

cat("Analyzed subtypes:\n")
for (st in names(summary_list)) {
  cat(sprintf("  OK %s\n", st))
}

cat("\nOutput directories:\n")
cat("  08.fibro_subtype_analysis/CF/ - Individual CF subtype results\n")
cat("  08.fibro_subtype_analysis/Normal/ - Individual Normal subtype results\n")
cat("  08.fibro_subtype_analysis/comparison/ - CF vs Normal comparisons\n")
cat("  08.fibro_subtype_analysis/inter_subtype/ - Inter-subtype communication\n")

cat("\nKey files:\n")
cat("  - subtype_summary_statistics.csv\n")
cat("  - inter_subtype/inter_subtype_LR_pairs.csv\n")
cat("  - inter_subtype/network_visualizations.pdf\n")
cat("  - inter_subtype/pathway_analysis.pdf\n")

cat("\nNEXT STEPS:\n")
cat("  1. Review subtype-specific L-R pairs in 08.fibro_subtype_analysis/CF/ and Normal/\n")
cat("  2. Examine inter-subtype communication network\n")
cat("  3. Integrate with DEG results (Script C)\n")
cat("  4. Integrate with SCENIC TF analysis (Script D)\n")

cat("\n============================================================================\n")
cat("Analysis completed successfully!\n")
cat("============================================================================\n")

N_CORES <- 80
MAX_RAM_GB <- 450
MEMORY_CHECK_INTERVAL <- 30

ENABLE_SPARSE <- TRUE
BATCH_SIZE <- 5000
AUTO_DOWNSAMPLE <- FALSE
MAX_CELLS_PER_TYPE <- 5000

CELLCHAT_DB <- "Secreted Signaling"
MIN_CELLS <- 10
NBOOT <- 100

SAVE_INCREMENTAL <- TRUE
RUN_FIBROBLAST_SUBTYPES <- TRUE

set.seed(123)

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("CELLCHAT ANALYSIS - CASE VS CONTROL\n")
cat(strrep("=", 80) %+% "\n\n")

if (basename(getwd()) != "13.cellchat") {
  if (dir.exists("13.cellchat")) {
    setwd("13.cellchat")
  } else {
    dir.create("13.cellchat", recursive = TRUE)
    setwd("13.cellchat")
  }
}

cat(sprintf("Working directory: %s\n", getwd()))

output_dirs <- c(
  "logs",
  "01.input_data",
  "02.CF_vs_Normal",
  "03.fibroblast_subtypes",
  "04.network_visualization",
  "05.pathway_analysis",
  "06.differential_analysis",
  "07.signaling_roles",
  "utils"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

cat("Output directories created.\n\n")

cat("Loading required libraries...\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(parallel)
  library(doParallel)
  library(future)
})

options(stringsAsFactors = FALSE)

cat("Libraries loaded successfully.\n")

cat("\nLoading utility modules...\n")

if (file.exists("utils/resource_monitor.R")) {
  source("utils/resource_monitor.R")
  cat("  OK Resource monitor loaded\n")
} else {
  stop("Resource monitor not found! Please ensure utils/resource_monitor.R exists.")
}

if (file.exists("utils/parallel_config.R")) {
  source("utils/parallel_config.R")
  cat("  OK Parallel config loaded\n")
} else {
  stop("Parallel config not found! Please ensure utils/parallel_config.R exists.")
}

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("SYSTEM RESOURCE SETUP\n")
cat(strrep("=", 80) %+% "\n")

sys_info <- get_system_ram_info()
cat(sprintf("System RAM: %.2f GB total, %.2f GB available\n",
           sys_info$total_gb, sys_info$available_gb))

max_cores <- detectCores()
if (N_CORES > max_cores) {
  warning(sprintf("Requested %d cores exceeds available %d. Adjusting to %d",
                 N_CORES, max_cores, max_cores - 2))
  N_CORES <- max_cores - 2
}

setup_parallel_computing(n_cores = N_CORES, backend = "multicore", verbose = TRUE)

cat("\n")
cat("Memory Monitoring:\n")
cat(sprintf("  RAM limit: %.1f GB (auto-kill if exceeded)\n", MAX_RAM_GB))
cat("  Current RAM: %.2f GB\n", get_ram_usage())
cat("  Check logs/memory_monitor.log for real-time RAM usage\n")
cat("\n")

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("STEP 1: LOADING AND PREPARING DATA\n")
cat(strrep("=", 80) %+% "\n\n")

start_time_load <- Sys.time()

cat("Loading annotated Seurat object...\n")
data_file <- "../4.cell.annotation/04_celltype/sce_annotated_merged.rds"

if (!file.exists(data_file)) {
  stop(sprintf("Data file not found: %s\n", data_file))
}

sce.all <- readRDS(data_file)
cat(sprintf("  OK Loaded: %d cells x %d genes\n", ncol(sce.all), nrow(sce.all)))

cat("\nCell type distribution:\n")
print(table(sce.all$celltype))

cat("\nGroup distribution:\n")
if ("group" %in% colnames(sce.all@meta.data)) {
  print(table(sce.all$group))
} else if ("orig.ident" %in% colnames(sce.all@meta.data)) {
  sce.all$group <- sce.all$orig.ident
  print(table(sce.all$group))
} else {
  stop("No group or orig.ident column found in metadata!")
}

cat("\nExtracting fibroblast cells...\n")
sce.fibro <- sce.all[, sce.all$celltype == "Fib"]
cat(sprintf("  OK %d fibroblast cells extracted\n", ncol(sce.fibro)))

fibro_anno_file <- "../7.0.fibro_extrction/07_fibroblast/fibroblast_annotated.rds"

if (file.exists(fibro_anno_file)) {
  cat("\nLoading fibroblast subtype annotations...\n")
  fibro.annotated <- readRDS(fibro_anno_file)
  cat(sprintf("  OK %d fibroblast cells with subtype annotations\n", ncol(fibro.annotated)))

  cat("\nFibroblast subtype distribution:\n")

  subtype_col <- if ("celltype" %in% colnames(fibro.annotated@meta.data)) {
    "celltype"
  } else if ("detailed_celltype" %in% colnames(fibro.annotated@meta.data)) {
    "detailed_celltype"
  } else if ("predicted_celltype" %in% colnames(fibro.annotated@meta.data)) {
    "predicted_celltype"
  } else {
    stop("Could not find cell type annotation column!")
  }

  cat(sprintf("Using column: %s\n", subtype_col))
  print(table(fibro.annotated@meta.data[[subtype_col]]))

  cat("\nSubtype by group:\n")
  print(table(fibro.annotated@meta.data[[subtype_col]], fibro.annotated$group))

  fibro.annotated$detailed_celltype <- fibro.annotated@meta.data[[subtype_col]]

} else {
  warning("Fibroblast annotation file not found. Subtype analysis will be skipped.")
  fibro.annotated <- NULL
}

cat("\nPreparing data for CellChat...\n")
if (length(sce.all@assays$RNA@layers) > 1) {
  sce.all <- JoinLayers(sce.all)
  cat("  OK Layers joined\n")
}

if (ENABLE_SPARSE) {
  sce.all <- enable_memory_efficient_mode(sce.all,
                                         sparse = TRUE,
                                         downsample = AUTO_DOWNSAMPLE,
                                         max_cells = MAX_CELLS_PER_TYPE)
}

CF.data <- sce.all[, sce.all$group == "CF"]
Normal.data <- sce.all[, sce.all$group == "Normal"]

cat(sprintf("\nData split:\n"))
cat(sprintf("  CF: %d cells\n", ncol(CF.data)))
cat(sprintf("  Normal: %d cells\n", ncol(Normal.data)))

cat("\n")
estimate_memory_need(CF.data, n_cores = N_CORES)

if (SAVE_INCREMENTAL) {
  cat("\nSaving input data...\n")
  save(sce.all, CF.data, Normal.data, fibro.annotated,
       file = "01.input_data/input.Rdata")
  cat("  OK Input data saved\n")
}

log_performance("Data Loading", start_time_load)

cellchat_analysis <- function(seurat_obj,
                              sample_name,
                              db_type = CELLCHAT_DB,
                              n_cores = N_CORES) {

  cat("\n" %+% strrep("-", 80) %+% "\n")
  cat(sprintf("Processing: %s\n", sample_name))
  cat(strrep("-", 80) %+% "\n")

  start_time <- Sys.time()

  cat("Creating CellChat object...\n")
  data_input <- GetAssayData(seurat_obj, layer = "data", assay = "RNA")
  cellchat <- createCellChat(object = data_input,
                             meta = seurat_obj@meta.data,
                             group.by = "celltype")

  cat(sprintf("  Cell groups: %d\n", length(levels(cellchat@idents))))
  cat(sprintf("  Total cells: %d\n", sum(table(cellchat@idents))))

  cat(sprintf("Loading CellChat database: %s...\n", db_type))
  CellChatDB <- CellChatDB.human

  if (db_type == "all") {
    CellChatDB.use <- CellChatDB
  } else {
    CellChatDB.use <- subsetDB(CellChatDB, search = db_type)
  }
  cellchat@DB <- CellChatDB.use
  cat(sprintf("  L-R pairs in database: %d\n", nrow(CellChatDB.use$interaction)))

  cat("Preprocessing expression data...\n")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cat(sprintf("Computing communication probability (%d cores)...\n", n_cores))
  cellchat <- computeCommunProb(cellchat,
                                type = "triMean",
                                nboot = NBOOT,
                                seed = 123)

  if (MIN_CELLS > 0) {
    cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
  }

  cat("Extracting communication networks...\n")
  df.net <- subsetCommunication(cellchat, slot.name = 'net')
  cat(sprintf("  Total L-R interactions: %d\n", nrow(df.net)))

  cat("Computing pathway-level communication...\n")
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = 'netP')
  cat(sprintf("  Active pathways: %d\n", nrow(df.netp)))

  cat("Aggregating communication network...\n")
  cellchat <- aggregateNet(cellchat)

  cat("Computing network centrality...\n")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  cat("Saving results...\n")
  output_dir <- "02.CF_vs_Normal"
  write.csv(df.net,
           file = file.path(output_dir, paste0(sample_name, "_LR_interactions.csv")),
           row.names = FALSE)
  write.csv(df.netp,
           file = file.path(output_dir, paste0(sample_name, "_pathway_communications.csv")),
           row.names = FALSE)
  saveRDS(cellchat,
         file = file.path(output_dir, paste0(sample_name, "_cellchat.rds")))

  log_performance(paste("CellChat Analysis -", sample_name), start_time)

  return(cellchat)
}

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("STEP 2: CELLCHAT ANALYSIS - CF VS NORMAL\n")
cat(strrep("=", 80) %+% "\n")

CF.cellchat <- cellchat_analysis(CF.data, "CF", db_type = CELLCHAT_DB, n_cores = N_CORES)
Normal.cellchat <- cellchat_analysis(Normal.data, "Normal", db_type = CELLCHAT_DB, n_cores = N_CORES)

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("STEP 3: COMPARATIVE ANALYSIS - CF VS NORMAL\n")
cat(strrep("=", 80) %+% "\n")

start_time_comp <- Sys.time()

object.list <- list(CF = CF.cellchat, Normal = Normal.cellchat)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat.merged, file = "02.CF_vs_Normal/merged_cellchat.rds")

cat("\nComparing interaction numbers and strengths...\n")
pdf("02.CF_vs_Normal/01_overview_interactions.pdf", width = 10, height = 5)
gg1 <- compareInteractions(cellchat.merged, show.legend = TRUE,
                          measure = "count",
                          color.use = c("#E64B35", "#4DBBD5"))
gg2 <- compareInteractions(cellchat.merged, show.legend = TRUE,
                          measure = "weight",
                          color.use = c("#E64B35", "#4DBBD5"))
print(gg1 + gg2)
dev.off()

cat("Generating differential interaction heatmaps...\n")
pdf("02.CF_vs_Normal/02_differential_interactions.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE, measure = "count")
netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE, measure = "weight")
dev.off()

cat("Generating network circle plots...\n")
pdf("04.network_visualization/01_circle_plots_CF_vs_Normal.pdf", width = 14, height = 7)
par(mfrow = c(1,2))
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                  weight.scale = TRUE,
                  label.edge = FALSE,
                  edge.weight.max = weight.max[2],
                  edge.width.max = 12,
                  title.name = paste0("Interaction count - ", names(object.list)[i]))
}
dev.off()

log_performance("Comparative Analysis", start_time_comp)

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("STEP 4: PATHWAY-LEVEL ANALYSIS\n")
cat(strrep("=", 80) %+% "\n")

start_time_pathway <- Sys.time()

cat("Ranking signaling pathways...\n")
pdf("05.pathway_analysis/01_pathway_ranking.pdf", width = 12, height = 8)
gg1 <- rankNet(cellchat.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE,
              color.use = c("#E64B35", "#4DBBD5"))
gg2 <- rankNet(cellchat.merged, mode = "comparison", stacked = FALSE, do.stat = TRUE,
              color.use = c("#E64B35", "#4DBBD5"))
print(gg1 + gg2)
dev.off()

all.pathways <- unique(unlist(lapply(object.list, function(x) x@netP$pathways)))
cat(sprintf("Total pathways detected: %d\n", length(all.pathways)))

key.pathways <- c("TGFb", "PDGF", "FGF", "VEGF", "CTGF", "IL6", "TNF",
                  "CCL", "CXCL", "SPP1", "FN1", "COLLAGEN", "LAMININ")
key.pathways <- intersect(key.pathways, all.pathways)

cat(sprintf("Key pathways for visualization: %s\n", paste(key.pathways, collapse = ", ")))

for (pathway in key.pathways) {
  cat(sprintf("Processing pathway: %s...\n", pathway))

  pathway.exists <- sapply(object.list, function(x) pathway %in% x@netP$pathways)

  if (sum(pathway.exists) > 0) {

    weight.max <- tryCatch({
      getMaxWeight(object.list, slot.name = "netP", attribute = pathway)
    }, error = function(e) {
      cat(sprintf("  Warning: Could not compute max weight for %s, using default\n", pathway))
      return(c(1, 1))
    })

    pdf(file.path("05.pathway_analysis", paste0(pathway, "_circle.pdf")),
        width = 14, height = 7)
    par(mfrow = c(1,2), xpd = TRUE)

    for (i in 1:length(object.list)) {
      if (pathway.exists[i]) {
        tryCatch({
          netVisual_aggregate(object.list[[i]],
                            signaling = pathway,
                            layout = "circle",
                            edge.weight.max = weight.max[1],
                            edge.width.max = 10,
                            signaling.name = paste(pathway, "-", names(object.list)[i]))
        }, error = function(e) {
          plot.new()
          text(0.5, 0.5, paste("Error visualizing", pathway, "\n", names(object.list)[i]), cex = 1.2)
        })
      } else {
        plot.new()
        text(0.5, 0.5, paste("No", pathway, "signaling"), cex = 1.5)
      }
    }
    dev.off()

    tryCatch({
      pdf(file.path("05.pathway_analysis", paste0(pathway, "_heatmap.pdf")),
          width = 14, height = 7)
      par(mfrow = c(1,2), xpd = TRUE)
      ht <- list()
      for (i in 1:length(object.list)) {
        if (pathway.exists[i]) {
          tryCatch({
            ht[[i]] <- netVisual_heatmap(object.list[[i]],
                                         signaling = pathway,
                                         color.heatmap = "Reds",
                                         title.name = paste(pathway, "-", names(object.list)[i]))
          }, error = function(e) {
            cat(sprintf("  Warning: Could not create heatmap for %s in %s\n", pathway, names(object.list)[i]))
          })
        }
      }
      if (length(ht) == 2) {
        ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
      } else if (length(ht) == 1) {
        ComplexHeatmap::draw(ht[[1]])
      }
      dev.off()
    }, error = function(e) {
      cat(sprintf("  Warning: Heatmap visualization failed for %s\n", pathway))
      dev.off()
    })
  }
}

log_performance("Pathway Analysis", start_time_pathway)

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("STEP 5: SIGNALING ROLES ANALYSIS\n")
cat(strrep("=", 80) %+% "\n")

start_time_roles <- Sys.time()

for (i in 1:length(object.list)) {
  sample_name <- names(object.list)[i]
  cat(sprintf("Processing signaling roles for: %s...\n", sample_name))

  pdf(file.path("07.signaling_roles", paste0(sample_name, "_signaling_scatter.pdf")),
      width = 10, height = 8)
  gg <- netAnalysis_signalingRole_scatter(object.list[[i]], label.size = 3)
  print(gg)
  dev.off()

  pdf(file.path("07.signaling_roles", paste0(sample_name, "_signaling_heatmap.pdf")),
      width = 12, height = 10)
  ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                          pattern = "outgoing",
                                          font.size = 8)
  ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                          pattern = "incoming",
                                          font.size = 8)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
}

log_performance("Signaling Roles Analysis", start_time_roles)

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("STEP 6: SUMMARY STATISTICS\n")
cat(strrep("=", 80) %+% "\n")

summary_stats <- data.frame(
  Sample = names(object.list),
  Total_Cells = sapply(object.list, function(x) sum(table(x@idents))),
  Cell_Types = sapply(object.list, function(x) length(levels(x@idents))),
  Total_Interactions = sapply(object.list, function(x) sum(x@net$count)),
  Significant_LR_Pairs = sapply(object.list, function(x) {
    nrow(subsetCommunication(x, slot.name = 'net'))
  }),
  Active_Pathways = sapply(object.list, function(x) length(x@netP$pathways)),
  Avg_Interaction_Strength = sapply(object.list, function(x) {
    mean(x@net$weight[x@net$weight > 0], na.rm = TRUE)
  })
)

write.csv(summary_stats, file = "02.CF_vs_Normal/summary_statistics.csv", row.names = FALSE)

cat("\nSummary Statistics:\n")
print(summary_stats)

pathway_info <- data.frame(
  Pathway = all.pathways,
  Present_in_CF = all.pathways %in% object.list$CF@netP$pathways,
  Present_in_Normal = all.pathways %in% object.list$Normal@netP$pathways
)

write.csv(pathway_info, file = "05.pathway_analysis/pathway_presence.csv", row.names = FALSE)

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 80) %+% "\n")

cleanup_parallel()

final_ram <- get_ram_usage()
final_sys <- get_system_ram_info()

cat(sprintf("\nFinal RAM usage: %.2f GB (%.1f%% of system)\n",
           final_ram, (final_ram / final_sys$total_gb) * 100))

cat("\nAll results saved to organized subdirectories:\n")
for (dir in output_dirs[2:8]) {
  cat(sprintf("  - %s/\n", dir))
}

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("NEXT STEPS:\n")
cat("  1. Review results in 02.CF_vs_Normal/\n")
cat("  2. Check pathway analysis in 05.pathway_analysis/\n")
cat("  3. Examine signaling roles in 07.signaling_roles/\n")
if (!is.null(fibro.annotated)) {
  cat("  4. Run fibroblast subtype analysis (see B.cellchat_fibro_subtypes.R)\n")
}
cat(strrep("=", 80) %+% "\n\n")

message("CellChat analysis completed successfully!")

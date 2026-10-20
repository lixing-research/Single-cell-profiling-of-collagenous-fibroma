setup_parallel_computing <- function(n_cores = 80,
                                    backend = "multicore",
                                    verbose = TRUE) {

  max_cores <- parallel::detectCores()

  if (verbose) {
    cat("\n" %+% strrep("=", 80) %+% "\n")
    cat("PARALLEL COMPUTING SETUP\n")
    cat(strrep("=", 80) %+% "\n")
    cat(sprintf("System cores detected: %d\n", max_cores))
  }

  if (n_cores > 100) {
    warning(sprintf("Requested %d cores exceeds maximum (100). Setting to 100.", n_cores))
    n_cores <- 100
  }

  if (n_cores > max_cores) {
    warning(sprintf("Requested %d cores exceeds available %d cores. Adjusting to %d.",
                   n_cores, max_cores, max_cores - 2))
    n_cores <- max_cores - 2
  }

  if (n_cores < 1) {
    stop("n_cores must be at least 1")
  }

  if (backend == "multicore") {

    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' required. Install with: install.packages('doParallel')")
    }
    library(doParallel)

    cl <- makeCluster(n_cores, type = "FORK")
    registerDoParallel(cl)

    if (verbose) {
      cat(sprintf("Backend: doParallel (multicore)\n"))
    }

  } else if (backend == "multisession") {

    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' required. Install with: install.packages('future')")
    }
    library(future)
    library(future.apply)

    plan(multisession, workers = n_cores)

    if (verbose) {
      cat(sprintf("Backend: future (multisession)\n"))
    }

  } else if (backend == "cluster") {

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    if (verbose) {
      cat(sprintf("Backend: parallel cluster\n"))
    }
  } else {
    stop("Invalid backend. Choose: 'multicore', 'multisession', or 'cluster'")
  }

  Sys.setenv(OMP_NUM_THREADS = n_cores)
  Sys.setenv(MKL_NUM_THREADS = n_cores)
  Sys.setenv(OPENBLAS_NUM_THREADS = n_cores)

  if (verbose) {
    cat(sprintf("Active cores: %d (%.1f%% of system)\n",
               n_cores, (n_cores / max_cores) * 100))
    cat(strrep("=", 80) %+% "\n\n")
  }

  config <- list(
    n_cores = n_cores,
    max_cores = max_cores,
    backend = backend,
    utilization = (n_cores / max_cores) * 100
  )

  return(invisible(config))
}

optimize_cellchat_parallel <- function(cellchat,
                                      n_cores = 80,
                                      nboot = 100,
                                      seed = 123) {

  set.seed(seed)

  message(sprintf("Optimizing CellChat with %d cores", n_cores))

  cellchat <- computeCommunProb(cellchat,
                                type = "triMean",
                                nboot = nboot,
                                cores = n_cores,
                                seed = seed)

  message("Communication probability computed with parallel backend")

  return(cellchat)
}

batch_process <- function(seurat_obj,
                         batch_size = 5000,
                         func,
                         n_cores = 80) {

  n_cells <- ncol(seurat_obj)
  n_batches <- ceiling(n_cells / batch_size)

  message(sprintf("Processing %d cells in %d batches (size: %d)",
                 n_cells, n_batches, batch_size))

  cell_indices <- 1:n_cells
  batches <- split(cell_indices, ceiling(seq_along(cell_indices) / batch_size))

  library(parallel)
  results <- mclapply(batches, function(idx) {
    batch_obj <- seurat_obj[, idx]
    func(batch_obj)
  }, mc.cores = n_cores)

  message("Batch processing complete")

  return(results)
}

auto_tune_cores <- function(seurat_obj,
                           max_ram_gb = 450,
                           initial_cores = 80) {

  source("utils/resource_monitor.R")

  est_mem <- estimate_memory_need(seurat_obj, n_cores = initial_cores)

  if (est_mem <= max_ram_gb) {
    message(sprintf("Initial cores (%d) is safe. Estimated RAM: %.2f GB",
                   initial_cores, est_mem))
    return(initial_cores)
  }

  low <- 10
  high <- 100
  optimal_cores <- initial_cores

  while (low <= high) {
    mid <- floor((low + high) / 2)
    est_mem <- estimate_memory_need(seurat_obj, n_cores = mid)

    if (est_mem <= max_ram_gb * 0.9) {
      optimal_cores <- mid
      low <- mid + 1
    } else {
      high <- mid - 1
    }
  }

  message(sprintf("Auto-tuned cores: %d (estimated RAM: %.2f GB)",
                 optimal_cores,
                 estimate_memory_need(seurat_obj, n_cores = optimal_cores)))

  return(optimal_cores)
}

enable_memory_efficient_mode <- function(seurat_obj,
                                        sparse = TRUE,
                                        downsample = FALSE,
                                        max_cells = 5000) {

  message("Enabling memory-efficient mode...")

  if (sparse) {
    if (!inherits(GetAssayData(seurat_obj, layer = "counts"), "dgCMatrix")) {
      message("Converting to sparse matrix...")
      counts <- as(GetAssayData(seurat_obj, layer = "counts"), "dgCMatrix")
      seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = counts)

      data <- as(GetAssayData(seurat_obj, layer = "data"), "dgCMatrix")
      seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = data)
      message("  Sparse matrix conversion complete")
    }
  }

  if (downsample) {
    cell_counts <- table(seurat_obj$celltype)

    if (any(cell_counts > max_cells)) {
      message("Downsampling large cell types...")

      cells_to_keep <- c()
      for (ct in names(cell_counts)) {
        ct_cells <- colnames(seurat_obj)[seurat_obj$celltype == ct]

        if (length(ct_cells) > max_cells) {
          set.seed(123)
          sampled_cells <- sample(ct_cells, max_cells)
          message(sprintf("  %s: %d -> %d cells",
                         ct, length(ct_cells), max_cells))
          cells_to_keep <- c(cells_to_keep, sampled_cells)
        } else {
          cells_to_keep <- c(cells_to_keep, ct_cells)
        }
      }

      seurat_obj <- seurat_obj[, cells_to_keep]
      message(sprintf("  Total cells after downsampling: %d", ncol(seurat_obj)))
    }
  }

  gc(verbose = FALSE, full = TRUE)

  return(seurat_obj)
}

create_progress_bar <- function(n, title = "Processing") {
  if (requireNamespace("progress", quietly = TRUE)) {
    library(progress)
    pb <- progress_bar$new(
      format = paste0(title, " [:bar] :percent | ETA: :eta | Elapsed: :elapsed"),
      total = n,
      clear = FALSE,
      width = 80
    )
    return(pb)
  } else {
    message("Install 'progress' package for progress bar: install.packages('progress')")
    return(NULL)
  }
}

cleanup_parallel <- function() {

  if (exists("cl", envir = .GlobalEnv)) {
    parallel::stopCluster(get("cl", envir = .GlobalEnv))
    message("Parallel cluster stopped")
  }

  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(sequential)
    message("Future plan reset to sequential")
  }

  gc(verbose = FALSE, full = TRUE)
}

get_optimal_threads <- function(task_type = "cellchat") {
  max_cores <- parallel::detectCores()

  config <- switch(task_type,
    "cellchat" = list(
      cores = min(80, max_cores - 2),
      description = "CellChat communication probability",
      ram_per_core = 5
    ),
    "seurat" = list(
      cores = min(60, max_cores - 2),
      description = "Seurat normalization/scaling",
      ram_per_core = 3
    ),
    "general" = list(
      cores = min(40, max_cores - 2),
      description = "General parallel tasks",
      ram_per_core = 2
    ),
    list(cores = 10, description = "Default", ram_per_core = 2)
  )

  cat("\n")
  cat("Recommended Configuration:\n")
  cat(sprintf("  Task: %s\n", config$description))
  cat(sprintf("  Cores: %d (%.1f%% of %d)\n",
             config$cores,
             (config$cores / max_cores) * 100,
             max_cores))
  cat(sprintf("  Est. RAM: %.1f GB\n", config$cores * config$ram_per_core))
  cat("\n")

  return(config)
}

message("Parallel Configuration Module Loaded Successfully")
message("Available functions:")
message("  - setup_parallel_computing()")
message("  - optimize_cellchat_parallel()")
message("  - batch_process()")
message("  - auto_tune_cores()")
message("  - enable_memory_efficient_mode()")
message("  - cleanup_parallel()")
message("  - get_optimal_threads()")

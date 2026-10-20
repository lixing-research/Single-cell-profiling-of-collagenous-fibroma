get_ram_usage <- function() {

  pid <- Sys.getpid()
  cmd <- sprintf("ps -p %d -o rss= | awk '{print $1/1024/1024}'", pid)
  ram_gb <- as.numeric(system(cmd, intern = TRUE))

  return(ram_gb)
}

get_system_ram_info <- function() {

  meminfo <- readLines("/proc/meminfo")

  extract_value <- function(pattern) {
    line <- grep(pattern, meminfo, value = TRUE)
    as.numeric(gsub("[^0-9]", "", line)) / 1024 / 1024
  }

  total <- extract_value("^MemTotal:")
  free <- extract_value("^MemFree:")
  available <- extract_value("^MemAvailable:")
  used <- total - available

  list(
    total_gb = total,
    used_gb = used,
    free_gb = free,
    available_gb = available,
    used_percent = (used / total) * 100
  )
}

monitor_memory <- function(max_ram_gb = 450,
                          check_interval = 30,
                          log_file = "logs/memory_monitor.log",
                          auto_kill = TRUE) {

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  log_msg <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    full_msg <- sprintf("%s [%s] %s", timestamp, level, msg)
    message(full_msg)
    cat(full_msg, "\n", file = log_file, append = TRUE)
  }

  log_msg(sprintf("Memory monitoring started - Max: %.1f GB, Interval: %d sec",
                  max_ram_gb, check_interval))

  monitor_func <- function() {
    while (TRUE) {
      ram_usage <- get_ram_usage()
      sys_info <- get_system_ram_info()

      log_msg(sprintf("Current RAM: %.2f GB / %.2f GB (%.1f%% of system)",
                     ram_usage, sys_info$total_gb,
                     (ram_usage / sys_info$total_gb) * 100))

      if (ram_usage > max_ram_gb) {
        log_msg(sprintf("CRITICAL: RAM usage (%.2f GB) exceeded limit (%.2f GB)!",
                       ram_usage, max_ram_gb), level = "CRITICAL")

        if (auto_kill) {

          save_file <- sprintf("emergency_save_%s.RData",
                              format(Sys.time(), "%Y%m%d_%H%M%S"))
          log_msg(sprintf("Saving workspace to: %s", save_file), level = "WARN")
          save.image(file = save_file)

          report_file <- "logs/termination_report.txt"
          sink(report_file)
          cat("=" * 80, "\n")
          cat("CELLCHAT ANALYSIS TERMINATED - RAM EXCEEDED\n")
          cat("=" * 80, "\n\n")
          cat("Termination Time:", format(Sys.time()), "\n")
          cat("RAM Usage:", sprintf("%.2f GB", ram_usage), "\n")
          cat("RAM Limit:", sprintf("%.2f GB", max_ram_gb), "\n")
          cat("System RAM:", sprintf("%.2f GB total, %.2f GB available",
                                     sys_info$total_gb, sys_info$available_gb), "\n\n")
          cat("RECOMMENDATION:\n")
          cat("  Your analysis exceeded the RAM limit.\n")
          cat("  Please try one of the following:\n\n")
          cat("  1. Reduce thread count:\n")
          cat("     Current: 80-100 threads suggested\n")
          cat("     Try: 60 threads (N_CORES <- 60)\n\n")
          cat("  2. Enable memory-efficient mode:\n")
          cat("     ENABLE_SPARSE <- TRUE\n")
          cat("     BATCH_SIZE <- 3000\n\n")
          cat("  3. Downsample data:\n")
          cat("     AUTO_DOWNSAMPLE <- TRUE\n")
          cat("     MAX_CELLS_PER_TYPE <- 5000\n\n")
          cat("Workspace saved to:", save_file, "\n")
          cat("=" * 80, "\n")
          sink()

          log_msg(sprintf("Termination report saved to: %s", report_file), level = "CRITICAL")

          cat("\n\n")
          cat("" %+% strrep("=", 78) %+% "\n")
          cat("|  RAM LIMIT EXCEEDED - ANALYSIS TERMINATED                                     |\n")
          cat("" %+% strrep("=", 78) %+% "\n")
          cat(sprintf("|  Current RAM: %.2f GB / Limit: %.2f GB                                |\n",
                     ram_usage, max_ram_gb))
          cat("|                                                                                |\n")
          cat("|  RECOMMENDATION: Reduce thread count and rerun                                |\n")
          cat("|                                                                                |\n")
          cat("|  Try: N_CORES <- 60  (instead of current 80-100)                              |\n")
          cat("|                                                                                |\n")
          cat(sprintf("|  Workspace saved: %s                                    |\n",
                     basename(save_file)))
          cat(sprintf("|  Report: %s                                   |\n",
                     report_file))
          cat("" %+% strrep("=", 78) %+% "\n\n")

          log_msg("Terminating R process...", level = "CRITICAL")
          quit(save = "no", status = 1)
        } else {
          log_msg("Auto-kill disabled. Continuing with warning.", level = "WARN")
        }
      }

      if (ram_usage > max_ram_gb * 0.9 && ram_usage <= max_ram_gb) {
        log_msg(sprintf("WARNING: RAM usage at %.1f%% of limit (%.2f / %.2f GB)",
                       (ram_usage / max_ram_gb) * 100, ram_usage, max_ram_gb),
                level = "WARN")
      }

      Sys.sleep(check_interval)
    }
  }

  return(monitor_func)
}

start_memory_monitor <- function(max_ram_gb = 450, check_interval = 30) {

  if (!requireNamespace("future", quietly = TRUE)) {
    warning("Package 'future' not installed. Install with: install.packages('future')")
    return(NULL)
  }

  library(future)
  plan(multisession, workers = 1)

  monitor_func <- monitor_memory(max_ram_gb = max_ram_gb,
                                 check_interval = check_interval)

  f <- future({
    monitor_func()
  })

  message(sprintf("Memory monitor started (Max: %.1f GB, Interval: %d sec)",
                 max_ram_gb, check_interval))

  return(f)
}

estimate_memory_need <- function(seurat_obj, n_cores = 1) {
  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj)
  n_celltypes <- length(unique(seurat_obj$celltype))

  base_mem_gb <- (n_cells * 1e-6 * 1024)
  lr_mem_gb <- (n_celltypes^2 * 200 * 1e-6 * 1024)
  parallel_mem_gb <- n_cores * 2

  total_est <- base_mem_gb + lr_mem_gb + parallel_mem_gb + 10

  cat("\n=== Memory Requirement Estimation ===\n")
  cat(sprintf("  Cells: %d\n", n_cells))
  cat(sprintf("  Genes: %d\n", n_genes))
  cat(sprintf("  Cell types: %d\n", n_celltypes))
  cat(sprintf("  Cores: %d\n", n_cores))
  cat(sprintf("  \n"))
  cat(sprintf("  Estimated memory:\n"))
  cat(sprintf("    Base (data): %.2f GB\n", base_mem_gb))
  cat(sprintf("    L-R network: %.2f GB\n", lr_mem_gb))
  cat(sprintf("    Parallel overhead: %.2f GB\n", parallel_mem_gb))
  cat(sprintf("    Buffer: 10.00 GB\n"))
  cat(sprintf("  \n"))
  cat(sprintf("  TOTAL ESTIMATED: %.2f GB\n", total_est))
  cat("=====================================\n\n")

  if (total_est > 450) {
    warning(sprintf("Estimated memory (%.2f GB) exceeds limit (450 GB)!", total_est))
    cat("\nRECOMMENDATION:\n")

    safe_cores <- floor((450 - base_mem_gb - lr_mem_gb - 10) / 2)
    safe_cores <- max(10, min(safe_cores, 100))
    cat(sprintf("  - Reduce cores to: %d (currently: %d)\n", safe_cores, n_cores))

    if (n_cells > 50000) {
      safe_cells <- floor(50000 / n_celltypes) * n_celltypes
      cat(sprintf("  - Downsample to: %d cells (currently: %d)\n", safe_cells, n_cells))
    }

    cat("\n")
  } else {
    message(sprintf("Estimated memory (%.2f GB) is within safe limits.", total_est))
  }

  return(total_est)
}

clean_memory <- function(keep_objects = character(0)) {
  all_objs <- ls(envir = .GlobalEnv)
  to_remove <- setdiff(all_objs, keep_objects)

  if (length(to_remove) > 0) {
    rm(list = to_remove, envir = .GlobalEnv)
    gc(verbose = FALSE, full = TRUE)
    message(sprintf("Cleaned %d objects from memory", length(to_remove)))
  }
}

log_performance <- function(task_name, start_time, log_file = "logs/performance.log") {
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  ram_usage <- get_ram_usage()
  sys_info <- get_system_ram_info()

  log_entry <- sprintf(
    "[%s] Task: %s | Duration: %.2f sec | RAM: %.2f GB (%.1f%%)\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    task_name,
    duration,
    ram_usage,
    (ram_usage / sys_info$total_gb) * 100
  )

  cat(log_entry, file = log_file, append = TRUE)

  message(sprintf("%s completed in %.2f seconds (RAM: %.2f GB)",
                 task_name, duration, ram_usage))
}

message("Resource Monitor Module Loaded Successfully")
message("Available functions:")
message("  - get_ram_usage()")
message("  - get_system_ram_info()")
message("  - monitor_memory()")
message("  - start_memory_monitor()")
message("  - estimate_memory_need()")
message("  - clean_memory()")
message("  - log_performance()")

load_packages <- function(packages, quietly = TRUE) {
  not_installed <- packages[!packages %in% installed.packages()[,"Package"]]

  if(length(not_installed) > 0) {
    cat("========================================\n")
    cat("WARNING Warning: The following packages are not installed:\n")
    cat(paste("  -", not_installed, "\n"))
    cat("\nPlease run: source('install_all_packages.R')\n")
    cat("Or install manually:\n")

    bioc_packages <- c("SingleR", "celldex", "scater", "scuttle", "SingleCellExperiment",
                       "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot",
                       "ComplexHeatmap", "limma", "edgeR", "DESeq2", "GSVA", "GSEABase",
                       "AUCell", "RcisTarget", "GENIE3", "BiocParallel", "ExperimentHub")

    github_packages <- c("COSG", "CellChat", "nichenetr", "DoubletFinder", "monocle3")

    for(pkg in not_installed) {
      if(pkg %in% bioc_packages) {
        cat(paste("  BiocManager::install('", pkg, "')\n", sep=""))
      } else if(pkg %in% github_packages) {
        cat(paste("  # See install_all_packages.R for", pkg, "installation\n"))
      } else {
        cat(paste("  install.packages('", pkg, "')\n", sep=""))
      }
    }
    cat("========================================\n")
    stop("Required packages not installed. Please install them first.")
  }

  if(quietly) {
    suppressPackageStartupMessages({
      for(pkg in packages) {
        library(pkg, character.only = TRUE)
      }
    })
  } else {
    for(pkg in packages) {
      library(pkg, character.only = TRUE)
      cat(paste("OK Loaded:", pkg, "\n"))
    }
  }

  invisible(TRUE)
}

load_basic_packages <- function(quietly = TRUE) {
  basic_packages <- c("Seurat", "dplyr", "ggplot2", "cowplot", "patchwork")
  load_packages(basic_packages, quietly = quietly)
}

load_singleR_packages <- function(quietly = TRUE) {
  singleR_packages <- c("SingleR", "celldex", "scater", "scuttle",
                        "SingleCellExperiment", "BiocParallel", "ExperimentHub")
  load_packages(singleR_packages, quietly = quietly)
}

load_enrichment_packages <- function(quietly = TRUE) {
  enrichment_packages <- c("clusterProfiler", "org.Hs.eg.db", "DOSE",
                           "enrichplot", "GSVA", "GSEABase")
  load_packages(enrichment_packages, quietly = quietly)
}

load_visualization_packages <- function(quietly = TRUE) {
  viz_packages <- c("ggplot2", "cowplot", "patchwork", "pheatmap",
                    "RColorBrewer", "viridis", "ComplexHeatmap")
  load_packages(viz_packages, quietly = quietly)
}

load_trajectory_packages <- function(quietly = TRUE) {
  trajectory_packages <- c("monocle3", "princurve", "slingshot")
  if(!"monocle3" %in% installed.packages()[,"Package"]) {
    cat("Note: monocle3 not installed. Some trajectory analyses may not be available.\n")
    trajectory_packages <- trajectory_packages[trajectory_packages != "monocle3"]
  }
  if(length(trajectory_packages) > 0) {
    load_packages(trajectory_packages, quietly = quietly)
  }
}

check_package_versions <- function() {
  cat("========================================\n")
  cat("Package Versions:\n")
  cat("========================================\n")

  key_packages <- c("Seurat", "harmony", "SingleR", "clusterProfiler")

  for(pkg in key_packages) {
    if(pkg %in% installed.packages()[,"Package"]) {
      version <- packageVersion(pkg)
      cat(sprintf("%-20s: %s\n", pkg, version))
    } else {
      cat(sprintf("%-20s: Not installed\n", pkg))
    }
  }
  cat("========================================\n")
}

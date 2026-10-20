if(require("COSG", quietly = TRUE)) {
  library(COSG)

} else {
  cat("Note: COSG not available, will use Seurat methods for marker identification\n")
}

if(require("harmony", quietly = TRUE)) {
  library(harmony)
} else {
  cat("Note: harmony not available, skipping (needed for batch correction)\n")
}

library(ggsci)
library(dplyr)
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)

if(require("celldex", quietly = TRUE)) {
  library(celldex)
} else {
  cat("Note: celldex not available (needed for cell type annotation)\n")
}

if(require("singleseqgset", quietly = TRUE)) {
  library(singleseqgset)
} else {
  cat("Note: singleseqgset not available, continuing without it\n")
}

library(devtools)

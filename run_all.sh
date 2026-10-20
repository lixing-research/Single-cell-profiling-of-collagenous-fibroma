#!/usr/bin/env bash

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"

run_r () {
  local dir="$1"; shift
  for f in "$@"; do
    echo "== ${dir}/${f} =="
    ( cd "$ROOT/$dir" && Rscript "$f" )
  done
}

run_py () {
  local dir="$1"; shift
  for f in "$@"; do
    echo "== ${dir}/${f} =="
    ( cd "$ROOT/$dir" && python "$f" )
  done
}

run_r 1.QC 01_data_import_qc.R
run_r 2.harmony 02_batch_integration.R
run_r 3.clustering_analysis 03.a_clustering_analysis.R 03.c_merged_clustering.R
run_r 4.cell.annotation 04_cell_annotation.R
run_r 6.cell_proporion 06_cell_proportion.R
run_r 7.0.fibro_extrction 07_fibroblast_extraction.R
run_r 8.fibro_DEG 8.fibro_DEG_analysis_final.R
run_r 9.fibro_enrichment 9.fibro_enrichment_analysis.R
run_r 9.a_GSEA/01.0_data_preparation code/01_prepare_pseudobulk.R
run_r 9.a_GSEA/02.0_gsea_overall code/02_run_gsea.R
run_r 9.a_GSEA/04.0_leading_edge code/04_extract_leading_edge.R
run_r 10.1.fibro_trajectory code/01_trajectory_SPF_root.R

run_py 11.RNA_velocity 01_prepare_fibro_data.py 02_fibro_velocity_analysis.py

run_r 12.scenic A_scenic_publication_enhanced.R
run_r 13.cellchat A.cellchat_case_vs_control.R B.cellchat_fibro_subtypes.R
run_r . 16_infercnv_run.R
run_r 7.2.fibro_cnv scripts/generate_cnv_figures_v2.R

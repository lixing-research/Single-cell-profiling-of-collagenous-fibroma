# Single-cell RNA-seq analysis code for collagenous fibroma

Analysis code for a single-cell RNA-seq study of collagenous fibroma. The repository contains the
core R and Python scripts used for quality control, integration, clustering,
cell-type annotation, fibroblast subtyping, differential expression,
pathway analysis, trajectory analysis, RNA velocity, cell-cell communication,
regulon analysis, and expression-inferred copy-number analysis.

No data, patient information, manuscript text, author information, or journal
information are included.

---

## Input data

The analysis starts from a filtered count matrix in 10X-style format:

```text
data/
|-- CF_sample/output/filter_matrix/
|   |-- barcodes.tsv.gz
|   |-- features.tsv.gz
|   `-- matrix.mtx.gz
`-- Normal_sample/output/filter_matrix/
    |-- barcodes.tsv.gz
    |-- features.tsv.gz
    `-- matrix.mtx.gz
```

Data are not included in this repository. If starting from FASTQ, run
`00_upstream_dnbc4tools.sh` (dnbc4tools v2.1.3) to generate the matrices.

---

## Repository structure

```text
.
|-- 00_upstream_dnbc4tools.sh
|-- run_all.sh
|-- 1.QC/01_data_import_qc.R
|-- 2.harmony/02_batch_integration.R
|-- 3.clustering_analysis/
|   |-- 03.a_clustering_analysis.R
|   `-- 03.c_merged_clustering.R
|-- 4.cell.annotation/04_cell_annotation.R
|-- 6.cell_proporion/06_cell_proportion.R
|-- 7.0.fibro_extrction/07_fibroblast_extraction.R
|-- 8.fibro_DEG/8.fibro_DEG_analysis_final.R
|-- 9.fibro_enrichment/9.fibro_enrichment_analysis.R
|-- 9.a_GSEA/
|   |-- 01.0_data_preparation/code/01_prepare_pseudobulk.R
|   |-- 02.0_gsea_overall/code/02_run_gsea.R
|   `-- 04.0_leading_edge/code/04_extract_leading_edge.R
|-- 10.1.fibro_trajectory/code/01_trajectory_SPF_root.R
|-- 11.RNA_velocity/
|   |-- 01_prepare_fibro_data.py
|   |-- 02_fibro_velocity_analysis.py
|   |-- config.py
|   `-- resource_monitor.py
|-- 12.scenic/A_scenic_publication_enhanced.R
|-- 13.cellchat/
|   |-- A.cellchat_case_vs_control.R
|   |-- B.cellchat_fibro_subtypes.R
|   `-- utils/
|-- 7.2.fibro_cnv/scripts/generate_cnv_figures_v2.R
|-- 16_infercnv_run.R
|-- scRNA_scripts/          # shared helper functions
|-- environment/
`-- data/README.md
```

---

## Run order

Each script is executed from its own module directory. `run_all.sh` handles
the working directories:

```bash
bash run_all.sh
```

Manual order. Run each block from the repository root; every command uses a
subshell, so the working directory does not carry over between steps.

```bash
( cd 1.QC                      && Rscript 01_data_import_qc.R )
( cd 2.harmony                 && Rscript 02_batch_integration.R )
( cd 3.clustering_analysis     && Rscript 03.a_clustering_analysis.R )
( cd 3.clustering_analysis     && Rscript 03.c_merged_clustering.R )
( cd 4.cell.annotation         && Rscript 04_cell_annotation.R )
( cd 6.cell_proporion          && Rscript 06_cell_proportion.R )
( cd 7.0.fibro_extrction       && Rscript 07_fibroblast_extraction.R )
( cd 8.fibro_DEG               && Rscript 8.fibro_DEG_analysis_final.R )
( cd 9.fibro_enrichment        && Rscript 9.fibro_enrichment_analysis.R )
( cd 9.a_GSEA/01.0_data_preparation && Rscript code/01_prepare_pseudobulk.R )
( cd 9.a_GSEA/02.0_gsea_overall     && Rscript code/02_run_gsea.R )
( cd 9.a_GSEA/04.0_leading_edge     && Rscript code/04_extract_leading_edge.R )
( cd 10.1.fibro_trajectory     && Rscript code/01_trajectory_SPF_root.R )
```

RNA velocity:

```bash
( cd 11.RNA_velocity && python 01_prepare_fibro_data.py )
( cd 11.RNA_velocity && python 02_fibro_velocity_analysis.py )
```

Downstream analyses require intermediate files (see below):

```bash
( cd 12.scenic          && Rscript A_scenic_publication_enhanced.R )
( cd 13.cellchat        && Rscript A.cellchat_case_vs_control.R )
( cd 13.cellchat        && Rscript B.cellchat_fibro_subtypes.R )
( cd .                  && Rscript 16_infercnv_run.R )
( cd 7.2.fibro_cnv      && Rscript scripts/generate_cnv_figures_v2.R )
```

---

## Prerequisites

| Module | Needs | Produced by |
|---|---|---|
| 1.QC | `data/*/output/filter_matrix/` | dnbc4tools (`00_upstream_dnbc4tools.sh`) |
| 9.a GSEA / 01 | `9.a_GSEA/01.0_data_preparation/input.data.fibroblasts/fibroblast_annotated.rds` | `07_fibroblast_extraction.R` |
| 10.1 trajectory | `10.1.fibro_trajectory/input.data/fibro_sce.rds` | `07_fibroblast_extraction.R` |
| 16 inferCNV | `4.cell.annotation/04_celltype/sce_annotated_merged.rds` | `04_cell_annotation.R` |
| 7.2 CNV | `CNV_results/infercnv_output/*.infercnv_obj` | `16_infercnv_run.R` |
| 9.a GSEA / 02 | `02.0_gsea_overall/input.data.pseudobulk/ranked_genes_overall.rds` | `01_prepare_pseudobulk.R` |
| 11 velocity | `11.RNA_velocity/data/fibro_velocity_raw.h5ad` | `01_prepare_fibro_data.py` |
| 12 SCENIC | `int/scenicOptions.rds`, `int/3.4_regulonAUC.Rds`, `int/cellInfo.Rds`, `output/Step2_regulonTargetsInfo.tsv` | pySCENIC (GRNBoost2/GENIE3, cisTarget, AUCell) |
| 13 CellChat | `13.cellchat/01.input_data/input.Rdata` | `A.cellchat_case_vs_control.R` (first run) |

SCENIC upstream parameters used: GENIE3 on 200 representative cells (50 per
fibroblast subtype); hg38 RefSeq-r80 cisTarget databases; 412 initial
regulons, 290 non-redundant; AUCell on all fibroblasts.

---

## Analysis parameters

- **QC**: `min.cells = 3`, `min.features = 200`; no additional mitochondrial
  or feature filtering.
- **Integration**: Harmony (`theta = 2`).
- **Clustering**: major cell types at resolution 0.2; fibroblast
  subclustering at resolution 0.11 (first resolution giving four clusters).
- **Trajectory**: Monocle2 with the root assigned a priori to the
  SPF-enriched state. The ordering is exploratory and does not establish
  lineage direction.
- **SCENIC**: GENIE3 on 200 cells (50 per subtype); 412 -> 290 regulons;
  AUCell on all fibroblasts.
- **CellChat**: full human CellChatDB, `triMean`, `nboot = 100`,
  `filterCommunication(min.cells = 10)`.
- **Statistics**: with one scRNA-seq specimen, cell-level P values are
  descriptive; they are not between-patient inferences.

---

## Environment

| Component | Version |
|---|---|
| R | 4.5.1 |
| Seurat | 5.3.0 |
| Harmony | 1.0 |
| CellChat | 1.6.1 |
| Monocle2 | 2.24.0 |
| clusterProfiler | 4.12.0 |
| fgsea | 1.34.2 |
| GSVA | 2.2.0 |
| inferCNV | 1.23.0 |
| SCENIC / pySCENIC | 1.3.1 / 0.12.1 |
| scVelo | 0.2.5 |
| dnbc4tools | 2.1.3 |

See `environment/` for details. Newer versions may work with minor adjustments.

---

## License

MIT. See `LICENSE`.

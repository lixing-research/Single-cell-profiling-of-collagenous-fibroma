# Data

Input data are not included in this repository.

Place the filtered count matrix in 10X-style format:

```
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

If starting from FASTQ, run `00_upstream_dnbc4tools.sh` (dnbc4tools v2.1.3)
to generate these matrices. Data access is described in the associated
publication.

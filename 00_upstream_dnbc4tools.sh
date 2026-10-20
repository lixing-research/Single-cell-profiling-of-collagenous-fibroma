#!/usr/bin/env bash

set -euo pipefail

DNBC4="${DNBC4:-dnbc4tools}"
GENOME_DIR="${GENOME_DIR:-/path/to/Homo_sapiens}"
GTF="${GTF:-/path/to/gencode.v32.primary_assembly.annotation.gtf}"
RAW="${RAW:-/path/to/raw_fastq}"
OUTDIR="${OUTDIR:-$(pwd)/data}"

mkdir -p "$OUTDIR"

run_sample () {
  local name="$1"
  local prefix="$2"

  "$DNBC4" rna data \
    --cDNAfastq1  "$RAW/${prefix}_cDNA/${prefix}_cDNA_R1.fq.gz" \
    --cDNAfastq2  "$RAW/${prefix}_cDNA/${prefix}_cDNA_R2.fq.gz" \
    --oligofastq1 "$RAW/${prefix}_oligo/${prefix}_oligo_R1.fq.gz" \
    --oligofastq2 "$RAW/${prefix}_oligo/${prefix}_oligo_R2.fq.gz" \
    --threads 80 --name "$name" --chemistry auto --darkreaction auto \
    --outdir "$OUTDIR" --genomeDir "$GENOME_DIR" --gtf "$GTF" --chrMT chrM

  "$DNBC4" rna count \
    --name "$name" --calling_method emptydrops --expectcells 3000 \
    --threads 80 --outdir "$OUTDIR"

  "$DNBC4" rna analysis \
    --name "$name" --outdir "$OUTDIR" --species Homo_sapiens \
    --mtgenes "$GENOME_DIR/mtgene.list"

  "$DNBC4" rna report \
    --name "$name" --species Homo_sapiens --threads 80 --outdir "$OUTDIR"
}

run_sample CF_sample CF
run_sample Normal_sample Normal

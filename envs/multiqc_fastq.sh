#!/usr/bin/env bash
set -euo pipefail

# Dossiers du projet (adapte si besoin)
PROJ="/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
FASTQC_OUT="$PROJ/results/fastqc"
OUTDIR="$PROJ/results/multiqc"

mkdir -p "$OUTDIR"

echo "==> Génération rapport MultiQC (FASTQ/FASTQC)"
multiqc "$FASTQC_OUT" -o "$OUTDIR" -n "multiqc_fastqc.html"

echo "==> Rapport généré : $OUTDIR/multiqc_fastqc.html"

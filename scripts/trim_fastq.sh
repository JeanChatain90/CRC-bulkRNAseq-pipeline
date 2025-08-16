#!/usr/bin/env bash
set -euo pipefail

MANIFEST="data/meta/manifest.tsv"
IN="data/fastq"
OUT="data/trimmed"
REP="results/fastp"
THREADS=${THREADS:-8}   # change à volonté

mkdir -p "$OUT" "$REP"

# On boucle sur le manifest
awk -F'\t' 'NR>1{print $1,$2,$3,$4}' "$MANIFEST" | while read -r RUN LAYOUT R1 R2; do
  if [[ "$LAYOUT" == "PAIRED" || ( -n "$R1" && -n "$R2" ) ]]; then
    # Paired-end
    in1="$IN/$R1"
    in2="$IN/$R2"
    out1="$OUT/${RUN}_R1.trim.fastq.gz"
    out2="$OUT/${RUN}_R2.trim.fastq.gz"
    html="$REP/${RUN}.html"
    json="$REP/${RUN}.json"

    echo "[fastp] Paired: $RUN"
    fastp \
      -i "$in1" -I "$in2" \
      -o "$out1" -O "$out2" \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front --cut_tail --cut_mean_quality 20 \
      --length_required 30 \
      -q 20 -u 30 \
      -w "$THREADS" \
      -h "$html" -j "$json"

  else
    # Single-end
    in1="$IN/$R1"
    out1="$OUT/${RUN}.trim.fastq.gz"
    html="$REP/${RUN}.html"
    json="$REP/${RUN}.json"

    echo "[fastp] Single: $RUN"
    fastp \
      -i "$in1" \
      -o "$out1" \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front --cut_tail --cut_mean_quality 20 \
      --length_required 30 \
      -q 20 -u 30 \
      -w "$THREADS" \
      -h "$html" -j "$json"
  fi
done

echo "==> Terminé. FASTQ trim dans $OUT ; rapports fastp dans $REP"

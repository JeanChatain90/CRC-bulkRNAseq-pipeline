#!/usr/bin/env bash
#
set -euo pipefail

# Dossiers du projet (adapte si besoin)
PROJ="/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
STAR_OUT="$PROJ/results/star"
REPORT_DIR="$PROJ/results/star"
REPORT_TSV="$REPORT_DIR/star_alignment_summary.tsv"

mkdir -p "$REPORT_DIR"

# En-tête du rapport
echo -e "Run\tTotal_Reads\tUniquely_Mapped\tUniquely_Mapped_%\tMulti_Mapped_%\tUnmapped_%\tAvg_Read_Length" > "$REPORT_TSV"

shopt -s nullglob
for f in "$STAR_OUT"/*Log.final.out; do
  run="$(basename "$f" | cut -d. -f1)"

  total_reads=$(awk -F"|" '/Number of input reads/{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$f")
  uniq_reads=$(awk -F"|" '/Uniquely mapped reads number/{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$f")
  uniq_pct=$(awk -F"|" '/Uniquely mapped reads %/{gsub(/^[ \t]+|[ \t]+%?$/, "", $2); print $2}' "$f")
  multim_pct=$(awk -F"|" '/% of reads mapped to multiple loci/{gsub(/^[ \t]+|[ \t]+%?$/, "", $2); print $2}' "$f")

  # Les “unmapped” sont donnés en 3 lignes : too many mismatches / too short / other
  unm_tooshort=$(awk -F"|" '/% of reads unmapped: too short/{gsub(/^[ \t]+|[ \t]+%?$/, "", $2); print $2}' "$f")
  unm_mismatch=$(awk -F"|" '/% of reads unmapped: too many mismatches/{gsub(/^[ \t]+|[ \t]+%?$/, "", $2); print $2}' "$f")
  unm_other=$(awk -F"|" '/% of reads unmapped: other/{gsub(/^[ \t]+|[ \t]+%?$/, "", $2); print $2}' "$f")
  # Somme robuste (gère champs vides)
  unmapped_pct=$(awk -v a="${unm_tooshort:-0}" -v b="${unm_mismatch:-0}" -v c="${unm_other:-0}" 'BEGIN{print a+b+c}')

  avg_len=$(awk -F"|" '/Average input read length/{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$f")

  echo -e "${run}\t${total_reads}\t${uniq_reads}\t${uniq_pct}\t${multim_pct}\t${unmapped_pct}\t${avg_len}" >> "$REPORT_TSV"
done

echo "==> Rapport STAR compilé : $REPORT_TSV"

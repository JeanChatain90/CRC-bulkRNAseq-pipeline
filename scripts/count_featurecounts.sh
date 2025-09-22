#!/usr/bin/env bash
set -euo pipefail
PROJ="/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
STAR_OUT="$PROJ/results/star"
COUNTS_DIR="$PROJ/results/counts"
META="$PROJ/data/meta/manifest.tsv"
GTF="$HOME/ref/genome/gencode.v43.annotation.gtf"
THREADS="${THREADS:-8}"
mkdir -p "$COUNTS_DIR"/samples
: "${TMPDIR:=$HOME/work/tmp}"; mkdir -p "$TMPDIR"

# map run -> layout
declare -A LAYOUT
tail -n +2 "$META" | awk -F'\t' '{print $1"\t"$2}' | while IFS=$'\t' read -r RUN L; do
  LAYOUT["$RUN"]="$L"
done

# runs présents (BAM)
mapfile -t RUNS < <(find "$STAR_OUT" -maxdepth 1 -type f -name "*.Aligned.sortedByCoord.out.bam" \
  -printf "%f\n" | sed 's/\.Aligned\.sortedByCoord\.out\.bam$//' | sort)

# compter échantillon par échantillon (skip WGS)
for RUN in "${RUNS[@]}"; do
  [[ "$RUN" == "SRR748241" ]] && { echo "[INFO] Skip WGS $RUN"; continue; }
  BAM="$STAR_OUT/${RUN}.Aligned.sortedByCoord.out.bam"
  [[ ! -f "$BAM" ]] && { echo "[WARN] Missing $BAM"; continue; }
  layout="${LAYOUT[$RUN]:-PAIRED}"
  OUT_TXT="$COUNTS_DIR/samples/${RUN}.counts.txt"

  [[ -f "$OUT_TXT" ]] && { echo "[SKIP] $RUN done"; continue; }

  echo "==> featureCounts: $RUN (layout=$layout)"
  if [[ "$layout" == "SINGLE" ]]; then
    featureCounts -T "$THREADS" -a "$GTF" -o "$OUT_TXT" "$BAM" --tmpDir "$TMPDIR"
  else
    featureCounts -T "$THREADS" -a "$GTF" -o "$OUT_TXT" "$BAM" -p -B -C --tmpDir "$TMPDIR"
  fi
done

# fusion gène x échantillon
MERGED="$COUNTS_DIR/gene_counts.tsv"
FIRST=$(ls "$COUNTS_DIR"/samples/*.counts.txt 2>/dev/null | head -n1) || true
[[ -z "${FIRST:-}" ]] && { echo "[ERROR] aucun counts"; exit 1; }

awk 'BEGIN{FS=OFS="\t"} NR>=2{print $1,$7}' "$FIRST" | sort -k1,1 > "$COUNTS_DIR"/_base.tsv

for f in "$COUNTS_DIR"/samples/*.counts.txt; do
  [[ "$f" == "$FIRST" ]] && continue
  awk 'BEGIN{FS=OFS="\t"} NR>=2{print $1,$7}' "$f" | sort -k1,1 > "$COUNTS_DIR"/_one.tsv
  join -t $'\t' -1 1 -2 1 "$COUNTS_DIR"/_base.tsv "$COUNTS_DIR"/_one.tsv > "$COUNTS_DIR"/_join.tsv
  mv "$COUNTS_DIR"/_join.tsv "$COUNTS_DIR"/_base.tsv
done

{
  printf "Geneid"
  for f in "$COUNTS_DIR"/samples/*.counts.txt; do
    hdr=$(head -n1 "$f")
    sample=$(echo "$hdr" | awk 'BEGIN{FS="\t"} {print $7}')
    [[ -z "$sample" || "$sample" == "-" ]] && sample=$(basename "$f" .counts.txt)
    sample=$(basename "$sample" .Aligned.sortedByCoord.out.bam)
    printf "\t%s", sample
  done
  printf "\n"
  cat "$COUNTS_DIR"/_base.tsv
} > "$MERGED"

rm -f "$COUNTS_DIR"/_base.tsv "$COUNTS_DIR"/_one.tsv "$COUNTS_DIR"/_join.tsv
echo "==> Matrice : $MERGED"

#!/usr/bin/env bash
set -euo pipefail

# Répertoires (adapter si besoin)
PROJ="/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
TRIM="$PROJ/data/trimmed"
META="$PROJ/data/meta/manifest.tsv"
OUT="$PROJ/results/star"

# IMPORTANT: tmp sur partition Linux (home) -> évite l’erreur FIFO/NTFS
TMPSTAR="$HOME/work/_STARtmp"

STARIDX="$HOME/ref/STAR_index"
THREADS="${THREADS:-8}"

mkdir -p "$OUT" "$TMPSTAR"

echo "==> Alignement STAR en GRCh38 (index: $STARIDX)"
# Parcours du manifest (sauter l'en-tête)
tail -n +2 "$META" | while IFS=$'\t' read -r RUN LAYOUT R1 R2; do
  R1F="$TRIM/${RUN}_R1.trim.fastq.gz"
  R2F="$TRIM/${RUN}_R2.trim.fastq.gz"

  if [[ ! -f "$R1F" && ! -f "$R2F" ]]; then
    echo "[WARN] Pas de fastq trim pour $RUN → on passe."
    continue
  fi

  pref="$OUT/${RUN}."
  echo "----"
  echo "RUN=$RUN (layout=$LAYOUT)"

  if [[ -f "$R2F" ]]; then
    # Paired-end
    STAR --runThreadN "$THREADS" \
      --genomeDir "$STARIDX" \
      --readFilesIn "$R1F" "$R2F" \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "$pref" \
      --quantMode GeneCounts \
      --twopassMode Basic \
      --outSAMattrRGline ID:$RUN SM:$RUN PL:ILLUMINA \
      --outTmpDir "$TMPSTAR/${RUN}._STARtmp"
  else
    # Single-end
    STAR --runThreadN "$THREADS" \
      --genomeDir "$STARIDX" \
      --readFilesIn "$R1F" \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "$pref" \
      --quantMode GeneCounts \
      --twopassMode Basic \
      --outSAMattrRGline ID:$RUN SM:$RUN PL:ILLUMINA \
      --outTmpDir "$TMPSTAR/${RUN}._STARtmp"
  fi

  # Index BAM
  if [[ -f "${pref}Aligned.sortedByCoord.out.bam" ]]; then
    samtools index -@ "$THREADS" "${pref}Aligned.sortedByCoord.out.bam"
  fi
done

echo "==> Terminé. BAM triés + index dans $OUT"

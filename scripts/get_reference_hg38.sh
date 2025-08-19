#!/usr/bin/env bash
set -euo pipefail

# Dossiers de référence (dans ton $HOME côté WSL)
REF="$HOME/ref"
GENOME="$REF/genome"
STARIDX="$REF/STAR_index"
THREADS="${THREADS:-8}"

mkdir -p "$GENOME" "$STARIDX"

echo "==> (1/3) Téléchargement GENCODE v43 (GRCh38)"
cd "$GENOME"
# FASTA (primary assembly) + GTF
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

echo "==> (2/3) Décompression"
gunzip -f GRCh38.primary_assembly.genome.fa.gz
gunzip -f gencode.v43.annotation.gtf.gz

echo "==> (3/3) Construction index STAR"
STAR --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$STARIDX" \
  --genomeFastaFiles "$GENOME/GRCh38.primary_assembly.genome.fa" \
  --sjdbGTFfile "$GENOME/gencode.v43.annotation.gtf" \
  --sjdbOverhang 100   # adapté à des reads ~101bp

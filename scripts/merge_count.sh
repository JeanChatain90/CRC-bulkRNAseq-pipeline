#!/usr/bin/env bash
set -euo pipefail

# === Config projet ===
PROJ="/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
INDIR="$PROJ/results/counts/samples"      # *.counts.txt (1 fichier par run)
OUTDIR="$PROJ/results/counts"
OUT_TSV="$OUTDIR/gene_counts.tsv"         # matrice finale gènes x échantillons
EXCLUDE_REGEX="^SRR748241$"               # on exclut le WGS

mkdir -p "$OUTDIR"

# === Récupère la liste des fichiers counts (triés), en excluant SRR748241 ===
mapfile -t FILES < <(ls -1 "$INDIR"/*.counts.txt 2>/dev/null | sort || true)
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "[ERR] Aucun fichier *.counts.txt trouvé dans $INDIR"
  exit 1
fi

# Filtrer par nom d’échantillon (basename sans extension)
declare -a KEPT_FILES=()
declare -a SAMPLE_NAMES=()
for f in "${FILES[@]}"; do
  sample="$(basename "$f" .counts.txt)"
  if [[ "$sample" =~ $EXCLUDE_REGEX ]]; then
    echo "[INFO] Exclu: $sample"
    continue
  fi
  KEPT_FILES+=("$f")
  SAMPLE_NAMES+=("$sample")
done

if [[ ${#KEPT_FILES[@]} -eq 0 ]]; then
  echo "[ERR] Tous les fichiers ont été exclus. Vérifie EXCLUDE_REGEX."
  exit 1
fi

# === Construction de la matrice ===
TMPDIR="$(mktemp -d)"
GENES="$TMPDIR/genes.txt"
MATRIX="$TMPDIR/matrix.tsv"

# Fonction utilitaire : extraire Geneid + dernière colonne (comptes)
extract_counts () {
  local in="$1"
  local out="$2"
  # On saute les lignes commençant par '#' et la ligne d'entête "Geneid ..."
  # Puis on imprime: col1 (Geneid) et la dernière colonne (les comptes du sample)
  awk 'BEGIN{OFS="\t"}
       $1 !~ /^#/ && $1 != "Geneid" { print $1, $NF }' "$in" > "$out"
}

# 1er fichier : initialise la liste des gènes + 1ère colonne de la matrice
extract_counts "${KEPT_FILES[0]}" "$TMPDIR/col_0.tsv"
cut -f1 "$TMPDIR/col_0.tsv" > "$GENES"             # Geneid (ordre de référence)
cut -f2 "$TMPDIR/col_0.tsv" > "$MATRIX"            # 1ère colonne = counts 1er sample

# Fichiers suivants : on vérifie la cohérence des gènes puis on ajoute la colonne
for i in $(seq 1 $(( ${#KEPT_FILES[@]} - 1 ))); do
  extract_counts "${KEPT_FILES[$i]}" "$TMPDIR/col_${i}.tsv"

  # Vérifier que l’ordre des gènes est le même (sinon, on aligne explicitement)
  if ! paste "$GENES" <(cut -f1 "$TMPDIR/col_${i}.tsv") \
        | awk -F'\t' 'BEGIN{ok=1} {if($1!=$2){ok=0; exit}} END{exit ok?0:1}'; then
    # Ordre différent : on réaligne par Geneid
    join -t $'\t' -1 1 -2 1 \
      <(sort -k1,1 "$GENES") \
      <(sort -k1,1 "$TMPDIR/col_${i}.tsv") \
      | awk -F'\t' '{print $2}' > "$TMPDIR/col_${i}.only.tsv"

    # Réordonner aussi la matrice courante selon $GENES si besoin
    # (ici on suppose que la matrice suit déjà l’ordre de $GENES)
    paste "$MATRIX" "$TMPDIR/col_${i}.only.tsv" > "$TMPDIR/m_new.tsv"
    mv "$TMPDIR/m_new.tsv" "$MATRIX"
  else
    # Même ordre : on peut juste piocher la 2e colonne et coller
    cut -f2 "$TMPDIR/col_${i}.tsv" > "$TMPDIR/col_${i}.only.tsv"
    paste "$MATRIX" "$TMPDIR/col_${i}.only.tsv" > "$TMPDIR/m_new.tsv"
    mv "$TMPDIR/m_new.tsv" "$MATRIX"
  fi
done

# Entête : gene_id + noms d’échantillons
{
  printf "gene_id"
  for s in "${SAMPLE_NAMES[@]}"; do printf "\t%s" "$s"; done
  printf "\n"
} > "$OUT_TSV"

# Corps : gènes + matrice de comptes
paste "$GENES" "$MATRIX" >> "$OUT_TSV"

echo "==> Matrice fusionnée : $OUT_TSV"
echo "    N échantillons retenus : ${#SAMPLE_NAMES[@]}"
echo "    Exemple d’aperçu :"
head -n 5 "$OUT_TSV"

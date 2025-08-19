#!/usr/bin/env bash
# ==============================================================================
# build_meta.sh
# Récupère les métadonnées ENA pour une étude SRA (par défaut SRP010181),
# infère un statut Tumor/Normal/Unknown pour chaque RUN (SRR) et SAMPLE (biosample),
# et écrit des tables prêtes à l’emploi dans data/meta/.
#
# Usage:
#   scripts/build_meta.sh [SRP_ACCESSION]
#
# Exemples:
#   scripts/build_meta.sh
#   scripts/build_meta.sh SRP010181
#
# Variables d’environnement utiles:
#   ENA_FIELDS   : champs ENA à récupérer (optionnel)
#   TUMOR_RE     : regex (insensible à la casse) pour détecter "Tumor" (optionnel)
#   NORMAL_RE    : regex (insensible à la casse) pour détecter "Normal" (optionnel)
# ==============================================================================

set -euo pipefail

# ---------- Config ----------
ACCESSION="${1:-SRP010181}"
PROJ_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${PROJ_ROOT}/data/meta"

# Champs par défaut (riches et stables)
ENA_FIELDS_DEFAULT="run_accession,sample_accession,sample_alias,sample_title,experiment_title,library_layout,study_accession,study_title,center_name,fastq_ftp,collection_date"
ENA_FIELDS="${ENA_FIELDS:-$ENA_FIELDS_DEFAULT}"

# Heuristiques par défaut (tu peux les surcharger via env si besoin)
TUMOR_RE="${TUMOR_RE:-tumou?r|cancer|carcinoma|adenocarcinoma|lesion|malignan}"
NORMAL_RE="${NORMAL_RE:-normal|adjacent|healthy|control|non[- ]tumou?r}"

# Fichiers de sortie
ENA_TSV="${OUT_DIR}/ena_read_run_full.tsv"
RUN_STATUS_TSV="${OUT_DIR}/run_status.tsv"
SAMPLE_STATUS_TSV="${OUT_DIR}/sample_status.tsv"

# ---------- Pré-checks ----------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' introuvable."; exit 1; }; }
need curl
need awk
need sort

mkdir -p "${OUT_DIR}"

echo "==> Étude : ${ACCESSION}"
echo "==> Dossier sortie : ${OUT_DIR}"
echo "==> Champs ENA    : ${ENA_FIELDS}"

# ---------- (1) Télécharger la table ENA ----------
URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&download=true&fields=${ENA_FIELDS}"
echo "==> (1/5) Téléchargement ENA (${URL})"
# Tentative avec 3 retries
curl -fsSL --retry 3 --retry-delay 2 "${URL}" > "${ENA_TSV}"

# Sanity check
if [[ ! -s "${ENA_TSV}" ]]; then
  echo "ERROR: ${ENA_TSV} est vide. Arrêt."
  exit 1
fi

# ---------- (2) Inférence Tumor/Normal/Unknown par RUN ----------
echo "==> (2/5) Détection du statut (RUN)"
awk -F'\t' -v OFS='\t' \
    -v tumor_re="${TUMOR_RE}" \
    -v normal_re="${NORMAL_RE}" \
    'BEGIN {
       IGNORECASE=1;
     }
     NR==1 {
       for (i=1;i<=NF;i++) { h[$i]=i }
       print "run_accession","sample_accession","sample_alias","status","sample_title","experiment_title";
       next
     }
     {
       # Concatène texte libre (titre échantillon + libellé expérience)
       t = tolower($h["sample_title"] " " $h["experiment_title"]);
       status = "Unknown";
       if (t ~ normal_re) { status = "Normal" }
       if (t ~ tumor_re)  { status = "Tumor" }  # priorité tumor si match mixte
       print $h["run_accession"], $h["sample_accession"], $h["sample_alias"], status, $h["sample_title"], $h["experiment_title"];
     }' "${ENA_TSV}" > "${RUN_STATUS_TSV}"

# ---------- (3) Résumé: comptage par RUN ----------
echo "==> (3/5) Résumé par RUN :"
awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' "${RUN_STATUS_TSV}" | sort

# ---------- (4) Résumé: comptage par SAMPLE (biosample) ----------
echo "==> (4/5) Résumé par SAMPLE :"
awk -F'\t' 'NR>1 && !seen[$2]++{c[$4]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' "${RUN_STATUS_TSV}" | sort

# ---------- (5) Table finale par SAMPLE ----------
echo "==> (5/5) Écriture table par SAMPLE : ${SAMPLE_STATUS_TSV}"
awk -F'\t' -v OFS='\t' '
  NR==1 { next }
  {
    sa = $2; st = $4; al = $3; ti = $5
    if (!(sa in status)) { status[sa]=st; alias[sa]=al }
    title[sa]=ti
  }
  END {
    print "sample_accession","status","sample_alias","example_title";
    for (s in status) { print s, status[s], alias[s], title[s] }
  }' "${RUN_STATUS_TSV}" \
  | sort > "${SAMPLE_STATUS_TSV}"

# ---------- Résumé final ----------
echo "----"
echo "Fichiers créés :"
echo " - ${ENA_TSV}"
echo " - ${RUN_STATUS_TSV}"
echo " - ${SAMPLE_STATUS_TSV}"
echo "----"
echo "Astuce :"
echo "  # Comptes par RUN    : awk -F\"\t\" 'NR>1{c[\$4]++} END{for(k in c) print k, c[k]}' ${RUN_STATUS_TSV}"
echo "  # Comptes par SAMPLE : awk -F\"\t\" 'NR>1 && !seen[\$2]++{c[\$4]++} END{for(k in c) print k, c[k]}' ${RUN_STATUS_TSV}"

#!/usr/bin/env bash
# multiqc_star.sh — Génère un rapport MultiQC à partir des logs STAR
# Usage:
#   bash scripts/multiqc_star.sh [--logs DIR_LOGS] [--out OUTDIR] [--name FICHIER.html]
# Par défaut:
#   --logs  = <racine_du_projet>/results/star
#   --out   = <racine_du_projet>/results/multiqc
#   --name  = multiqc_star.html
#

set -euo pipefail

# Déduire la racine du projet depuis l’emplacement du script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LOGS_DIR="${PROJ_ROOT}/results/star"
OUT_DIR="${PROJ_ROOT}/results/multiqc"
REPORT_NAME="multiqc_star.html"

# Parse arguments simples
while [[ $# -gt 0 ]]; do
  case "$1" in
    --logs)
      LOGS_DIR="$2"; shift 2;;
    --out)
      OUT_DIR="$2"; shift 2;;
    --name)
      REPORT_NAME="$2"; shift 2;;
    -h|--help)
      grep -E '^# ' "$0" | sed 's/^# //'; exit 0;;
    *)
      echo "[ERR] Argument inconnu: $1" >&2; exit 1;;
  esac
done

# Créer dossier de sortie
mkdir -p "$OUT_DIR"

# Assurer un TMPDIR sur partition Linux (utile sous WSL)
if grep -qi microsoft /proc/version 2>/dev/null; then
  # On est probablement sous WSL
  mkdir -p "$HOME/work/tmp"
  export TMPDIR="$HOME/work/tmp"
fi

# Vérifier la présence de logs STAR
N_LOGS=$(find "$LOGS_DIR" -type f -name "Log.final.out" | wc -l | tr -d ' ')
if [[ "$N_LOGS" -eq 0 ]]; then
  echo "[ERR] Aucun fichier Log.final.out trouvé dans: $LOGS_DIR" >&2
  echo "      Lance d’abord l’alignement STAR ou corrige --logs." >&2
  exit 1
fi

echo "==> MultiQC sur ${N_LOGS} logs STAR"
echo "    - Logs : $LOGS_DIR"
echo "    - Sortie : $OUT_DIR/$REPORT_NAME"
echo "    - TMPDIR : ${TMPDIR:-<non défini>}"

multiqc "$LOGS_DIR" \
  -o "$OUT_DIR" \
  -n "$REPORT_NAME" \
  --no-data-dir \
  --force

echo "==> Rapport MultiQC prêt : $OUT_DIR/$REPORT_NAME"

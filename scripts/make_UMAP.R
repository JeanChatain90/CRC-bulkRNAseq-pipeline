#!/usr/bin/env Rscript
# ============================================================
# UMAP pro pour RNA-seq (VST DESeq2)
# - Lit results/deseq2/deseq2_objects.rds (vst ou dds)
# - Sélectionne les gènes les + variables (ntop)
# - UMAP (uwot) avec métrique au choix (cosine par défaut)
# - Exporte: PNG + TSV (coordonnées, métadonnées)
# Usage:
#   Rscript scripts/umap_pro.R \
#     --rds results/deseq2/deseq2_objects.rds \
#     --out-tab results/deseq2 \
#     --out-plot results/plots \
#     --group condition \
#     --ntop 1000 \
#     --neighbors 15 \
#     --min-dist 0.2 \
#     --metric cosine
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(matrixStats)
  library(data.table)
})

# Vérifier uwot (UMAP)
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("Le package 'uwot' n'est pas installé. Dans R: install.packages('uwot')")
}

logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ------------------ #
#   CLI options      #
# ------------------ #
option_list <- list(
  make_option("--rds",      type="character", default="results/deseq2/deseq2_objects.rds",
              help="Chemin du RDS (dds/vst/meta)"),
  make_option("--out-tab",  type="character", default="results/deseq2",
              help="Dossier des TSV (coordonnées)"),
  make_option("--out-plot", type="character", default="results/plots",
              help="Dossier des figures"),
  make_option("--group",    type="character", default="condition",
              help="Colonne de colData utilisée pour la couleur/shape"),
  make_option("--ntop",     type="integer",   default=1000,
              help="Nb de gènes les plus variables pour l’UMAP"),
  make_option("--neighbors",type="integer",   default=15,
              help="n_neighbors UMAP"),
  make_option("--min-dist", type="double",    default=0.2,
              help="min_dist UMAP"),
  make_option("--metric",   type="character", default="cosine",
              help="Métrique UMAP (cosine, euclidean, manhattan, ... )"),
  make_option("--seed",     type="integer",   default=1,
              help="Seed pour la reproductibilité")
)
opt <- parse_args(OptionParser(option_list=option_list))

rds_fp   <- opt$rds
out_tab  <- opt$`out-tab`
out_plot <- opt$`out-plot`
groupvar <- opt$group
ntop     <- opt$ntop
n_neighbors <- opt$neighbors
min_dist <- opt$`min-dist`
metric   <- opt$metric
seed     <- opt$seed

dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# ------------------ #
#   Chargement       #
# ------------------ #
logmsg("Lecture RDS:", rds_fp)
if (!file.exists(rds_fp)) stop("RDS introuvable: ", rds_fp)
obj <- readRDS(rds_fp)

# VST (ou recalcul)
if (!is.null(obj$vst) && inherits(obj$vst, "DESeqTransform")) {
  vst_obj <- obj$vst
  logmsg("VST trouvé dans le RDS.")
} else if (!is.null(obj$dds) && inherits(obj$dds, "DESeqDataSet")) {
  logmsg("VST absent -> recalcul depuis dds (blind=TRUE).")
  vst_obj <- vst(obj$dds, blind = TRUE)
} else stop("Ni 'vst' ni 'dds' disponibles dans le RDS.")

mat <- assay(vst_obj)               # gènes x samples
cd  <- as.data.frame(colData(vst_obj))

if (!groupvar %in% colnames(cd)) {
  stop("La colonne '", groupvar, "' n'existe pas dans colData. Colonnes: ",
       paste(colnames(cd), collapse=", "))
}
cd[[groupvar]] <- droplevels(as.factor(cd[[groupvar]]))

# ------------------ #
#   Sélection gènes  #
# ------------------ #
ntop <- min(ntop, nrow(mat))
rv   <- matrixStats::rowVars(mat)
top  <- order(rv, decreasing = TRUE)[seq_len(ntop)]
mat_top <- t(mat[top, , drop = FALSE])   # samples x gènes

logmsg("UMAP sur", nrow(mat_top), "échantillons x", ncol(mat_top), "gènes (top variance)")

# ------------------ #
#        UMAP        #
# ------------------ #
emb <- uwot::umap(
  X = mat_top,
  n_neighbors = n_neighbors,
  min_dist    = min_dist,
  metric      = metric,
  n_components = 2,
  n_epochs     = 0,            # laisser uwot choisir
  verbose      = TRUE,
  ret_model    = FALSE,
  init         = "spectral",
  fast_sgd     = TRUE
)

colnames(emb) <- c("UMAP1","UMAP2")
df <- data.frame(emb, cd, Sample = rownames(mat_top), check.names = FALSE)

# ------------------ #
#     Sorties        #
# ------------------ #
# TSV coordonnées + métadonnées
f_tsv <- file.path(out_tab, sprintf("UMAP_%s_ntop%d_nn%d_mindist%s_%s.tsv",
                                    groupvar, ntop, n_neighbors,
                                    gsub("\\.", "_", as.character(min_dist)), metric))
data.table::fwrite(df, f_tsv, sep="\t", quote=FALSE)
logmsg("Coordonnées UMAP ->", f_tsv)

# Figure PNG
p <- ggplot(df, aes(UMAP1, UMAP2, color = .data[[groupvar]], shape = .data[[groupvar]])) +
  geom_point(size = 2.2, alpha = 0.9) +
  labs(title = sprintf("UMAP (VST, top %d gènes) — %s", ntop, groupvar),
       color = groupvar, shape = groupvar) +
  theme_bw(base_size = 12)
f_png <- file.path(out_plot, sprintf("UMAP_%s_ntop%d_nn%d_mindist%s_%s.png",
                                     groupvar, ntop, n_neighbors,
                                     gsub("\\.", "_", as.character(min_dist)), metric))
ggsave(f_png, p, width = 7, height = 5.5, dpi = 150)
logmsg("Figure UMAP ->", f_png)

logmsg("✅ Terminé.")

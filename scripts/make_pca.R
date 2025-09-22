#!/usr/bin/env Rscript

# ============================================================
# PCA pour RNA-seq (DESeq2/VST)
# - Lit un RDS produit par le pipeline (contient dds, vst, meta)
# - Recalcule VST si besoin
# - PCA (DESeq2::plotPCA) + PCA "maison" (top gènes variables)
# - Heatmap des distances
# - Exporte tables (scores/loadings/variance) + figures
# Usage :
#   Rscript scripts/pca_pro.R \
#     --rds results/deseq2/deseq2_objects.rds \
#     --out-tab results/deseq2 \
#     --out-plot results/plots \
#     --group condition \
#     --ntop 500
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
  library(data.table)
})

logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ------------------ #
#   CLI options      #
# ------------------ #
option_list <- list(
  make_option("--rds",      type="character", default="results/deseq2/deseq2_objects.rds",
              help="Chemin du RDS (avec dds/vst/meta)"),
  make_option("--out-tab",  type="character", default="results/deseq2",
              help="Dossier des sorties tabulaires (TSV)"),
  make_option("--out-plot", type="character", default="results/plots",
              help="Dossier des figures"),
  make_option("--group",    type="character", default="condition",
              help="Nom de la colonne d’annotation échantillon pour color/shape"),
  make_option("--ntop",     type="integer",   default=500,
              help="Nb de gènes les plus variables pour la PCA 'maison'")
)
opt <- parse_args(OptionParser(option_list=option_list))

rds_fp   <- opt$rds
out_tab  <- opt$`out-tab`
out_plot <- opt$`out-plot`
groupvar <- opt$group
ntop     <- opt$ntop

dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot, recursive = TRUE, showWarnings = FALSE)

# ------------------ #
#   Chargement       #
# ------------------ #
logmsg("Lecture RDS:", rds_fp)
if (!file.exists(rds_fp)) stop("RDS introuvable: ", rds_fp)
obj <- readRDS(rds_fp)

# Récupérer VST ou le recalculer
vst_obj <- NULL
if (!is.null(obj$vst) && inherits(obj$vst, "DESeqTransform")) {
  vst_obj <- obj$vst
  logmsg("VST trouvé dans le RDS (obj$vst).")
} else if (!is.null(obj$dds) && inherits(obj$dds, "DESeqDataSet")) {
  logmsg("VST absent -> recalcul depuis dds (blind=TRUE).")
  vst_obj <- vst(obj$dds, blind = TRUE)
} else {
  stop("Ni 'vst' ni 'dds' disponibles dans le RDS.")
}

# Vérifier le groupvar
cd <- as.data.frame(colData(vst_obj))
if (!groupvar %in% colnames(cd)) {
  stop("La colonne '", groupvar, "' n'existe pas dans colData. Colonnes dispo: ",
       paste(colnames(cd), collapse=", "))
}
cd[[groupvar]] <- droplevels(as.factor(cd[[groupvar]]))

# ------------------ #
#   PCA (DESeq2)     #
# ------------------ #
logmsg("PCA via DESeq2::plotPCA(", groupvar, ")")
p1 <- plotPCA(vst_obj, intgroup = groupvar) +
  ggtitle(paste0("PCA (VST) — ", groupvar)) +
  theme_bw(base_size = 12)
ggsave(file.path(out_plot, paste0("PCA_", groupvar, ".png")),
       p1, width = 7, height = 5.5, dpi = 150)

# ------------------ #
# PCA "maison" top g #
# ------------------ #
mat <- assay(vst_obj)
ntop <- min(ntop, nrow(mat))
rv   <- matrixStats::rowVars(mat)
top  <- order(rv, decreasing = TRUE)[seq_len(ntop)]

logmsg("PCA 'maison' sur", ntop, "gènes les + variables")
pc <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
expl <- 100 * (pc$sdev^2 / sum(pc$sdev^2))
expl2 <- round(expl[1:2], 1)

scores <- data.frame(Sample = rownames(pc$x), pc$x, cd, check.names = FALSE)
loadings <- data.frame(Gene = rownames(pc$rotation), pc$rotation, check.names = FALSE)
var_expl <- data.frame(PC = paste0("PC", seq_along(expl)),
                       Variance_explained_percent = expl, check.names = FALSE)

# Export tables
f_scores   <- file.path(out_tab, "PCA_scores.tsv")
f_loadings <- file.path(out_tab, "PCA_loadings.tsv")
f_var      <- file.path(out_tab, "PCA_variance_explained.tsv")
fwrite(scores,   f_scores, sep="\t", quote=FALSE)
fwrite(loadings, f_loadings, sep="\t", quote=FALSE)
fwrite(var_expl, f_var,      sep="\t", quote=FALSE)
logmsg("Tables écrites:", basename(f_scores), ",", basename(f_loadings), ",", basename(f_var))

# Plot PCA ggplot
p2 <- ggplot(scores, aes(PC1, PC2, color = .data[[groupvar]], shape = .data[[groupvar]])) +
  geom_point(size = 2.5, alpha = 0.9) +
  labs(title = paste0("PCA (VST, top ", ntop, " gènes)"),
       x = paste0("PC1 (", expl2[1], "%)"),
       y = paste0("PC2 (", expl2[2], "%)"),
       color = groupvar, shape = groupvar) +
  theme_bw(base_size = 12)
ggsave(file.path(out_plot, paste0("PCA_top", ntop, "_", groupvar, ".png")),
       p2, width = 7, height = 5.5, dpi = 150)

# ------------------ #
#  Heatmap distances #
# ------------------ #
logmsg("Heatmap des distances (échantillons)")
dist_mat <- dist(t(mat))
dist_mat <- as.matrix(dist_mat)
rownames(dist_mat) <- colnames(mat)
colnames(dist_mat) <- colnames(mat)
pheatmap(dist_mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = cd[, groupvar, drop=FALSE],
         main = paste0("Sample-to-sample distances (VST) — ", groupvar),
         filename = file.path(out_plot, "sample_distance_vst.png"),
         width = 7.5, height = 6)

# ------------------ #
#     Fin            #
# ------------------ #
logmsg("✅ Terminé. Figures ->", out_plot, "| Tables ->", out_tab)

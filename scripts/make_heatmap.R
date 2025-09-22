#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(data.table)
})

rds_file   <- "results/deseq2/deseq2_objects.rds"
out_tab    <- "results/deseq2"
out_plot   <- "results/plots"
groupvar   <- "condition"   # colonne pour l’annotation
use_spearman <- TRUE        # TRUE = 1 - corr Spearman ; FALSE = euclidienne
cutree_k  <- NA             # ex. 4 pour couper en 4 clusters

dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot, recursive = TRUE, showWarnings = FALSE)

cat("[dist_heatmap] Lecture:", rds_file, "\n")
obj <- readRDS(rds_file)
vst_obj <- if (!is.null(obj$vst)) obj$vst else vst(obj$dds, blind = TRUE)

mat <- assay(vst_obj)
cd  <- as.data.frame(colData(vst_obj))

if (!groupvar %in% names(cd)) stop("Colonne '", groupvar, "' absente de colData.")
annot <- cd[, groupvar, drop = FALSE]     # <-- corrigé (pas de virgule finale)
rownames(annot) <- colnames(mat)          # indispensable pour pheatmap

# ---- distance : corr vs euclidienne ----
if (use_spearman) {
  cat("[dist_heatmap] Distance = 1 - Spearman correlation\n")
  cor_mat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  d <- as.dist(1 - cor_mat)
  suffix <- "cor"
} else {
  cat("[dist_heatmap] Distance = euclidienne sur VST\n")
  d <- dist(t(mat), method = "euclidean")
  suffix <- "euclid"
}

# ---- clustering & ordre ----
hc <- hclust(d, method = "ward.D2")
ord <- hc$labels[hc$order]
f_order <- file.path(out_tab, paste0("sample_order_", suffix, ".tsv"))
fwrite(data.frame(sample = ord), f_order, sep = "\t")
cat("[dist_heatmap] Ordre des échantillons -> ", f_order, "\n", sep="")

# ---- Aperçu PNG (sans noms) ----
png_small <- file.path(out_plot, paste0("distance_heatmap_", suffix, "_thumbnail.png"))
pheatmap(as.matrix(d),
         clustering_distance_rows = d,
         clustering_distance_cols = d,
         clustering_method = "ward.D2",
         annotation_col = annot,
         show_rownames = FALSE, show_colnames = FALSE,
         legend = TRUE, border_color = NA,
         main = paste0("Sample-to-sample distances (VST) — ", groupvar),
         filename = png_small, width = 2000/96, height = 1800/96)
cat("[dist_heatmap] Aperçu -> ", png_small, "\n", sep="")

# ---- PDF grand format ----
pdf_big <- file.path(out_plot, paste0("distance_heatmap_", suffix, "_large.pdf"))
pheatmap(as.matrix(d),
         clustering_distance_rows = d,
         clustering_distance_cols = d,
         clustering_method = "ward.D2",
         annotation_col = annot,
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA, legend = TRUE,
         cutree_rows = if (is.na(cutree_k)) NA else cutree_k,
         cutree_cols = if (is.na(cutree_k)) NA else cutree_k,
         main = paste0("Sample-to-sample distances (VST) — ", groupvar),
         filename = pdf_big, width = 12, height = 10)
cat("[dist_heatmap] PDF -> ", pdf_big, "\n", sep="")

cat("[dist_heatmap] ✅ Terminé\n")

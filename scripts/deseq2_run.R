#!/usr/bin/env Rscript

# ============================================================
#  DESeq2 - Pipeline pro (CRC / SRP010181)
#  - Lecture counts featureCounts (gene_id + SRR*)
#  - Lecture méta: run_status.tsv (run_accession, status) ou coldata_runs.tsv
#  - QC: library size, PCA (VST), distances
#  - DE: Tumor vs Normal, shrinkage (apeglm), MA, volcano
#  - Annotation: org.Hs.eg.db (ENSEMBL -> SYMBOL, ENTREZID, GENENAME)
#  - Sauvegardes: tables, figures, RDS, sessionInfo
#  - Design auto: ~ patient_id + condition | ~ batch + condition | ~ condition
#  Usage :
#    Rscript scripts/deseq2_pro.R \
#      --proj "/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline" \
#      --counts "results/counts/gene_counts.tsv" \
#      --meta "data/meta/run_status.tsv" \
#      --out "results/deseq2" \
#      --plots "results/plots" \
#      --annot TRUE
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(apeglm)         # shrinkage
  library(matrixStats)    # rowVars
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ------------------ #
#   CLI options      #
# ------------------ #
option_list <- list(
  make_option("--proj",   type="character", default=getwd(), help="Dossier projet (racine)"),
  make_option("--counts", type="character", default="results/counts/gene_counts.tsv",
              help="Chemin vers le TSV de counts (featureCounts)"),
  make_option("--meta",   type="character", default=NULL,
              help="Chemin vers la méta. Par défaut cherche data/meta/run_status.tsv puis results/counts/coldata_runs.tsv"),
  make_option("--out",    type="character", default="results/deseq2", help="Dossier résultats tabulaires"),
  make_option("--plots",  type="character", default="results/plots",  help="Dossier figures"),
  make_option("--annot",  type="logical",   default=TRUE, help="Annoter les gènes (TRUE/FALSE)"),
  make_option("--minCount", type="integer", default=10, help="Seuil de count par échantillon"),
  make_option("--minProp",  type="double",  default=0.20, help="Proportion d'échantillons pour filtrage (0-1)")
)
opt <- parse_args(OptionParser(option_list=option_list))

proj_dir   <- normalizePath(opt$proj, mustWork = FALSE)
counts_fp  <- file.path(proj_dir, opt$counts)
out_dir    <- file.path(proj_dir, opt$out)
plot_dir   <- file.path(proj_dir, opt$plots)
meta_fp    <- if (!is.null(opt$meta)) file.path(proj_dir, opt$meta) else NA_character_
do_annot   <- isTRUE(opt$annot)
min_count  <- opt$minCount
min_prop   <- opt$minProp

dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
contrast_name   <- "condition"
contrast_levels <- c("Normal","Tumor")  # référence -> comparé

# ------------------ #
#    Helpers         #
# ------------------ #
sanitize_names <- function(x) {
  x <- sub("^.*/", "", x)
  x <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", x)
  x <- sub("\\.Aligned\\.out\\.bam$", "", x)
  x <- sub("\\.bam$", "", x)
  x
}

read_meta_smart <- function(meta_fp, proj_dir, counts_colnames) {
  # priorité: meta_fp si fourni ; sinon essaie run_status.tsv ; sinon coldata_runs.tsv
  if (is.na(meta_fp) || !file.exists(meta_fp)) {
    cand1 <- file.path(proj_dir, "data/meta/run_status.tsv")
    cand2 <- file.path(proj_dir, "results/counts/coldata_runs.tsv")
    meta_fp <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else NA_character_
  }
  if (is.na(meta_fp)) stop("Aucun fichier méta trouvé. Fournis --meta ou crée data/meta/run_status.tsv.")

  meta <- fread(meta_fp)
  # Détection de la colonne identifiant (SRR…)
  sample_id_candidates <- c(
    "sample_id","run","Run","RUN","SRR","srr","accession","Accession",
    "run_accession","Run_accession","RUN_ACCESSION",
    "sample","sample_accession","ENA_RUN","SRA_RUN"
  )
  sid <- intersect(names(meta), sample_id_candidates)
  if (length(sid) == 0) {
    inter_sizes <- sapply(names(meta), function(cn) sum(as.character(meta[[cn]]) %in% counts_colnames))
    if (max(inter_sizes) == 0) {
      stop("Impossible d'identifier la colonne échantillon (SRR). Vérifie ta méta (run_accession attendu).")
    }
    sid <- names(meta)[which.max(inter_sizes)]
    logmsg("Colonne échantillon déduite automatiquement:", sid)
  }
  setnames(meta, sid[1], "sample_id")

  # Condition
  cond_candidates <- c("condition","status","group","phenotype")
  cc <- intersect(names(meta), cond_candidates)
  if (length(cc) == 0) stop("Aucune colonne condition/status/group/phenotype dans la méta (Tumor/Normal attendu).")
  setnames(meta, cc[1], "condition")

  # Optionnels
  if ("batch" %in% names(meta))  meta[, batch := factor(batch)]
  if ("patient_id" %in% names(meta)) meta[, patient_id := factor(patient_id)]

  # Normaliser labels
  meta[, condition := factor(condition)]
  levels(meta$condition) <- sub("^tumou?r$","Tumor", tolower(levels(meta$condition)))
  levels(meta$condition) <- sub("^normal$","Normal", levels(meta$condition))
  # Exclure valeurs inconnues si présentes
  if (any(levels(meta$condition) %in% c("unknown","Unknown"))) {
    meta <- meta[condition %in% c("Normal","Tumor")]
    meta[, condition := droplevels(condition)]
  }

  if (!all(contrast_levels %in% levels(meta$condition))) {
    stop("Niveaux attendus absents (Normal/Tumor). Niveaux trouvés: ",
         paste(levels(meta$condition), collapse=", "))
  }
  meta[, condition := factor(as.character(condition), levels = contrast_levels)]
  return(meta[])
}

coef_index <- function(dds, contrast_name, contrast_levels) {
  rn <- resultsNames(dds)
  idx <- grep(paste0("^", contrast_name, "_", contrast_levels[2], "_vs_", contrast_levels[1], "$"), rn)
  if (length(idx) == 0) idx <- grep(paste0("^", contrast_name, ".?\\Q", contrast_levels[2], "\\E"), rn)
  if (length(idx) == 0) stop("Coef introuvable pour lfcShrink. resultsNames: ", paste(rn, collapse=", "))
  idx[1]
}

annotate_results <- function(dt) {
  # Ensembl sans version
  ens <- sub("\\.\\d+$","", dt$gene_id)
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = ens,
                               keytype = "ENSEMBL",
                               columns = c("SYMBOL","ENTREZID","GENENAME"))
  map <- unique(map)
  dt$ENSEMBL <- ens
  merge(dt, map, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
}

# ------------------ #
#  1) Lectures       #
# ------------------ #
logmsg("Lecture counts:", counts_fp)
if (!file.exists(counts_fp)) stop("Counts introuvables: ", counts_fp)
ct <- fread(counts_fp)

gene_col_candidates <- c("Geneid","gene_id","GeneID","ENSEMBL","gene")
gene_col <- intersect(names(ct), gene_col_candidates)
gene_col <- if (length(gene_col) == 0) names(ct)[1] else gene_col[1]
logmsg("Colonne gène détectée:", gene_col)

annot_cols <- c("Chr","Start","End","Strand","Length","chr","start","end","strand","length")
keep_cols  <- c(gene_col, setdiff(names(ct), annot_cols))
ct <- ct[, ..keep_cols]

gene_ids    <- as.character(ct[[gene_col]])
sample_cols <- setdiff(names(ct), gene_col)
counts      <- as.data.frame(ct[, ..sample_cols])
rownames(counts) <- gene_ids
colnames(counts) <- sanitize_names(colnames(counts))

# Vérif entiers
is_integer_like <- function(v) all(grepl("^[0-9]+$", as.character(v)))
non_numeric_cols <- names(counts)[vapply(counts, function(z) !is_integer_like(z), logical(1))]
if (length(non_numeric_cols) > 0) {
  stop("Colonnes non numériques détectées (extrait): ",
       paste(head(non_numeric_cols, 10), collapse=", "),
       "\nVérifie ton 'gene_counts.tsv' (comptes entiers après la colonne gène).")
}
counts[] <- lapply(counts, function(x) as.integer(as.character(x)))

# Méta
logmsg("Recherche/lecture méta…")
meta <- read_meta_smart(meta_fp, proj_dir, colnames(counts))

# Intersection
common <- intersect(colnames(counts), meta$sample_id)
if (length(common) < 4)
  stop("Trop peu d’échantillons communs counts/meta: ", length(common))
counts <- counts[, common, drop=FALSE]
meta   <- meta[match(common, meta$sample_id)]

# ------------------ #
#  2) Filtrage       #
# ------------------ #
nsamps <- ncol(counts)
min_samples <- max(1, floor(min_prop * nsamps))
keep <- rowSums(counts >= min_count) >= min_samples
logmsg("Pré-filtrage:", sum(keep), "gènes conservés /", nrow(counts),
       sprintf("(%.1f%%)", 100*mean(keep)))
counts <- counts[keep, , drop=FALSE]

# ------------------ #
#  3) Design         #
# ------------------ #
has_batch   <- "batch" %in% names(meta)
has_patient <- "patient_id" %in% names(meta)

if (has_patient) {
  design_formula <- ~ patient_id + condition
} else if (has_batch) {
  design_formula <- ~ batch + condition
} else {
  design_formula <- ~ condition
}
logmsg("Design utilisé:", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = as.data.frame(meta),
                              design    = design_formula)

# ------------------ #
#  4) DESeq2 & QC    #
# ------------------ #
dds <- DESeq(dds, parallel = FALSE)

# Library size boxplot
libsize <- colSums(counts(dds))
df_lib <- data.frame(sample=names(libsize), libsize=as.numeric(libsize),
                     condition=colData(dds)[,contrast_name])
p_lib <- ggplot(df_lib, aes(x=condition, y=libsize)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_y_log10() + theme_bw(base_size = 12) +
  labs(title="Library sizes (post-filter)", x="", y="Reads (log10)")
ggsave(file.path(plot_dir, "libsize_boxplot.png"), p_lib, width=6, height=4, dpi=150)

# VST (blind TRUE si confondeur patient/batch)
vst_assay <- vst(dds, blind = (has_batch || has_patient))
vst_mat <- assay(vst_assay)

# PCA (top 500)
ntop <- min(500, nrow(vst_mat))
rv <- matrixStats::rowVars(vst_mat)
top <- order(rv, decreasing = TRUE)[seq_len(ntop)]
pca <- prcomp(t(vst_mat[top, ]))
pca_df <- data.frame(pca$x[,1:2], colData(dds))
expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
p_pca <- ggplot(pca_df, aes(PC1, PC2, shape = !!as.name(contrast_name), color = !!as.name(contrast_name))) +
  geom_point(size=2.5, alpha=0.8) +
  labs(title = "PCA (VST, top500)",
       x = paste0("PC1 (", expl[1], "%)"),
       y = paste0("PC2 (", expl[2], "%)")) +
  theme_bw(base_size = 12)
ggsave(file.path(plot_dir, "pca_vst_top500.png"), p_pca, width=6, height=5, dpi=150)

# Distances inter-échantillons
dist_mat <- dist(t(vst_mat))
dist_mat <- as.matrix(dist_mat)
rownames(dist_mat) <- colnames(vst_mat)
colnames(dist_mat) <- colnames(vst_mat)
pheatmap(dist_mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Sample-to-sample distances (VST)",
         filename = file.path(plot_dir, "sample_distance_vst.png"),
         width=7, height=6)

# MA plot (sur résultats bruts)
png(file.path(plot_dir, "MAplot_deseq2.png"), width=800, height=600)
plotMA(dds, main="MA-plot (DESeq2)", ylim=c(-5,5))
dev.off()

# ------------------ #
#  5) Résultats      #
# ------------------ #
# Résultats bruts (avec pvalue/padj/stat)
res_raw <- results(dds, contrast = c(contrast_name, contrast_levels[2], contrast_levels[1]))

# Shrinkage LFC
coef_idx <- coef_index(dds, contrast_name, contrast_levels)
res_shrunk <- lfcShrink(dds, coef = coef_idx, type = "apeglm")

# Fusion : pvalue/padj/stat depuis res_raw, LFC/lfcSE depuis res_shrunk
res <- as.data.frame(res_raw)
res$log2FoldChange <- as.numeric(res_shrunk$log2FoldChange)
if ("lfcSE" %in% colnames(res_shrunk)) res$lfcSE <- as.numeric(res_shrunk$lfcSE)
res$gene_id <- rownames(res)

final_cols <- intersect(c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"),
                        colnames(res))
ord <- order(is.na(res$padj), res$padj, res$pvalue, na.last = TRUE)
res <- res[ord, final_cols, drop=FALSE]

# Export tables brutes
f_all <- file.path(out_dir, "deseq2_results_all.tsv")
fwrite(res, f_all, sep="\t", quote=FALSE, na="NA")

deg <- subset(res, ("padj" %in% names(res)) & !is.na(padj) & padj < 0.05 &
                    ("log2FoldChange" %in% names(res)) & abs(log2FoldChange) >= 1)
f_deg <- file.path(out_dir, "deseq2_DEGs_filtered.tsv")
fwrite(deg, f_deg, sep="\t", quote=FALSE, na="NA")

logmsg("Résultats écrits:")
logmsg(" -", f_all)
logmsg(" -", f_deg, paste0("(", nrow(deg), " DEGs)"))

# Volcano (si pvalue dispo)
if (all(c("log2FoldChange","pvalue") %in% names(res))) {
  volc_df <- within(res, {
    sig <- ifelse(("padj" %in% names(res)) & !is.na(padj) &
                    padj < 0.05 & abs(log2FoldChange) >= 1, "DEG", "NS")
    neglog10p <- -log10(pvalue)
  })
  p_volc <- ggplot(volc_df, aes(log2FoldChange, neglog10p)) +
    geom_point(aes(alpha = sig == "DEG"), size = 1.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = "Volcano plot (Tumor vs Normal)",
         x = "log2 Fold Change",
         y = "-log10(p-value)") +
    theme_bw(base_size = 12)
  ggsave(file.path(plot_dir, "volcano_tumor_vs_normal.png"), p_volc, width=6, height=5, dpi=150)
} else {
  logmsg("(!) Volcano non tracé: colonnes 'pvalue' ou 'log2FoldChange' manquantes.")
}

# ------------------ #
#  6) Annotation     #
# ------------------ #
if (do_annot) {
  logmsg("Annotation (org.Hs.eg.db)…")
  res_annot <- annotate_results(as.data.frame(res))
  deg_annot <- annotate_results(as.data.frame(deg))
  fwrite(res_annot, file.path(out_dir, "deseq2_results_all_annot.tsv"), sep="\t", quote=FALSE, na="NA")
  fwrite(deg_annot, file.path(out_dir, "deseq2_DEGs_filtered_annot.tsv"), sep="\t", quote=FALSE, na="NA")
  logmsg(" -> écrit:", file.path(out_dir, "deseq2_results_all_annot.tsv"))
  logmsg(" -> écrit:", file.path(out_dir, "deseq2_DEGs_filtered_annot.tsv"))
}

# ------------------ #
#  7) Sauvegardes    #
# ------------------ #
saveRDS(list(dds=dds, vst=vst_assay, res=res, deg=deg, meta=meta),
        file = file.path(out_dir, "deseq2_objects.rds"))

sink(file.path(out_dir, "sessionInfo.txt")); sessionInfo(); sink()
logmsg("✅ Terminé. Figures:", plot_dir, " | Tables:", out_dir)

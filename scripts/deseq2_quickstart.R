#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

proj <- "/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
counts_file <- file.path(proj, "results/counts/gene_counts.tsv")
coldata_file <- file.path(proj, "results/counts/coldata_runs.tsv")
outdir <- file.path(proj, "results/deseq2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1) Charger counts
cat("Reading counts...\n")
cts <- read_tsv(counts_file, show_col_types = FALSE)
stopifnot("gene_id" %in% names(cts))
rownames_mat <- cts$gene_id
cts <- as.data.frame(cts[,-1])
rownames(cts) <- rownames_mat

# 2) Charger colData
cat("Reading colData...\n")
coldata <- read_tsv(coldata_file, show_col_types = FALSE) %>%
  distinct(run_accession, .keep_all = TRUE)
rownames(coldata) <- coldata$run_accession

# Aligner colonnes et colData (et enlever SRR manquants si besoin)
common <- intersect(colnames(cts), rownames(coldata))
cts <- cts[, common, drop=FALSE]
coldata <- coldata[common, , drop=FALSE]

# 3) Filtre de base (peu exprimés)
keep <- rowSums(cts) >= 10   # gènes avec au moins 10 reads au total
cts_f <- cts[keep, , drop=FALSE]
cat("Genes kept after filtering:", nrow(cts_f), "of", nrow(cts), "\n")

# 4) Design simple ~ status (on ne garde que Tumor/Normal)
coldata$status <- factor(coldata$status, levels = c("Normal","Tumor","Unknown"))
keep_samples <- coldata$status %in% c("Normal","Tumor")
cts_f <- cts_f[, keep_samples, drop=FALSE]
coldata2 <- droplevels(coldata[keep_samples, , drop=FALSE])

# Vérifier qu’on a encore des colonnes
stopifnot(ncol(cts_f) > 1)

dds <- DESeqDataSetFromMatrix(countData = round(cts_f),
                              colData = coldata2,
                              design = ~ status)

# 5) QC rapide : tailles de bibliothèques
libsize <- colSums(counts(dds))
df_lib <- tibble(sample = names(libsize), libsize = as.numeric(libsize), status = coldata2$status)
p_lib <- ggplot(df_lib, aes(x=status, y=libsize)) + geom_boxplot() + scale_y_log10() +
  ggtitle("Library sizes (post-filter)") + theme_bw()
ggsave(file.path(outdir, "libsize_boxplot.png"), p_lib, width = 6, height = 4, dpi = 150)

# 6) Variance stabilisation + PCA
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "status")
ggsave(file.path(outdir, "PCA_status.png"), pca, width = 6, height = 5, dpi = 150)

# 7) DE simple Tumor vs Normal (si 2 groupes)
if (nlevels(droplevels(coldata2$status)) == 2) {
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("status", "Tumor", "Normal"))
  res <- lfcShrink(dds, coef="status_Tumor_vs_Normal", res=res, type="apeglm")
  res_tbl <- as_tibble(res, rownames = "gene_id")
  write_tsv(res_tbl, file.path(outdir, "deseq2_Tumor_vs_Normal.tsv"))
  png(file.path(outdir, "MAplot_Tumor_vs_Normal.png"), width=800, height=600)
  DESeq2::plotMA(res, ylim=c(-5,5))
  dev.off()
  cat("DE results written to results/deseq2/\n")
} else {
  cat("Skipping DE: need exactly two groups (Normal vs Tumor).\n")
}

# 8) Sauvegarder objets utiles
saveRDS(dds, file.path(outdir, "dds.rds"))
saveRDS(vsd, file.path(outdir, "vsd.rds"))

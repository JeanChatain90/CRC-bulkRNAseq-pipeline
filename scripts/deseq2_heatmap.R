#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(DESeq2); library(pheatmap); library(RColorBrewer)
  library(tibble); library(matrixStats)
})

proj <- "/mnt/c/Users/chata/Documents/GitHub/CRC-bulkRNAseq-pipeline"
counts_file <- file.path(proj, "results/counts/gene_counts.tsv")
meta_file   <- file.path(proj, "data/meta/run_status.tsv")

message("Lecture des matrices et métadonnées…")
cts  <- read_tsv(counts_file, show_col_types = FALSE)
meta <- read_tsv(meta_file,   show_col_types = FALSE)

stopifnot("gene_id" %in% colnames(cts))
stopifnot(all(c("run_accession","status") %in% colnames(meta)))

runs_in_counts <- setdiff(colnames(cts), "gene_id")

meta_keep <- meta %>%
  filter(status %in% c("Normal","Tumor")) %>%
  select(run_accession, status)

keep_runs <- intersect(runs_in_counts, meta_keep$run_accession)
if(length(keep_runs) < 3){
  stop("Trop peu d’échantillons Tumor/Normal présents dans gene_counts.tsv.")
}

meta_keep <- meta_keep %>%
  filter(run_accession %in% keep_runs) %>%
  distinct(run_accession, .keep_all = TRUE) %>%
  arrange(run_accession)

cts_mat <- cts %>%
  select(gene_id, all_of(keep_runs)) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

coldata <- meta_keep %>% tibble::column_to_rownames("run_accession")
coldata$status <- factor(coldata$status, levels = c("Normal","Tumor"))

dds <- DESeqDataSetFromMatrix(countData = round(cts_mat),
                              colData = coldata,
                              design = ~ status)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet = TRUE)

vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

rv <- rowVars(mat)
sel <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
mat_sel <- mat[sel, ]

ann_col <- data.frame(Status = coldata$status)
rownames(ann_col) <- rownames(coldata)

paletteLength <- 50
myColor <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(paletteLength)

out_png <- file.path(proj, "results/plots/heatmap_topvar_vst.png")
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

png(out_png, width = 1400, height = 1200, res = 120)
pheatmap(mat_sel,
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_col,
         color = myColor,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Top variable genes (VST) - Tumor vs Normal (run-level)")
dev.off()

message("✅ Heatmap écrite : ", out_png)

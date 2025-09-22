#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

out_dir  <- "results/deseq2"
plot_dir <- "results/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

res <- fread(file.path(out_dir, "deseq2_results_all.tsv"))

stopifnot(all(c("log2FoldChange","pvalue","padj") %in% names(res)))

volc_df <- within(res, {
  sig <- ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, "DEG", "NS")
  neglog10p <- -log10(pvalue)
})

p_volc <- ggplot(volc_df, aes(log2FoldChange, neglog10p, color = sig)) +
  geom_point(size = 1.4) +
  scale_color_manual(values = c("NS" = "grey50", "DEG" = "red3")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano plot (Tumor vs Normal)",
       x = "log2 Fold Change", y = "-log10(p-value)", color = NULL) +
  theme_bw(base_size = 12)

ggsave(file.path(plot_dir, "volcano_tumor_vs_normal.png"),
       p_volc, width = 6, height = 5, dpi = 150)

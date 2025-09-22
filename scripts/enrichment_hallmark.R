#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# -------- CLI --------
option_list <- list(
  make_option("--all",   type="character", default="results/deseq2/deseq2_results_all_annot.tsv",
              help="TSV résultats DESeq2 complets annotés (avec SYMBOL, log2FC, pvalue, stat?)"),
  make_option("--out",   type="character", default="results/deseq2", help="Dossier TSV"),
  make_option("--plots", type="character", default="results/plots",  help="Dossier figures"),
  make_option("--minSize", type="integer", default=15,  help="taille min des pathways"),
  make_option("--maxSize", type="integer", default=500, help="taille max des pathways"),
  make_option("--topN",    type="integer", default=20,  help="Nb pathways affichés dans la figure")
)
opt <- parse_args(OptionParser(option_list = option_list))

f_all   <- opt$all
out_dir <- opt$out
plt_dir <- opt$plots
minSize <- opt$minSize
maxSize <- opt$maxSize
topN    <- opt$topN

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plt_dir, recursive = TRUE, showWarnings = FALSE)

# -------- Lecture + Ranking --------
logmsg("[H] Lecture:", f_all)
res <- fread(f_all)
need <- c("SYMBOL","log2FoldChange","pvalue")
miss <- setdiff(need, names(res))
if (length(miss)) stop("[H] Colonnes manquantes: ", paste(miss, collapse=", "))

if ("stat" %in% names(res) && !all(is.na(res$stat))) {
  res$rank <- res$stat
  logmsg("[H] Ranking = 'stat'")
} else {
  res$pvalue <- pmax(res$pvalue, 1e-300)
  res$rank <- sign(res$log2FoldChange) * -log10(res$pvalue)
  logmsg("[H] Ranking = sign(LFC) * -log10(pval)")
}
res <- res[!is.na(SYMBOL) & !is.na(rank)]
# gérer doublons de SYMBOL : on garde le rang le plus extrême
res <- res %>% group_by(SYMBOL) %>% summarise(rank = rank[which.max(abs(rank))], .groups="drop")
ranks <- res$rank; names(ranks) <- res$SYMBOL
ranks <- sort(ranks, decreasing = TRUE)
logmsg("[H] Gènes classés:", length(ranks))

# -------- Hallmarks --------
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
          split(.$gene_symbol, .$gs_name)

logmsg("[H] fgsea (minSize=", minSize, ", maxSize=", maxSize, ")")
fg <- fgsea(pathways = msig_h, stats = ranks, minSize = minSize, maxSize = maxSize)
fg <- fg[order(padj, -abs(NES)), ]
f_tsv <- file.path(out_dir, "fgsea_hallmark.tsv")
fwrite(fg, f_tsv, sep = "\t")
logmsg("[H] TSV ->", f_tsv, sprintf("(%d pathways)", nrow(fg)))

# Leading-edge en long
if (nrow(fg)) {
  le <- fg[, .(pathway, leadingEdge)]
  le_long <- le[, .(gene = unlist(leadingEdge)), by = pathway]
  fwrite(le_long, file.path(out_dir, "fgsea_hallmark_leading_edge.tsv"), sep = "\t")
}

# -------- Figure topN --------
if (nrow(fg)) {
  topn <- head(fg, topN)
  topn$pathway <- factor(topn$pathway, levels = rev(topn$pathway))
  p <- ggplot(topn, aes(x = pathway, y = NES, fill = padj < 0.05)) +
    geom_col() + coord_flip() +
    labs(title = paste0("Hallmarks GSEA (top ", topN, ")"),
         x = "", y = "NES", fill = "padj<0.05") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face="bold"))
  ggsave(file.path(plt_dir, "fgsea_hallmark_top.png"), p, width = 8, height = 6, dpi = 150)
}

# -------- Enrichment plots (top 5 up & 5 down) --------
if (nrow(fg)) {
  up  <- fg[NES > 0][order(padj)][1:min(5, .N)]$pathway
  dn  <- fg[NES < 0][order(padj)][1:min(5, .N)]$pathway
  sel <- c(up, dn)
  pdf(file.path(plt_dir, "fgsea_hallmark_enrichment_top10.pdf"), width = 7, height = 5)
  for (pw in sel) {
    try(print(plotEnrichment(msig_h[[pw]], ranks) + ggtitle(pw)), silent = TRUE)
  }
  dev.off()
  logmsg("[H] PDF des courbes -> fgsea_hallmark_enrichment_top10.pdf")
}

# -------- session --------
sink(file.path(out_dir, "fgsea_hallmark_sessionInfo.txt")); sessionInfo(); sink()
logmsg("✅ Terminé. Figures ->", plt_dir, "| Tables ->", out_dir)

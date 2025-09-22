#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
})

logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# -------- CLI --------
option_list <- list(
  make_option("--deg",  type="character", default="results/deseq2/deseq2_DEGs_filtered_annot.tsv",
              help="TSV des DEGs annotés (avec SYMBOL)"),
  make_option("--all",  type="character", default="results/deseq2/deseq2_results_all_annot.tsv",
              help="TSV des résultats complets annotés (univers)"),
  make_option("--out",  type="character", default="results/deseq2", help="Dossier TSV"),
  make_option("--plots",type="character", default="results/plots",  help="Dossier figures"),
  make_option("--ont",  type="character", default="BP", help="Ontologie GO: BP | MF | CC"),
  make_option("--qcut", type="double",    default=0.05, help="q-value cutoff (BH)"),
  make_option("--top",  type="integer",   default=30,   help="Nb de termes à afficher dans les plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

f_deg   <- opt$deg
f_all   <- opt$all
out_dir <- opt$out
plt_dir <- opt$plots
ont     <- toupper(opt$ont)
qcut    <- opt$qcut
topn    <- opt$top

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plt_dir, recursive = TRUE, showWarnings = FALSE)

# -------- Lecture --------
logmsg("[GO] Lecture:", f_deg, "et", f_all)
if (!file.exists(f_deg) || !file.exists(f_all))
  stop("[GO] Fichiers introuvables (deg/all).")

deg  <- fread(f_deg)
allg <- fread(f_all)

deg_symbols <- na.omit(unique(deg$SYMBOL))
bg_symbols  <- na.omit(unique(allg$SYMBOL))
logmsg("[GO] Nb DEGs (SYMBOL):", length(deg_symbols))
logmsg("[GO] Taille univers:", length(bg_symbols))
if (length(deg_symbols) < 5) stop("[GO] Trop peu de gènes pour enrichissement.")

# -------- Enrichissement --------
logmsg("[GO] enrichGO ont =", ont, "| qcut =", qcut)
ego <- enrichGO(
  gene          = deg_symbols,
  universe      = bg_symbols,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = ont,
  pAdjustMethod = "BH",
  qvalueCutoff  = qcut,
  readable      = TRUE
)

ego_df <- as.data.frame(ego)
f_go   <- file.path(out_dir, paste0("GO_", ont, "_enrichment.tsv"))
fwrite(ego_df, f_go, sep = "\t")
logmsg("[GO] TSV ->", f_go, sprintf("(%d termes)", nrow(ego_df)))

# -------- Simplify (si GOSemSim dispo) --------
ego_s <- NULL
if (nrow(ego_df) > 0 && requireNamespace("GOSemSim", quietly = TRUE)) {
  logmsg("[GO] simplify() avec GOSemSim (Wang, cutoff=0.7)")
  ego_s <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
  ego_s_df <- as.data.frame(ego_s)
  f_go_s <- file.path(out_dir, paste0("GO_", ont, "_enrichment_simplified.tsv"))
  fwrite(ego_s_df, f_go_s, sep = "\t")
  logmsg("[GO] TSV simplifié ->", f_go_s, sprintf("(%d termes)", nrow(ego_s_df)))
} else if (nrow(ego_df) > 0) {
  logmsg("[GO] GOSemSim indisponible -> pas de simplify().")
}

# -------- Plots --------
if (nrow(ego_df) > 0) {
  # Dotplot
  p_dot <- dotplot(ego, showCategory = topn,
                   title = paste0("GO ", ont, " (DEGs)")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(plt_dir, paste0("GO_", ont, "_dotplot.png")), p_dot, width = 8, height = 6, dpi = 150)

  # Barplot
  p_bar <- barplot(ego, showCategory = topn, title = paste0("GO ", ont, " (DEGs)")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(plt_dir, paste0("GO_", ont, "_barplot.png")), p_bar, width = 8, height = 6, dpi = 150)

  # Dotplot simplified (si dispo)
  if (!is.null(ego_s) && nrow(as.data.frame(ego_s)) > 0) {
    p_dot_s <- dotplot(ego_s, showCategory = topn,
                       title = paste0("GO ", ont, " (simplified)")) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    ggsave(file.path(plt_dir, paste0("GO_", ont, "_dotplot_simplified.png")),
           p_dot_s, width = 8, height = 6, dpi = 150)
  }
} else {
  logmsg("[GO] Aucun terme enrichi à q<", qcut)
}

# -------- Métadonnées run --------
params <- data.frame(
  timestamp = as.character(Sys.time()),
  deg_file  = f_deg,
  all_file  = f_all,
  ontology  = ont,
  q_cutoff  = qcut,
  top_shown = topn,
  n_deg     = length(deg_symbols),
  n_bg      = length(bg_symbols)
)
f_meta <- file.path(out_dir, paste0("GO_", ont, "_run_metadata.tsv"))
fwrite(params, f_meta, sep = "\t")

sink(file.path(out_dir, paste0("GO_", ont, "_sessionInfo.txt"))); sessionInfo(); sink()
logmsg("✅ Terminé. Figures ->", plt_dir, "| Tables ->", out_dir)

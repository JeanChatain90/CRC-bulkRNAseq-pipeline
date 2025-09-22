#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(clusterProfiler); library(enrichplot)
  library(org.Hs.eg.db); library(ggplot2); library(stringr)
})

# --- Thème blanc propre pour TOUTES les figures ---
white_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_line(color = "grey90"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

out_dir  <- "results/deseq2"
plot_dir <- "results/plots"
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

f_deg <- file.path(out_dir, "deseq2_DEGs_filtered_annot.tsv")
f_all <- file.path(out_dir, "deseq2_results_all_annot.tsv")

cat("[GO dir] Lecture:\n -", f_deg, "\n -", f_all, "\n")
deg  <- fread(f_deg)
allg <- fread(f_all)

# Univers
bg_symbols <- na.omit(unique(allg$SYMBOL))

# UP / DOWN (mêmes seuils que le volcano)
deg_up   <- deg %>% filter(!is.na(SYMBOL), log2FoldChange >=  1, padj < 0.05) %>% pull(SYMBOL) %>% unique()
deg_down <- deg %>% filter(!is.na(SYMBOL), log2FoldChange <= -1, padj < 0.05) %>% pull(SYMBOL) %>% unique()
cat(sprintf("[GO dir] UP: %d | DOWN: %d | BG: %d\n",
            length(deg_up), length(deg_down), length(bg_symbols)))

# ---- enrichGO + dot/bar (fond blanc) ----
run_go <- function(genes, ont = "BP", label = "UP", top_show = 25) {
  if (length(genes) < 5) {
    message("[GO dir] ", ont, " ", label, ": <5 gènes -> skip")
    return(NULL)
  }
  ego <- enrichGO(
    gene          = genes,
    universe      = bg_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  egodf <- as.data.frame(ego)
  if (!nrow(egodf)) {
    message("[GO dir] ", ont, " ", label, ": aucun terme significatif.")
    return(NULL)
  }

  # Wrap pour étiquettes lisibles
  ego@result$Description <- stringr::str_wrap(ego@result$Description, width = 40)

  # TSV
  fwrite(egodf,
         file.path(out_dir, sprintf("GO_%s_%s_enrichment.tsv", ont, label)),
         sep = "\t", quote = FALSE)

  # Dot & bar (thème blanc)
  p_dot <- dotplot(ego, showCategory = top_show,
                   title = sprintf("GO %s (%s DEGs)", ont, label)) + white_theme
  ggsave(file.path(plot_dir, sprintf("GO_%s_%s_dotplot.png", ont, label)),
         p_dot, width = 8, height = 6, dpi = 150)

  p_bar <- barplot(ego, showCategory = top_show,
                   title = sprintf("GO %s (%s DEGs)", ont, label)) + white_theme
  ggsave(file.path(plot_dir, sprintf("GO_%s_%s_barplot.png", ont, label)),
         p_bar, width = 8, height = 6, dpi = 150)

  invisible(ego)
}

# ---- lancer enrichissements & stocker objets ----
gos <- list()  # clé = "<ONT>_<UP/DOWN>"
for (ont in c("BP","MF","CC")) {
  message("[GO dir] Ont: ", ont)
  gos[[paste0(ont, "_UP")]]   <- run_go(deg_up,   ont, "UP",   top_show = 25)
  gos[[paste0(ont, "_DOWN")]] <- run_go(deg_down, ont, "DOWN", top_show = 25)
}

# ---- Préparer couleurs gènes = log2FC (pour cnetplot) ----
fc <- allg$log2FoldChange
names(fc) <- allg$SYMBOL
fc <- fc[!is.na(names(fc)) & !is.na(fc)]

# ---- Réseaux terme–terme (emap) et terme–gène (cnet) en fond blanc ----
if (!requireNamespace("GOSemSim", quietly = TRUE)) {
  message("[GO dir] ⚠ GOSemSim non installé -> emap/cnet sautés. Installe: bioconductor-gosemsim")
} else {
  make_network_figs <- function(ego, ont, label) {
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(invisible(NULL))

    # emapplot: similarité sémantique + thème blanc
    e2 <- tryCatch(enrichplot::pairwise_termsim(ego), error = function(e) NULL)
    if (!is.null(e2)) {
      showCat <- min(30, nrow(as.data.frame(e2)))
      p_emap <- enrichplot::emapplot(e2, showCategory = showCat) + white_theme
      ggsave(file.path(plot_dir, sprintf("GO_%s_%s_emap.png", ont, label)),
             p_emap, width = 10, height = 7, dpi = 150)
    }

    # cnetplot: liens terme–gène, gènes colorés par LFC, fond blanc
    showCat2 <- min(10, nrow(as.data.frame(ego)))
    if (showCat2 >= 1) {
      p_cnet <- enrichplot::cnetplot(
        ego, showCategory = showCat2, circular = FALSE,
        foldChange = fc, colorEdge = TRUE
      ) + white_theme
      ggsave(file.path(plot_dir, sprintf("GO_%s_%s_cnet.png", ont, label)),
             p_cnet, width = 10, height = 7, dpi = 150)
    }

    # Export long "terme ↔ gènes" (+ LFC/padj) pour les top catégories
    res_df <- as.data.frame(ego@result)
    if (nrow(res_df)) {
      long <- res_df |>
        dplyr::slice_head(n = showCat2) |>
        dplyr::select(ID, Description, p.adjust, Count, geneID) |>
        tidyr::separate_rows(geneID, sep = "/") |>
        dplyr::rename(SYMBOL = geneID) |>
        dplyr::left_join(allg[, c("SYMBOL","log2FoldChange","padj")], by = "SYMBOL")
      data.table::fwrite(
        long,
        file.path(out_dir, sprintf("GO_%s_%s_top%d_genes.tsv", ont, label, showCat2)),
        sep = "\t", quote = FALSE, na = "NA"
      )
    }
    invisible(NULL)
  }

  for (ont in c("BP","MF","CC")) {
    for (lab in c("UP","DOWN")) {
      key <- paste0(ont, "_", lab)
      make_network_figs(gos[[key]], ont, lab)
    }
  }
}

cat("[GO dir] ✅ Terminé. Figures -> results/plots | Tables -> results/deseq2\n")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

res <- fread(file.path(out_dir, "deseq2_results_all.tsv"))
ens <- sub("\\.\\d+$","", res$gene_id) # retire version ENSG.x
map <- AnnotationDbi::select(org.Hs.eg.db, keys=ens,
                             keytype="ENSEMBL", columns=c("SYMBOL","ENTREZID","GENENAME"))
map <- unique(map)
res$ENSEMBL <- ens
res_anno <- merge(res, map, by.x="ENSEMBL", by.y="ENSEMBL", all.x=TRUE)
fwrite(res_anno, file.path(out_dir, "deseq2_results_all_annot.tsv"), sep="\t", quote=FALSE, na="NA")

deg <- fread(file.path(out_dir, "deseq2_DEGs_filtered.tsv"))
ens2 <- sub("\\.\\d+$","", deg$gene_id)
map2 <- AnnotationDbi::select(org.Hs.eg.db, keys=ens2, keytype="ENSEMBL",
                              columns=c("SYMBOL","ENTREZID","GENENAME"))
deg$ENSEMBL <- ens2
deg_anno <- merge(deg, map2, by.x="ENSEMBL", by.y="ENSEMBL", all.x=TRUE)
fwrite(deg_anno, file.path(out_dir, "deseq2_DEGs_filtered_annot.tsv"), sep="\t", quote=FALSE, na="NA")
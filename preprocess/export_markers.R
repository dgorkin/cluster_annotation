#!/usr/bin/env Rscript
# Export marker RDS -> gzipped TSV + reconstruct the feature-plot page index.
#
# Usage:
#   Rscript export_markers.R <genes_rds> <motif_rds> <out_dir> [gene_topn=20] [motif_topn=30] [gene_separator=Gapdh]
#
# Outputs in <out_dir>:
#   gene_markers.tsv.gz    feature,p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster,delta_pct
#   motif_markers.tsv.gz   feature,p_val,avg_diff,pct.1,pct.2,p_val_adj,cluster,delta_pct
#   feature_plot_index.tsv kind,cluster,page,feature   (page 1 = __UMAP__ highlight)
#
# The feature-plot ordering reproduces the plotting pipeline exactly:
#   GENE : page1=UMAP; unique(c(head[order(p_val_adj,-avg_log2FC)],topn), separator, head[order(p_val_adj,-pct.1)],topn))
#   MOTIF: page1=UMAP; head[order(p_val_adj,-avg_diff)], topn   (single list, no separator)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("need: <genes_rds> <motif_rds> <out_dir> [gene_topn] [motif_topn] [gene_sep]")
genes_rds <- args[1]; motif_rds <- args[2]; out_dir <- args[3]
gene_topn  <- if (length(args) >= 4) as.integer(args[4]) else 20L
motif_topn <- if (length(args) >= 5) as.integer(args[5]) else 30L
gene_sep   <- if (length(args) >= 6) args[6] else "Gapdh"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

unwrap <- function(x) if (is.list(x) && !is.data.frame(x)) x[[1]] else x

# ---- load ----
genes <- unwrap(readRDS(genes_rds))
motifs <- unwrap(readRDS(motif_rds))

# feature column: genes live in rownames; motifs already have a `gene` column
if (!"gene" %in% colnames(genes)) genes$gene <- rownames(genes)
if (!"gene" %in% colnames(motifs)) motifs$gene <- rownames(motifs)

add_delta <- function(df) { df$delta_pct <- df[["pct.1"]] - df[["pct.2"]]; df }
genes  <- add_delta(genes)
motifs <- add_delta(motifs)

write_tsv_gz <- function(df, path, feat_cols) {
  df <- df[, feat_cols, drop = FALSE]
  colnames(df)[colnames(df) == "gene"] <- "feature"
  gz <- gzfile(path, "w")
  write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
}
write_tsv_gz(genes,  file.path(out_dir, "gene_markers.tsv.gz"),
             c("gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","delta_pct"))
write_tsv_gz(motifs, file.path(out_dir, "motif_markers.tsv.gz"),
             c("gene","p_val","avg_diff","pct.1","pct.2","p_val_adj","cluster","delta_pct"))

# ---- reconstruct feature-plot page index ----
gene_features <- function(t) {
  # the separator (e.g. Gapdh) is always plotted, even when a cluster has zero real markers
  unique(c(head(t$gene[order(t$p_val_adj, (t$avg_log2FC * -1))], gene_topn),
           gene_sep,
           head(t$gene[order(t$p_val_adj, (t[["pct.1"]] * -1))], gene_topn)))
}
motif_features <- function(t) {
  if (nrow(t) == 0) return(character(0))
  head(t$gene[order(t$p_val_adj, (t$avg_diff * -1))], motif_topn)
}

index_rows <- function(df, kind, featfun) {
  # use factor levels when present so clusters with zero marker rows still get a UMAP page
  clusters <- if (is.factor(df$cluster)) sort(as.numeric(levels(df$cluster)))
              else sort(unique(as.numeric(as.character(df$cluster))))
  out <- list()
  for (cl in clusters) {
    feats <- featfun(df[as.character(df$cluster) == as.character(cl), ])
    pages <- c("__UMAP__", feats)          # page 1 is always the UMAP highlight
    out[[length(out) + 1]] <- data.frame(
      kind = kind, cluster = cl, page = seq_along(pages), feature = pages,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

idx <- rbind(index_rows(genes,  "gene",  gene_features),
             index_rows(motifs, "motif", motif_features))
write.table(idx, file.path(out_dir, "feature_plot_index.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote gene_markers (%d rows), motif_markers (%d rows), feature_plot_index (%d rows) to %s\n",
            nrow(genes), nrow(motifs), nrow(idx), out_dir))

#!/usr/bin/env Rscript
# Count cells per cluster in a Seurat object -> cluster_cells.tsv (+ a small .json sidecar).
#
# Usage:
#   Rscript export_cell_counts.R <object> <out_dir> [cluster_column] [object_var] [expected_clusters]
#
#     object            .rds (readRDS) or .rdata/.RData (load)
#     cluster_column    metadata column holding the cluster ids; "auto" (default) detects it
#     object_var        variable name inside an .rdata; "auto" (default) finds the Seurat object
#     expected_clusters comma-separated cluster ids from the marker tables, used both to pick the
#                       right metadata column in auto mode and to verify an explicit one
#
# Outputs in <out_dir>:
#   cluster_cells.tsv        cluster, ncells, pct_cells   (pct_cells = fraction of all cells)
#   cluster_cells.meta.json  which object/column was used, the total, and when
#
# Deliberately BASE R ONLY: no Seurat, no SeuratObject. Cell counts live in the object's
# `meta.data` slot, which is a plain data.frame reachable as an attribute — so this runs under the
# system Rscript rather than needing whichever Seurat env built the object. The consequence is
# that nothing here validates the object: it reads attributes and counts rows.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("need: <object> <out_dir> [cluster_column] [object_var] [expected_clusters]")
obj_path <- args[1]
out_dir  <- args[2]
want_col <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "auto"
want_var <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "auto"
expected <- if (length(args) >= 5 && nzchar(args[5])) trimws(strsplit(args[5], ",")[[1]]) else character(0)
if (!file.exists(obj_path)) stop(sprintf("object not found: %s", obj_path))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- metadata out of an object we cannot instantiate ------------------------------------------
# `@` validates the slot against a class definition that isn't loaded here, so read the attribute.
meta_of <- function(x) {
  m <- tryCatch(attr(x, "meta.data"), error = function(e) NULL)
  if (is.data.frame(m) && nrow(m) > 0) m else NULL
}
idents_of <- function(x) tryCatch(attr(x, "active.ident"), error = function(e) NULL)

# ---- load ------------------------------------------------------------------------------------
cat(sprintf("Reading %s (%.1f GB) …\n", obj_path, file.size(obj_path) / 1024^3))
t0 <- Sys.time()
# The extension is not reliable: this project's ".rdata" files are single serialized objects
# (readRDS), not save() archives (load). Try readRDS first, since a save() archive fails it fast
# and cheaply on the magic number, then fall back to load().
obj <- tryCatch(readRDS(obj_path), error = function(e) NULL)
if (!is.null(obj)) {
  obj_var <- basename(obj_path)
} else {
  env <- new.env(parent = emptyenv())
  loaded <- load(obj_path, envir = env)
  if (want_var != "auto") {
    if (!want_var %in% loaded)
      stop(sprintf("object_var '%s' is not in %s (found: %s)",
                   want_var, basename(obj_path), paste(loaded, collapse = ", ")))
    obj_var <- want_var
  } else {
    # The Seurat object is whichever loaded variable carries a populated meta.data.
    cands <- loaded[vapply(loaded, function(n) !is.null(meta_of(get(n, envir = env))), logical(1))]
    if (length(cands) == 0)
      stop(sprintf("no object with a meta.data slot in %s (found: %s). Set cell_counts.object_var.",
                   basename(obj_path), paste(loaded, collapse = ", ")))
    if (length(cands) > 1) {
      # Ambiguity is the user's call to make — guessing risks counting cells from the wrong object.
      sizes <- vapply(cands, function(n) nrow(meta_of(get(n, envir = env))), numeric(1))
      stop(sprintf("several objects have meta.data in %s (%s). Set cell_counts.object_var to one.",
                   basename(obj_path),
                   paste(sprintf("%s: %d cells", cands, sizes), collapse = "; ")))
    }
    obj_var <- cands[[1]]
  }
  obj <- get(obj_var, envir = env)
}
cat(sprintf("Loaded '%s' in %.0f s\n", obj_var, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

meta <- meta_of(obj)
if (is.null(meta)) stop("the object has no populated meta.data slot — is it a Seurat object?")
total <- nrow(meta)

# ---- pick the cluster column -----------------------------------------------------------------
as_ids <- function(v) as.character(if (is.factor(v)) as.character(v) else v)
# Plausible cluster columns only: few enough distinct values to be a clustering, and not unique
# per cell. `active.ident` is offered under a pseudo-name so an object whose idents are the
# clustering still works when no metadata column holds it.
cols <- c(as.list(meta), list(`__active.ident__` = idents_of(obj)))
cols <- cols[vapply(cols, function(v) !is.null(v) && length(v) == total, logical(1))]
cols <- cols[vapply(cols, function(v) {
  u <- unique(as_ids(v)); length(u) > 1 && length(u) <= 2000 && length(u) < total
}, logical(1))]

matches_expected <- function(v) {
  if (!length(expected)) return(FALSE)
  setequal(unique(as_ids(v)), expected)
}

if (want_col != "auto") {
  if (!want_col %in% names(cols))
    stop(sprintf("cluster_column '%s' is not a usable metadata column. Candidates: %s",
                 want_col, paste(names(cols), collapse = ", ")))
  chosen <- want_col
  if (length(expected) && !matches_expected(cols[[chosen]]))
    cat(sprintf(paste0("WARNING: column '%s' holds %d distinct values but the marker tables have ",
                       "%d clusters. The counts may not line up with the app's clusters.\n"),
                chosen, length(unique(as_ids(cols[[chosen]]))), length(expected)))
} else {
  # Auto-detect by agreement with the marker tables: the right column is the one whose value set
  # IS the cluster set. That is self-verifying, which guessing from column names is not.
  hits <- names(cols)[vapply(cols, matches_expected, logical(1))]
  if (!length(hits)) {
    diag <- vapply(names(cols), function(n)
      sprintf("  %-30s %d distinct", n, length(unique(as_ids(cols[[n]])))), character(1))
    stop(sprintf(paste0("could not identify the cluster column: none of the %d candidate columns ",
                        "matches the %d cluster ids in the marker tables%s.\nCandidates:\n%s\n",
                        "Set cell_counts.cluster_column in the dataset config."),
                 length(cols), length(expected),
                 if (length(expected)) "" else " (none were supplied)",
                 paste(diag, collapse = "\n")))
  }
  # Several columns can agree (a resolution sweep often leaves duplicates); they give identical
  # counts by definition of the match, so prefer the conventional name for the record.
  pref <- c(grep("seurat_clusters", hits, value = TRUE),
            grep("res", hits, value = TRUE), grep("clust", hits, value = TRUE), hits)
  chosen <- pref[[1]]
  if (length(hits) > 1)
    cat(sprintf("Cluster column: %d columns match the marker clusters (%s); using '%s'.\n",
                length(hits), paste(hits, collapse = ", "), chosen))
}
ids <- as_ids(cols[[chosen]])
cat(sprintf("Cluster column: %s · %d cells · %d clusters\n",
            chosen, total, length(unique(ids))))

# ---- write -----------------------------------------------------------------------------------
tab <- table(ids)
num <- suppressWarnings(as.numeric(names(tab)))
ord <- if (any(is.na(num))) order(names(tab)) else order(num)   # numeric ids sort numerically
out <- data.frame(cluster = names(tab)[ord],
                  ncells = as.integer(tab)[ord],
                  pct_cells = as.integer(tab)[ord] / total,
                  stringsAsFactors = FALSE)
write.table(out, file.path(out_dir, "cluster_cells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

json_str <- function(x) paste0('"', gsub('"', '\\\\"', x), '"')
writeLines(paste0("{\n",
  '  "object": ',         json_str(normalizePath(obj_path)),           ",\n",
  '  "object_var": ',     json_str(obj_var),                           ",\n",
  '  "cluster_column": ', json_str(chosen),                            ",\n",
  '  "total_cells": ',    total,                                       ",\n",
  '  "n_clusters": ',     length(unique(ids)),                         ",\n",
  '  "generated_at": ',   json_str(format(Sys.time(), "%Y-%m-%dT%H:%M:%S")), "\n}"),
  file.path(out_dir, "cluster_cells.meta.json"))

cat(sprintf("Wrote cluster_cells.tsv (%d clusters, %d cells) to %s\n", nrow(out), total, out_dir))

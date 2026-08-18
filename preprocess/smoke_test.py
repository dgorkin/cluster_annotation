"""Headless smoke test for the data/pdf layer (does not boot Streamlit's server)."""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "app"))
import data as D
import pdf as P

cfg = D.load_config("config/mydataset.yaml")
D.ensure_preprocessed(cfg)
g = D.load_markers(cfg, "gene")
m = D.load_markers(cfg, "motif")
print("genes", g.shape, "cols", list(g.columns))
print("motifs", m.shape, "cols", list(m.columns))
print("clusters", len(D.clusters(cfg)))
tm = D.tangram_index(cfg)
print("tangram mapped", len(tm), "| c0->", tm.get(0), "c10->", tm.get(10), "c2->", tm.get(2))
top = D.top_markers(g, 0, 5, "avg_log2FC", False)
print("top5 c0 genes by log2FC:", list(top.feature))
fp = D.feature_pages(cfg, "gene", 0)
print("gene c0 pages", len(fp), "| p2 =", fp[fp.page == 2].feature.iloc[0])
png = P.render_page(D.featureplot_pdf(cfg, "gene", 0), 1, 120, D.cache_dir(cfg) / "img")
print("rendered UMAP png:", os.path.exists(png), os.path.getsize(png), "bytes")
print("DATA_LAYER_OK")

"""PDF helpers: render a single page to a cached PNG, count pages, read page text."""
from __future__ import annotations
import hashlib
import os
from pathlib import Path

import fitz  # PyMuPDF


def page_count(pdf_path: str | Path) -> int:
    with fitz.open(pdf_path) as doc:
        return doc.page_count


def page_texts(pdf_path: str | Path) -> list[str]:
    out = []
    with fitz.open(pdf_path) as doc:
        for pg in doc:
            out.append(" ".join(pg.get_text().split()))
    return out


def render_page(pdf_path: str | Path, page_1based: int, dpi: int, cache_dir: str | Path) -> str:
    """Render one page (1-based) to a PNG, cached by (path, mtime, page, dpi). Returns PNG path."""
    pdf_path = Path(pdf_path)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    mtime = int(os.path.getmtime(pdf_path))
    key = f"{pdf_path.resolve()}|{mtime}|{page_1based}|{dpi}"
    h = hashlib.md5(key.encode()).hexdigest()[:16]
    out = cache_dir / f"{pdf_path.stem}_p{page_1based}_d{dpi}_{h}.png"
    if out.exists():
        return str(out)
    with fitz.open(pdf_path) as doc:
        page = doc[page_1based - 1]
        pix = page.get_pixmap(dpi=dpi)
        pix.save(out)
    return str(out)

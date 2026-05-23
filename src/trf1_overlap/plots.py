"""Matplotlib bar plots summarising matched gene counts."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

# matplotlib needs writable cache dirs; redirect to a tmp location so the
# package works in sandboxed environments where $HOME is read-only.
_CACHE_DIR = Path(tempfile.gettempdir()) / "trf1_promoter_overlap_cache"
_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from .promoter import PROMOTER_SIZES


def plot_bar(summaries: dict[int, dict[str, int]], outdir: Path) -> Path:
    sizes = list(PROMOTER_SIZES)
    counts = [summaries[size]["matched genes"] for size in sizes]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar([str(size) for size in sizes], counts, color="#4C78A8")
    ax.set_xlabel("Upstream window size (bp)")
    ax.set_ylabel("TRF1 matched genes count")
    ax.set_title("TRF1 promoter-upstream overlap")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path = outdir / "TRF1_matched_genes_barplot.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def plot_stacked_bar(summaries: dict[int, dict[str, int]], outdir: Path) -> Path:
    sizes = list(PROMOTER_SIZES)
    up_counts = [summaries[size]["UP matched"] for size in sizes]
    down_counts = [summaries[size]["DOWN matched"] for size in sizes]

    fig, ax = plt.subplots(figsize=(6, 4))
    labels = [str(size) for size in sizes]
    ax.bar(labels, up_counts, color="#D55E00", label="UP")
    ax.bar(labels, down_counts, bottom=up_counts, color="#0072B2", label="DOWN")
    ax.set_xlabel("Upstream window size (bp)")
    ax.set_ylabel("TRF1 matched genes count")
    ax.set_title("TRF1 promoter-upstream overlap by regulation")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path = outdir / "TRF1_matched_genes_stacked_barplot.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path

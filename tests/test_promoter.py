"""Promoter-window math: TSS placement and BED half-open intervals."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from trf1_overlap.promoter import build_promoters


def _make_target(strand, gene_start, gene_end):
    return pd.DataFrame(
        [
            {
                "gene_name": "G1",
                "gene_id": "ENSG0",
                "regulation": "UP",
                "chr": "chr1",
                "gene_start": gene_start,
                "gene_end": gene_end,
                "strand": strand,
                "_deg_key": "ENSG0",
            }
        ]
    )


def test_plus_strand_window_ends_at_tss():
    df = build_promoters(_make_target("+", 10_000, 12_000), size=1000)
    row = df.iloc[0]
    assert row["TSS"] == 10_000
    assert row["promoter_start"] == 9_000
    assert row["promoter_end"] == 10_000


def test_minus_strand_window_starts_at_tss():
    df = build_promoters(_make_target("-", 5_000, 6_000), size=1000)
    row = df.iloc[0]
    # TSS for minus strand is labelled as gene_end (BED-exclusive)
    assert row["TSS"] == 6_000
    assert row["promoter_start"] == 6_000
    assert row["promoter_end"] == 7_000


def test_plus_strand_clamps_at_zero():
    df = build_promoters(_make_target("+", 500, 1_500), size=1000)
    row = df.iloc[0]
    assert row["promoter_start"] == 0
    assert row["promoter_end"] == 500


def test_mt_nd1_window_matches_known_case():
    # MT-ND1: BED chrM:3306-4262 +; 1000bp upstream window must be [2306, 3306)
    df = build_promoters(_make_target("+", 3306, 4262), size=1000)
    row = df.iloc[0]
    assert (row["promoter_start"], row["promoter_end"]) == (2306, 3306)

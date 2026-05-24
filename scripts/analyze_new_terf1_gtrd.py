#!/usr/bin/env python3
"""Analyze newer GTRD/MSigDB TERF1 target genes against siTRF1 DEGs."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
GENESET_JSON = ROOT / "results/new_data/TERF1_TARGET_GENES.json"
DEG_REGULATION = ROOT / "data/siTRF1_DEG_list_gene_name.txt"
BETA_BSF = ROOT / "results/integration_analysis/beta/siTRF1_DEGs.bsf"
ANNOTATION = ROOT / "data/gene_annotation.tsv"
OUTDIR = ROOT / "results/new_data"


def hypergeom_sf(k: int, population: int, success: int, draws: int) -> float:
    """P[X >= k] for Hypergeometric(N=population, K=success, n=draws)."""
    if k <= 0:
        return 1.0
    upper = min(success, draws)
    denom = math.comb(population, draws)
    return sum(
        math.comb(success, i) * math.comb(population - success, draws - i)
        for i in range(k, upper + 1)
    ) / denom


def load_terf1_targets() -> tuple[set[str], dict[str, object]]:
    payload = json.loads(GENESET_JSON.read_text())
    record = payload["TERF1_TARGET_GENES"]
    genes = {str(g).strip() for g in record["geneSymbols"] if str(g).strip()}
    return genes, record


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    targets, record = load_terf1_targets()

    deg = pd.read_csv(DEG_REGULATION, sep="\t", dtype=str)
    deg["gene_name"] = deg["gene_name"].fillna("").str.strip()
    deg["regulation"] = deg["regulation"].fillna("").str.strip().str.upper()
    deg = deg.drop_duplicates(["gene_name", "regulation"])

    bsf = pd.read_csv(
        BETA_BSF,
        sep="\t",
        header=None,
        names=["gene_name", "log2FC", "pvalue"],
        dtype={"gene_name": str, "log2FC": float, "pvalue": float},
    )
    bsf = bsf.drop_duplicates("gene_name")
    deg = deg.merge(bsf, how="left", on="gene_name")

    annotation = pd.read_csv(ANNOTATION, sep="\t", dtype=str)
    universe = {
        str(g).strip()
        for g in annotation["gene_name"].dropna()
        if str(g).strip() and str(g).strip() != "nan"
    }
    deg_genes = set(deg["gene_name"])
    targets_in_universe = targets & universe
    deg_in_universe = deg_genes & universe
    overlap = targets & deg_genes
    overlap_in_universe = targets_in_universe & deg_in_universe

    overlap_df = (
        deg[deg["gene_name"].isin(overlap)]
        .sort_values(["regulation", "gene_name"])
        .reset_index(drop=True)
    )
    overlap_df.to_csv(OUTDIR / "terf1_gtrd_deg_overlap.tsv", sep="\t", index=False)

    summary_rows = [
        ("source", "MSigDB TERF1_TARGET_GENES"),
        ("systematic_name", record.get("systematicName", "")),
        ("pmid", record.get("pmid", "")),
        ("exact_source", record.get("exactSource", "")),
        ("collection", record.get("collection", "")),
        ("msigdb_url", record.get("msigdbURL", "")),
        ("external_url", ";".join(record.get("externalDetailsURL", []))),
        ("terf1_target_gene_symbols", len(targets)),
        ("deg_gene_symbols", len(deg_genes)),
        ("annotation_universe_gene_symbols", len(universe)),
        ("targets_in_annotation_universe", len(targets_in_universe)),
        ("degs_in_annotation_universe", len(deg_in_universe)),
        ("overlap_targets_vs_degs", len(overlap)),
        ("overlap_targets_vs_degs_in_universe", len(overlap_in_universe)),
        (
            "hypergeom_p_annotation_universe",
            hypergeom_sf(
                len(overlap_in_universe),
                len(universe),
                len(targets_in_universe),
                len(deg_in_universe),
            ),
        ),
        (
            "overlap_genes",
            ",".join(sorted(overlap)),
        ),
        (
            "overlap_up",
            ",".join(sorted(overlap_df.loc[overlap_df["regulation"] == "UP", "gene_name"])),
        ),
        (
            "overlap_down",
            ",".join(
                sorted(overlap_df.loc[overlap_df["regulation"] == "DOWN", "gene_name"])
            ),
        ),
    ]
    summary = pd.DataFrame(summary_rows, columns=["metric", "value"])
    summary.to_csv(OUTDIR / "terf1_gtrd_deg_overlap_summary.tsv", sep="\t", index=False)

    lines = [
        "# Newer TERF1 target gene set check",
        "",
        "Input: MSigDB `TERF1_TARGET_GENES`, derived from GTRD ChIP-seq harmonization.",
        "",
        f"- TERF1 target gene symbols: {len(targets)}",
        f"- siTRF1 DEG symbols: {len(deg_genes)}",
        f"- Overlap: {len(overlap)}",
        f"- Overlap genes: {', '.join(sorted(overlap)) if overlap else 'None'}",
        f"- Annotation-universe hypergeometric p-value: {summary_rows[14][1]:.6g}",
        "",
        "Interpretation rule: this is a harmonized public ChIP-seq target-gene set,",
        "not matched to the original siTRF1 experimental context. Treat overlap as",
        "independent database cross-check, not direct perturbation evidence.",
        "",
    ]
    (OUTDIR / "terf1_gtrd_deg_overlap_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

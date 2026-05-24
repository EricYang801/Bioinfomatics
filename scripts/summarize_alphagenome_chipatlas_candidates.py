#!/usr/bin/env python3
"""Create a compact per-gene AlphaGenome summary for ChIP-Atlas candidates."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/alphagenome/chipatlas_candidates"
TRACKS = OUTDIR / "alphagenome_bj_chipatlas_candidate_track_summary.tsv"
CHIPATLAS = ROOT / "results/new_data/chipatlas/chipatlas_terf1_deg_overlap.tsv"


def max_signal(table: pd.DataFrame, output_attr: str, **filters: str) -> pd.DataFrame:
    subset = table[table["output_attr"].eq(output_attr)].copy()
    for key, value in filters.items():
        subset = subset[subset[key].fillna("").eq(value)]
    if subset.empty:
        return pd.DataFrame(columns=["gene_name"])
    label = output_attr
    if "histone_mark" in filters:
        label = filters["histone_mark"].lower()
    if "transcription_factor" in filters:
        label = filters["transcription_factor"].lower()
    return (
        subset.groupby("gene_name", dropna=False)
        .agg(
            **{
                f"{label}_max_signal": ("max_signal", "max"),
                f"{label}_max_tss_bin_signal": ("tss_bin_signal", "max"),
            }
        )
        .reset_index()
    )


def main() -> None:
    tracks = pd.read_csv(TRACKS, sep="\t")
    chipatlas = pd.read_csv(CHIPATLAS, sep="\t")
    chipatlas = chipatlas[chipatlas["window"].eq("10kb")].copy()
    chipatlas = chipatlas[
        [
            "gene_name",
            "regulation",
            "log2FC",
            "pvalue",
            "average_score",
            "max_chipseq_score",
            "chipseq_experiments_with_signal",
        ]
    ].drop_duplicates("gene_name")

    summary = chipatlas.copy()
    for collapsed in [
        max_signal(tracks, "rna_seq"),
        max_signal(tracks, "dnase"),
        max_signal(tracks, "chip_histone", histone_mark="H3K4me3"),
        max_signal(tracks, "chip_histone", histone_mark="H3K36me3"),
        max_signal(tracks, "chip_histone", histone_mark="H3K27me3"),
        max_signal(tracks, "chip_tf", transcription_factor="CTCF"),
    ]:
        summary = summary.merge(collapsed, how="left", on="gene_name")

    summary = summary.sort_values(
        ["average_score", "dnase_max_tss_bin_signal", "h3k4me3_max_tss_bin_signal"],
        ascending=[False, False, False],
    )
    summary.to_csv(
        OUTDIR / "alphagenome_bj_chipatlas_candidate_gene_summary.tsv",
        sep="\t",
        index=False,
    )

    ranked_sections = [
        ("Top DNase at TSS", "dnase_max_tss_bin_signal"),
        ("Top H3K4me3 at TSS", "h3k4me3_max_tss_bin_signal"),
        ("Top RNA max signal", "rna_seq_max_signal"),
        ("Top CTCF at TSS", "ctcf_max_tss_bin_signal"),
    ]
    lines = [
        "# AlphaGenome BJ candidate gene compact summary",
        "",
        "Input: AlphaGenome BJ-context outputs for ChIP-Atlas TERF1-overlap DEGs.",
        "All values are model output aggregations over 131,072 bp TSS-centered hg38 intervals.",
        "",
        f"- Genes summarized: {summary['gene_name'].nunique()}",
        f"- Downregulated candidates: {(summary['regulation'] == 'DOWN').sum()}",
        f"- Upregulated candidates: {(summary['regulation'] == 'UP').sum()}",
        "",
    ]
    for title, column in ranked_sections:
        ranked = summary.sort_values(column, ascending=False).head(10)
        lines.append(f"## {title}")
        for _, row in ranked.iterrows():
            lines.append(
                f"- {row['gene_name']} ({row['regulation']}, log2FC={row['log2FC']:.3g}): "
                f"{column}={row[column]:.6g}; "
                f"ChIP-Atlas average_score={row['average_score']:.6g}"
            )
        lines.append("")
    lines.extend(
        [
            "Interpretation rule: these ranks only prioritize candidate loci for",
            "manual follow-up. They are not direct TERF1 binding evidence because",
            "the AlphaGenome metadata probe found no TERF1/TRF1 ChIP-TF track.",
            "",
        ]
    )
    (OUTDIR / "alphagenome_bj_chipatlas_candidate_gene_summary.md").write_text(
        "\n".join(lines)
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()

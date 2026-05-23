"""Result formatting, summary aggregation, and text report writing."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

from .promoter import PROMOTER_SIZES


OUTPUT_COLUMNS = [
    "gene_name",
    "gene_id",
    "regulation",
    "chr",
    "TSS",
    "strand",
    "promoter_start",
    "promoter_end",
    "peak_chr",
    "peak_start",
    "peak_end",
    "peak_fold",
    "peak_pvalue",
]


def format_overlap(overlap: pd.DataFrame) -> pd.DataFrame:
    if overlap.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    output = overlap[OUTPUT_COLUMNS + ["_deg_key"]]
    output = output.sort_values(
        ["gene_name", "gene_id", "chr", "promoter_start", "peak_start"]
    )
    return output.reset_index(drop=True)


def summarize_overlap(
    overlap: pd.DataFrame, target_genes: pd.DataFrame
) -> tuple[dict[str, int], pd.DataFrame]:
    total_genes = target_genes["_deg_key"].nunique()
    if overlap.empty:
        matched = pd.DataFrame(columns=["gene_name", "gene_id", "regulation", "_deg_key"])
    else:
        matched = (
            overlap[["gene_name", "gene_id", "regulation", "_deg_key"]]
            .drop_duplicates("_deg_key")
            .sort_values(["regulation", "gene_name", "gene_id"])
            .reset_index(drop=True)
        )

    summary = {
        "total genes": int(total_genes),
        "matched genes": int(matched["_deg_key"].nunique()),
        "UP matched": int((matched["regulation"] == "UP").sum()) if not matched.empty else 0,
        "DOWN matched": int((matched["regulation"] == "DOWN").sum()) if not matched.empty else 0,
    }
    return summary, matched


def write_summary_tsv(summaries: dict[int, dict[str, int]], outdir: Path) -> Path:
    rows = []
    for size in PROMOTER_SIZES:
        row = {"upstream_window_bp": size}
        row.update(summaries[size])
        rows.append(row)
    path = outdir / "TRF1_promoter_overlap_summary.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def write_report(
    outdir: Path,
    args: argparse.Namespace,
    peak_path: Path,
    cleaned_peak_path: Path,
    annotation_source: str,
    messages: list[str],
    peaks: pd.DataFrame,
    deg: pd.DataFrame,
    annotation: pd.DataFrame,
    target_genes: pd.DataFrame,
    summaries: dict[int, dict[str, int]],
    matched_genes: dict[int, pd.DataFrame],
    plot_paths: list[Path],
    intersect_mode: str,
) -> Path:
    report_path = outdir / "report.txt"
    missing_annotation_keys = sorted(set(deg["_deg_key"]) - set(target_genes["_deg_key"]))
    deg_lookup = deg.set_index("_deg_key")[["gene_name", "gene_id", "regulation"]].to_dict(
        orient="index"
    )

    lines: list[str] = []
    lines.append("TRF1 promoter overlap analysis")
    lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("Command")
    command = (
        "trf1-promoter-overlap "
        f"--peak {args.peak} --deg {args.deg} --annotation {args.annotation} "
        f"--outdir {args.outdir}"
    )
    if args.gtf:
        command += f" --gtf {args.gtf}"
    if args.exclude_chroms:
        command += f" --exclude-chroms {args.exclude_chroms}"
    lines.append(command)
    lines.append("")
    lines.append("Inputs")
    lines.append(f"Peak file used: {peak_path}")
    lines.append(f"DEG file: {args.deg}")
    lines.append(f"Annotation source: {annotation_source}")
    lines.append(f"Intersection engine: {intersect_mode}")
    lines.append("")
    lines.append("Window interpretation")
    lines.append(
        "The 1000bp, 2000bp, and 3000bp analyses are nested upstream "
        "search windows from each gene TSS, not three independent promoters."
    )
    lines.append(
        "A gene found in the 1000bp window is therefore also expected to "
        "remain matched in the 2000bp and 3000bp windows."
    )
    lines.append("")
    lines.append("Generated files")
    lines.append(f"Cleaned peak BED: {cleaned_peak_path}")
    for size in PROMOTER_SIZES:
        lines.append(f"Overlap {size}bp: {outdir / f'TRF1_promoter_overlap_{size}.tsv'}")
    for path in plot_paths:
        lines.append(f"Plot: {path}")
    lines.append("")

    if messages:
        lines.append("Notes")
        for message in messages:
            lines.append(f"- {message}")
        lines.append("")

    lines.append("Data counts")
    lines.append(f"TRF1 peaks parsed: {len(peaks)}")
    lines.append(f"DEG genes loaded: {deg['_deg_key'].nunique()}")
    lines.append(f"Annotation genes loaded: {len(annotation)}")
    lines.append(f"DEG genes with annotation: {target_genes['_deg_key'].nunique()}")
    lines.append(f"DEG genes missing annotation: {len(missing_annotation_keys)}")
    lines.append("")

    lines.append("Summary")
    for size in PROMOTER_SIZES:
        lines.append(f"{size}bp:")
        for key in ["total genes", "matched genes", "UP matched", "DOWN matched"]:
            lines.append(f"  {key}: {summaries[size][key]}")
        lines.append("")

    for size in PROMOTER_SIZES:
        lines.append(f"Matched genes ({size}bp)")
        matched = matched_genes[size]
        if matched.empty:
            lines.append("None")
        else:
            lines.append("gene_name\tgene_id\tregulation")
            for _, row in matched.iterrows():
                lines.append(f"{row['gene_name']}\t{row['gene_id']}\t{row['regulation']}")
        lines.append("")

    if missing_annotation_keys:
        lines.append("DEG genes missing annotation")
        lines.append("gene_name\tgene_id\tregulation")
        for key in missing_annotation_keys:
            row = deg_lookup[key]
            lines.append(f"{row['gene_name']}\t{row['gene_id']}\t{row['regulation']}")
        lines.append("")

    report_path.write_text("\n".join(lines) + "\n")
    return report_path

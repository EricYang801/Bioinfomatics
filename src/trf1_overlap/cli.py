"""Command-line entry point for the TRF1 promoter overlap pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .annotation import read_annotation
from .chroms import parse_chrom_list
from .deg import match_deg_to_annotation, read_deg
from .intersect import build_peak_bedtool, intersect_with_pybedtools, require_pybedtools
from .peaks import read_peak, write_cleaned_peak
from .plots import plot_bar, plot_stacked_bar
from .promoter import PROMOTER_SIZES, build_promoters, exclude_chromosomes
from .report import (
    OUTPUT_COLUMNS,
    format_overlap,
    summarize_overlap,
    write_report,
    write_summary_tsv,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find TRF1 ChIP-seq peak overlaps with DEG promoter regions "
            "defined from gene TSS."
        )
    )
    parser.add_argument("--peak", required=True, help="TRF1 peak file (.bed or .bed.gz)")
    parser.add_argument(
        "--deg", required=True, help="DEG gene list with gene_id/gene_name/regulation"
    )
    parser.add_argument(
        "--annotation",
        required=True,
        help=(
            "Gene annotation TSV with chr/gene_start/gene_end/gene_name/gene_id/strand. "
            "If this file is absent, provide --gtf to create it from gene features."
        ),
    )
    parser.add_argument(
        "--gtf",
        help=(
            "Optional GTF used to create --annotation when the annotation TSV is absent. "
            "GTF starts are converted from 1-based closed to BED-style 0-based starts."
        ),
    )
    parser.add_argument(
        "--exclude-chroms",
        default="",
        help=(
            "Optional comma-separated chromosome names to exclude after chromosome "
            "normalization, for example: chrM,chrY."
        ),
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    peaks, peak_messages, peak_path = read_peak(args.peak)
    deg = read_deg(args.deg)
    annotation, annotation_messages, annotation_source = read_annotation(
        args.annotation, args.gtf
    )
    excluded_chroms = parse_chrom_list(args.exclude_chroms)
    peaks, annotation, exclude_message = exclude_chromosomes(peaks, annotation, excluded_chroms)
    target_genes = match_deg_to_annotation(deg, annotation)
    if target_genes.empty:
        raise ValueError("No DEG genes matched the annotation by gene_id or gene_name.")

    cleaned_peak_path = write_cleaned_peak(peaks, outdir)
    intersect_mode = require_pybedtools()
    peak_bt = build_peak_bedtool(peaks)

    summaries: dict[int, dict[str, int]] = {}
    matched_genes: dict[int, "pd.DataFrame"] = {}  # noqa: F821 (lazy import)

    for size in PROMOTER_SIZES:
        promoters = build_promoters(target_genes, size)
        raw_overlap = intersect_with_pybedtools(promoters, peak_bt)

        overlap = format_overlap(raw_overlap)
        output_path = outdir / f"TRF1_promoter_overlap_{size}.tsv"
        overlap[OUTPUT_COLUMNS].to_csv(output_path, sep="\t", index=False)

        summaries[size], matched_genes[size] = summarize_overlap(raw_overlap, target_genes)

    summary_path = write_summary_tsv(summaries, outdir)
    plot_paths = [plot_bar(summaries, outdir), plot_stacked_bar(summaries, outdir)]
    report_path = write_report(
        outdir=outdir,
        args=args,
        peak_path=peak_path,
        cleaned_peak_path=cleaned_peak_path,
        annotation_source=annotation_source,
        messages=(
            peak_messages
            + annotation_messages
            + ([exclude_message] if exclude_message else [])
            + [f"Summary TSV: {summary_path}"]
        ),
        peaks=peaks,
        deg=deg,
        annotation=annotation,
        target_genes=target_genes,
        summaries=summaries,
        matched_genes=matched_genes,
        plot_paths=plot_paths,
        intersect_mode=intersect_mode,
    )

    print(f"Analysis complete. Report: {report_path}")
    for size in PROMOTER_SIZES:
        print(
            f"{size}bp: matched genes={summaries[size]['matched genes']} "
            f"(UP={summaries[size]['UP matched']}, DOWN={summaries[size]['DOWN matched']})"
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

#!/usr/bin/env python3
"""Extended QC/audit summaries for legacy TRF1 peaks and existing analyses."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from trf1_overlap.peaks import read_peak


ROOT = Path(__file__).resolve().parents[1]
PEAKS = ROOT / "data/GSM638201_TRF1-p0.05.bed.gz"
PROMOTER_SUMMARY = ROOT / "results/TRF1_promoter_overlap_summary.tsv"
MOTIF_SCAN = ROOT / "results/integration_analysis/motif_scan/ttaggg_motif_scan_full.tsv"
BETA_NEAREST = ROOT / "results/integration_analysis/beta/deg_nearest_peak.tsv"
OUTDIR = ROOT / "results/extended_audit"


def simple_table(df: pd.DataFrame) -> str:
    return df.to_csv(sep="\t", index=False).strip()


def numeric_peak_table() -> pd.DataFrame:
    peak, _, _ = read_peak(str(PEAKS))
    peak["peak_fold"] = pd.to_numeric(peak["peak_fold"], errors="coerce")
    peak["peak_pvalue"] = pd.to_numeric(peak["peak_pvalue"], errors="coerce")
    peak["peak_length"] = peak["peak_end"] - peak["peak_start"]
    peak["is_chrM"] = peak["peak_chr"].eq("chrM")
    return peak


def threshold_table(peak: pd.DataFrame) -> pd.DataFrame:
    rows = []
    p_thresholds = [0.05, 0.01, 0.001, 3.9e-4, 1e-4]
    fold_thresholds = [1, 2, 3, 5, 10, 20, 50]
    for p_cut in p_thresholds:
        for f_cut in fold_thresholds:
            subset = peak[(peak["peak_pvalue"] <= p_cut) & (peak["peak_fold"] >= f_cut)]
            rows.append(
                {
                    "pvalue_lte": p_cut,
                    "fold_gte": f_cut,
                    "peaks": len(subset),
                    "nuclear_peaks": int((~subset["is_chrM"]).sum()),
                    "chrM_peaks": int(subset["is_chrM"].sum()),
                    "median_fold": subset["peak_fold"].median(),
                    "max_fold": subset["peak_fold"].max(),
                }
            )
    return pd.DataFrame(rows)


def nearest_summary() -> tuple[pd.DataFrame, pd.DataFrame]:
    nearest = pd.read_csv(BETA_NEAREST, sep="\t")
    nearest["abs_distance"] = pd.to_numeric(
        nearest["distance(peak-TSS)"], errors="coerce"
    ).abs()
    bins = [0, 1_000, 3_000, 10_000, 50_000, 100_000, 500_000, 1_000_000, float("inf")]
    labels = ["<=1kb", "1-3kb", "3-10kb", "10-50kb", "50-100kb", "100-500kb", "500kb-1Mb", ">1Mb"]
    nearest["distance_bin"] = pd.cut(
        nearest["abs_distance"], bins=bins, labels=labels, include_lowest=True
    )
    summary = (
        nearest.groupby(["distance_bin", "direction"], observed=False)
        .size()
        .reset_index(name="n_degs")
    )
    return nearest, summary


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    peak = numeric_peak_table()
    peak.to_csv(OUTDIR / "legacy_trf1_peaks_parsed.tsv", sep="\t", index=False)

    by_chrom = (
        peak.groupby("peak_chr")
        .agg(
            peaks=("peak_chr", "size"),
            median_fold=("peak_fold", "median"),
            max_fold=("peak_fold", "max"),
            min_pvalue=("peak_pvalue", "min"),
            median_length=("peak_length", "median"),
        )
        .sort_values(["peaks", "max_fold"], ascending=False)
        .reset_index()
    )
    by_chrom.to_csv(OUTDIR / "legacy_trf1_peaks_by_chrom.tsv", sep="\t", index=False)

    thresh = threshold_table(peak)
    thresh.to_csv(OUTDIR / "legacy_trf1_peak_threshold_grid.tsv", sep="\t", index=False)

    top_fold = peak.sort_values("peak_fold", ascending=False).head(30)
    top_fold.to_csv(OUTDIR / "legacy_trf1_top_fold_peaks.tsv", sep="\t", index=False)

    promoter_summary = pd.read_csv(PROMOTER_SUMMARY, sep="\t")
    promoter_summary.to_csv(
        OUTDIR / "existing_promoter_overlap_summary.copy.tsv", sep="\t", index=False
    )

    motif = pd.read_csv(MOTIF_SCAN, sep="\t")
    motif_summary = pd.DataFrame(
        [
            {"metric": "motif_scan_degs", "value": len(motif)},
            {
                "metric": "degs_with_single_TTAGGG",
                "value": int((motif["n_TTAGGG_x1"] >= 1).sum()),
            },
            {
                "metric": "degs_with_tandem_TTAGGG_x2",
                "value": int((motif["n_TTAGGG_x2"] >= 1).sum()),
            },
            {
                "metric": "degs_with_tandem_TTAGGG_x3",
                "value": int((motif["n_TTAGGG_x3"] >= 1).sum()),
            },
        ]
    )
    motif_summary.to_csv(OUTDIR / "existing_ttaggg_motif_summary.tsv", sep="\t", index=False)

    nearest, nearest_bins = nearest_summary()
    nearest.to_csv(OUTDIR / "existing_beta_nearest_peak_with_bins.tsv", sep="\t", index=False)
    nearest_bins.to_csv(OUTDIR / "existing_beta_nearest_peak_bins.tsv", sep="\t", index=False)

    chrM = peak[peak["is_chrM"]]
    nuclear = peak[~peak["is_chrM"]]
    lines = [
        "# TRF1 extended audit",
        "",
        "Scope: legacy GSM638201 TRF1 peak QC and existing analysis summaries.",
        "No hg18-to-hg38 coordinate conversion is performed here.",
        "",
        f"- Parsed peaks: {len(peak)}",
        f"- Nuclear peaks: {len(nuclear)}",
        f"- chrM peaks: {len(chrM)}",
        f"- Median fold, all peaks: {peak['peak_fold'].median():.3g}",
        f"- Median fold, chrM peaks: {chrM['peak_fold'].median():.3g}" if len(chrM) else "- Median fold, chrM peaks: NA",
        f"- Median fold, nuclear peaks: {nuclear['peak_fold'].median():.3g}",
        f"- Max fold, all peaks: {peak['peak_fold'].max():.3g}",
        "",
        "Existing promoter-overlap summary:",
        simple_table(promoter_summary),
        "",
        "Motif scan summary:",
        simple_table(motif_summary),
        "",
        "Nearest peak bins from existing hg18/BETA run:",
        simple_table(nearest_bins),
        "",
    ]
    (OUTDIR / "trf1_extended_audit_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

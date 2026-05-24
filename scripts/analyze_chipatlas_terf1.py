#!/usr/bin/env python3
"""Analyze ChIP-Atlas hg38 TERF1 target genes against siTRF1 DEGs."""

from __future__ import annotations

import math
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CHIPATLAS_DIR = ROOT / "results/new_data/chipatlas"
DEG_REGULATION = ROOT / "data/siTRF1_DEG_list_gene_name.txt"
BETA_BSF = ROOT / "results/integration_analysis/beta/siTRF1_DEGs.bsf"
ANNOTATION = ROOT / "data/gene_annotation.tsv"
OUTDIR = ROOT / "results/new_data/chipatlas"


WINDOWS = {
    "1kb": CHIPATLAS_DIR / "TERF1.1.tsv",
    "5kb": CHIPATLAS_DIR / "TERF1.5.tsv",
    "10kb": CHIPATLAS_DIR / "TERF1.10.tsv",
}

EXPERIMENT_METADATA = [
    CHIPATLAS_DIR / "SRX4957970.metadata.json",
    CHIPATLAS_DIR / "SRX359962.metadata.json",
    CHIPATLAS_DIR / "SRX472275.metadata.json",
    CHIPATLAS_DIR / "SRX4957969.metadata.json",
]


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


def load_deg() -> pd.DataFrame:
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
    return deg.merge(bsf, how="left", on="gene_name")


def load_universe() -> set[str]:
    annotation = pd.read_csv(ANNOTATION, sep="\t", dtype=str)
    return {
        str(g).strip()
        for g in annotation["gene_name"].dropna()
        if str(g).strip() and str(g).strip() != "nan"
    }


def load_chipatlas(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    table = table.rename(columns={"Target_genes": "gene_name", "TERF1|Average": "average_score"})
    sample_cols = [c for c in table.columns if c.startswith("SRX")]
    for col in ["average_score", *sample_cols]:
        table[col] = pd.to_numeric(table[col], errors="coerce").fillna(0.0)
    table["max_chipseq_score"] = table[sample_cols].max(axis=1)
    table["n_chipseq_experiments_with_signal"] = (table[sample_cols] > 0).sum(axis=1)
    table["chipseq_experiments_with_signal"] = table.apply(
        lambda row: ",".join(col for col in sample_cols if row[col] > 0),
        axis=1,
    )
    # The downloaded target-gene tables contain a STRING-only row. For direct
    # TRF1 binding checks we keep only rows supported by ChIP-seq sample scores.
    return table[table["max_chipseq_score"] > 0].copy()


def summarize_experiment_metadata() -> pd.DataFrame:
    rows = []
    for path in EXPERIMENT_METADATA:
        if not path.exists():
            continue
        records = json.loads(path.read_text())
        hg38_records = [r for r in records if r.get("genome") == "hg38"]
        for record in hg38_records:
            rows.append(
                {
                    "expid": record.get("expid", ""),
                    "genome": record.get("genome", ""),
                    "agClass": record.get("agClass", ""),
                    "agSubClass": record.get("agSubClass", ""),
                    "clClass": record.get("clClass", ""),
                    "clSubClass": record.get("clSubClass", ""),
                    "title": record.get("title", ""),
                    "attributes": record.get("attributes", ""),
                    "readInfo": record.get("readInfo", ""),
                    "clSubClassInfo": record.get("clSubClassInfo", ""),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    deg = load_deg()
    deg_genes = set(deg["gene_name"])
    n_down_degs = int((deg["regulation"] == "DOWN").sum())
    n_up_degs = int((deg["regulation"] == "UP").sum())
    universe = load_universe()
    deg_in_universe = deg_genes & universe

    summary_rows = []
    overlap_frames = []
    score_cutoff_rows = []
    overlap_experiment_rows = []

    for window, path in WINDOWS.items():
        targets = load_chipatlas(path)
        target_genes = set(targets["gene_name"])
        targets_in_universe = target_genes & universe
        overlap_genes = target_genes & deg_genes
        overlap_in_universe = targets_in_universe & deg_in_universe

        overlap = (
            deg[deg["gene_name"].isin(overlap_genes)]
            .merge(targets, how="left", on="gene_name")
            .sort_values(["average_score", "gene_name"], ascending=[False, True])
            .reset_index(drop=True)
        )
        overlap.insert(0, "window", window)
        overlap_frames.append(overlap)

        pvalue = hypergeom_sf(
            len(overlap_in_universe),
            len(universe),
            len(targets_in_universe),
            len(deg_in_universe),
        )
        overlap_down = int((overlap["regulation"] == "DOWN").sum())
        down_bias_pvalue = hypergeom_sf(
            overlap_down,
            len(deg_genes),
            n_down_degs,
            len(overlap_genes),
        )
        summary_rows.append(
            {
                "source": "ChIP-Atlas hg38 TERF1 target genes",
                "window": window,
                "target_genes_chipseq_supported": len(target_genes),
                "target_genes_in_annotation_universe": len(targets_in_universe),
                "deg_genes": len(deg_genes),
                "deg_genes_in_annotation_universe": len(deg_in_universe),
                "overlap_targets_vs_degs": len(overlap_genes),
                "overlap_targets_vs_degs_in_universe": len(overlap_in_universe),
                "overlap_up": int((overlap["regulation"] == "UP").sum()),
                "overlap_down": overlap_down,
                "hypergeom_p_annotation_universe": pvalue,
                "down_bias_hypergeom_p_within_degs": down_bias_pvalue,
                "overlap_genes": ",".join(sorted(overlap_genes)),
            }
        )
        for experiment in [c for c in targets.columns if c.startswith("SRX")]:
            hit_genes = set(overlap.loc[overlap[experiment] > 0, "gene_name"])
            overlap_experiment_rows.append(
                {
                    "window": window,
                    "experiment": experiment,
                    "overlap_degs": len(hit_genes),
                    "overlap_genes": ",".join(sorted(hit_genes)),
                }
            )

        for cutoff in [25, 50, 100, 250, 500, 1000]:
            filtered = targets[targets["average_score"] >= cutoff]
            filtered_genes = set(filtered["gene_name"])
            filtered_in_universe = filtered_genes & universe
            filtered_overlap = filtered_genes & deg_genes
            filtered_overlap_in_universe = filtered_in_universe & deg_in_universe
            score_cutoff_rows.append(
                {
                    "window": window,
                    "average_score_cutoff": cutoff,
                    "target_genes": len(filtered_genes),
                    "targets_in_annotation_universe": len(filtered_in_universe),
                    "overlap_degs": len(filtered_overlap),
                    "overlap_degs_in_universe": len(filtered_overlap_in_universe),
                    "hypergeom_p_annotation_universe": hypergeom_sf(
                        len(filtered_overlap_in_universe),
                        len(universe),
                        len(filtered_in_universe),
                        len(deg_in_universe),
                    ),
                    "overlap_genes": ",".join(sorted(filtered_overlap)),
                }
            )

    summary = pd.DataFrame(summary_rows)
    overlap_all = pd.concat(overlap_frames, ignore_index=True)
    score_cutoffs = pd.DataFrame(score_cutoff_rows)
    overlap_experiments = pd.DataFrame(overlap_experiment_rows)
    experiment_metadata = summarize_experiment_metadata()

    summary.to_csv(OUTDIR / "chipatlas_terf1_deg_overlap_summary.tsv", sep="\t", index=False)
    overlap_all.to_csv(OUTDIR / "chipatlas_terf1_deg_overlap.tsv", sep="\t", index=False)
    score_cutoffs.to_csv(OUTDIR / "chipatlas_terf1_score_cutoff_summary.tsv", sep="\t", index=False)
    overlap_experiments.to_csv(
        OUTDIR / "chipatlas_terf1_overlap_by_experiment.tsv", sep="\t", index=False
    )
    experiment_metadata.to_csv(
        OUTDIR / "chipatlas_terf1_experiment_metadata.tsv", sep="\t", index=False
    )

    lines = [
        "# ChIP-Atlas hg38 TERF1 target-gene check",
        "",
        "Input: ChIP-Atlas precomputed hg38 TERF1 target-gene tables.",
        "Rows supported only by STRING, not by a ChIP-seq sample score, are excluded.",
        "",
        "Experiment columns present in the TERF1 table:",
        "- SRX4957970: 293",
        "- SRX359962: LoVo",
        "- SRX472275: lymphoblastoid cell line",
        "- SRX4957969: smooth muscle",
        "",
        f"DEG direction background: {n_up_degs} up, {n_down_degs} down.",
        "",
        "Summary:",
    ]
    for row in summary_rows:
        genes = row["overlap_genes"] or "None"
        cutoff = score_cutoffs[
            (score_cutoffs["window"] == row["window"])
            & (score_cutoffs["average_score_cutoff"] == 100)
        ].iloc[0]
        by_experiment = overlap_experiments[overlap_experiments["window"] == row["window"]]
        dominant = by_experiment.sort_values("overlap_degs", ascending=False).iloc[0]
        lines.extend(
            [
                f"- TSS +/- {row['window']}: {row['target_genes_chipseq_supported']} target genes; "
                f"{row['overlap_targets_vs_degs']} overlap DEGs "
                f"({row['overlap_up']} up, {row['overlap_down']} down); "
                f"target-set hypergeometric p={row['hypergeom_p_annotation_universe']:.6g}; "
                f"down-skew p={row['down_bias_hypergeom_p_within_degs']:.6g}",
                f"  dominant overlap source: {dominant['experiment']} "
                f"({dominant['overlap_degs']} overlap DEGs)",
                f"  average_score >=100: {cutoff['overlap_degs']} overlap DEGs",
                f"  overlap: {genes}",
            ]
        )
    lines.extend(
        [
            "",
            "Interpretation rule: ChIP-Atlas is a newer hg38 public-data cross-check,",
            "but these experiments are not matched to the original siTRF1 context.",
            "The broad low/moderate-score target sets overlap the DEG list, but the",
            "strongest average-score tier has no DEG overlap and most hits come from",
            "one LoVo experiment. This is evidence to flag candidates, not enough to",
            "claim direct TRF1 control of the siTRF1 DEG program.",
            "",
        ]
    )
    (OUTDIR / "chipatlas_terf1_deg_overlap_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

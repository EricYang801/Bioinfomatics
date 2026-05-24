#!/usr/bin/env python3
"""AlphaGenome BJ-cell pilot for selected siTRF1 DEG loci."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client


ROOT = Path(__file__).resolve().parents[1]
ANNOTATION = ROOT / "data/gene_annotation.tsv"
BETA_BSF = ROOT / "results/integration_analysis/beta/siTRF1_DEGs.bsf"
OUTDIR = ROOT / "results/alphagenome/pilot"

PILOT_GENES = [
    "TERF1",
    "DDIT4",
    "ASNS",
    "CHAC1",
    "DDB2",
    "PHLDA3",
    "IFITM1",
    "RNASEL",
    "SLFN5",
]

OUTPUT_ATTRS = [
    ("rna_seq", dna_client.OutputType.RNA_SEQ),
    ("dnase", dna_client.OutputType.DNASE),
    ("chip_histone", dna_client.OutputType.CHIP_HISTONE),
    ("chip_tf", dna_client.OutputType.CHIP_TF),
]


def load_loci() -> pd.DataFrame:
    ann = pd.read_csv(ANNOTATION, sep="\t", dtype={"chr": str, "gene_name": str})
    bsf = pd.read_csv(
        BETA_BSF,
        sep="\t",
        header=None,
        names=["gene_name", "log2FC", "pvalue"],
        dtype={"gene_name": str, "log2FC": float, "pvalue": float},
    ).drop_duplicates("gene_name")
    ann = ann[ann["gene_name"].isin(PILOT_GENES)].copy()
    ann["gene_start"] = pd.to_numeric(ann["gene_start"], errors="raise").astype(int)
    ann["gene_end"] = pd.to_numeric(ann["gene_end"], errors="raise").astype(int)
    ann["tss"] = np.where(ann["strand"].eq("+"), ann["gene_start"], ann["gene_end"])
    ann["requested_order"] = ann["gene_name"].map({g: i for i, g in enumerate(PILOT_GENES)})
    ann = (
        ann.sort_values(["requested_order", "chr", "gene_start", "gene_end"])
        .drop_duplicates("gene_name")
        .merge(bsf, how="left", on="gene_name")
    )
    return ann


def summarize_track(gene_row: pd.Series, interval: genome.Interval, attr: str, track) -> list[dict[str, object]]:
    if track is None or getattr(track, "num_tracks", 0) == 0:
        return []
    values = np.asarray(track.values)
    metadata = track.metadata.reset_index(drop=True)
    if values.ndim == 1:
        values = values[:, None]
    tss_bin = int((int(gene_row["tss"]) - interval.start) // int(track.resolution))
    rows = []
    for idx, meta in metadata.iterrows():
        column = values[:, idx]
        tss_value = column[tss_bin] if 0 <= tss_bin < len(column) else np.nan
        rows.append(
            {
                "gene_name": gene_row["gene_name"],
                "log2FC": gene_row.get("log2FC", np.nan),
                "pvalue": gene_row.get("pvalue", np.nan),
                "gene_chr": gene_row["chr"],
                "gene_start": gene_row["gene_start"],
                "gene_end": gene_row["gene_end"],
                "strand": gene_row["strand"],
                "tss": gene_row["tss"],
                "query_interval": f"{interval.chromosome}:{interval.start}-{interval.end}",
                "output_attr": attr,
                "resolution": track.resolution,
                "track_name": meta.get("name", ""),
                "assay_title": meta.get("Assay title", ""),
                "ontology_curie": meta.get("ontology_curie", ""),
                "biosample_name": meta.get("biosample_name", ""),
                "histone_mark": meta.get("histone_mark", ""),
                "transcription_factor": meta.get("transcription_factor", ""),
                "track_strand": meta.get("strand", ""),
                "mean_signal": float(np.nanmean(column)),
                "max_signal": float(np.nanmax(column)),
                "sum_signal": float(np.nansum(column)),
                "tss_bin_signal": float(tss_value) if not np.isnan(tss_value) else np.nan,
            }
        )
    return rows


def main() -> None:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        raise SystemExit("ALPHAGENOME_API_KEY is required.")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    loci = load_loci()
    loci.to_csv(OUTDIR / "alphagenome_bj_pilot_loci.tsv", sep="\t", index=False)

    model = dna_client.create(api_key)
    rows: list[dict[str, object]] = []
    failures: list[str] = []
    requested_outputs = [item[1] for item in OUTPUT_ATTRS]

    for _, gene_row in loci.iterrows():
        interval = genome.Interval(
            chromosome=str(gene_row["chr"]),
            start=int(gene_row["tss"]),
            end=int(gene_row["tss"]),
        ).resize(2**17)
        try:
            output = model.predict_interval(
                interval=interval,
                requested_outputs=requested_outputs,
                ontology_terms=["EFO:0002779"],
            )
        except Exception as exc:  # API errors should be recorded, not hidden.
            failures.append(f"{gene_row['gene_name']}: {type(exc).__name__}: {exc}")
            continue

        for attr, _output_type in OUTPUT_ATTRS:
            rows.extend(summarize_track(gene_row, interval, attr, getattr(output, attr)))

    result = pd.DataFrame(rows)
    result.to_csv(OUTDIR / "alphagenome_bj_pilot_track_summary.tsv", sep="\t", index=False)

    if not result.empty:
        compact = (
            result.groupby(["gene_name", "output_attr"], dropna=False)
            .agg(
                n_tracks=("track_name", "size"),
                max_signal=("max_signal", "max"),
                max_tss_bin_signal=("tss_bin_signal", "max"),
                tracks=("track_name", lambda x: ";".join(map(str, x))),
            )
            .reset_index()
        )
    else:
        compact = pd.DataFrame(
            columns=["gene_name", "output_attr", "n_tracks", "max_signal", "max_tss_bin_signal", "tracks"]
        )
    compact.to_csv(OUTDIR / "alphagenome_bj_pilot_compact_summary.tsv", sep="\t", index=False)

    lines = [
        "# AlphaGenome BJ pilot",
        "",
        "Context: BJ is the closest AlphaGenome metadata match to GSM638201",
        "(`TRF1 monoclonal ChIP BJ fibroblasts`). This pilot uses 131,072 bp",
        "TSS-centered hg38 intervals and ontology `EFO:0002779`.",
        "",
        f"- Requested genes: {', '.join(PILOT_GENES)}",
        f"- Loci found in annotation: {len(loci)}",
        f"- Track summary rows: {len(result)}",
        f"- Failed genes: {len(failures)}",
        "",
    ]
    if failures:
        lines.append("Failures:")
        lines.extend(f"- {item}" for item in failures)
        lines.append("")
    lines.extend(
        [
            "Interpretation rule: this describes predicted basal regulatory context",
            "in BJ-like metadata only. It is not a siTRF1 perturbation model and",
            "cannot establish upstream/downstream causality.",
            "",
        ]
    )
    (OUTDIR / "alphagenome_bj_pilot_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

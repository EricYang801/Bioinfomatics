#!/usr/bin/env python3
"""Probe AlphaGenome metadata without writing the API key to disk."""

from __future__ import annotations

import sys
import os
from pathlib import Path

import pandas as pd
from alphagenome.models import dna_client


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/alphagenome/metadata"


def simple_table(df: pd.DataFrame) -> str:
    return df.to_csv(sep="\t", index=False).strip()


def main() -> None:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        api_key = sys.stdin.read().strip()
    if not api_key:
        raise SystemExit("AlphaGenome API key was not provided on stdin.")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    model = dna_client.create(api_key)
    metadata = model.output_metadata(dna_client.Organism.HOMO_SAPIENS).concatenate()
    metadata.to_csv(OUTDIR / "alphagenome_human_output_metadata.tsv", sep="\t", index=False)

    counts = (
        metadata.groupby("output_type", dropna=False)
        .size()
        .reset_index(name="n_tracks")
        .sort_values("output_type")
    )
    counts.to_csv(OUTDIR / "alphagenome_human_output_type_counts.tsv", sep="\t", index=False)

    text_cols = [
        col
        for col in ["name", "Assay title", "biosample_name", "transcription_factor"]
        if col in metadata.columns
    ]
    haystack = metadata[text_cols].fillna("").astype(str).agg(" ".join, axis=1)
    trf_mask = haystack.str.contains("TERF1|TRF1|telomeric repeat", case=False, regex=True)
    trf = metadata[trf_mask].copy()
    trf.to_csv(OUTDIR / "alphagenome_terf1_trf1_candidate_tracks.tsv", sep="\t", index=False)

    keywords = [
        "fibroblast",
        "epithelial",
        "HeLa",
        "HEK",
        "K562",
        "HepG2",
        "A549",
        "HCT116",
        "U2OS",
        "MCF",
        "T cell",
        "B cell",
    ]
    rows = []
    for keyword in keywords:
        mask = metadata.get("biosample_name", pd.Series("", index=metadata.index)).fillna("").astype(str).str.contains(
            keyword, case=False, regex=False
        )
        subset = metadata[mask]
        rows.append(
            {
                "keyword": keyword,
                "n_tracks": len(subset),
                "n_biosamples": subset.get("biosample_name", pd.Series(dtype=str)).nunique(),
                "ontology_terms": ",".join(
                    sorted(
                        str(x)
                        for x in subset.get("ontology_curie", pd.Series(dtype=str)).dropna().unique()
                    )[:20]
                ),
                "example_biosamples": ",".join(
                    sorted(
                        str(x)
                        for x in subset.get("biosample_name", pd.Series(dtype=str)).dropna().unique()
                    )[:20]
                ),
            }
        )
    biosamples = pd.DataFrame(rows)
    biosamples.to_csv(
        OUTDIR / "alphagenome_biosample_keyword_candidates.tsv", sep="\t", index=False
    )

    chip_tf = metadata[metadata["output_type"].astype(str).str.contains("CHIP_TF")]
    tf_counts = (
        chip_tf.groupby("transcription_factor", dropna=False)
        .size()
        .reset_index(name="n_tracks")
        .sort_values(["n_tracks", "transcription_factor"], ascending=[False, True])
    )
    tf_counts.to_csv(OUTDIR / "alphagenome_chip_tf_factor_counts.tsv", sep="\t", index=False)

    lines = [
        "# AlphaGenome metadata probe",
        "",
        f"- Total human metadata rows: {len(metadata)}",
        f"- Output types: {', '.join(counts['output_type'].astype(str))}",
        f"- CHIP_TF rows: {len(chip_tf)}",
        f"- Unique CHIP_TF transcription factors: {chip_tf['transcription_factor'].nunique(dropna=True)}",
        f"- TERF1/TRF1 candidate metadata rows: {len(trf)}",
        "",
        "Output type counts:",
        simple_table(counts),
        "",
        "Biosample keyword candidates:",
        simple_table(biosamples),
        "",
    ]
    if len(trf):
        lines.extend(["TERF1/TRF1 candidate rows:", simple_table(trf.head(20)), ""])
    else:
        lines.extend(
            [
                "TERF1/TRF1 candidate rows: none found in metadata text fields.",
                "Interpretation: AlphaGenome cannot currently be used here as direct",
                "TERF1/TRF1 ChIP-TF evidence unless a hidden/renamed track is later",
                "identified manually.",
                "",
            ]
        )
    (OUTDIR / "alphagenome_metadata_probe_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

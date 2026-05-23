"""DEG list reading and matching to annotation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .chroms import gene_id_base


def read_deg(deg_text: str) -> pd.DataFrame:
    deg_path = Path(deg_text)
    deg = pd.read_csv(deg_path, sep="\t", dtype=str)
    deg.columns = [str(col).strip() for col in deg.columns]
    required = {"gene_id", "gene_name", "regulation"}
    missing = required.difference(deg.columns)
    if missing:
        raise ValueError("DEG file is missing required columns: " + ", ".join(sorted(missing)))

    deg = deg[["gene_id", "gene_name", "regulation"]].copy()
    deg["gene_id"] = deg["gene_id"].fillna("").astype(str).str.strip()
    deg["gene_name"] = deg["gene_name"].fillna("").astype(str).str.strip()
    deg["regulation"] = deg["regulation"].fillna("").astype(str).str.strip().str.upper()
    invalid = sorted(set(deg["regulation"]) - {"UP", "DOWN"})
    if invalid:
        raise ValueError(f"Unsupported DEG regulation value(s): {', '.join(invalid)}")

    deg["_gene_id_base"] = deg["gene_id"].map(gene_id_base)
    deg["_deg_key"] = deg.apply(
        lambda row: row["_gene_id_base"] or row["gene_name"], axis=1
    )
    deg = deg.drop_duplicates(subset=["_deg_key", "regulation"]).reset_index(drop=True)
    return deg


def match_deg_to_annotation(deg: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    by_gene_id: dict[str, pd.DataFrame] = {
        gene_id: group for gene_id, group in annotation.groupby("_gene_id_base", sort=False)
    }
    by_gene_name: dict[str, pd.DataFrame] = {
        gene_name: group for gene_name, group in annotation.groupby("gene_name", sort=False)
    }
    records: list[dict[str, object]] = []

    for _, deg_row in deg.iterrows():
        matches = pd.DataFrame()
        if deg_row["_gene_id_base"]:
            matches = by_gene_id.get(deg_row["_gene_id_base"], pd.DataFrame())
        if matches.empty and deg_row["gene_name"]:
            matches = by_gene_name.get(deg_row["gene_name"], pd.DataFrame())
        if matches.empty:
            continue

        for _, ann_row in matches.iterrows():
            records.append(
                {
                    "gene_name": deg_row["gene_name"] or ann_row["gene_name"],
                    "gene_id": deg_row["gene_id"] or ann_row["gene_id"],
                    "regulation": deg_row["regulation"],
                    "chr": ann_row["chr"],
                    "gene_start": ann_row["gene_start"],
                    "gene_end": ann_row["gene_end"],
                    "strand": ann_row["strand"],
                    "_deg_key": deg_row["_deg_key"],
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "gene_name",
                "gene_id",
                "regulation",
                "chr",
                "gene_start",
                "gene_end",
                "strand",
                "_deg_key",
            ]
        )

    matched = pd.DataFrame.from_records(records)
    matched = matched.drop_duplicates(
        subset=["_deg_key", "chr", "gene_start", "gene_end", "strand"]
    )
    return matched.reset_index(drop=True)

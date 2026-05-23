"""SISSRs ChIP-seq peak parsing and cleaned BED output."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from .chroms import normalize_chr


PEAK_BED_COLUMNS = ["peak_chr", "peak_start", "peak_end", "peak_fold", "peak_pvalue"]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def resolve_existing_path(path_text: str, kind: str) -> tuple[Path, str | None]:
    path = Path(path_text)
    if path.exists():
        return path, None

    if path.suffix == ".gz":
        uncompressed = path.with_suffix("")
        if uncompressed.exists():
            return uncompressed, (
                f"{kind} file {path} was not found; using uncompressed file "
                f"{uncompressed}."
            )

    raise FileNotFoundError(f"{kind} file not found: {path}")


def read_peak(peak_text: str) -> tuple[pd.DataFrame, list[str], Path]:
    peak_path, message = resolve_existing_path(peak_text, "Peak")
    messages = [message] if message else []
    records: list[dict[str, object]] = []
    header_seen = False

    with open_text(peak_path) as peak_file:
        for line in peak_file:
            fields = line.strip().split()
            if not fields:
                continue
            if fields[:3] == ["Chr", "cStart", "cEnd"]:
                header_seen = True
                continue
            if not header_seen:
                continue
            if len(fields) < 6:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue

            records.append(
                {
                    "peak_chr": normalize_chr(fields[0]),
                    "peak_start": start,
                    "peak_end": end,
                    "peak_fold": fields[4],
                    "peak_pvalue": fields[5],
                }
            )

    if not records:
        raise ValueError(f"No peaks were parsed from {peak_path}")

    peak = pd.DataFrame.from_records(records)
    peak = peak[peak["peak_end"] > peak["peak_start"]].copy()
    peak = peak.sort_values(["peak_chr", "peak_start", "peak_end"]).reset_index(drop=True)
    return peak, messages, peak_path


def write_cleaned_peak(peaks: pd.DataFrame, outdir: Path) -> Path:
    path = outdir / "TRF1_cleaned_peaks.bed"
    peaks[["peak_chr", "peak_start", "peak_end"]].to_csv(
        path, sep="\t", header=False, index=False
    )
    return path

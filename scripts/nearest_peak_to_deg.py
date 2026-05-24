#!/usr/bin/env python3
"""For each DEG, find the nearest TRF1 peak to its TSS (hg18)."""
import bisect
from collections import defaultdict

REF = "/Users/ericyang/Documents/Codex/2026-05-24/base-ericyang-erics-macbook-air-beta/.venv-beta/lib/python3.13/site-packages/beta/references/hg18.refseq"
PEAKS = "/Users/ericyang/Downloads/BETA_run/TRF1_peaks.bed"
DEGS = "/Users/ericyang/Downloads/BETA_run/siTRF1_DEGs.bsf"
OUT = "/Users/ericyang/Downloads/BETA_run/deg_nearest_peak.tsv"

# 1. gene symbol -> list of TSS (one per transcript)
gene_tss = defaultdict(list)
with open(REF) as f:
    next(f)  # header
    for ln in f:
        nm, chrom, strand, tx_start, tx_end, sym = ln.rstrip("\n").split("\t")[:6]
        tss = int(tx_start) if strand == "+" else int(tx_end)
        gene_tss[sym].append((chrom, tss))

# 2. peaks per chr, sorted by center
peak_centers = defaultdict(list)
with open(PEAKS) as f:
    for ln in f:
        c, s, e = ln.split("\t")[:3]
        peak_centers[c].append((int(s) + int(e)) // 2)
for c in peak_centers:
    peak_centers[c].sort()

def nearest_peak(chrom, pos):
    arr = peak_centers.get(chrom)
    if not arr:
        return None, None
    i = bisect.bisect_left(arr, pos)
    cands = []
    if i < len(arr):
        cands.append(arr[i])
    if i > 0:
        cands.append(arr[i - 1])
    best = min(cands, key=lambda p: abs(p - pos))
    return best, best - pos  # signed distance peak - TSS

# 3. iterate DEGs
rows = []
with open(DEGS) as f:
    for ln in f:
        sym, lfc, p = ln.rstrip("\n").split("\t")
        lfc = float(lfc)
        tss_list = gene_tss.get(sym, [])
        if not tss_list:
            rows.append((sym, lfc, "NA", "NA", "NA", "not_in_hg18_refseq"))
            continue
        best = None
        for chrom, tss in tss_list:
            peak, dist = nearest_peak(chrom, tss)
            if peak is None:
                continue
            if best is None or abs(dist) < abs(best[2]):
                best = (chrom, tss, dist, peak)
        if best is None:
            rows.append((sym, lfc, "NA", "NA", "NA", "no_peak_on_chr"))
        else:
            chrom, tss, dist, peak = best
            rows.append((sym, lfc, chrom, tss, dist, peak))

rows.sort(key=lambda r: (abs(r[4]) if isinstance(r[4], int) else 10**12))

with open(OUT, "w") as f:
    f.write("symbol\tlog2FC\tdirection\tchrom\tTSS\tnearest_peak_center\tdistance(peak-TSS)\n")
    for sym, lfc, chrom, tss, dist, peak in rows:
        direction = "up" if lfc > 0 else "down"
        if isinstance(dist, int):
            f.write(f"{sym}\t{lfc}\t{direction}\t{chrom}\t{tss}\t{peak}\t{dist}\n")
        else:
            f.write(f"{sym}\t{lfc}\t{direction}\tNA\tNA\tNA\t{peak}\n")

within_100k = [r for r in rows if isinstance(r[4], int) and abs(r[4]) <= 100_000]
print(f"DEGs total: {len(rows)}")
print(f"DEGs with peak within 100 kb of TSS: {len(within_100k)}")
print()
print(f"{'symbol':<15}{'direction':<10}{'log2FC':>10}{'chrom':>8}{'distance':>15}")
for sym, lfc, chrom, tss, dist, peak in within_100k:
    direction = "up" if lfc > 0 else "down"
    print(f"{sym:<15}{direction:<10}{lfc:>10.3f}{chrom:>8}{dist:>15,}")
print()
print(f"Wrote ranked table to {OUT}")

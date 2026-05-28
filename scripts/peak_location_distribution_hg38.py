"""TRF1 peak location distribution (hg38) + DEG overlap matrix.

Input:  ENCODE ENCSR031WWM conservative IDR peaks (hg38, HepG2 CRISPR-tagged TERF1)
        274 siTRF1 DEGs with hg38 TSS from data/gene_annotation.tsv
Output: results/integration_analysis/peak_location_analysis_hg38/
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BED_GZ = REPO / "data/ENCFF540ATM_TERF1_hg38_conservative_IDR.bed.gz"
GENE_ANN = REPO / "data/gene_annotation.tsv"
DEG_LIST = REPO / "data/siTRF1_DEG_list_gene_name.txt"
OUTDIR = REPO / "results/integration_analysis/peak_location_analysis_hg38"
OUTDIR.mkdir(parents=True, exist_ok=True)

# hg38 chromosome lengths
HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,
}

SUBTEL = 100_000
NEAR_TEL = 1_000_000
CLUSTER_RADIUS = 1_000_000
CLUSTER_MIN = 5


def load_peaks(path: Path):
    import gzip, io
    opener = gzip.open if str(path).endswith(".gz") else open
    peaks = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            center = (start + end) // 2
            # col7 = fold enrichment in narrowPeak
            try:
                fold = float(parts[6]) if len(parts) >= 7 and parts[6] not in {"", ".", "-1"} else 0.0
            except ValueError:
                fold = 0.0
            peaks.append({"chrom": chrom, "start": start, "end": end,
                          "center": center, "fold": fold})
    return peaks


def classify(peak):
    chrom = peak["chrom"]
    center = peak["center"]
    if chrom == "chrM":
        return "chrM"
    L = HG38.get(chrom)
    if L is None:
        return "other"
    d5 = center
    d3 = L - center
    if d5 < SUBTEL:
        return "5p_subtelomere"
    if d3 < SUBTEL:
        return "3p_subtelomere"
    if d5 < NEAR_TEL:
        return "5p_near_telomere"
    if d3 < NEAR_TEL:
        return "3p_near_telomere"
    return "interstitial"


def detect_clusters(peaks):
    by_chr = defaultdict(list)
    for p in peaks:
        by_chr[p["chrom"]].append(p)
    clusters = []
    for chrom, plist in by_chr.items():
        plist.sort(key=lambda x: x["center"])
        n = len(plist)
        i = 0
        while i < n:
            j = i
            while j + 1 < n and (plist[j + 1]["center"] - plist[i]["center"]) <= CLUSTER_RADIUS:
                j += 1
            count = j - i + 1
            if count >= CLUSTER_MIN:
                clusters.append({
                    "chrom": chrom,
                    "start": plist[i]["center"],
                    "end": plist[j]["center"],
                    "n_peaks": count,
                    "max_fold": max(p["fold"] for p in plist[i:j + 1]),
                })
                i = j + 1
            else:
                i += 1
    clusters.sort(key=lambda c: -c["n_peaks"])
    return clusters


def load_degs_with_tss(gene_ann: Path, deg_list: Path):
    # read DEG symbols — file may have header (gene_id\tgene_name\tregulation)
    degs = set()
    with deg_list.open() as fh:
        first = fh.readline().strip()
        if "\t" in first:
            # TSV with header; gene_name is col index 1
            parts = first.split("\t")
            name_idx = parts.index("gene_name") if "gene_name" in parts else 1
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) > name_idx:
                    degs.add(cols[name_idx])
        else:
            degs.add(first)
            for line in fh:
                if line.strip():
                    degs.add(line.strip())

    # read gene annotation — find TSS for each DEG
    tss_map = {}
    with gene_ann.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sym = row.get("gene_name") or row.get("symbol") or row.get("name", "")
            if sym not in degs:
                continue
            chrom = row.get("chr") or row.get("chrom") or row.get("chromosome", "")
            strand = row.get("strand", "+")
            try:
                start = int(row.get("gene_start") or row.get("start") or row.get("txStart") or 0)
                end = int(row.get("gene_end") or row.get("end") or row.get("txEnd") or 0)
            except ValueError:
                continue
            tss = start if strand == "+" else end
            if sym not in tss_map:
                tss_map[sym] = {"chrom": chrom, "tss": tss}
    return tss_map


def nearest_peak_to_tss(tss_map, peaks):
    # index peaks by chrom
    by_chr = defaultdict(list)
    for p in peaks:
        by_chr[p["chrom"]].append(p)
    results = []
    for sym, info in tss_map.items():
        chrom = info["chrom"]
        tss = info["tss"]
        plist = by_chr.get(chrom, [])
        if not plist:
            results.append({"gene": sym, "chrom": chrom, "tss": tss,
                             "nearest_dist": None, "nearest_fold": None, "peak_chrom": chrom})
            continue
        best = min(plist, key=lambda p: abs(p["center"] - tss))
        dist = best["center"] - tss
        results.append({"gene": sym, "chrom": chrom, "tss": tss,
                         "nearest_dist": dist, "nearest_fold": best["fold"],
                         "peak_chrom": chrom})
    return results


def main():
    peaks = load_peaks(BED_GZ)
    total = len(peaks)
    print(f"Loaded {total} peaks from {BED_GZ.name}")

    # 1. distribution
    bucket_counts = defaultdict(int)
    bucket_strong = defaultdict(int)
    for p in peaks:
        b = classify(p)
        bucket_counts[b] += 1
        if p["fold"] >= 10.0:
            bucket_strong[b] += 1

    dist_path = OUTDIR / "peak_distribution_by_region.tsv"
    order = ["5p_subtelomere", "3p_subtelomere", "5p_near_telomere",
             "3p_near_telomere", "interstitial", "chrM", "other"]
    with dist_path.open("w") as fh:
        fh.write("region\tn_peaks\tn_strong_fold_ge_10\tpct_of_total\n")
        for b in order:
            n = bucket_counts[b]
            if n == 0:
                continue
            fh.write(f"{b}\t{n}\t{bucket_strong[b]}\t{n/total*100:.1f}\n")

    # 2. clusters
    clusters = detect_clusters(peaks)
    total_in_clusters = sum(c["n_peaks"] for c in clusters)
    cluster_path = OUTDIR / "peak_clusters_1Mb_min5.tsv"
    with cluster_path.open("w") as fh:
        fh.write("chrom\tcluster_start\tcluster_end\tspan_bp\tn_peaks\tmax_fold\tregion_class\n")
        for c in clusters:
            span = c["end"] - c["start"]
            mid = (c["start"] + c["end"]) // 2
            cls = classify({"chrom": c["chrom"], "center": mid})
            fh.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{span}\t"
                     f"{c['n_peaks']}\t{c['max_fold']:.2f}\t{cls}\n")

    # 3. top 20 strongest
    strongest = sorted(peaks, key=lambda p: -p["fold"])[:20]
    strong_path = OUTDIR / "top20_strongest_peaks.tsv"
    with strong_path.open("w") as fh:
        fh.write("chrom\tstart\tend\tcenter\tfold\tregion_class\tdist_from_nearest_end_bp\n")
        for p in strongest:
            cls = classify(p)
            L = HG38.get(p["chrom"], 0)
            d_end = min(p["center"], L - p["center"]) if L else -1
            fh.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t{p['center']}\t"
                     f"{p['fold']:.2f}\t{cls}\t{d_end}\n")

    # 4. DEG overlap matrix
    tss_map = load_degs_with_tss(GENE_ANN, DEG_LIST)
    print(f"DEGs with hg38 TSS: {len(tss_map)}")
    deg_results = nearest_peak_to_tss(tss_map, peaks)

    windows = [10_000, 100_000, 1_000_000, 10_000_000]
    overlap_path = OUTDIR / "deg_peak_overlap_matrix.tsv"
    with overlap_path.open("w") as fh:
        fh.write("window_bp\tdeg_within_window\tpct_of_degs\n")
        n_deg = len(tss_map)
        for w in windows:
            n = sum(1 for r in deg_results if r["nearest_dist"] is not None and abs(r["nearest_dist"]) <= w)
            fh.write(f"{w}\t{n}\t{n/n_deg*100:.2f}\n")

    # DEGs within 100kb
    within_100k = [r for r in deg_results if r["nearest_dist"] is not None and abs(r["nearest_dist"]) <= 100_000]
    near_path = OUTDIR / "deg_within_100kb.tsv"
    with near_path.open("w") as fh:
        fh.write("gene\tchrom\ttss\tnearest_dist\tnearest_fold\tpeak_region\n")
        for r in sorted(within_100k, key=lambda x: abs(x["nearest_dist"])):
            cls = classify({"chrom": r["chrom"], "center": r["tss"]})
            fh.write(f"{r['gene']}\t{r['chrom']}\t{r['tss']}\t{r['nearest_dist']}\t"
                     f"{r['nearest_fold']:.2f}\t{cls}\n")

    # 5. summary
    summary = OUTDIR / "peak_location_summary_hg38.md"
    with summary.open("w") as fh:
        fh.write("# TRF1 hg38 peak location distribution and DEG overlap\n\n")
        fh.write(f"Input: ENCODE ENCSR031WWM conservative IDR peaks ({total} peaks, hg38, HepG2 CRISPR-tagged TERF1).\n\n")
        fh.write("## 1. Peak distribution by chromosomal region\n\n")
        fh.write("| region | n_peaks | strong (fold>=10) | % of total |\n")
        fh.write("|---|---:|---:|---:|\n")
        for b in order:
            n = bucket_counts[b]
            if n == 0:
                continue
            fh.write(f"| {b} | {n} | {bucket_strong[b]} | {n/total*100:.1f} |\n")
        fh.write("\n")
        fh.write(f"## 2. Cluster detection (>={CLUSTER_MIN} peaks within 1 Mb)\n\n")
        fh.write(f"- Total clusters: {len(clusters)}\n")
        fh.write(f"- Peaks in clusters: {total_in_clusters} ({total_in_clusters/total*100:.1f}%)\n\n")
        fh.write("Top clusters:\n\n")
        fh.write("| chrom | center range | n_peaks | max_fold | region_class |\n")
        fh.write("|---|---|---:|---:|---|\n")
        for c in clusters[:12]:
            mid = (c["start"] + c["end"]) // 2
            cls = classify({"chrom": c["chrom"], "center": mid})
            fh.write(f"| {c['chrom']} | {c['start']:,}–{c['end']:,} | "
                     f"{c['n_peaks']} | {c['max_fold']:.2f} | {cls} |\n")
        n_deg = len(tss_map)
        fh.write(f"\n## 3. DEG overlap matrix ({n_deg} siTRF1 DEGs with hg38 TSS)\n\n")
        fh.write("| window (bp) | DEG within window | % of DEGs |\n")
        fh.write("|---:|---:|---:|\n")
        for w in windows:
            n = sum(1 for r in deg_results if r["nearest_dist"] is not None and abs(r["nearest_dist"]) <= w)
            fh.write(f"| {w:,} | {n} | {n/n_deg*100:.2f}% |\n")
        fh.write(f"\n## 4. DEGs within 100 kb of any peak ({len(within_100k)} genes)\n\n")
        if within_100k:
            fh.write("| gene | chrom | dist_bp | fold | region |\n")
            fh.write("|---|---|---:|---:|---|\n")
            for r in sorted(within_100k, key=lambda x: abs(x["nearest_dist"]))[:20]:
                cls = classify({"chrom": r["chrom"], "center": r["tss"]})
                fh.write(f"| {r['gene']} | {r['chrom']} | {r['nearest_dist']:,} | "
                         f"{r['nearest_fold']:.2f} | {cls} |\n")
        else:
            fh.write("None.\n")

    print("Wrote:")
    for p in [dist_path, cluster_path, strong_path, overlap_path, near_path, summary]:
        print(" ", p.relative_to(REPO))


if __name__ == "__main__":
    main()

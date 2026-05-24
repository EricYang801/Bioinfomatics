#!/usr/bin/env python3
"""TTAGGG promoter motif scan — batched Ensembl REST POST (50 regions/call)."""
import json, re, time, urllib.request, urllib.error, sys
from collections import defaultdict

DEG = "/Users/ericyang/Github/Bioinfomatics/results/integration_analysis/beta/siTRF1_DEGs.bsf"
ANNO = "/Users/ericyang/Github/Bioinfomatics/data/gene_annotation.tsv"
OUT_TSV = "/Users/ericyang/Github/Bioinfomatics/results/integration_analysis/motif_scan/ttaggg_motif_scan_full.tsv"
OUT_SUM = "/Users/ericyang/Github/Bioinfomatics/results/integration_analysis/motif_scan/ttaggg_motif_summary_full.txt"

ENSEMBL = "https://rest.ensembl.org/sequence/region/human"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json",
           "User-Agent": "Mozilla/5.0 (BETA_run motif scan batched)"}
BATCH = 50

def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

# Load DEGs
degs = {}
with open(DEG) as f:
    for ln in f:
        s, l, p = ln.rstrip("\n").split("\t")
        degs[s] = float(l)

# Load annotation, build promoter windows per gene
gene_windows = defaultdict(list)  # gene -> [(chrom, start, end, strand)]
with open(ANNO) as f:
    next(f)
    for ln in f:
        chrom, gs, ge, gn, gid, strand = ln.rstrip("\n").split("\t")[:6]
        if gn not in degs: continue
        gs, ge = int(gs), int(ge)
        chrom = chrom.lstrip("chr")
        if chrom == "M": chrom = "MT"
        if chrom == "MT": continue  # skip chrM (analyzed separately)
        if strand == "+":
            ws, we = max(0, gs - 1000), gs
        else:
            ws, we = ge, ge + 1000
        # Ensembl 1-based inclusive: 0-based [ws,we) -> 1-based ws+1..we
        gene_windows[gn].append((chrom, ws + 1, we, strand))

all_windows = []  # (gene, idx_in_gene, chrom, e_start, e_end, strand)
for gn, ws in gene_windows.items():
    for i, (c, s, e, st) in enumerate(ws):
        all_windows.append((gn, i, c, s, e, st))

print(f"DEGs total: {len(degs)}")
print(f"DEGs in annotation (non-chrM): {len(gene_windows)}")
print(f"Promoter windows to fetch: {len(all_windows)}")
print(f"Batches (size={BATCH}): {(len(all_windows) + BATCH - 1) // BATCH}")

# Fetch in batches
seqs = {}  # (gene, idx) -> sequence (already reverse-complemented if minus)
fail = 0
for i in range(0, len(all_windows), BATCH):
    batch = all_windows[i:i + BATCH]
    regions = [f"{c}:{s}..{e}:1" for _, _, c, s, e, _ in batch]
    body = json.dumps({"regions": regions}).encode()
    req = urllib.request.Request(ENSEMBL, data=body, headers=HEADERS, method="POST")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            for (gn, idx, c, s, e, st), item in zip(batch, data):
                seq = (item.get("seq") or "").upper().replace("\n", "")
                if st == "-":
                    seq = revcomp(seq)
                seqs[(gn, idx)] = seq
            break
        except urllib.error.HTTPError as exc:
            err_body = exc.read()[:200].decode("utf-8", errors="replace")
            print(f"  batch {i // BATCH}: HTTP {exc.code} (attempt {attempt+1}): {err_body}")
            time.sleep((2 ** attempt) * 1.0)
        except Exception as exc:
            print(f"  batch {i // BATCH}: {type(exc).__name__} (attempt {attempt+1}): {exc}")
            time.sleep((2 ** attempt) * 1.0)
    else:
        fail += len(batch)
        for tup in batch:
            seqs[(tup[0], tup[1])] = ""
    if (i // BATCH) % 5 == 4:
        print(f"  done batch {i//BATCH + 1} / {(len(all_windows)+BATCH-1)//BATCH}, fetched={len(seqs)}")
    time.sleep(0.1)  # be polite between batches

print(f"Total sequences fetched: {len([s for s in seqs.values() if s])}, failed/empty: {sum(1 for s in seqs.values() if not s)}")

# Scan motifs, take MAX per gene
pat1 = re.compile(r"TTAGGG")
pat_tandem = re.compile(r"(?:TTAGGG){2,}")
pat_tandem3 = re.compile(r"(?:TTAGGG){3,}")
pat_ccctaa = re.compile(r"CCCTAA")

rows = []
for gn, ws in gene_windows.items():
    best = {"x1": 0, "x2": 0, "x3": 0, "cc": 0, "win": None}
    for i, (c, s, e, st) in enumerate(ws):
        seq = seqs.get((gn, i), "")
        if not seq or set(seq) <= {"N"}: continue
        x1 = len(pat1.findall(seq))
        x2 = len(pat_tandem.findall(seq))
        x3 = len(pat_tandem3.findall(seq))
        cc = len(pat_ccctaa.findall(seq))
        score = (x2, x1, x3)
        if score > (best["x2"], best["x1"], best["x3"]):
            best.update(x1=x1, x2=x2, x3=x3, cc=cc, win=(c, s, e, st, seq))
    rows.append((gn, degs[gn], best))

rows.sort(key=lambda r: (-r[2]["x2"], -r[2]["x1"], r[0]))

with open(OUT_TSV, "w") as f:
    f.write("gene_name\tlog2FC\tregulation\tn_TTAGGG_x1\tn_TTAGGG_x2\tn_TTAGGG_x3\tn_CCCTAA\tbest_chr\tbest_start\tbest_end\tstrand\n")
    for gn, lfc, b in rows:
        reg = "up" if lfc > 0 else "down"
        if b["win"]:
            c, s, e, st, _ = b["win"]
            f.write(f"{gn}\t{lfc}\t{reg}\t{b['x1']}\t{b['x2']}\t{b['x3']}\t{b['cc']}\t{c}\t{s}\t{e}\t{st}\n")
        else:
            f.write(f"{gn}\t{lfc}\t{reg}\t0\t0\t0\t0\tNA\tNA\tNA\tNA\n")

x1c = sum(1 for r in rows if r[2]["x1"] >= 1)
x2c = sum(1 for r in rows if r[2]["x2"] >= 1)
x3c = sum(1 for r in rows if r[2]["x3"] >= 1)
scanned = sum(1 for r in rows if r[2]["win"])

summary = [
    "TTAGGG promoter motif scan — full (Ensembl REST batched POST, 50 region/call)",
    f"Total DEGs in input: {len(degs)}",
    f"DEGs with annotation (chrM excluded — handled separately): {len(gene_windows)}",
    f"DEGs WITHOUT annotation / chrM-only: {len(degs) - len(gene_windows)}",
    f"Promoter windows fetched: {len(seqs)} / {len(all_windows)}  (failed/empty: {fail})",
    f"DEGs with ≥1 scanned non-empty window: {scanned}",
    f"DEGs with n_TTAGGG_x1 >= 1: {x1c}",
    f"DEGs with n_TTAGGG_x2 >= 1 (TRF1 minimum): {x2c}",
    f"DEGs with n_TTAGGG_x3 >= 1 (TRF1 strong):  {x3c}",
    "",
    "Top 10 by motif count:",
]
for gn, lfc, b in rows[:10]:
    reg = "up" if lfc > 0 else "down"
    win = f"{b['win'][0]}:{b['win'][1]}-{b['win'][2]}({b['win'][3]})" if b["win"] else "NA"
    summary.append(f"  {gn:<15} log2FC={lfc:+.3f} {reg:>4}  x1={b['x1']:>2} x2={b['x2']:>2} x3={b['x3']:>2} CCCTAA={b['cc']:>2}  win={win}")

summary.append("")
summary.append(f"Key question: Does any DEG promoter contain ≥2 tandem TTAGGG repeats for stable TRF1 binding? {'YES (' + str(x2c) + ' DEGs)' if x2c else 'NO'}")
out_text = "\n".join(summary) + "\n"

with open(OUT_SUM, "w") as f:
    f.write(out_text)
print()
print(out_text)

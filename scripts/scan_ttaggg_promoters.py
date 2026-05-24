#!/usr/bin/env python3
import csv
import re
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta


DEG_PATH = Path("/Users/ericyang/Downloads/BETA_run/siTRF1_DEGs.bsf")
ANNOT_PATH = Path("/Users/ericyang/Github/Bioinfomatics/data/gene_annotation.tsv")
FASTA_PATH = Path("/Users/ericyang/Github/Bioinfomatics/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")
OUT_TSV = Path("/Users/ericyang/Downloads/BETA_run/ttaggg_motif_scan.tsv")
OUT_SUMMARY = Path("/Users/ericyang/Downloads/BETA_run/ttaggg_motif_summary.txt")
FALLBACK_OUT_TSV = Path("/private/tmp/ttaggg_motif_scan.tsv")
FALLBACK_OUT_SUMMARY = Path("/private/tmp/ttaggg_motif_summary.txt")
TMP_FAI = Path("/private/tmp/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fasta_key(chrom, keys):
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append("chr" + chrom)
    for key in candidates:
        if key in keys:
            return key
    return None


def promoter_coords(gene_start, gene_end, strand, convention):
    if convention == "bed0":
        if strand == "+":
            return max(0, gene_start - 1000), gene_start
        return gene_end, gene_end + 1000

    # Interpret annotation positions as 1-based closed and convert to
    # 0-based half-open slicing. This differs by one base from BED-style.
    if strand == "+":
        end = max(0, gene_start - 1)
        return max(0, end - 1000), end
    start = gene_end
    return start, start + 1000


def extract_promoter(fasta, keys, chrom, gene_start, gene_end, strand, convention):
    key = fasta_key(chrom, keys)
    if key is None:
        return None, None, None, None
    start, end = promoter_coords(gene_start, gene_end, strand, convention)
    if end - start != 1000:
        return key, start, end, ""
    seq = str(fasta[key][start:end]).upper()
    if strand == "-":
        seq = revcomp(seq).upper()
    return key, start, end, seq


def motif_counts(seq):
    return {
        "n_TTAGGG_x1": len(re.findall("TTAGGG", seq)),
        "n_TTAGGG_x2": len(re.findall(r"(?:TTAGGG){2,}", seq)),
        "n_TTAGGG_x3": len(re.findall(r"(?:TTAGGG){3,}", seq)),
        "n_CCCTAA": len(re.findall("CCCTAA", seq)),
    }


def load_degs():
    degs = []
    with DEG_PATH.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            gene, log2fc, pvalue = row[:3]
            degs.append((gene, log2fc, pvalue))
    return degs


def load_annotation():
    ann = defaultdict(list)
    with ANNOT_PATH.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ann[row["gene_name"]].append(
                (
                    row["chr"],
                    int(row["gene_start"]),
                    int(row["gene_end"]),
                    row["strand"],
                )
            )
    return ann


def convention_probe(fasta, keys, degs, ann):
    probes = []
    for gene, _, _ in degs:
        for entry in ann.get(gene, []):
            if fasta_key(entry[0], keys) is None:
                continue
            probes.append((gene, entry))
            if len(probes) == 5:
                break
        if len(probes) == 5:
            break

    stats = {}
    examples = []
    for convention in ("bed0", "one_based"):
        exact_len = 0
        non_all_n = 0
        for gene, (chrom, start, end, strand) in probes:
            key, pstart, pend, seq = extract_promoter(
                fasta, keys, chrom, start, end, strand, convention
            )
            if seq is None:
                examples.append(
                    f"{convention}\t{gene}\t{chrom}:NA-NA({strand})\tlen=0\tall_N=False\tmissing_chrom=True"
                )
                continue
            ok_len = len(seq) == 1000
            exact_len += int(ok_len)
            non_all_n += int(ok_len and set(seq) != {"N"})
            examples.append(
                f"{convention}\t{gene}\t{chrom}:{pstart}-{pend}({strand})\tlen={len(seq)}\tall_N={ok_len and set(seq) == {'N'}}"
            )
        stats[convention] = (non_all_n, exact_len)

    chosen = max(stats, key=lambda c: (stats[c][0], stats[c][1], c == "bed0"))
    return chosen, stats, examples


def main():
    degs = load_degs()
    ann = load_annotation()

    fasta = Fasta(
        str(FASTA_PATH),
        indexname=str(TMP_FAI),
        sequence_always_upper=True,
        rebuild=True,
        build_index=True,
    )
    keys = set(fasta.keys())

    convention, convention_stats, convention_examples = convention_probe(
        fasta, keys, degs, ann
    )

    rows = []
    not_found = []
    all_n_windows = []
    missing_chrom = []
    scanned_genes = set()

    for gene, log2fc, _pvalue in degs:
        entries = ann.get(gene, [])
        if not entries:
            not_found.append(gene)
            row = {
                "gene_name": gene,
                "log2FC": log2fc,
                "regulation": "up" if float(log2fc) > 0 else "down",
                "n_TTAGGG_x1": 0,
                "n_TTAGGG_x2": 0,
                "n_TTAGGG_x3": 0,
                "n_CCCTAA": 0,
                "best_window_chrom": "NA",
                "best_window_start": "NA",
                "best_window_end": "NA",
                "strand": "NA",
                "seq_preview": "NA",
            }
            rows.append(row)
            continue

        best = None
        for chrom, start, end, strand in entries:
            key, pstart, pend, seq = extract_promoter(
                fasta, keys, chrom, start, end, strand, convention
            )
            if key is None:
                missing_chrom.append(f"{gene}:{chrom}")
                continue
            if len(seq) != 1000:
                continue
            if set(seq) == {"N"}:
                all_n_windows.append(f"{gene}:{key}:{pstart}-{pend}:{strand}")
                continue
            scanned_genes.add(gene)
            counts = motif_counts(seq)
            candidate = {
                "gene_name": gene,
                "log2FC": log2fc,
                "regulation": "up" if float(log2fc) > 0 else "down",
                **counts,
                "best_window_chrom": key,
                "best_window_start": pstart,
                "best_window_end": pend,
                "strand": strand,
                "seq_preview": seq[:80],
                "full_seq": seq,
            }
            score = (
                candidate["n_TTAGGG_x2"],
                candidate["n_TTAGGG_x3"],
                candidate["n_TTAGGG_x1"],
                candidate["n_CCCTAA"],
            )
            if best is None or score > best[0]:
                best = (score, candidate)

        if best is None:
            row = {
                "gene_name": gene,
                "log2FC": log2fc,
                "regulation": "up" if float(log2fc) > 0 else "down",
                "n_TTAGGG_x1": 0,
                "n_TTAGGG_x2": 0,
                "n_TTAGGG_x3": 0,
                "n_CCCTAA": 0,
                "best_window_chrom": "NA",
                "best_window_start": "NA",
                "best_window_end": "NA",
                "strand": "NA",
                "seq_preview": "NA",
                "full_seq": "",
            }
        else:
            row = best[1]
        rows.append(row)

    rows.sort(
        key=lambda r: (
            int(r["n_TTAGGG_x2"]),
            int(r["n_TTAGGG_x1"]),
            int(r["n_TTAGGG_x3"]),
            r["gene_name"],
        ),
        reverse=True,
    )

    fieldnames = [
        "gene_name",
        "log2FC",
        "regulation",
        "n_TTAGGG_x1",
        "n_TTAGGG_x2",
        "n_TTAGGG_x3",
        "n_CCCTAA",
        "best_window_chrom",
        "best_window_start",
        "best_window_end",
        "strand",
        "seq_preview",
    ]
    actual_out_tsv = OUT_TSV
    actual_out_summary = OUT_SUMMARY
    try:
        fh = OUT_TSV.open("w", newline="")
    except PermissionError:
        actual_out_tsv = FALLBACK_OUT_TSV
        actual_out_summary = FALLBACK_OUT_SUMMARY
        fh = actual_out_tsv.open("w", newline="")

    with fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})

    found_count = sum(1 for gene, _, _ in degs if gene in ann)
    genes_without_scanned_window = len(degs) - len(scanned_genes)
    x1_count = sum(1 for r in rows if int(r["n_TTAGGG_x1"]) >= 1)
    x2_count = sum(1 for r in rows if int(r["n_TTAGGG_x2"]) >= 1)
    x3_count = sum(1 for r in rows if int(r["n_TTAGGG_x3"]) >= 1)
    top10 = rows[:10]

    lines = []
    lines.append("TTAGGG promoter motif scan summary")
    if actual_out_tsv != OUT_TSV:
        lines.append(
            f"requested output directory was not writable in this sandbox; wrote TSV to: {actual_out_tsv}"
        )
        lines.append(
            f"requested output directory was not writable in this sandbox; wrote summary to: {actual_out_summary}"
        )
    lines.append(f"total DEGs in input: {len(degs)}")
    lines.append(f"FASTA contigs available: {len(keys)} ({', '.join(list(fasta.keys())[:10])})")
    lines.append(f"DEGs found in annotation: {found_count}")
    lines.append(f"DEGs not found in annotation: {len(not_found)}")
    lines.append(
        "first 10 not-found names: "
        + (", ".join(not_found[:10]) if not_found else "None")
    )
    lines.append(f"coordinate convention selected: {convention}")
    lines.append(f"coordinate convention probe stats: {convention_stats}")
    lines.append("coordinate convention probe examples:")
    lines.extend("  " + item for item in convention_examples)
    lines.append(f"all-N promoter transcript windows skipped: {len(all_n_windows)}")
    lines.append(
        "first 10 all-N skipped windows: "
        + (", ".join(all_n_windows[:10]) if all_n_windows else "None")
    )
    lines.append(f"missing chromosome windows skipped: {len(missing_chrom)}")
    lines.append(
        "first 10 missing chromosome windows: "
        + (", ".join(missing_chrom[:10]) if missing_chrom else "None")
    )
    lines.append(f"DEGs with at least one scanned non-all-N promoter window: {len(scanned_genes)}")
    lines.append(f"DEGs without a scanned non-all-N promoter window: {genes_without_scanned_window}")
    lines.append(f"DEGs with n_TTAGGG_x1 >= 1: {x1_count}")
    lines.append(f"DEGs with n_TTAGGG_x2 >= 1: {x2_count}")
    lines.append(f"DEGs with n_TTAGGG_x3 >= 1: {x3_count}")
    lines.append("")
    lines.append("Top 10 DEGs by motif count:")
    for row in top10:
        lines.append(
            f">{row['gene_name']} log2FC={row['log2FC']} regulation={row['regulation']} "
            f"TTAGGG_x1={row['n_TTAGGG_x1']} TTAGGG_x2={row['n_TTAGGG_x2']} "
            f"TTAGGG_x3={row['n_TTAGGG_x3']} CCCTAA={row['n_CCCTAA']} "
            f"window={row['best_window_chrom']}:{row['best_window_start']}-{row['best_window_end']} "
            f"strand={row['strand']}"
        )
        lines.append(row.get("full_seq", ""))

    direct = (
        "YES"
        if x2_count > 0
        else "NO"
    )
    lines.append("")
    lines.append(
        f"Key question: Does any DEG promoter contain at least two tandem TTAGGG repeats for stable TRF1 binding? {direct}"
    )

    summary = "\n".join(lines) + "\n"
    actual_out_summary.write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()

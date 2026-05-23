import argparse
import sys
from pathlib import Path


GTF_PATH = Path("Homo_sapiens.GRCh38.108.gtf")
GTF_COLUMNS = {
    "seqname": 0,
    "chrom": 0,
    "source": 1,
    "feature": 2,
    "start": 3,
    "end": 4,
    "score": 5,
    "strand": 6,
    "frame": 7,
}


def parse_key_value(items, option_name):
    filters = []
    for item in items:
        if "=" not in item:
            raise SystemExit(f"{option_name} must use KEY=VALUE, got: {item}")
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key or not value:
            raise SystemExit(f"{option_name} must use KEY=VALUE, got: {item}")
        filters.append((key, value))
    return filters


def parse_attributes(attributes_text):
    attributes = {}
    for part in attributes_text.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        attributes[key] = value.strip().strip('"')
    return attributes


def matches_terms(line, terms, ignore_case):
    if ignore_case:
        line = line.lower()
        terms = [term.lower() for term in terms]
    return all(term in line for term in terms)


def matches_columns(fields, column_filters):
    for key, expected in column_filters:
        column_index = GTF_COLUMNS.get(key)
        if column_index is None:
            valid = ", ".join(sorted(GTF_COLUMNS))
            raise SystemExit(f"Unknown column '{key}'. Valid columns: {valid}")
        if fields[column_index] != expected:
            return False
    return True


def matches_attributes(fields, attr_filters):
    attributes = parse_attributes(fields[8])
    return all(attributes.get(key) == expected for key, expected in attr_filters)


def search_gtf(path, terms, column_filters, attr_filters, ignore_case, limit):
    found = 0

    with path.open() as gtf:
        for line_number, line in enumerate(gtf, 1):
            if line.startswith("#") or not matches_terms(line, terms, ignore_case):
                continue

            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) != 9:
                continue

            if not matches_columns(fields, column_filters):
                continue
            if attr_filters and not matches_attributes(fields, attr_filters):
                continue

            found += 1
            print(f"{line_number}:{line}", end="")

            if limit and found >= limit:
                break

    return found


def ask_limit(default=20):
    value = input(f"最多顯示幾筆？直接 Enter = {default}，0 = 不限制：").strip()
    if not value:
        return default
    try:
        limit = int(value)
    except ValueError:
        raise SystemExit("顯示筆數必須是整數。")
    if limit < 0:
        raise SystemExit("顯示筆數不能小於 0。")
    return limit


def run_text_search(prompt, default_limit=20):
    text = input(prompt).strip()
    if not text:
        raise SystemExit("沒有輸入搜尋內容。")
    return search_gtf(
        GTF_PATH,
        terms=text.split(),
        column_filters=[],
        attr_filters=[],
        ignore_case=True,
        limit=ask_limit(default=default_limit),
    )


def run_attr_search(attr_name, prompt, feature=None, default_limit=20):
    value = input(prompt).strip()
    if not value:
        raise SystemExit(f"沒有輸入 {attr_name}。")

    column_filters = []
    if feature:
        column_filters.append(("feature", feature))

    return search_gtf(
        GTF_PATH,
        terms=[],
        column_filters=column_filters,
        attr_filters=[(attr_name, value)],
        ignore_case=False,
        limit=ask_limit(default=default_limit),
    )


def run_column_search(column_name, prompt, default_limit=20):
    value = input(prompt).strip()
    if not value:
        raise SystemExit(f"沒有輸入 {column_name}。")
    return search_gtf(
        GTF_PATH,
        terms=[],
        column_filters=[(column_name, value)],
        attr_filters=[],
        ignore_case=False,
        limit=ask_limit(default=default_limit),
    )


def interactive_search():
    print("搜尋 Homo_sapiens.GRCh38.108.gtf")
    print("1. 任意關鍵字，例如 ADA、protein_coding、ENSG00000196839、transcript_id")
    print("2. gene name，例如 ADA")
    print("3. gene id，例如 ENSG00000196839")
    print("4. transcript id，例如 ENST00000673477")
    print("5. transcript name，例如 ATAD3B-206")
    print("6. gene biotype，例如 protein_coding、lncRNA")
    print("7. feature，例如 gene、transcript、exon、CDS")
    print("8. chromosome / seqname，例如 1、20、X、MT")
    print("9. 自訂 attribute，例如 gene_source = ensembl_havana")
    print("10. 自訂 GTF 欄位，例如 source = ensembl_havana、strand = +")

    mode = input("請選擇搜尋方式 [1-10，直接 Enter = 1]：").strip() or "1"

    if mode == "1":
        found = run_text_search("你要找什麼：")
    elif mode == "2":
        found = run_attr_search("gene_name", "你要找哪個基因：", feature="gene")
    elif mode == "3":
        found = run_attr_search("gene_id", "你要找哪個 gene_id：", feature="gene")
    elif mode == "4":
        found = run_attr_search("transcript_id", "你要找哪個 transcript_id：")
    elif mode == "5":
        found = run_attr_search("transcript_name", "你要找哪個 transcript_name：")
    elif mode == "6":
        found = run_attr_search("gene_biotype", "你要找哪個 gene_biotype：")
    elif mode == "7":
        found = run_column_search("feature", "你要找哪個 feature：")
    elif mode == "8":
        found = run_column_search("seqname", "你要找哪條 chromosome / seqname：")
    elif mode == "9":
        key = input("attribute 名稱，例如 gene_name、gene_id、transcript_id：").strip()
        value = input("attribute 值：").strip()
        if not key or not value:
            raise SystemExit("attribute 名稱和值都要輸入。")
        found = search_gtf(
            GTF_PATH,
            terms=[],
            column_filters=[],
            attr_filters=[(key, value)],
            ignore_case=False,
            limit=ask_limit(default=20),
        )
    elif mode == "10":
        key = input("欄位名稱，例如 seqname、feature、start、end、strand：").strip()
        value = input("欄位值：").strip()
        if not key or not value:
            raise SystemExit("欄位名稱和值都要輸入。")
        found = search_gtf(
            GTF_PATH,
            terms=[],
            column_filters=[(key, value)],
            attr_filters=[],
            ignore_case=False,
            limit=ask_limit(default=20),
        )
    else:
        raise SystemExit("請選 1 到 10。")

    if found == 0:
        raise SystemExit("No matches found.")


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Search Homo_sapiens.GRCh38.108.gtf by text, GTF column, or attribute. "
            "Run without arguments for interactive mode."
        )
    )
    parser.add_argument(
        "terms",
        nargs="*",
        help="Text that must appear in the GTF line. Multiple terms mean AND search.",
    )
    parser.add_argument(
        "--attr",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help='Match a GTF attribute exactly, e.g. --attr gene_name=ATAD3B.',
    )
    parser.add_argument(
        "--where",
        action="append",
        default=[],
        metavar="COLUMN=VALUE",
        help="Match a GTF column exactly, e.g. --where feature=gene --where seqname=1.",
    )
    parser.add_argument(
        "-i",
        "--ignore-case",
        action="store_true",
        help="Ignore case when matching text terms.",
    )
    parser.add_argument(
        "-n",
        "--limit",
        type=int,
        default=20,
        help="Maximum matches to print. Use 0 for no limit. Default: 20.",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=GTF_PATH,
        help=f"GTF file path. Default: {GTF_PATH}",
    )
    return parser


def main():
    if len(sys.argv) == 1:
        interactive_search()
        return

    parser = build_parser()
    args = parser.parse_args()

    if not args.terms and not args.attr and not args.where:
        parser.print_help()
        return

    column_filters = parse_key_value(args.where, "--where")
    attr_filters = parse_key_value(args.attr, "--attr")
    found = search_gtf(
        args.gtf,
        args.terms,
        column_filters,
        attr_filters,
        args.ignore_case,
        args.limit,
    )

    if found == 0:
        raise SystemExit("No matches found.")


if __name__ == "__main__":
    main()

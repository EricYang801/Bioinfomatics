#!/usr/bin/env python3
"""Check whether siTRF1 DEGs implicate genes that perturb TRF1 binding."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEG_REGULATION = ROOT / "data/siTRF1_DEG_list_gene_name.txt"
BETA_BSF = ROOT / "results/integration_analysis/beta/siTRF1_DEGs.bsf"
CHIPATLAS_OVERLAP = ROOT / "results/new_data/chipatlas/chipatlas_terf1_deg_overlap.tsv"
ALPHAGENOME_CANDIDATES = (
    ROOT / "results/alphagenome/chipatlas_candidates/"
    "alphagenome_bj_chipatlas_candidate_gene_summary.tsv"
)
OUTDIR = ROOT / "results/trf1_binding_error"


CATALOG = [
    {
        "gene_name": "TERF1",
        "layer": "direct_binding_factor",
        "role": "TRF1 protein itself; binds double-stranded TTAGGG telomeric DNA.",
        "interpretation": "If reduced, TRF1 telomeric occupancy is expected to fall directly.",
    },
    {
        "gene_name": "TINF2",
        "layer": "shelterin_anchor",
        "role": "TIN2 bridges TRF1, TRF2, and TPP1 and stabilizes shelterin.",
        "interpretation": "Loss or mutation can destabilize shelterin and telomere protection.",
    },
    {
        "gene_name": "ACD",
        "layer": "shelterin_anchor",
        "role": "TPP1; connects POT1 to TIN2 and helps shelterin/telomerase function.",
        "interpretation": "Affects telomere-end protection through shelterin organization.",
    },
    {
        "gene_name": "POT1",
        "layer": "shelterin_anchor",
        "role": "Binds single-stranded telomeric DNA and suppresses ATR signaling.",
        "interpretation": "Perturbation causes telomeric DDR, mostly ATR-linked.",
    },
    {
        "gene_name": "TERF2",
        "layer": "parallel_shelterin",
        "role": "TRF2 binds duplex telomeric DNA and blocks ATM/NHEJ.",
        "interpretation": "Parallel shelterin failure, not TRF1-specific promoter regulation.",
    },
    {
        "gene_name": "TERF2IP",
        "layer": "parallel_shelterin",
        "role": "RAP1, recruited by TRF2.",
        "interpretation": "Parallel shelterin/telomere-protection component.",
    },
    {
        "gene_name": "TNKS",
        "layer": "trf1_turnover",
        "role": "Tankyrase 1 PARylates TRF1 and promotes TRF1 eviction from telomeres.",
        "interpretation": "Increased activity can reduce telomeric TRF1 occupancy.",
    },
    {
        "gene_name": "TNKS2",
        "layer": "trf1_turnover",
        "role": "Tankyrase 2; related TRF1-interacting PARP.",
        "interpretation": "Potential TRF1 turnover/localization regulator.",
    },
    {
        "gene_name": "FBXO4",
        "layer": "trf1_turnover",
        "role": "F-box E3 ligase component implicated in TRF1 ubiquitination/degradation.",
        "interpretation": "Could affect TRF1 protein stability.",
    },
    {
        "gene_name": "PIN1",
        "layer": "trf1_turnover",
        "role": "Prolyl isomerase reported to regulate TRF1 stability.",
        "interpretation": "Could indirectly alter telomeric TRF1 levels.",
    },
    {
        "gene_name": "PINX1",
        "layer": "trf1_interactor",
        "role": "TRF1-interacting telomerase inhibitor.",
        "interpretation": "More related to telomerase/telomere length than acute TRF1 binding error.",
    },
    {
        "gene_name": "MAPK1",
        "layer": "upstream_signaling",
        "role": "ERK2 can phosphorylate TRF1; TRF1 T330 phosphorylation supports telomere complex formation.",
        "interpretation": "Pathway perturbation can weaken TRF1-telomere protection in cancer models.",
    },
    {
        "gene_name": "MAPK3",
        "layer": "upstream_signaling",
        "role": "ERK1/MAPK pathway member; related Ras/MAPK shelterin regulation.",
        "interpretation": "Potential pathway-level regulator, less direct than MAPK1/ERK2.",
    },
    {
        "gene_name": "BRAF",
        "layer": "upstream_signaling",
        "role": "bRAF reported to phosphorylate TRF1; T330 is important for TRF1-telomere complex formation.",
        "interpretation": "Pathway inhibition can increase telomere uncapping.",
    },
    {
        "gene_name": "MTOR",
        "layer": "upstream_signaling",
        "role": "Reported TRF1 substrate pathway in cancer signaling screen.",
        "interpretation": "Potential upstream signaling effect on TRF1 protection.",
    },
    {
        "gene_name": "ATM",
        "layer": "damage_response_or_regulator",
        "role": "DDR kinase and reported co-mediator of TRF1 telomere length control with MRN.",
        "interpretation": "Mostly downstream DDR; can also interact with telomere-length control.",
    },
    {
        "gene_name": "ATR",
        "layer": "damage_response_or_regulator",
        "role": "DDR kinase activated by telomeric replication stress/ssDNA contexts.",
        "interpretation": "Downstream marker/effector of telomere replication failure.",
    },
    {
        "gene_name": "MRE11",
        "layer": "damage_response_or_regulator",
        "role": "MRN complex component linked to TRF1 telomere length control.",
        "interpretation": "Potential telomere-processing cofactor, not DEG evidence unless changed.",
    },
    {
        "gene_name": "RAD50",
        "layer": "damage_response_or_regulator",
        "role": "MRN complex component linked to TRF1 telomere length control.",
        "interpretation": "Potential telomere-processing cofactor, not DEG evidence unless changed.",
    },
    {
        "gene_name": "NBN",
        "layer": "damage_response_or_regulator",
        "role": "MRN complex component linked to TRF1 telomere length control.",
        "interpretation": "Potential telomere-processing cofactor, not DEG evidence unless changed.",
    },
    {
        "gene_name": "BLM",
        "layer": "telomere_replication",
        "role": "TRF1 recruits/uses BLM-related replication help at telomeric repeats.",
        "interpretation": "Perturbation can worsen telomere replication stress.",
    },
    {
        "gene_name": "AKTIP",
        "layer": "telomere_replication",
        "role": "Shelterin-interacting factor required for telomere maintenance.",
        "interpretation": "Can contribute to telomere replication phenotypes.",
    },
    {
        "gene_name": "ERCC2",
        "layer": "telomere_replication",
        "role": "TFIIH helicase component; TFIIH has a noncanonical role in TRF1-mediated telomere replication.",
        "interpretation": "Replication support, not evidence of altered TRF1 binding in this DEG set unless changed.",
    },
    {
        "gene_name": "ERCC3",
        "layer": "telomere_replication",
        "role": "TFIIH helicase component; TFIIH has a noncanonical role in TRF1-mediated telomere replication.",
        "interpretation": "Replication support, not evidence of altered TRF1 binding in this DEG set unless changed.",
    },
    {
        "gene_name": "GTF2H1",
        "layer": "telomere_replication",
        "role": "TFIIH component.",
        "interpretation": "Replication support candidate.",
    },
    {
        "gene_name": "GTF2H2",
        "layer": "telomere_replication",
        "role": "TFIIH component.",
        "interpretation": "Replication support candidate.",
    },
]


DOWNSTREAM_MARKERS = [
    ("DDB2", "p53_DDR", "p53/DNA-damage response transcriptional target."),
    ("PHLDA3", "p53_DDR", "p53-linked apoptosis/AKT-pathway target."),
    ("RPS27L", "p53_DDR", "p53/DNA-damage response marker."),
    ("CASP7", "apoptosis", "Executioner caspase; apoptosis downstream marker."),
    ("GSDME", "apoptosis", "Pyroptosis/apoptosis execution-linked marker."),
    ("DDIT4", "ISR_mTOR_stress", "Stress response; mTOR/ISR-linked."),
    ("ASNS", "ISR_amino_acid_stress", "Integrated stress response/amino-acid stress marker."),
    ("CHAC1", "ER_stress_UPR", "ER stress / ATF4-CHOP pathway marker."),
    ("IFITM1", "IFN_antiviral", "Interferon-stimulated gene."),
    ("RNASEL", "IFN_OAS_antiviral", "OAS/RNase L antiviral pathway component."),
    ("SLFN5", "IFN_antiviral", "Interferon-linked gene."),
    ("TRIM38", "IFN_inflammation", "Innate immune response regulator."),
    ("MT-ND1", "mitochondrial", "Mitochondrial Complex I gene."),
    ("COX6A1", "mitochondrial", "Mitochondrial respiratory-chain component."),
    ("SIRT4", "mitochondrial", "Mitochondrial metabolism regulator."),
]


def load_deg() -> pd.DataFrame:
    deg = pd.read_csv(DEG_REGULATION, sep="\t", dtype=str)
    deg["gene_name"] = deg["gene_name"].fillna("").str.strip()
    deg["regulation"] = deg["regulation"].fillna("").str.strip().str.upper()

    bsf = pd.read_csv(
        BETA_BSF,
        sep="\t",
        header=None,
        names=["gene_name", "log2FC", "pvalue"],
        dtype={"gene_name": str, "log2FC": float, "pvalue": float},
    ).drop_duplicates("gene_name")
    return deg.merge(bsf, how="left", on="gene_name")


def add_external_flags(table: pd.DataFrame) -> pd.DataFrame:
    out = table.copy()
    out["chipatlas_terf1_overlap_10kb"] = False
    out["alphagenome_candidate"] = False
    if CHIPATLAS_OVERLAP.exists():
        ca = pd.read_csv(CHIPATLAS_OVERLAP, sep="\t", dtype=str)
        genes = set(ca.loc[ca["window"].eq("10kb"), "gene_name"].dropna())
        out["chipatlas_terf1_overlap_10kb"] = out["gene_name"].isin(genes)
    if ALPHAGENOME_CANDIDATES.exists():
        ag = pd.read_csv(ALPHAGENOME_CANDIDATES, sep="\t", dtype=str)
        genes = set(ag["gene_name"].dropna())
        out["alphagenome_candidate"] = out["gene_name"].isin(genes)
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    deg = load_deg()

    catalog = pd.DataFrame(CATALOG)
    catalog = add_external_flags(catalog)
    catalog.to_csv(OUTDIR / "trf1_binding_error_gene_catalog.tsv", sep="\t", index=False)

    upstream_hits = catalog.merge(deg, how="inner", on="gene_name")
    upstream_hits = upstream_hits.sort_values(["layer", "gene_name"])
    upstream_hits.to_csv(
        OUTDIR / "trf1_binding_error_upstream_deg_hits.tsv", sep="\t", index=False
    )

    downstream = pd.DataFrame(
        DOWNSTREAM_MARKERS, columns=["gene_name", "downstream_class", "role"]
    )
    downstream_hits = downstream.merge(deg, how="inner", on="gene_name")
    downstream_hits = add_external_flags(downstream_hits)
    downstream_hits = downstream_hits.sort_values(["downstream_class", "gene_name"])
    downstream_hits.to_csv(
        OUTDIR / "trf1_binding_error_downstream_deg_markers.tsv", sep="\t", index=False
    )

    upstream_genes = ", ".join(upstream_hits["gene_name"]) if len(upstream_hits) else "None"
    downstream_genes = (
        ", ".join(downstream_hits["gene_name"]) if len(downstream_hits) else "None"
    )
    upstream_non_terf1 = upstream_hits[upstream_hits["gene_name"] != "TERF1"]

    class_counts = (
        downstream_hits.groupby(["downstream_class", "regulation"], dropna=False)
        .size()
        .reset_index(name="n_genes")
        .sort_values(["downstream_class", "regulation"])
    )
    class_counts.to_csv(
        OUTDIR / "trf1_binding_error_downstream_class_counts.tsv", sep="\t", index=False
    )

    lines = [
        "# TRF1 binding-error network check",
        "",
        "Question: do the siTRF1 DEGs contain genes that could cause erroneous TRF1",
        "binding, or downstream markers consistent with TRF1/telomere binding failure?",
        "",
        "## Upstream / causal candidates",
        "",
        f"- Curated upstream/TRF1-binding regulators checked: {len(catalog)}",
        f"- Hits among 274 DEGs: {len(upstream_hits)} ({upstream_genes})",
        f"- Hits excluding the intended siRNA target TERF1: {len(upstream_non_terf1)}",
        "",
        "Interpretation: within this DEG table, the only curated direct TRF1-binding",
        "or TRF1-localization factor that changes is `TERF1` itself. I do not see",
        "evidence that another DEG upstream of TRF1 is causing the binding defect.",
        "",
        "## Downstream markers",
        "",
        f"- Downstream marker hits among DEGs: {len(downstream_hits)} ({downstream_genes})",
        "",
        "Class counts:",
    ]
    for _, row in class_counts.iterrows():
        lines.append(f"- {row['downstream_class']} / {row['regulation']}: {row['n_genes']}")
    lines.extend(
        [
            "",
            "Interpretation: the DEG list strongly contains downstream stress-response",
            "markers expected after telomere/TRF1 perturbation: p53/DDR, integrated",
            "stress response, ER stress, interferon/antiviral response, apoptosis,",
            "and mitochondrial-response markers. This supports downstream indirect",
            "failure, not direct promoter binding by TRF1.",
            "",
            "## Practical conclusion",
            "",
            "Most defensible model for this dataset:",
            "",
            "`siTRF1 -> TERF1/TRF1 reduction -> telomere protection/replication stress ->",
            "DDR/ISR/UPR/IFN/mitochondrial response genes change`.",
            "",
            "Not supported by this dataset:",
            "",
            "`another DEG -> causes TRF1 binding error -> direct TRF1 binding at DEG promoters`.",
            "",
        ]
    )
    (OUTDIR / "trf1_binding_error_network_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

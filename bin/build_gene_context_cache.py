#!/usr/bin/env python3
"""
build_gene_context_cache.py
===========================
One-time helper to build or refresh:
    assets/gene_annotations/somatic_gene_context.tsv

Sources:
  1. COSMIC Cancer Gene Census (CGC) TSV  — optional, pass via --cosmic_cgc
     Fields used: Gene Symbol, Role in Cancer, Tier
  2. NCBI gene summaries via MyGene.info  — free REST API, no auth required
     Fields used: summary

Usage
-----
    python3 bin/build_gene_context_cache.py \
        --genes  assets/gene_panels/somatic_snv_panel.txt \
        [--cosmic_cgc  /path/to/Census_all.tsv] \
        [--out   assets/gene_annotations/somatic_gene_context.tsv]

If --cosmic_cgc is not provided the script uses a built-in fallback dict for
the 60 standard panel genes and warns for any unknown genes.

The output TSV is checked into the repository pre-populated, so end-users do
not need to run this script unless they want to refresh gene summaries or
update the panel.
"""
import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.parse

# ── built-in fallback for CGC role / tier  ──────────────────────────────────
# Format: GENE -> (role, tier, pathway)
# role values used in consequence logic: TSG | Oncogene | Both | Unknown
FALLBACK_CGC: dict[str, tuple[str, str, str]] = {
    "APC":      ("TSG",      "1", "WNT Signaling"),
    "AR":       ("Oncogene", "1", "Hormone Signaling"),
    "ARAF":     ("Oncogene", "2", "RTK/RAS/MAPK"),
    "ATM":      ("TSG",      "1", "DNA Damage Response / Homologous Recombination"),
    "ATRX":     ("TSG",      "1", "Chromatin Remodeling"),
    "AXIN1":    ("TSG",      "1", "WNT Signaling"),
    "BAP1":     ("TSG",      "1", "Chromatin Remodeling"),
    "BRAF":     ("Oncogene", "1", "RTK/RAS/MAPK"),
    "BRCA1":    ("TSG",      "1", "DNA Damage Response / Homologous Recombination"),
    "BRCA2":    ("TSG",      "1", "DNA Damage Response / Homologous Recombination"),
    "CDH1":     ("TSG",      "1", "Cell Adhesion"),
    "CDKN2A":   ("TSG",      "1", "Cell Cycle"),
    "CTNNB1":   ("Oncogene", "1", "WNT Signaling"),
    "DICER1":   ("TSG",      "1", "RNA Processing"),
    "DNMT3A":   ("TSG",      "1", "Epigenetic Regulation"),
    "EGFR":     ("Oncogene", "1", "RTK/RAS/MAPK"),
    "ERBB2":    ("Oncogene", "1", "RTK/RAS/MAPK"),
    "EZH2":     ("Both",     "1", "Epigenetic Regulation"),
    "FBXW7":    ("TSG",      "1", "Cell Cycle"),
    "FGFR1":    ("Oncogene", "1", "RTK/RAS/MAPK"),
    "FGFR2":    ("Oncogene", "1", "RTK/RAS/MAPK"),
    "FGFR3":    ("Oncogene", "1", "RTK/RAS/MAPK"),
    "FLT3":     ("Oncogene", "1", "RTK/RAS/MAPK"),
    "GNAS":     ("Oncogene", "1", "G Protein Signaling"),
    "HNF1A":    ("TSG",      "2", "Transcription Factor"),
    "HRAS":     ("Oncogene", "1", "RTK/RAS/MAPK"),
    "IDH1":     ("Both",     "1", "Metabolic / Epigenetic"),
    "IDH2":     ("Both",     "1", "Metabolic / Epigenetic"),
    "JAK2":     ("Oncogene", "1", "JAK/STAT Signaling"),
    "KDR":      ("Oncogene", "2", "RTK/RAS/MAPK"),
    "KIT":      ("Oncogene", "1", "RTK/RAS/MAPK"),
    "KRAS":     ("Oncogene", "1", "RTK/RAS/MAPK"),
    "MAP2K1":   ("Oncogene", "2", "RTK/RAS/MAPK"),
    "MET":      ("Oncogene", "1", "RTK/RAS/MAPK"),
    "MLH1":     ("TSG",      "1", "Mismatch Repair"),
    "MSH2":     ("TSG",      "1", "Mismatch Repair"),
    "MSH6":     ("TSG",      "1", "Mismatch Repair"),
    "MTOR":     ("Oncogene", "1", "PI3K/AKT/mTOR"),
    "MYC":      ("Oncogene", "1", "Transcription Factor / Cell Cycle"),
    "NF1":      ("TSG",      "1", "RTK/RAS/MAPK"),
    "NOTCH1":   ("Both",     "1", "NOTCH Signaling"),
    "NRAS":     ("Oncogene", "1", "RTK/RAS/MAPK"),
    "PALB2":    ("TSG",      "1", "DNA Damage Response / Homologous Recombination"),
    "PDGFRA":   ("Oncogene", "1", "RTK/RAS/MAPK"),
    "PIK3CA":   ("Oncogene", "1", "PI3K/AKT/mTOR"),
    "PIK3R1":   ("TSG",      "1", "PI3K/AKT/mTOR"),
    "PMS2":     ("TSG",      "1", "Mismatch Repair"),
    "PTEN":     ("TSG",      "1", "PI3K/AKT/mTOR"),
    "RB1":      ("TSG",      "1", "Cell Cycle"),
    "RET":      ("Oncogene", "1", "RTK/RAS/MAPK"),
    "RHOA":     ("Both",     "2", "RHO/RAC Signaling"),
    "RNF43":    ("TSG",      "1", "WNT Signaling"),
    "SMAD4":    ("TSG",      "1", "TGF-β Signaling"),
    "SMARCB1":  ("TSG",      "1", "Chromatin Remodeling"),
    "STK11":    ("TSG",      "1", "LKB1/AMPK Signaling"),
    "TERT":     ("Oncogene", "1", "Telomere Maintenance"),
    "TP53":     ("TSG",      "1", "Cell Cycle / Apoptosis"),
    "TSC1":     ("TSG",      "1", "PI3K/AKT/mTOR"),
    "TSC2":     ("TSG",      "1", "PI3K/AKT/mTOR"),
    "VHL":      ("TSG",      "1", "Hypoxia / HIF Signaling"),
}


def parse_args():
    p = argparse.ArgumentParser(description="Build somatic_gene_context.tsv")
    p.add_argument("--genes",      required=True,
                   help="Gene list file (one gene per line, # comments allowed)")
    p.add_argument("--cosmic_cgc", default=None,
                   help="COSMIC Cancer Gene Census TSV (Census_all.tsv). "
                        "If not provided, fallback dict is used.")
    p.add_argument("--out",        default="assets/gene_annotations/somatic_gene_context.tsv",
                   help="Output TSV path")
    p.add_argument("--delay",      type=float, default=0.3,
                   help="Seconds to wait between MyGene.info requests (default 0.3)")
    return p.parse_args()


def load_gene_list(path: str) -> list[str]:
    genes = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.append(line.split()[0])
    return genes


def load_cosmic_cgc(path: str) -> dict[str, tuple[str, str]]:
    """Return {gene_symbol: (role, tier)} from COSMIC CGC TSV."""
    cgc = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("Gene Symbol") or row.get("Gene name") or "").strip()
            role = (row.get("Role in Cancer") or "Unknown").strip()
            tier = str(row.get("Tier") or "").strip()
            if gene:
                # Normalise role to TSG / Oncogene / Both
                role_norm = role
                if "TSG" in role.upper() and "oncogene" in role.lower():
                    role_norm = "Both"
                elif "TSG" in role.upper():
                    role_norm = "TSG"
                elif "oncogene" in role.lower():
                    role_norm = "Oncogene"
                cgc[gene] = (role_norm, tier)
    return cgc


def fetch_ncbi_summary(gene_symbol: str) -> str:
    """Query MyGene.info for the NCBI gene summary."""
    url = (
        "https://mygene.info/v3/query?"
        + urllib.parse.urlencode({
            "q": f"symbol:{gene_symbol}",
            "species": "human",
            "fields": "summary",
            "size": "1",
        })
    )
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
        hits = data.get("hits", [])
        if hits:
            return hits[0].get("summary", "").replace("\n", " ").replace("\t", " ")
    except Exception as exc:
        print(f"  WARNING: MyGene.info lookup failed for {gene_symbol}: {exc}", file=sys.stderr)
    return ""


def main():
    args = parse_args()
    genes = load_gene_list(args.genes)
    print(f"Loaded {len(genes)} genes from {args.genes}")

    # CGC data
    cgc_data: dict[str, tuple[str, str]] = {}
    if args.cosmic_cgc:
        cgc_data = load_cosmic_cgc(args.cosmic_cgc)
        print(f"Loaded COSMIC CGC for {len(cgc_data)} genes from {args.cosmic_cgc}")
    else:
        print("No --cosmic_cgc provided; using built-in fallback dict.")

    rows = []
    for gene in genes:
        # Role + tier
        if gene in cgc_data:
            role, tier = cgc_data[gene]
            pathway = FALLBACK_CGC.get(gene, ("Unknown", "", "Unknown"))[2]
        elif gene in FALLBACK_CGC:
            role, tier, pathway = FALLBACK_CGC[gene]
        else:
            print(f"  WARNING: {gene} not in CGC or fallback — marking Unknown", file=sys.stderr)
            role, tier, pathway = "Unknown", "", "Unknown"

        # NCBI summary
        print(f"  Fetching summary for {gene} ...", end=" ", flush=True)
        summary = fetch_ncbi_summary(gene)
        if summary:
            print("OK")
        else:
            print("(empty)")
        time.sleep(args.delay)

        rows.append({
            "GENE": gene,
            "ROLE": role,
            "CGC_TIER": tier,
            "PATHWAY": pathway,
            "GENE_SUMMARY": summary,
        })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["GENE", "ROLE", "CGC_TIER", "PATHWAY", "GENE_SUMMARY"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} genes to {args.out}")


if __name__ == "__main__":
    main()

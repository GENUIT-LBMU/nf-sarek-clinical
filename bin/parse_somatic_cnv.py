#!/usr/bin/env python3
"""
parse_somatic_cnv.py
====================
Interprets somatic CNV events from CNVkit and ASCAT outputs (post nf-core/sarek).

Outputs:
  <sample>_somatic_cnv_report.tsv  – master TSV with all events
  <sample>_somatic_cnv_report.txt  – human-readable summary for clinical report

CNVkit events (filtered to somatic gene panel):
  HOMDEL  CN=0
  HETDEL  CN=1
  AMP     CN >= amp_threshold (default 5)

ASCAT events (filtered to germline gene panel – LOH as second hit):
  cnLOH          nMinor=0, totCN=2
  HETDEL+LOH     nMinor=0, totCN=1
  LOH_major_amp  nMinor=0, totCN>=3
"""

import argparse
import os
import sys
import re
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_any(path: str) -> pd.DataFrame | None:
    """Tab-first, CSV fallback reader. Skips comment lines."""
    if not path or not os.path.isfile(path):
        return None
    for sep in ("\t", ","):
        try:
            return pd.read_csv(path, sep=sep, comment="#")
        except Exception:
            continue
    return None


def find_col(cols: list[str], keys: list[str]) -> str | None:
    low = {c.lower(): c for c in cols}
    for k in keys:
        if k.lower() in low:
            return low[k.lower()]
    return None


def load_gene_list(path: str) -> set[str]:
    genes: set[str] = set()
    if not path or not os.path.isfile(path):
        return genes
    with open(path) as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.add(g.upper())
    return genes


def norm_chr(s: str) -> str:
    s = str(s).replace("chr", "").replace("CHR", "")
    if s in ("23", "X", "x"): return "X"
    if s in ("24", "Y", "y"): return "Y"
    if s.upper() in ("M", "MT", "Mt", "m"): return "MT"
    return s


# ---------------------------------------------------------------------------
# Gene-field parsing
# ---------------------------------------------------------------------------

def split_genes_field(series: pd.Series) -> list[list[str]]:
    out = []
    for v in series.fillna("").astype(str):
        if not v or v in (".", "NA", "None"):
            out.append([])
            continue
        s = re.sub(r"[;|+]", ",", v)
        lst = [p.strip().upper() for p in s.split(",") if p.strip()]
        lst = [g for g in lst if g not in (".", "NA")]
        out.append(lst)
    return out


def read_gene_coords(path: str) -> dict[str, pd.DataFrame] | None:
    """
    Reads a BED/TSV with columns: chrom, start, end, gene/name/symbol.
    Returns a dict keyed by (norm) chrom → DataFrame with [chrom, start, end, gene].
    Handles multi-gene fields and explodes to one gene per row.
    """
    if not path or not os.path.isfile(path):
        return None

    # A plain BED has no header: if the first data line has numeric start/end,
    # it is an interval, not column names (header=0 would drop it).
    first = None
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                first = line.rstrip("\n").split("\t")
                break
    if first is None:
        return None
    has_header = not (len(first) > 2 and first[1].isdigit() and first[2].isdigit())

    try:
        df = pd.read_csv(path, sep="\t", comment="#",
                         header=0 if has_header else None)
    except Exception:
        return None

    cols = list(df.columns)
    low = {str(c).lower(): c for c in cols}
    chrom = (low.get("chromosome") or low.get("chrom") or low.get("#chrom")
             or cols[0] if cols else None)
    start = (low.get("start") or low.get("tx_start")
             or cols[1] if len(cols) > 1 else None)
    end   = (low.get("end") or low.get("tx_end")
             or cols[2] if len(cols) > 2 else None)
    gene  = (low.get("gene") or low.get("name") or low.get("symbol")
             or cols[3] if len(cols) > 3 else None)

    if any(c is None for c in (chrom, start, end, gene)):
        return None

    sub = df[[chrom, start, end, gene]].copy()
    sub.columns = ["chrom", "start", "end", "gene_raw"]

    def norm_tokens(s: str) -> list[str]:
        parts = re.split(r"[;,|/+]", str(s))
        toks = []
        for p in parts:
            t = p.strip().upper()
            if not t or t in (".", "NA", "NONE"):
                continue
            if "-" in t:
                left = t.split("-")[0]
                if left and left not in (".", "NA"):
                    toks.append(left)
            toks.append(t)
        return sorted(set(toks))

    sub["gene_tokens"] = sub["gene_raw"].apply(norm_tokens)
    sub = sub.explode("gene_tokens").rename(columns={"gene_tokens": "gene"})
    sub = sub.drop(columns=["gene_raw"])
    sub["chrom"] = sub["chrom"].map(norm_chr)
    sub["start"] = pd.to_numeric(sub["start"], errors="coerce")
    sub["end"] = pd.to_numeric(sub["end"], errors="coerce")
    sub = sub.dropna(subset=["start", "end"])
    return {k: v for k, v in sub.groupby("chrom")}


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def parse_cnvkit(cns_path: str, somatic_set: set[str], amp_threshold: int) -> list[dict]:
    cns = read_any(cns_path)
    if cns is None:
        sys.exit(f"ERROR: Cannot read CNVkit file: {cns_path}")

    cols = list(cns.columns)
    chrom_c = find_col(cols, ["chromosome", "chrom", "chr"])
    start_c = find_col(cols, ["start", "startpos"])
    end_c   = find_col(cols, ["end", "endpos"])
    log2_c  = find_col(cols, ["log2"])
    cn_c    = find_col(cols, ["cn", "copy_number"])
    gene_c  = find_col(cols, ["gene", "genes", "name", "symbol"])

    if not all([chrom_c, start_c, end_c, log2_c]) or cn_c is None:
        sys.exit("ERROR: Missing required columns in CNVkit .cns file")

    keep = [chrom_c, start_c, end_c, log2_c, cn_c] + ([gene_c] if gene_c else [])
    df = cns[keep].copy()
    df.columns = ["chrom", "start", "end", "log2", "cn"] + (["genes"] if gene_c else [])
    df["cn"]   = pd.to_numeric(df["cn"],   errors="coerce")
    df["log2"] = pd.to_numeric(df["log2"], errors="coerce")
    if "genes" not in df.columns:
        df["genes"] = ""
    df["chrom"] = df["chrom"].map(norm_chr)
    df["genes_list"] = split_genes_field(df["genes"])

    rows = []

    def add(row, etype, cn_val):
        for g in row["genes_list"]:
            if g in somatic_set:
                rows.append({
                    "event_type": etype,
                    "gene": g,
                    "source": "CNVkit",
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "cn": cn_val,
                    "log2": round(row["log2"], 3) if pd.notna(row["log2"]) else "",
                    "minor": "",
                    "loh_type": "",
                })

    for _, r in df[df["cn"] == 0].iterrows():
        add(r, "HOMDEL", 0)
    for _, r in df[df["cn"] == 1].iterrows():
        add(r, "HETDEL", 1)
    for _, r in df[df["cn"] >= amp_threshold].iterrows():
        add(r, "AMP", int(r["cn"]))

    return rows


def parse_ascat_loh(
    segments_path: str,
    metrics_path: str,
    germline_set: set[str],
    gene_coords: dict[str, pd.DataFrame] | None,
) -> tuple[list[dict], dict]:
    rows = []
    metrics = {}

    if metrics_path and os.path.isfile(metrics_path):
        met = read_any(metrics_path)
        if met is not None and not met.empty:
            metrics = met.iloc[0].to_dict()

    ascat_seg = read_any(segments_path) if segments_path else None
    if ascat_seg is None or ascat_seg.empty or gene_coords is None:
        return rows, metrics

    a_cols = list(ascat_seg.columns)
    achr   = find_col(a_cols, ["chromosome", "chrom", "chr"])
    astart = find_col(a_cols, ["start", "startpos", "startposition"])
    aend   = find_col(a_cols, ["end", "endpos", "endposition"])
    nmaj   = find_col(a_cols, ["nmajor", "nmaj", "nmaj1_a", "nmajor1"])
    nmin   = find_col(a_cols, ["nminor", "nmin", "nmin1_a", "nminor1"])

    if not all([achr, astart, aend, nmaj, nmin]):
        return rows, metrics

    seg = ascat_seg[[achr, astart, aend, nmaj, nmin]].copy()
    seg.columns = ["chrom", "start", "end", "nMajor", "nMinor"]
    seg["chrom"]  = seg["chrom"].map(norm_chr)
    seg["nMajor"] = pd.to_numeric(seg["nMajor"], errors="coerce")
    seg["nMinor"] = pd.to_numeric(seg["nMinor"], errors="coerce")
    seg["totCN"]  = seg["nMajor"] + seg["nMinor"]

    # Male tumors (ASCAT sex XY): chrX / chrY are hemizygous, nMinor = 0 is the
    # normal state and must not read as LOH (AR was reported as a second hit)
    sex = str(metrics.get("sex", "")).strip().upper()
    if not sex:
        print("WARNING: ASCAT metrics without sex; chrX/chrY LOH kept (may be hemizygous)", file=sys.stderr)

    for chrom, segs in seg.groupby("chrom"):
        if chrom not in gene_coords:
            continue
        if chrom in ("X", "Y") and sex == "XY":
            continue
        gsub = gene_coords[chrom]
        for _, ar in segs.iterrows():
            if ar["nMinor"] != 0:
                continue
            a_s, a_e = int(ar["start"]), int(ar["end"])
            ghit = gsub[(gsub["end"] >= a_s) & (gsub["start"] <= a_e)]
            if ghit.empty:
                continue
            genes = sorted(set(ghit["gene"].astype(str).str.upper()) & germline_set)
            if not genes:
                continue

            a_tot, a_nmin = int(ar["totCN"]), int(ar["nMinor"])
            if a_tot == 1:
                loh_type = "HETDEL+LOH"
            elif a_tot == 2:
                loh_type = "cnLOH"
            elif a_tot >= 3:
                loh_type = f"LOH_major_amp(CN={a_tot})"
            else:
                loh_type = "LOH_minor0"

            for g in genes:
                rows.append({
                    "event_type": "2nd_hit_LOH",
                    "gene": g,
                    "source": "ASCAT",
                    "chrom": str(chrom),
                    "start": a_s,
                    "end": a_e,
                    "cn": a_tot,
                    "log2": "",
                    "minor": a_nmin,
                    "loh_type": loh_type,
                })

    return rows, metrics


def build_report_txt(sample: str, rows: list[dict], metrics: dict, amp_threshold: int) -> list[str]:
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["event_type", "gene", "source", "cn", "loh_type"])

    homdel = sorted(df.loc[df["event_type"] == "HOMDEL", "gene"].unique())
    hetdel = sorted(df.loc[df["event_type"] == "HETDEL", "gene"].unique())
    amp    = sorted(df.loc[df["event_type"] == "AMP",    "gene"].unique())
    loh    = sorted(df.loc[df["event_type"] == "2nd_hit_LOH", "gene"].unique())

    amp_detail = df[df["event_type"] == "AMP"][["gene", "cn"]].drop_duplicates()
    amp_detail = amp_detail.sort_values(["cn", "gene"], ascending=[False, True])
    amp_cn = [f"{r['gene']}(CN={r['cn']})" for _, r in amp_detail.iterrows()]

    loh_cnloh  = sorted(df[(df["event_type"] == "2nd_hit_LOH") & (df["loh_type"] == "cnLOH")]["gene"].unique())
    loh_hetdel = sorted(df[(df["event_type"] == "2nd_hit_LOH") & (df["loh_type"] == "HETDEL+LOH")]["gene"].unique())
    loh_amp    = sorted(df[(df["event_type"] == "2nd_hit_LOH") & (df["loh_type"].str.contains("major_amp", na=False))]["gene"].unique())

    lines = [f"Sample: {sample}", "=" * 60]

    if metrics:
        lines += [
            "",
            "ASCAT METRICS:",
            f"  Tumor purity:          {metrics.get('purity', 'N/D')}",
            f"  Ploidy:                {metrics.get('ploidy', 'N/D')}",
            f"  Whole Genome Doubling: {metrics.get('WGD', 'N/D')}",
            f"  Goodness of fit:       {metrics.get('goodness_of_fit', 'N/D')}",
            f"  LOH fraction:          {metrics.get('LOH', 'N/D')}",
            f"  Genomic Instability:   {metrics.get('GI', 'N/D')}",
        ]

    lines += [
        "", "=" * 60,
        f"CNVkit EVENTS (somatic 38-gene panel)",
        "=" * 60,
        "",
        f"Amplifications (CN>={amp_threshold}):",
        "  " + (", ".join(amp_cn) if amp_cn else "N/D"),
        "",
        "Homozygous deletions (CN=0):",
        "  " + (", ".join(homdel) if homdel else "N/D"),
        "",
        "Heterozygous deletions (CN=1):",
        "  " + (", ".join(hetdel) if hetdel else "N/D"),
        "", "=" * 60,
        "ASCAT EVENTS – LOH (germline 90-gene panel)",
        "=" * 60,
        "",
        "Second hit (LOH) – all types:",
        "  " + (", ".join(loh) if loh else "N/D"),
        "",
        "  LOH breakdown by type:",
        f"    cnLOH (CN=2, nMinor=0):       {', '.join(loh_cnloh) if loh_cnloh else 'N/D'}",
        f"    HETDEL+LOH (CN=1):            {', '.join(loh_hetdel) if loh_hetdel else 'N/D'}",
        f"    LOH with amplification:       {', '.join(loh_amp) if loh_amp else 'N/D'}",
        "",
    ]
    return lines


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Interpret somatic CNV events from CNVkit and ASCAT outputs."
    )
    ap.add_argument("--cns",            required=True,  help="CNVkit .somatic.call.cns or .call.cns file")
    ap.add_argument("--ascat_segments", default="",     help="ASCAT segments file (*segments.txt)")
    ap.add_argument("--ascat_metrics",  default="",     help="ASCAT metrics file (*metrics.txt)")
    ap.add_argument("--genecoords",     default="",     help="BED/TSV with gene coordinates (chrom,start,end,gene)")
    ap.add_argument("--germline",       required=True,  help="Germline gene list (90 genes, for LOH)")
    ap.add_argument("--somatic_genes",  required=True,  help="Somatic gene list (38 genes, for CNVkit events)")
    ap.add_argument("--sample",         required=True,  help="Sample name")
    ap.add_argument("--out_tsv",        required=True,  help="Output TSV file")
    ap.add_argument("--out_txt",        required=True,  help="Output TXT summary file")
    ap.add_argument("--amp_threshold",  type=int, default=5, help="Amplification CN threshold (default: 5)")
    args = ap.parse_args()

    somatic_set  = load_gene_list(args.somatic_genes)
    germline_set = load_gene_list(args.germline)
    gene_coords  = read_gene_coords(args.genecoords) if args.genecoords else None

    cnvkit_rows = parse_cnvkit(args.cns, somatic_set, args.amp_threshold)
    ascat_rows, metrics = parse_ascat_loh(
        args.ascat_segments, args.ascat_metrics, germline_set, gene_coords
    )

    all_rows = cnvkit_rows + ascat_rows
    for r in all_rows:
        r["sample"] = args.sample

    out_df = pd.DataFrame(all_rows, columns=[
        "sample", "event_type", "gene", "source", "chrom", "start", "end",
        "cn", "log2", "minor", "loh_type",
    ]).drop_duplicates()

    # Append ASCAT summary metrics as columns
    for k in ["purity", "ploidy", "goodness_of_fit", "tumour_mapd", "normal_mapd", "LOH", "WGD", "GI"]:
        out_df[k] = metrics.get(k, np.nan)

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out_df.to_csv(args.out_tsv, sep="\t", index=False)

    txt_lines = build_report_txt(args.sample, all_rows, metrics, args.amp_threshold)
    with open(args.out_txt, "w") as fh:
        fh.write("\n".join(txt_lines))

    print(f"[OK] TSV  → {args.out_tsv}")
    print(f"[OK] TXT  → {args.out_txt}")
    print(f"[OK] CNVkit events: {len(cnvkit_rows)}  |  ASCAT LOH events: {len(ascat_rows)}")


if __name__ == "__main__":
    main()

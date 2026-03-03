#!/usr/bin/env python3
"""
generate_signature_plots.py
============================
Generates two sets of publication-quality PNG plots for the clinical report:

  1. SBS96 matrix profile  — 96-channel mutational spectrum (reconstructed from
                              Activity × COSMIC reference signature profiles)
  2. ID83  matrix profile  — 83-channel indel spectrum (same approach)
  3. Etiology activity bar chart for SBS96 — one bar per active signature,
     coloured and labelled by COSMIC v3.4 aetiology
  4. Etiology activity bar chart for ID     — same

All PNGs are saved inside the respective Assignment_Solution/Activities/
directory, next to the existing Activities.txt files, so the final_report
module can locate them by the same act_dir path it already knows.

Usage
-----
    python3 generate_signature_plots.py \\
        --sbs96_dir  <path_to_{run}_SBS96>   \\
        --id_dir     <path_to_{run}_ID>       \\
        --sample     <sample_name>
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── COSMIC v3.4 signature aetiology ────────────────────────────────────────
# Each entry: (aetiology_label, bar_colour)
SBS_ETIOLOGY = {
    "SBS1":  ("Clock-like",          "#4e79a7"),
    "SBS2":  ("APOBEC",              "#f28e2b"),
    "SBS3":  ("HRD (BRCA1/2)",       "#e15759"),
    "SBS4":  ("Tobacco smoking",     "#9c755f"),
    "SBS5":  ("Clock-like",          "#4e79a7"),
    "SBS6":  ("MMR deficiency",      "#b07aa1"),
    "SBS7a": ("UV radiation",        "#edc948"),
    "SBS7b": ("UV radiation",        "#edc948"),
    "SBS7c": ("UV radiation",        "#edc948"),
    "SBS7d": ("UV radiation",        "#edc948"),
    "SBS8":  ("HRD-associated",      "#c0392b"),
    "SBS9":  ("Somatic hypermut.",   "#bab0ac"),
    "SBS10a":("POLE/POLD1",          "#ff9da7"),
    "SBS10b":("POLE/POLD1",          "#ff9da7"),
    "SBS10c":("POLE/POLD1",          "#ff9da7"),
    "SBS10d":("POLE/POLD1",          "#ff9da7"),
    "SBS11": ("Alkylating agents",   "#59a14f"),
    "SBS12": ("Unknown",             "#c7c7c7"),
    "SBS13": ("APOBEC",              "#f28e2b"),
    "SBS14": ("MMR def. + POLE",     "#7b2d8b"),
    "SBS15": ("MMR deficiency",      "#b07aa1"),
    "SBS16": ("Alcohol",             "#ff7043"),
    "SBS17a":("Unknown (ROS?)",      "#c7c7c7"),
    "SBS17b":("Unknown (ROS?)",      "#c7c7c7"),
    "SBS18": ("ROS / MUTYH",         "#76b7b2"),
    "SBS19": ("Unknown",             "#c7c7c7"),
    "SBS20": ("MMR deficiency",      "#b07aa1"),
    "SBS21": ("MMR deficiency",      "#b07aa1"),
    "SBS22a":("Aristolochic acid",   "#8b4513"),
    "SBS22b":("Aristolochic acid",   "#8b4513"),
    "SBS23": ("Unknown",             "#c7c7c7"),
    "SBS24": ("Aflatoxin",           "#a0522d"),
    "SBS25": ("Chemotherapy",        "#2e8b57"),
    "SBS26": ("MMR deficiency",      "#b07aa1"),
    "SBS27": ("Possible artifact",   "#eeeeee"),
    "SBS28": ("POLE/POLD1",          "#ff9da7"),
    "SBS29": ("Tobacco chewing",     "#6d4c41"),
    "SBS30": ("BER def. (NTHL1)",    "#00acc1"),
    "SBS31": ("Platinum agents",     "#1b5e20"),
    "SBS32": ("Azathioprine",        "#558b2f"),
    "SBS33": ("Unknown",             "#c7c7c7"),
    "SBS34": ("Unknown",             "#c7c7c7"),
    "SBS35": ("Platinum agents",     "#1b5e20"),
    "SBS36": ("BER def. (MUTYH)",    "#0097a7"),
    "SBS37": ("Unknown",             "#c7c7c7"),
    "SBS38": ("UV (indirect)",       "#f9a825"),
    "SBS39": ("Unknown",             "#c7c7c7"),
    "SBS40a":("Unknown",             "#c7c7c7"),
    "SBS40b":("Unknown",             "#c7c7c7"),
    "SBS40c":("Unknown",             "#c7c7c7"),
    "SBS41": ("Unknown",             "#c7c7c7"),
    "SBS42": ("Haloalkane",          "#a1887f"),
    "SBS44": ("MMR deficiency",      "#b07aa1"),
    "SBS84": ("Somatic hypermut.",   "#bab0ac"),
    "SBS85": ("Somatic hypermut.",   "#bab0ac"),
    "SBS86": ("Chemotherapy",        "#2e8b57"),
    "SBS87": ("Thiopurine",          "#388e3c"),
    "SBS88": ("Colibactin (E. coli)","#8bc34a"),
    "SBS89": ("Unknown",             "#c7c7c7"),
    "SBS90": ("Duocarmycin",         "#43a047"),
    "SBS91": ("Unknown",             "#c7c7c7"),
    "SBS92": ("Tobacco (bladder)",   "#5d4037"),
    "SBS93": ("Unknown",             "#c7c7c7"),
    "SBS94": ("Unknown",             "#c7c7c7"),
    "SBS95": ("Unknown",             "#c7c7c7"),
    "SBS96": ("Unknown",             "#c7c7c7"),
    "SBS97": ("Unknown",             "#c7c7c7"),
    "SBS98": ("Unknown",             "#c7c7c7"),
    "SBS99": ("Unknown",             "#c7c7c7"),
}

ID_ETIOLOGY = {
    "ID1":  ("Clock-like / Replication slippage", "#4e79a7"),
    "ID2":  ("Clock-like / Replication slippage", "#4e79a7"),
    "ID3":  ("Tobacco smoking",                   "#9c755f"),
    "ID4":  ("Unknown",                           "#c7c7c7"),
    "ID5":  ("Unknown",                           "#c7c7c7"),
    "ID6":  ("MMR deficiency",                    "#b07aa1"),
    "ID7":  ("MMR deficiency",                    "#b07aa1"),
    "ID8":  ("Unknown",                           "#c7c7c7"),
    "ID9":  ("Unknown",                           "#c7c7c7"),
    "ID10": ("Unknown",                           "#c7c7c7"),
    "ID11": ("Platinum agents",                   "#1b5e20"),
    "ID12": ("Unknown",                           "#c7c7c7"),
    "ID13": ("APOBEC",                            "#f28e2b"),
    "ID14": ("POLE/POLD1",                        "#ff9da7"),
    "ID15": ("Unknown",                           "#c7c7c7"),
    "ID16": ("Unknown",                           "#c7c7c7"),
    "ID17": ("Unknown",                           "#c7c7c7"),
    "ID18": ("Unknown",                           "#c7c7c7"),
    "ID19": ("Unknown",                           "#c7c7c7"),
    "ID20": ("Unknown",                           "#c7c7c7"),
    "ID21": ("Unknown",                           "#c7c7c7"),
    "ID22": ("Unknown",                           "#c7c7c7"),
    "ID23": ("Unknown",                           "#c7c7c7"),
}

# ── SBS96 COSMIC colour scheme ──────────────────────────────────────────────
SBS_COLORS = {
    "C>A": "#1f97f5",
    "C>G": "#231f20",
    "C>T": "#e52221",
    "T>A": "#ccc9c9",
    "T>C": "#97d54c",
    "T>G": "#e8aac2",
}

# ── ID83 colour scheme (by indel category) ──────────────────────────────────
# Pattern: {length}:{Del|Ins}:{C|T|R|M}:{n}
# Categories:  1bp C Del, 1bp T Del, 1bp C Ins, 1bp T Ins,
#              >1bp Del (repeat), >1bp Ins (repeat), Microhomology Del, Complex
ID_CAT_COLORS = {
    ("Del", "C"):  "#fcbb8a",   # 1bp C deletion
    ("Del", "T"):  "#f98465",   # 1bp T deletion
    ("Del", "R"):  "#e85d04",   # repeat-mediated deletion (>1bp)
    ("Del", "M"):  "#9b2226",   # microhomology deletion
    ("Ins", "C"):  "#9ecae1",   # 1bp C insertion
    ("Ins", "T"):  "#6baed6",   # 1bp T insertion
    ("Ins", "R"):  "#2171b5",   # repeat-mediated insertion
    ("Ins", "M"):  "#084594",   # complex insertion
}
ID_CAT_DEFAULT = "#aaaaaa"

# ── Helpers ─────────────────────────────────────────────────────────────────

def _read_tsv(path):
    """Return (header_list, list_of_row_lists) for a tab-separated file."""
    rows = []
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    for line in lines[1:]:
        rows.append(line.split("\t"))
    return header, rows


def _find_sample_row(rows, sample_name):
    for row in rows:
        if row and (row[0] == sample_name or sample_name in row[0] or row[0] in sample_name):
            return row
    return None


def _reconstruct_spectrum(act_file, sig_file, sample_name):
    """
    Compute the reconstructed mutation spectrum (mut_type → count) by
    multiplying COSMIC reference signature profiles by the sample's activities.

    Returns: (mut_types: list[str], counts: list[float])
             or ([], []) on error.
    """
    if not os.path.isfile(act_file) or not os.path.isfile(sig_file):
        return [], []

    act_hdr, act_rows = _read_tsv(act_file)
    sig_hdr, sig_rows = _read_tsv(sig_file)

    sample_row = _find_sample_row(act_rows, sample_name)
    if sample_row is None:
        print(f"  [WARN] sample '{sample_name}' not found in {act_file}", file=sys.stderr)
        return [], []

    # activities dict: sig_name → count
    activities = {}
    for i, col in enumerate(act_hdr[1:], 1):
        if i < len(sample_row):
            try:
                v = float(sample_row[i])
                if v > 0:
                    activities[col] = v
            except ValueError:
                pass

    if not activities:
        return [], []

    # sig_rows: each row is [mut_type, val_sig1, val_sig2, ...]
    # sig_hdr: [MutationType, sig1, sig2, ...]
    mut_types = []
    counts = []
    for row in sig_rows:
        if not row:
            continue
        mut_type = row[0]
        total = 0.0
        for i, sig in enumerate(sig_hdr[1:], 1):
            if sig in activities and i < len(row):
                try:
                    total += activities[sig] * float(row[i])
                except ValueError:
                    pass
        mut_types.append(mut_type)
        counts.append(total)

    return mut_types, counts


# ── Plot: SBS96 matrix profile ───────────────────────────────────────────────

def _sbs96_color(mut_type):
    """Return COSMIC colour for a trinucleotide mutation type like 'A[C>A]A'."""
    m = re.search(r"\[([ACGT]>[ACGT])\]", mut_type)
    if m:
        return SBS_COLORS.get(m.group(1), "#888888")
    return "#888888"


def _read_matrix_file(matrix_path, sample_name):
    """Read observed mutation counts from a SigProfilerMatrixGenerator matrix file.
    Returns (mut_types, counts) or ([], []) on failure.
    """
    if not matrix_path or not os.path.isfile(matrix_path):
        return [], []
    hdr, rows = _read_tsv(matrix_path)
    if not hdr or len(hdr) < 2:
        return [], []
    # Columns: MutationType | sample1 | sample2 | ...
    # Find the column matching sample_name
    sample_col = None
    for i, col in enumerate(hdr[1:], 1):
        if col == sample_name or sample_name in col or col in sample_name:
            sample_col = i
            break
    if sample_col is None:
        return [], []
    mut_types = []
    counts = []
    for row in rows:
        if not row:
            continue
        mut_types.append(row[0])
        try:
            counts.append(float(row[sample_col]) if sample_col < len(row) else 0.0)
        except ValueError:
            counts.append(0.0)
    return mut_types, counts


def plot_sbs96_profile(sbs96_dir, sample_name, matrix_dir=""):
    sol_dir  = os.path.join(sbs96_dir, "Assignment_Solution")
    act_file = os.path.join(sol_dir, "Activities", "Assignment_Solution_Activities.txt")
    sig_file = os.path.join(sol_dir, "Signatures",  "Assignment_Solution_Signatures.txt")
    out_png  = os.path.join(sol_dir, "Activities",  "SBS96_matrix_profile.png")

    # Prefer the real observed matrix from MatrixGenerator output
    mut_types, counts, is_observed = [], [], False
    if matrix_dir:
        sbs_dir = os.path.join(matrix_dir, "SBS")
        if os.path.isdir(sbs_dir):
            # Sort to prefer project-named files over Input_vcffiles.SBS96.exome
            fnames = sorted(os.listdir(sbs_dir),
                            key=lambda f: (f.startswith("Input_vcffiles"), f))
            for fname in fnames:
                if "SBS96" in fname and (fname.endswith(".exome") or fname.endswith(".all")):
                    mt, ct = _read_matrix_file(os.path.join(sbs_dir, fname), sample_name)
                    if mt:
                        mut_types, counts, is_observed = mt, ct, True
                        print(f"  [INFO] Using observed SBS96 matrix: {fname}")
                        break

    # Fallback: reconstruct from Activity × Signatures
    if not mut_types:
        mut_types, counts = _reconstruct_spectrum(act_file, sig_file, sample_name)
    if not mut_types:
        print("  [SKIP] SBS96 matrix profile — no data", file=sys.stderr)
        return

    colors = [_sbs96_color(mt) for mt in mut_types]
    n = len(mut_types)

    fig, ax = plt.subplots(figsize=(22, 4.5), facecolor="white")
    x = np.arange(n)
    ax.bar(x, counts, color=colors, width=0.9, linewidth=0)

    # Substitution type header boxes + dividers
    group_size = 16  # 16 trinucleotide contexts per substitution type
    max_y = max(counts) if max(counts) > 0 else 1
    for gi, (sub, sub_col) in enumerate(SBS_COLORS.items()):
        start = gi * group_size
        end   = start + group_size
        mid   = (start + end) / 2 - 0.5
        ax.axvspan(start - 0.45, end - 0.55, ymin=0, ymax=1,
                   alpha=0.10, color=sub_col, zorder=0)
        ax.text(mid, max_y * 1.07, sub,
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=sub_col)
        if gi > 0:
            ax.axvline(start - 0.5, color="#cccccc", linewidth=0.6, zorder=1)

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, max_y * 1.20)
    ax.set_xticks([])
    ax.set_ylabel("Mutations (observed)" if is_observed else "Mutations (reconstructed)", fontsize=10)
    ax.set_title(f"SBS96 Mutational Spectrum — {sample_name}", fontsize=13, fontweight="bold", pad=18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_tick_params(labelsize=9)

    legend_handles = [mpatches.Patch(facecolor=c, label=s) for s, c in SBS_COLORS.items()]
    ax.legend(handles=legend_handles, ncol=6, loc="upper right",
              fontsize=9, framealpha=0.7, edgecolor="#cccccc")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  [OK] SBS96 matrix profile → {out_png}")


# ── Plot: ID83 matrix profile ────────────────────────────────────────────────

def _id83_color(mut_type):
    """Return colour for ID83 mutation type like '1:Del:C:0' or '5:Ins:R:2'."""
    parts = mut_type.split(":")
    if len(parts) < 3:
        return ID_CAT_DEFAULT
    indel = parts[1]   # Del | Ins
    base  = parts[2]   # C | T | R | M
    return ID_CAT_COLORS.get((indel, base), ID_CAT_DEFAULT)


def _id83_category_label(mut_type):
    """Return a human-readable category label for legend grouping."""
    parts = mut_type.split(":")
    if len(parts) < 3:
        return "Other"
    length = parts[0]
    indel  = parts[1]
    base   = parts[2]
    if base == "M":
        return f"Microhomology {indel}"
    if base == "R":
        return f"Repeat-mediated {indel} (>{1 if int(length) <= 1 else length}bp)" if length.isdigit() else f"Repeat {indel}"
    return f"1bp {base} {indel}"


def plot_id83_profile(id_dir, sample_name, matrix_dir=""):
    sol_dir  = os.path.join(id_dir, "Assignment_Solution")
    act_file = os.path.join(sol_dir, "Activities", "Assignment_Solution_Activities.txt")
    sig_file = os.path.join(sol_dir, "Signatures",  "Assignment_Solution_Signatures.txt")
    out_png  = os.path.join(sol_dir, "Activities",  "ID83_matrix_profile.png")

    # Prefer the real observed matrix from MatrixGenerator output
    mut_types, counts, is_observed = [], [], False
    if matrix_dir:
        id_dir_mg = os.path.join(matrix_dir, "ID")
        if os.path.isdir(id_dir_mg):
            # Sort to prefer project-named files (e.g. 20260225_871_MG.ID83.exome)
            # over generic Input_vcffiles.ID83.exome
            fnames = sorted(os.listdir(id_dir_mg),
                            key=lambda f: (f.startswith("Input_vcffiles"), f))
            for fname in fnames:
                if "ID83" in fname and (fname.endswith(".exome") or fname.endswith(".all")):
                    mt, ct = _read_matrix_file(os.path.join(id_dir_mg, fname), sample_name)
                    if mt:
                        mut_types, counts, is_observed = mt, ct, True
                        print(f"  [INFO] Using observed ID83 matrix: {fname}")
                        break

    # Fallback: reconstruct from Activity × Signatures
    if not mut_types:
        mut_types, counts = _reconstruct_spectrum(act_file, sig_file, sample_name)
    if not mut_types:
        print("  [SKIP] ID83 matrix profile — no data", file=sys.stderr)
        return

    colors = [_id83_color(mt) for mt in mut_types]
    n = len(mut_types)

    fig, ax = plt.subplots(figsize=(22, 4.5), facecolor="white")
    x = np.arange(n)
    ax.bar(x, counts, color=colors, width=0.9, linewidth=0)

    # Build legend from unique category → colour pairs
    seen = {}
    for mt, col in zip(mut_types, colors):
        lbl = _id83_category_label(mt)
        seen[lbl] = col
    legend_handles = [mpatches.Patch(facecolor=c, label=l) for l, c in seen.items()]
    ax.legend(handles=legend_handles, ncol=4, loc="upper right",
              fontsize=8, framealpha=0.7, edgecolor="#cccccc")

    max_y = max(counts) if max(counts) > 0 else 1
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, max_y * 1.25)
    ax.set_xticks([])
    ax.set_ylabel("Mutations (observed)" if is_observed else "Mutations (reconstructed)", fontsize=10)
    ax.set_title(f"ID83 Indel Spectrum — {sample_name}", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_tick_params(labelsize=9)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  [OK] ID83 matrix profile → {out_png}")


# ── Plot: Etiology-annotated activity bar chart ─────────────────────────────

def plot_etiology_activity(sig_dir, sample_name, sig_type):
    """
    Horizontal bar chart showing mutation count per active signature,
    colour-coded and labelled by COSMIC aetiology.

    sig_type: 'SBS96' or 'ID'
    """
    etiology_map = SBS_ETIOLOGY if sig_type == "SBS96" else ID_ETIOLOGY
    sol_dir  = os.path.join(sig_dir, "Assignment_Solution")
    act_file = os.path.join(sol_dir, "Activities", "Assignment_Solution_Activities.txt")
    out_png  = os.path.join(sol_dir, "Activities", "etiology_activity_plot.png")

    if not os.path.isfile(act_file):
        return

    act_hdr, act_rows = _read_tsv(act_file)
    sample_row = _find_sample_row(act_rows, sample_name)
    if sample_row is None:
        return

    # Collect non-zero signatures
    sigs = []
    for i, col in enumerate(act_hdr[1:], 1):
        if i < len(sample_row):
            try:
                v = float(sample_row[i])
                if v > 0:
                    sigs.append((col, v))
            except ValueError:
                pass

    if not sigs:
        print(f"  [SKIP] {sig_type} etiology plot — no active signatures", file=sys.stderr)
        return

    # Sort descending by count
    sigs.sort(key=lambda x: x[1], reverse=True)
    total = sum(v for _, v in sigs)

    sig_names = [s for s, _ in sigs]
    sig_counts = [c for _, c in sigs]
    sig_pcts   = [100.0 * c / total for c in sig_counts]

    # Colours and aetiology labels
    bar_colors  = []
    etiol_labels = []
    for sig, _ in sigs:
        info = etiology_map.get(sig, ("Unknown", "#c7c7c7"))
        etiol_labels.append(info[0])
        bar_colors.append(info[1])

    n = len(sigs)
    fig_h = max(3.5, 0.55 * n + 1.5)
    fig, ax = plt.subplots(figsize=(11, fig_h), facecolor="white")

    y = np.arange(n)
    bars = ax.barh(y, sig_counts, color=bar_colors, height=0.65,
                   edgecolor="white", linewidth=0.5)

    # Annotation: count + % + aetiology on the right of each bar
    max_c = max(sig_counts) if sig_counts else 1
    for i, (bar, lbl, pct) in enumerate(zip(bars, etiol_labels, sig_pcts)):
        w = bar.get_width()
        # Text inside/outside bar depending on length
        x_text = w + max_c * 0.02
        ax.text(x_text, bar.get_y() + bar.get_height() / 2,
                f"{int(round(w))} mut ({pct:.1f}%)  •  {lbl}",
                va="center", ha="left", fontsize=9, color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels(sig_names, fontsize=11, fontweight="bold")
    ax.set_xlabel("Mutations assigned", fontsize=10)
    ax.set_title(f"COSMIC {sig_type} Signature Assignment — {sample_name}\n"
                 f"({int(round(total))} total mutations)",
                 fontsize=12, fontweight="bold")
    ax.set_xlim(0, max_c * 1.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    ax.xaxis.set_tick_params(labelsize=9)

    # Unique-aetiology legend
    seen = {}
    for lbl, col in zip(etiol_labels, bar_colors):
        seen[lbl] = col
    legend_handles = [mpatches.Patch(facecolor=c, label=l) for l, c in seen.items()]
    ax.legend(handles=legend_handles, loc="lower right",
              fontsize=8, framealpha=0.8, edgecolor="#cccccc", title="Aetiology")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  [OK] {sig_type} etiology activity plot → {out_png}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sbs96_dir", required=False, default="",
                        help="Path to the SBS96 assignment directory (containing Assignment_Solution/)")
    parser.add_argument("--id_dir", required=False, default="",
                        help="Path to the ID assignment directory (containing Assignment_Solution/)")
    parser.add_argument("--sample", required=True,
                        help="Sample name as it appears in the Activities.txt Samples column")
    parser.add_argument("--matrix_dir", required=False, default="",
                        help="Path to the matrix_profiles/output/ directory (SigProfilerMatrixGenerator output). "
                             "If provided, real observed spectra are used instead of Activity×Signatures reconstruction.")
    args = parser.parse_args()

    print(f"[generate_signature_plots] sample={args.sample}")
    mdir = args.matrix_dir if args.matrix_dir else ""

    if args.sbs96_dir and os.path.isdir(args.sbs96_dir):
        print(f"  Processing SBS96: {args.sbs96_dir}")
        try:
            plot_sbs96_profile(args.sbs96_dir, args.sample, matrix_dir=mdir)
        except Exception as e:
            print(f"  [ERROR] SBS96 profile: {e}", file=sys.stderr)
        try:
            plot_etiology_activity(args.sbs96_dir, args.sample, "SBS96")
        except Exception as e:
            print(f"  [ERROR] SBS96 etiology plot: {e}", file=sys.stderr)
    else:
        print(f"  [SKIP] SBS96 dir not found: {args.sbs96_dir}", file=sys.stderr)

    if args.id_dir and os.path.isdir(args.id_dir):
        print(f"  Processing ID: {args.id_dir}")
        try:
            plot_id83_profile(args.id_dir, args.sample, matrix_dir=mdir)
        except Exception as e:
            print(f"  [ERROR] ID83 profile: {e}", file=sys.stderr)
        try:
            plot_etiology_activity(args.id_dir, args.sample, "ID")
        except Exception as e:
            print(f"  [ERROR] ID etiology plot: {e}", file=sys.stderr)
    else:
        print(f"  [SKIP] ID dir not found: {args.id_dir}", file=sys.stderr)

    print("[generate_signature_plots] done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
harmonize_cnvs.py
=================
Merges CNVkit+ASCAT (parse_somatic_cnv.py output) with GATK somatic CNV
(MANE_ANNOTATE_SOMATIC_CNV output) and produces a unified, clinically
actionable TSV.

Reportable events
-----------------
  HOMDEL  CN = 0               somatic 38-gene panel
  AMP     CN >= 5              somatic 38-gene panel
  LOH     cnLOH / HETDEL+LOH / LOH_AMP    germline panel

Confidence
----------
  HIGH      event confirmed by both CNVkit/ASCAT AND GATK
  MODERATE  event from a single source

Output columns
--------------
  gene | event_class | cn_primary | cn_secondary | loh_type
  tools | confidence | chrom | log2cr
"""

import argparse
import csv
import os
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_panel(path):
    genes = set()
    if not path or not os.path.isfile(path):
        return genes
    with open(path) as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith('#'):
                genes.add(g.upper())
    return genes


def load_ascat_cnvkit(path, somatic_panel, germline_panel, amp_threshold):
    homdel = {}   # gene -> {cn, chrom, log2}
    amp    = {}   # gene -> {cn, chrom, log2}
    loh    = {}   # gene -> {cn, loh_type, chrom}

    if not path or not os.path.isfile(path):
        return homdel, amp, loh

    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gene  = row.get('gene', '').strip().upper()
            etype = row.get('event_type', '').strip()
            try:
                cn = int(float(row.get('cn', -1)))
            except (ValueError, TypeError):
                cn = -1
            chrom    = row.get('chrom', '')
            log2     = row.get('log2', '')
            loh_type = row.get('loh_type', '').strip()

            if etype == 'HOMDEL' and gene in somatic_panel:
                homdel[gene] = {'cn': 0, 'chrom': chrom, 'log2': log2}

            elif etype == 'AMP' and gene in somatic_panel and cn >= amp_threshold:
                if gene not in amp or cn > amp[gene]['cn']:
                    amp[gene] = {'cn': cn, 'chrom': chrom, 'log2': log2}

            elif etype == '2nd_hit_LOH' and gene in germline_panel:
                loh[gene] = {'cn': cn, 'loh_type': loh_type, 'chrom': chrom}

    return homdel, amp, loh


def load_gatk_somatic(path, somatic_panel, germline_panel, amp_threshold):
    homdel = {}   # gene -> {cn, chrom}
    amp    = {}   # gene -> {cn, chrom}
    loh    = {}   # gene -> {cn, event_type, chrom}

    if not path or not os.path.isfile(path):
        return homdel, amp, loh

    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gene  = row.get('GENE', '').strip().upper()
            etype = row.get('EVENT_TYPE', '').strip()
            try:
                cn = int(float(row.get('CN_ABSOLUTE', -1)))
            except (ValueError, TypeError):
                cn = -1
            chrom = row.get('CONTIG', '')

            if etype == 'DELETION' and cn == 0 and gene in somatic_panel:
                homdel[gene] = {'cn': 0, 'chrom': chrom}

            elif etype == 'AMPLIFICATION' and cn >= amp_threshold and gene in somatic_panel:
                if gene not in amp or cn > amp[gene]['cn']:
                    amp[gene] = {'cn': cn, 'chrom': chrom}

            elif etype in ('CN-LOH', 'DEL+LOH') and gene in germline_panel:
                # Map GATK event type to ASCAT-style loh_type
                loh_type = 'cnLOH' if etype == 'CN-LOH' else 'HETDEL+LOH'
                loh[gene] = {'cn': cn, 'loh_type': loh_type, 'chrom': chrom}

    return homdel, amp, loh


def merge_events(a_homdel, a_amp, a_loh, g_homdel, g_amp, g_loh):
    rows = []

    # HOMDEL
    for gene in sorted(set(a_homdel) | set(g_homdel)):
        a = a_homdel.get(gene)
        g = g_homdel.get(gene)
        tools      = '+'.join(t for t in ['CNVkit+ASCAT' if a else '', 'GATK' if g else ''] if t)
        confidence = 'HIGH' if (a and g) else 'MODERATE'
        rows.append({
            'gene':         gene,
            'event_class':  'HOMDEL',
            'cn_primary':   0,
            'cn_secondary': g['cn'] if g else '',
            'loh_type':     '',
            'tools':        tools,
            'confidence':   confidence,
            'chrom':        (a or g)['chrom'],
            'log2cr':       a['log2'] if a else '',
        })

    # AMP
    for gene in sorted(set(a_amp) | set(g_amp)):
        a = a_amp.get(gene)
        g = g_amp.get(gene)
        tools      = '+'.join(t for t in ['CNVkit+ASCAT' if a else '', 'GATK' if g else ''] if t)
        confidence = 'HIGH' if (a and g) else 'MODERATE'
        rows.append({
            'gene':         gene,
            'event_class':  'AMP',
            'cn_primary':   a['cn'] if a else (g['cn'] if g else ''),
            'cn_secondary': g['cn'] if g else '',
            'loh_type':     '',
            'tools':        tools,
            'confidence':   confidence,
            'chrom':        (a or g)['chrom'],
            'log2cr':       a['log2'] if a else '',
        })

    # LOH
    for gene in sorted(set(a_loh) | set(g_loh)):
        a = a_loh.get(gene)
        g = g_loh.get(gene)
        tools      = '+'.join(t for t in ['ASCAT' if a else '', 'GATK' if g else ''] if t)
        confidence = 'HIGH' if (a and g) else 'MODERATE'
        loh_type   = a['loh_type'] if a else (g['loh_type'] if g else '')
        rows.append({
            'gene':         gene,
            'event_class':  loh_type or 'LOH',
            'cn_primary':   a['cn'] if a else '',
            'cn_secondary': g['cn'] if g else '',
            'loh_type':     loh_type,
            'tools':        tools,
            'confidence':   confidence,
            'chrom':        (a or g)['chrom'],
            'log2cr':       '',
        })

    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='Harmonise CNVkit/ASCAT and GATK somatic CNV outputs.')
    ap.add_argument('--ascat-cnvkit',    required=True,  help='parse_somatic_cnv output TSV')
    ap.add_argument('--gatk-somatic',    default='',     help='MANE-annotated GATK somatic CNV TSV')
    ap.add_argument('--somatic-panel',   required=True,  help='Somatic 38-gene panel file')
    ap.add_argument('--germline-panel',  required=True,  help='Germline panel file')
    ap.add_argument('--amp-threshold',   type=int, default=5)
    ap.add_argument('--out',             required=True,  help='Output TSV path')
    args = ap.parse_args()

    somatic_panel  = load_panel(args.somatic_panel)
    germline_panel = load_panel(args.germline_panel)

    a_homdel, a_amp, a_loh = load_ascat_cnvkit(
        args.ascat_cnvkit, somatic_panel, germline_panel, args.amp_threshold)
    g_homdel, g_amp, g_loh = load_gatk_somatic(
        args.gatk_somatic, somatic_panel, germline_panel, args.amp_threshold)

    rows = merge_events(a_homdel, a_amp, a_loh, g_homdel, g_amp, g_loh)

    COLS = ['gene', 'event_class', 'cn_primary', 'cn_secondary',
            'loh_type', 'tools', 'confidence', 'chrom', 'log2cr']

    with open(args.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

    n_homdel = sum(1 for r in rows if r['event_class'] == 'HOMDEL')
    n_amp    = sum(1 for r in rows if r['event_class'] == 'AMP')
    n_loh    = sum(1 for r in rows if r['event_class'] not in ('HOMDEL', 'AMP'))
    print(f'[OK] {len(rows)} events — HOMDEL:{n_homdel}  AMP:{n_amp}  LOH:{n_loh}')
    print(f'[OK] Written → {args.out}')


if __name__ == '__main__':
    main()

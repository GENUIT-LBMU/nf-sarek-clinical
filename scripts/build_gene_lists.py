#!/usr/bin/env python3
"""
build_gene_lists.py — derive every gene list of the pipeline from the single
gene table assets/gene_panels/genes.tsv (one row per gene, 0/1 columns).

  genes.tsv column   generated file                                   used by
  germline           assets/gene_panels/germline_cancer_predisposition.txt   germline SNV / CNV, LOH
  second_hit         assets/gene_panels/germline_second_hit.txt              somatic SNV second hits
  somatic_core       assets/gene_panels/somatic_snv_panel.txt                somatic SNV report
  tso500             assets/gene_panels/somatic_tso500_panel.txt             TSO500 definition (523)
  second_hit - core  assets/gene_panels/somatic_track_second_hit.txt         second-hit SNV track
  tso500 - core      assets/gene_panels/somatic_track_tso500.txt             pan-tumor SNV track
  somatic_cnv        assets/gene_panels/somatic_38_oncogenes.txt             somatic CNV report
  mmr, hrr           bin/gene_sets.py (MMR_GENES, HRR_GENES)                 correlation, reports

Edit genes.tsv, then run (from the repository root):
  python3 scripts/build_gene_lists.py           # write the files
  python3 scripts/build_gene_lists.py --check   # exit 1 if any file is out of date

tso500: Illumina TruSight Oncology 500 DNA content (523 genes), data sheet
M-GL-00173 v9.0 Table 3, symbols updated to HGNC (tso500_symbol keeps the
data-sheet spelling); see assets/gene_panels/genes.manifest.txt.
second_hit: germline genes where a somatic loss-of-function event completes
biallelic inactivation (activating / risk-allele / non-tumor germline genes
are 0).
"""
import argparse
import csv
import os
import sys

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE = os.path.join(ROOT, 'assets', 'gene_panels', 'genes.tsv')

LISTS = [
    ('germline', 'germline_cancer_predisposition.txt', [
        '# Germline cancer predisposition panel',
        '# Used for:',
        '#   - Germline SNV filtering (HIGH/MODERATE impact variants)',
        '#   - Somatic LOH second-hit detection (ASCAT)',
        '# References: ACMG SF v3.2, ESMO guidelines, published WES cancer panels']),
    ('second_hit', 'germline_second_hit.txt', [
        '# Germline genes for somatic SNV second hits (biallelic inactivation)',
        '# Used for the somatic SNV second-hit track (germline-somatic correlation only)']),
    ('somatic_core', 'somatic_snv_panel.txt', [
        '# Somatic SNV gene panel',
        '# Used for Mutect2 somatic SNV filtering (HIGH/MODERATE impact variants)',
        '# Covers actionable oncogenes, tumor suppressors, and DNA repair genes']),
    ('tso500', 'somatic_tso500_panel.txt', [
        '# Pan-tumor somatic panel: Illumina TruSight Oncology 500 DNA content (523 genes)',
        '# Data sheet M-GL-00173 v9.0 Table 3, HGNC symbols (see genes.manifest.txt)']),
    ('track_second_hit', 'somatic_track_second_hit.txt', [
        '# Somatic SNV second-hit track: second_hit genes outside the core somatic panel',
        '# Variants feed the germline-somatic correlation only (not the somatic SNV table)']),
    ('track_tso500', 'somatic_track_tso500.txt', [
        '# Somatic SNV pan-tumor track: TSO500 genes outside the core somatic panel',
        '# Reported apart from the core panel (Tier I/II with the actionable findings)']),
    ('somatic_cnv', 'somatic_38_oncogenes.txt', [
        '# 38-gene somatic oncogene/TSG panel',
        '# Used for CNVkit somatic event reporting:',
        '#   - HOMDEL (CN=0): homozygous deletion',
        '#   - HETDEL (CN=1): heterozygous deletion',
        '#   - AMP (CN>=5): amplification',
        '# Includes key oncogenes, tumor suppressor genes, and actionable targets']),
]
FLAGS = ['germline', 'second_hit', 'somatic_core', 'tso500', 'somatic_cnv', 'mmr', 'hrr']


def read_table(path):
    with open(path) as fh:
        rows = list(csv.DictReader((l[1:] if l.startswith('#gene') else l for l in fh), delimiter='\t'))
    seen = set()
    for r in rows:
        g = r['gene']
        if g in seen:
            sys.exit(f'ERROR: {g} listed twice in {path}')
        seen.add(g)
        for k in FLAGS:
            if r[k] not in ('0', '1'):
                sys.exit(f'ERROR: {g}: column {k} must be 0 or 1, got {r[k]!r}')
        if r['second_hit'] == '1' and r['germline'] != '1':
            sys.exit(f'ERROR: {g}: second_hit=1 needs germline=1')
    return rows


def render(rows):
    out = {}
    derived = {
        'track_second_hit': lambda r: r['second_hit'] == '1' and r['somatic_core'] == '0',
        'track_tso500':     lambda r: r['tso500'] == '1' and r['somatic_core'] == '0',
    }
    for col, name, header in LISTS:
        keep = derived.get(col, lambda r, c=col: r[c] == '1')
        genes = sorted(r['gene'] for r in rows if keep(r))
        out[os.path.join('assets', 'gene_panels', name)] = '\n'.join(header + genes) + '\n'
    sets = []
    for col, var in (('mmr', 'MMR_GENES'), ('hrr', 'HRR_GENES')):
        genes = sorted(r['gene'] for r in rows if r[col] == '1')
        body = ''.join(f"    '{g}',\n" for g in genes)
        sets.append(f'{var} = frozenset({{\n{body}}})\n')
    out[os.path.join('bin', 'gene_sets.py')] = (
        '"""Gene sets generated by scripts/build_gene_lists.py from\n'
        'assets/gene_panels/genes.tsv. Do not edit: edit genes.tsv and rebuild."""\n\n'
        + '\n'.join(sets))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true', help='only check the generated files are up to date')
    a = ap.parse_args()
    stale = []
    for rel, text in render(read_table(TABLE)).items():
        path = os.path.join(ROOT, rel)
        cur = open(path).read() if os.path.isfile(path) else None
        if cur == text:
            continue
        stale.append(rel)
        if not a.check:
            with open(path, 'w') as fh:
                fh.write(text)
    if a.check and stale:
        sys.exit('ERROR: out of date (run scripts/build_gene_lists.py): ' + ', '.join(stale))
    print(('stale: ' if a.check else 'written: ') + (', '.join(stale) or 'none'))


if __name__ == '__main__':
    main()

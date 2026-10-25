#!/usr/bin/env python3
"""
build_gene_context.py — add the genes missing from
assets/gene_annotations/somatic_gene_context.tsv (pan-tumor TSO500 track),
from sources without commercial-use restrictions. Existing rows (the curated
core panel) are kept unchanged.

  ROLE          CancerMine collated (Lever et al., Nat Methods 2019; CC0, Zenodo
                10.5281/zenodo.1156241): citations summed over cancer types;
                TSG / Oncogene when one has >= --ratio (5) times the other's
                citations and both together >= --min-citations (20); otherwise
                empty ("Driver" citations give no direction and are not used).
                Text mining is noisy: against the 60 curated core genes, a 2x / 3
                citations rule gave 4 opposite roles (MSH6, ATR, BLM, RECQL4 as
                Oncogene); 5x / 20 gives none (21 equal, 39 empty), so a role is
                only set when the literature is one-sided. Never "Both".
  GENE_SUMMARY  CIViC gene description when present (CC0, from the civicpy cache
                of the run), else the NCBI Gene RefSeq summary (E-utilities).
  CGC_TIER, PATHWAY  left empty (COSMIC CGC needs a licence).
  SOURCE        provenance of ROLE and GENE_SUMMARY ("curated" for existing rows).

  python3 scripts/build_gene_context.py \\
      --genes assets/gene_panels/somatic_tso500_panel.txt \\
      --cancermine cancermine_collated.tsv --hgnc hgnc_complete_set.txt \\
      --civic-genes civic_gene_descriptions.tsv \\
      --out assets/gene_annotations/somatic_gene_context.tsv
"""
import argparse
import csv
import json
import sys
import time
import urllib.request

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'


def read_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith('#')]


def cancermine_roles(path, hgnc_symbol, min_cit, ratio):
    cit = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            sym = hgnc_symbol.get(r['gene_hugo_id'], r['gene_normalized'])
            c = cit.setdefault(sym, {'Oncogene': 0, 'Tumor_Suppressor': 0})
            if r['role'] in c:
                c[r['role']] += int(r['citation_count'] or 0)
    roles = {}
    for g, c in cit.items():
        onc, tsg = c['Oncogene'], c['Tumor_Suppressor']
        if onc + tsg < min_cit:
            continue
        if tsg >= ratio * onc:
            roles[g] = ('TSG', f'CancerMine TSG {tsg} / Oncogene {onc} citations')
        elif onc >= ratio * tsg:
            roles[g] = ('Oncogene', f'CancerMine Oncogene {onc} / TSG {tsg} citations')
    return roles


def ncbi_summaries(entrez_ids):
    out, ids = {}, sorted(set(i for i in entrez_ids if i))
    for k in range(0, len(ids), 200):
        url = f'{EUTILS}?db=gene&retmode=json&tool=nf-sarek-clinical&id=' + ','.join(ids[k:k + 200])
        with urllib.request.urlopen(url, timeout=60) as resp:
            res = json.load(resp)['result']
        for u in res.get('uids', []):
            out[u] = ' '.join((res[u].get('summary') or '').split())
        time.sleep(0.4)          # NCBI: <= 3 requests/s without an API key
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genes', required=True)
    ap.add_argument('--cancermine', required=True)
    ap.add_argument('--hgnc', required=True)
    ap.add_argument('--civic-genes', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-citations', type=int, default=20)
    ap.add_argument('--ratio', type=float, default=5.0)
    a = ap.parse_args()

    hgnc_symbol, entrez = {}, {}
    with open(a.hgnc) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r.get('status') == 'Approved':
                hgnc_symbol[r['hgnc_id']] = r['symbol']
                entrez[r['symbol']] = r['entrez_id']
    roles = cancermine_roles(a.cancermine, hgnc_symbol, a.min_citations, a.ratio)
    with open(a.civic_genes) as fh:
        civic = {r['gene']: r['description'] for r in csv.DictReader(fh, delimiter='\t')}

    with open(a.out) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        fields = list(rd.fieldnames)
        rows = list(rd)
    if 'SOURCE' not in fields:
        fields.append('SOURCE')
        for r in rows:
            r['SOURCE'] = 'curated'
    have = {r['GENE'] for r in rows}
    new = [g for g in read_list(a.genes) if g not in have]
    summ = ncbi_summaries([entrez.get(g, '') for g in new if g not in civic])

    n_role = n_sum = 0
    for g in new:
        role, role_src = roles.get(g, ('', ''))
        if g in civic:
            text, sum_src = civic[g], 'CIViC gene description (CC0)'
        else:
            text = summ.get(entrez.get(g, ''), '')
            sum_src = 'NCBI Gene RefSeq summary' if text else ''
        n_role += bool(role)
        n_sum += bool(text)
        rows.append({'GENE': g, 'ROLE': role, 'CGC_TIER': '', 'PATHWAY': '',
                     'GENE_SUMMARY': text.replace('\t', ' '),
                     'SOURCE': '; '.join(x for x in (role_src, sum_src) if x)})
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(new)} genes added ({n_role} with ROLE, {n_sum} with GENE_SUMMARY); {len(rows)} rows', file=sys.stderr)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
germline_somatic_correlate.py
==============================
Germline-Somatic Correlation Analysis for nf-sarek-clinical.

For each germline P/LP/VUS variant:
  - Identifies the functional pathway (MMR, HRR, Other TSG)
  - Detects somatic second hits (SNV or CNV/LOH)
  - Correlates with pathway biomarkers (TMB/MSI for MMR; HRD for HRR)
  - Correlates with mutational signatures (SBS6/15/20/21/26/44 for MMR; SBS3 for HRR)
  - Produces a clinical interpretation (CONCORDANT, PARTIAL, etc.)

Outputs:
  --html-out   HTML fragment (table) for embedding in the HTML report
  --json-out   JSON with structured data for the PDF report
"""

import argparse
import csv
import json
import os
import sys
import html as html_lib


# ---------------------------------------------------------------------------
# Gene pathway definitions, signatures, and ACMG — single source of truth in
# bin/acmg_shared.py. Imported here for pathway classification + germline P/LP
# filtering of correlation input.
# ---------------------------------------------------------------------------
from acmg_shared import (
    MMR_GENES,
    HRR_GENES,
    MMR_SBS_SIGNATURES as MMR_SIGNATURES,
    HRD_SBS_SIGNATURES as HRD_SIGNATURES,
    MMR_ID_SIGNATURES,
    acmg_class as _acmg_class_shared,
    getter_for_dict,
    signature_calls,
)


def _hgvs(row, key):
    """HGVS field without the '.' placeholder of SnpSift extractFields."""
    v = (row.get(key) or '').strip()
    return '' if v == '.' else v

# Biomarker thresholds (specific to germline-somatic correlation, not ACMG)
_E = os.environ.get
TMB_HIGH       = float(_E('BIO_TMB_HIGH', '10.0'))
TMB_BORDER     = float(_E('BIO_TMB_BORDERLINE', '6.0'))
MSI_HIGH_PCT   = float(_E('BIO_MSI_HIGH', '3.5'))
MSI_BORDER     = float(_E('BIO_MSI_BORDERLINE', '1.5'))
HRD_HIGH       = float(_E('BIO_HRD_HIGH', '0.45'))
HRD_BORDER     = float(_E('BIO_HRD_BORDERLINE', '0.3'))


def acmg_class(row_dict):
    """Returns ACMG class (P/LP/VUS/LB/B) for a row dict."""
    return _acmg_class_shared(getter_for_dict(row_dict))


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------
def load_germline_variants(path):
    variants = []
    if not path or not os.path.isfile(path):
        return variants
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            acmg = acmg_class(row)
            if acmg not in ('P', 'LP'):
                continue
            gene = (row.get('ANN[*].GENE') or '').strip().upper()
            if not gene or gene in ('.', 'NA', ''):
                continue
            variants.append({
                'gene':    gene,
                'hgvs_c':  _hgvs(row, 'ANN[*].HGVS_C'),
                'hgvs_p':  _hgvs(row, 'ANN[*].HGVS_P'),
                'impact':  (row.get('ANN[*].IMPACT') or '').strip(),
                'effect':  (row.get('ANN[*].EFFECT') or '').strip(),
                'clnsig':  (row.get('ClinVar_sig') or row.get('dbNSFP_clinvar_clnsig') or '').strip(),
                'revel':   (row.get('REVEL_score') or row.get('dbNSFP_REVEL_score') or '').strip(),
                'gnomad':  (row.get('gnomAD_AF') or row.get('dbNSFP_gnomAD_exomes_AF') or '').strip(),
                'acmg':    acmg,
            })
    return variants


def _int_or_none(value):
    """'4' and '4.0' -> 4; empty or non-numeric -> None."""
    try:
        return int(float((value or '').strip()))
    except (ValueError, OverflowError):
        return None


def _float_or_none(value):
    try:
        return float((value or '').strip())
    except ValueError:
        return None


def load_germline_cnv_variants(path, panel_genes=None):
    """Load GATK gCNV germline deletions (ALT <DEL>, QA>=4, not low-mappability) in panel genes.

    Returns one record per gene per segment with affected exons listed.
    Duplications (<DUP>) are excluded — they are not classic 1st TSG hits.
    Only genes in the germline cancer-predisposition panel are reported.
    QA >= 4 : ALT allele Phred quality >= 4 (~60% confidence in the DEL/DUP call)
    """
    MIN_QA = int(os.environ.get('FILT_GERMLINE_CNV_MIN_QA', '4'))
    LOWMAP_FRAC = float(os.environ.get('FILT_GERMLINE_CNV_LOWMAP_FRAC', '0.5'))
    variants = []
    if not path or not os.path.isfile(path):
        return variants

    # Group rows by segment (CHROM, POS, END) → collect genes/exons
    segments = {}
    seg_order = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            try:
                cn_int = int(row.get('CN', '2').strip())
            except ValueError:
                cn_int = 2
            # gCNV sets ALT against the contig ploidy (male chrX: CN=1 is
            # reference); old tables without ALT fall back to CN < 2.
            if 'ALT' in row:
                if (row.get('ALT') or '').strip() != '<DEL>':
                    continue
            elif cn_int >= 2:
                continue  # diploid or duplication — skip
            key = (row.get('CHROM', ''), row.get('POS', ''), row.get('END', ''))
            if key not in segments:
                np_val = _int_or_none(row.get('NP'))
                if np_val is None:
                    np_val = 0
                segments[key] = {
                    'chrom': row.get('CHROM', ''),
                    'pos':   row.get('POS', ''),
                    'end':   row.get('END', ''),
                    'cn':    cn_int,
                    'np':    np_val,
                    'qa':    _int_or_none(row.get('QA')),
                    'qs':    row.get('QS', '.'),
                    # fraction of panel exon bp in low-mappability regions
                    # (pseudogenes / paralogs); None for old tables
                    'lowmap': _float_or_none(row.get('LOWMAP_SEG')),
                    'genes': {},   # gene -> set of exon ints
                }
                seg_order.append(key)
            gene = row.get('gen', '').strip().upper()
            exon = row.get('exon', '').strip()
            if gene:
                segments[key]['genes'].setdefault(gene, set())
                if exon:
                    try:
                        segments[key]['genes'][gene].add(int(exon))
                    except ValueError:
                        pass

    for key in seg_order:
        seg = segments[key]
        if seg['qa'] is None:
            print('WARNING: germline CNV ' + seg['chrom'] + ':' + seg['pos'] + '-'
                  + seg['end'] + ' has non-numeric QA, skipped', file=sys.stderr)
            continue
        if seg['qa'] < MIN_QA:
            continue
        if seg['lowmap'] is not None and seg['lowmap'] >= LOWMAP_FRAC:
            continue  # shown apart in the reports, never a first hit
        cn  = seg['cn']
        locus = seg['chrom'] + ':' + seg['pos'] + '-' + seg['end']
        for gene, exon_set in seg['genes'].items():
            if not gene or gene in ('.', 'NA'):
                continue
            # Only report panel genes
            if panel_genes and gene.upper() not in panel_genes:
                continue
            # Format exon list
            exons_sorted = sorted(exon_set)
            if exons_sorted:
                if len(exons_sorted) == 1:
                    exon_str = 'exon ' + str(exons_sorted[0])
                elif exons_sorted == list(range(min(exons_sorted), max(exons_sorted) + 1)):
                    exon_str = 'exons ' + str(exons_sorted[0]) + '-' + str(exons_sorted[-1])
                else:
                    exon_str = 'exons ' + ', '.join(str(e) for e in exons_sorted)
                cnv_desc = exon_str + (' homozygous deletion' if cn == 0 else ' deletion')
            else:
                exon_str = ''
                cnv_desc = 'heterozygous deletion' if cn == 1 else 'homozygous deletion'
            acmg   = 'P' if cn == 0 else 'LP'
            impact = 'HIGH' if cn == 0 else 'MODERATE'
            variants.append({
                'gene':         gene,
                'hgvs_c':       locus + ' (' + cnv_desc + ')',
                'hgvs_p':       '',
                'impact':       impact,
                'effect':       'copy_number_variation',
                'clnsig':       '',
                'revel':        '',
                'gnomad':       '',
                'acmg':         acmg,
                'variant_type': 'CNV',
                'exons':        exon_str,
                'cn':           cn,
                'qs':           seg['qs'],
                'np':           seg['np'],
            })
    return variants


def load_somatic_snv_by_gene(path):
    idx = {}
    if not path or not os.path.isfile(path):
        return idx
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            gene   = (row.get('ANN[*].GENE') or '').strip().upper()
            impact = (row.get('ANN[*].IMPACT') or '').strip().upper()
            if not gene or impact not in ('HIGH', 'MODERATE'):
                continue
            # Common population polymorphism (AMP Tier IV): not a somatic hit
            if (row.get('POP_COMMON') or '').strip() == '1':
                continue
            rec = {
                'hgvs_c': _hgvs(row, 'ANN[*].HGVS_C'),
                'hgvs_p': _hgvs(row, 'ANN[*].HGVS_P'),
                'effect': (row.get('ANN[*].EFFECT') or '').strip(),
                'impact': impact,
                'vaf':    (row.get('GEN[1].AF') or '').strip(),
            }
            idx.setdefault(gene, []).append(rec)
    return idx


def load_somatic_cnv_by_gene(path):
    idx = {}
    if not path or not os.path.isfile(path):
        return idx
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            gene  = (row.get('gene') or '').strip().upper()
            etype = (row.get('event_type') or '').strip()
            ltype = (row.get('loh_type') or '').strip()
            cn    = (row.get('cn') or '').strip()
            if not gene:
                continue
            idx.setdefault(gene, []).append({
                'event_type': etype,
                'loh_type':   ltype,
                'cn':         cn,
            })
    return idx


def find_somatic_biallelic(somatic_snv_by_gene, somatic_cnv_by_gene,
                           already_reported_genes=None):
    """Find genes with TWO independent somatic hits (somatic SNV + somatic CNV/LOH).

    These represent somatic biallelic inactivation without a germline first hit.
    Genes already covered in the germline–somatic correlation table are excluded.
    Returns a list of dicts describing each event pair.
    """
    if already_reported_genes is None:
        already_reported_genes = set()

    biallelic = []
    snv_genes = set(somatic_snv_by_gene.keys())
    cnv_genes = set(somatic_cnv_by_gene.keys())

    for gene in sorted(snv_genes & cnv_genes):
        if gene in already_reported_genes:
            continue  # already covered in germline-somatic correlation

        cnvs = somatic_cnv_by_gene[gene]
        # Only loss-of-function CNV events qualify as a second hit (not AMP)
        lof_cnvs = [c for c in cnvs
                    if c['event_type'] in ('HOMDEL', 'HETDEL', '2nd_hit_LOH')]
        if not lof_cnvs:
            continue

        snvs = somatic_snv_by_gene[gene]

        # Pathway classification
        if gene in MMR_GENES:
            pathway = 'MMR'
        elif gene in HRR_GENES:
            pathway = 'HRR'
        else:
            pathway = 'Other TSG'

        # Prioritise HIGH-impact SNV; fall back to first available
        best_snv = sorted(snvs, key=lambda x: 0 if x['impact'] == 'HIGH' else 1)[0]
        best_cnv = lof_cnvs[0]

        snv_disp = best_snv['hgvs_p'] or best_snv['hgvs_c'] or best_snv['effect'] or '?'
        vaf_s    = (' VAF=' + best_snv['vaf']) if best_snv['vaf'] and best_snv['vaf'] not in ('.', '') else ''

        cnv_etype = best_cnv['event_type']
        cnv_ltype = best_cnv.get('loh_type', '')
        cnv_cn    = best_cnv.get('cn', '')

        if cnv_etype == 'HOMDEL':
            cnv_disp  = 'homozygous deletion (CN=0)'
            mechanism = ('Two somatic hits resulting in complete gene loss: somatic '
                         + best_snv['impact'].lower() + ' variant ' + snv_disp + vaf_s
                         + ' followed by homozygous deletion (CN=0). Both alleles are inactivated.')
        elif cnv_etype == 'HETDEL':
            cnv_disp = 'heterozygous deletion (CN=' + str(cnv_cn) + ')'
            try:
                vaf_f = float(best_snv['vaf'])
                if vaf_f > 0.5:
                    mechanism = (
                        'Biallelic somatic inactivation: somatic '
                        + best_snv['impact'].lower() + ' variant '
                        + snv_disp + ' (VAF=' + best_snv['vaf'] + ') '
                        'combined with heterozygous deletion of the locus (CN=1). '
                        'The elevated VAF (>' + str(round(vaf_f, 2)) + ') is consistent '
                        'with preferential loss of the wild-type allele.'
                    )
                else:
                    mechanism = (
                        'Two somatic hits: '
                        + best_snv['impact'].lower() + ' variant '
                        + snv_disp + vaf_s
                        + ' and heterozygous deletion (CN=1). '
                        'May represent biallelic inactivation depending on clonal architecture.'
                    )
            except (ValueError, TypeError):
                mechanism = (
                    'Two somatic hits: '
                    + best_snv['impact'].lower() + ' variant '
                    + snv_disp + vaf_s
                    + ' and heterozygous deletion (CN=1).'
                )
        elif cnv_etype == '2nd_hit_LOH':
            if cnv_ltype == 'cnLOH':
                cnv_disp  = 'cnLOH (CN=2, minor allele=0)'
                mechanism = (
                    'Biallelic somatic inactivation via copy-neutral LOH: somatic '
                    + best_snv['impact'].lower() + ' variant ' + snv_disp + vaf_s
                    + ' followed by duplication of the mutant allele (cnLOH, CN=2 minor=0), '
                    'eliminating the wild-type copy.'
                )
            else:
                cnv_disp  = 'LOH (' + cnv_ltype + ')'
                mechanism = (
                    'Two somatic hits: '
                    + best_snv['impact'].lower() + ' variant '
                    + snv_disp + vaf_s
                    + ' and LOH (' + cnv_ltype + ').'
                )
        else:
            cnv_disp  = cnv_etype
            mechanism = (
                'Two somatic hits: '
                + best_snv['impact'].lower() + ' variant '
                + snv_disp + vaf_s
                + ' and somatic CNV (' + cnv_etype + ').'
            )

        biallelic.append({
            'gene':      gene,
            'pathway':   pathway,
            'snv_label': 'Somatic ' + best_snv['impact'] + ': ' + snv_disp + vaf_s,
            'cnv_label': 'Somatic ' + cnv_disp,
            'mechanism': mechanism,
            'all_snvs':  snvs,
            'all_cnvs':  lof_cnvs,
        })

    return biallelic


# ---------------------------------------------------------------------------
# Correlation engine
# ---------------------------------------------------------------------------
def _fmt_second_hits(gene, somatic_snv_by_gene, somatic_cnv_by_gene):
    hits = []
    # Somatic SNV hits
    for s in somatic_snv_by_gene.get(gene, []):
        label = s['hgvs_p'] or s['hgvs_c'] or s['effect'] or '?'
        vaf_s = (' VAF=' + s['vaf']) if s['vaf'] and s['vaf'] not in ('.', '') else ''
        hits.append({'type': 'SNV', 'label': 'Somatic SNV: ' + label + vaf_s, 'impact': s['impact']})
    # CNV / LOH hits
    for c in somatic_cnv_by_gene.get(gene, []):
        etype = c['event_type']
        ltype = c['loh_type']
        cn    = c['cn']
        if etype == '2nd_hit_LOH':
            if ltype == 'cnLOH':
                hits.append({'type': 'cnLOH',      'label': 'cnLOH (CN=2, minor allele=0)'})
            elif ltype == 'HETDEL+LOH':
                hits.append({'type': 'LOH',         'label': 'LOH + copy loss (CN=1)'})
            elif 'major_amp' in ltype:
                hits.append({'type': 'LOH_amp',     'label': 'LOH + amplification (' + ltype + ')'})
            else:
                hits.append({'type': 'LOH_other',   'label': 'LOH (' + ltype + ')'})
        elif etype == 'HOMDEL':
            hits.append({'type': 'HOMDEL', 'label': 'Homozygous deletion (CN=0)'})
        elif etype == 'HETDEL':
            hits.append({'type': 'HETDEL', 'label': 'Heterozygous deletion (CN=1)'})
    return hits


def _sig_status(calls, wanted, kind):
    """('detected: SBS6' | 'not detected' | 'not assessable (n=12 indels)', found list, status)."""
    if not calls['assessable']:
        unit = 'indels' if kind == 'ID' else 'SBS'
        return f"not assessable (n={calls['total']:.0f} {unit})", [], 'n/a'
    found = sorted(wanted & calls['detected'])
    return (('detected: ' + ', '.join(found)) if found else 'not detected'), found, ('HIGH' if found else 'absent')


def _pathway_biomarkers(pathway, tmb, msi, hrd_prob, sbs_calls, id_calls):
    """MMR: MSI-high, or TMB-high explained by an MMR signature (SBS or ID).
    HRR: HRD probability, or an HRD SBS signature. Signatures count only when
    detected by acmg_shared.signature_calls (enough mutations)."""
    items = []
    match = False

    if pathway == 'MMR':
        tmb_s  = 'HIGH' if tmb  >= TMB_HIGH    else ('borderline' if tmb  >= TMB_BORDER  else 'low')
        msi_s  = 'HIGH' if msi  >= MSI_HIGH_PCT else ('borderline' if msi  >= MSI_BORDER  else 'low')
        items.append({'key': 'TMB',  'val': str(round(tmb, 1)) + ' mut/Mb', 'status': tmb_s,
                      'threshold': '>= ' + str(TMB_HIGH)})
        items.append({'key': 'MSI',  'val': str(round(msi, 1)) + '%',        'status': msi_s,
                      'threshold': '>= ' + str(MSI_HIGH_PCT) + '%'})
        sbs_v, sbs_found, sbs_s = _sig_status(sbs_calls, MMR_SIGNATURES, 'SBS')
        items.append({'key': 'MMR SBS signatures', 'val': sbs_v, 'status': sbs_s,
                      'threshold': '/'.join(sorted(MMR_SIGNATURES))})
        id_v, id_found, id_s = _sig_status(id_calls, MMR_ID_SIGNATURES, 'ID')
        items.append({'key': 'MMR ID signatures', 'val': id_v, 'status': id_s,
                      'threshold': '/'.join(sorted(MMR_ID_SIGNATURES))})
        match = msi >= MSI_HIGH_PCT or (tmb >= TMB_HIGH and bool(sbs_found or id_found))

    elif pathway == 'HRR':
        hrd_s = 'HIGH' if hrd_prob >= HRD_HIGH else ('borderline' if hrd_prob >= HRD_BORDER else 'low')
        items.append({'key': 'HRD probability', 'val': str(round(hrd_prob, 4)),
                      'status': hrd_s, 'threshold': '>= ' + str(HRD_HIGH)})
        sig_v, found_hrd, sig_s = _sig_status(sbs_calls, HRD_SIGNATURES, 'SBS')
        items.append({'key': 'HRD signatures', 'val': sig_v, 'status': sig_s,
                      'threshold': '/'.join(sorted(HRD_SIGNATURES))})
        match = hrd_prob >= HRD_HIGH or bool(found_hrd)

    else:
        items.append({'key': 'Pathway', 'val': 'No specific biomarker panel for this gene',
                      'status': 'n/a', 'threshold': 'n/a'})

    return items, match


def _interpretation(acmg, has_second_hit, biomarker_match, pathway):
    if acmg == 'VUS':
        if has_second_hit and biomarker_match:
            return 'VUS_UPGRADED', 'badge-LP', 'VUS — Evidence supports pathogenicity'
        return 'VUS_UNCERTAIN', 'badge-VUS', 'VUS — Insufficient somatic evidence'

    if has_second_hit and biomarker_match:
        return 'CONCORDANT', 'badge-P', 'CONCORDANT'
    if has_second_hit and not biomarker_match and pathway in ('MMR', 'HRR'):
        return 'BIALLELIC_NO_BIOMARKER', 'badge-LP', 'Biallelic — pathway biomarkers absent'
    if not has_second_hit and biomarker_match:
        return 'BIOMARKER_NO_2ND_HIT', 'badge-LP', 'Pathway active — 2nd hit undetected'
    if has_second_hit:
        return 'BIALLELIC_OTHER', 'badge-LP', 'Biallelic — Other TSG (no pathway biomarkers)'
    return 'GERMLINE_ONLY', 'badge-III', 'Germline only — no somatic confirmation'


def run_correlation(args):
    # Parse metric float values safely
    def safe_float(s, default=0.0):
        try:
            return float(str(s).strip())
        except (ValueError, TypeError):
            return default

    tmb      = safe_float(args.tmb)
    msi      = safe_float(args.msi)
    hrd_prob = safe_float(args.hrd)

    # Load all data
    germ_snv_variants   = load_germline_variants(args.germline)
    for v in germ_snv_variants:
        v.setdefault('variant_type', 'SNV')   # tag source type
    # Load panel genes for CNV filtering
    _panel_genes = set()
    _panel_path = getattr(args, 'germline_panel', '') or ''
    if _panel_path and os.path.isfile(_panel_path):
        with open(_panel_path) as _pf:
            for _ln in _pf:
                _g = _ln.strip().split('\t')[0].strip().upper()
                if _g and not _g.startswith('#'):
                    _panel_genes.add(_g)
    germ_cnv_variants   = load_germline_cnv_variants(
        getattr(args, 'cnv_germline', None) or '',
        panel_genes=_panel_genes or None)
    germ_variants       = germ_snv_variants + germ_cnv_variants

    somatic_snv_by_gene = load_somatic_snv_by_gene(args.somatic)
    # Second-hit track: germline genes outside the core somatic panel (no gene overlap)
    for _g, _recs in load_somatic_snv_by_gene(args.somatic_second_hit).items():
        somatic_snv_by_gene.setdefault(_g, []).extend(_recs)
    somatic_cnv_by_gene = load_somatic_cnv_by_gene(args.cnv)
    sbs_calls           = signature_calls(args.sig, args.tumor, 'SBS')
    id_calls            = signature_calls(getattr(args, 'sig_id', None) or '', args.tumor, 'ID')

    results = []
    for v in germ_variants:
        gene = v['gene']

        # Pathway
        if gene in MMR_GENES:
            pathway = 'MMR'
        elif gene in HRR_GENES:
            pathway = 'HRR'
        else:
            pathway = 'Other TSG'

        # Second hits
        second_hits = _fmt_second_hits(gene, somatic_snv_by_gene, somatic_cnv_by_gene)
        has_second_hit = bool(second_hits)

        # Biomarkers
        biomarker_items, biomarker_match = _pathway_biomarkers(
            pathway, tmb, msi, hrd_prob, sbs_calls, id_calls)

        # Interpretation
        interp_key, interp_class, interp_label = _interpretation(
            v['acmg'], has_second_hit, biomarker_match, pathway)

        # Skip entries with no somatic evidence whatsoever — they have nothing
        # actionable to show in a germline–somatic correlation table.
        # (A P/LP-only finding without any somatic second hit or pathway
        #  biomarker concordance belongs in the germline section, not here.)
        if not has_second_hit and not biomarker_match:
            continue

        results.append({
            'gene':            gene,
            'hgvs_c':          v['hgvs_c'],
            'hgvs_p':          v['hgvs_p'],
            'impact':          v['impact'],
            'acmg':            v['acmg'],
            'clnsig':          v['clnsig'],
            'pathway':         pathway,
            'second_hits':     second_hits,
            'biomarkers':      biomarker_items,
            'biomarker_match': biomarker_match,
            'interp_key':      interp_key,
            'interp_class':    interp_class,
            'interp_label':    interp_label,
            'variant_type':    v.get('variant_type', 'SNV'),
            'exons':           v.get('exons', ''),
        })

    # Track genes already in germline correlation so we don't double-report them
    reported_genes = {r['gene'] for r in results}

    # Find purely somatic biallelic events (SNV + CNV/LOH, no germline anchor)
    biallelic_somatic = find_somatic_biallelic(
        somatic_snv_by_gene, somatic_cnv_by_gene,
        already_reported_genes=reported_genes,
    )

    return results, biallelic_somatic, sbs_calls, id_calls


# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------
BADGE_CSS = {
    'P':   'background:#ef4444;color:#fff',
    'LP':  'background:#f97316;color:#fff',
    'VUS': 'background:#eab308;color:#111',
    'LB':  'background:#3b82f6;color:#fff',
    'B':   'background:#22c55e;color:#fff',
    'badge-P':   'background:#ef4444;color:#fff',
    'badge-LP':  'background:#f97316;color:#fff',
    'badge-VUS': 'background:#eab308;color:#111',
    'badge-III': 'background:#6b7280;color:#fff',
    'badge-IV':  'background:#374151;color:#888',
}

STATUS_ICON = {
    'HIGH':      '<span style="color:#22c55e">&#10003;</span>',
    'borderline':'<span style="color:#eab308">&#9679;</span>',
    'low':       '<span style="color:#ef4444">&#10007;</span>',
    'absent':    '<span style="color:#ef4444">&#10007;</span>',
    'detected':  '<span style="color:#22c55e">&#10003;</span>',
    'n/a':       '<span style="color:#94a3b8">&#8212;</span>',
}

INTERP_MEANING = {
    'CONCORDANT':
        'Germline variant confirmed biallelic with a somatic second hit. '
        'Pathway functional biomarkers are consistent. High confidence for this pathway.',
    'VUS_UPGRADED':
        'Variant of Uncertain Significance with supporting somatic evidence '
        '(second hit + pathway biomarkers). Consider reclassification.',
    'BIOMARKER_NO_2ND_HIT':
        'Pathway biomarkers are active (high TMB/MSI or HRD), but no second '
        'somatic hit was detected in the target gene. Possible explanations: '
        'somatic hit outside coding regions, epigenetic silencing, or phenocopy.',
    'BIALLELIC_NO_BIOMARKER':
        'Biallelic inactivation detected (germline + somatic second hit), but '
        'pathway biomarkers are not elevated. Consider tissue specificity or '
        'incomplete penetrance.',
    'BIALLELIC_OTHER':
        'Biallelic inactivation detected in a gene without a specific biomarker '
        'panel. Loss of function confirmed.',
    'VUS_UNCERTAIN':
        'Variant of Uncertain Significance. No somatic second hit detected and '
        'pathway biomarkers are not conclusive.',
    'GERMLINE_ONLY':
        'Germline variant without detected somatic second hit and without '
        'pathway biomarker concordance. Cannot confirm functional impact at '
        'the tumor level from available data.',
}

PATHWAY_GENES_DISPLAYED = {
    'MMR': 'MLH1, MSH2, MSH6, PMS2, EPCAM, MSH3',
    'HRR': 'BRCA1, BRCA2, PALB2, RAD51C, RAD51D, ATM, CHEK2, BARD1, BRIP1 …',
    'Other TSG': 'Tumor suppressor gene — no specific biomarker panel',
}

SECOND_HIT_ICON = {
    'SNV':       '&#9632;',   # filled square
    'cnLOH':     '&#9670;',   # diamond
    'LOH':       '&#9670;',
    'LOH_amp':   '&#9670;',
    'LOH_other': '&#9670;',
    'HOMDEL':    '&#9650;',   # triangle
    'HETDEL':    '&#9651;',
}


def _badge(text, style_key):
    style = BADGE_CSS.get(style_key, 'background:#334155;color:#fff')
    return ('<span style="display:inline-block;padding:2px 9px;border-radius:12px;'
            'font-size:.72rem;font-weight:700;' + style + '">'
            + html_lib.escape(text) + '</span>')


def build_html(results, sbs_calls, tmb, msi, hrd_prob, biallelic_results=None):
    detected_sigs = sbs_calls['detected']
    if biallelic_results is None:
        biallelic_results = []
    if not results and not biallelic_results:
        return ('<p style="color:#fed7aa;background:rgba(124,45,18,.3);'
                'border-radius:8px;padding:8px 12px;font-size:.85rem;">'
                'No germline P/LP/VUS variants or somatic biallelic events available '
                'for correlation analysis.'
                '</p>')

    # --- biomarker summary bar ---
    mmr_sigs_found = sorted(MMR_SIGNATURES & detected_sigs)
    hrd_sigs_found = sorted(HRD_SIGNATURES & detected_sigs)

    tmb_ok  = tmb  >= TMB_HIGH
    msi_ok  = msi  >= MSI_HIGH_PCT
    hrd_ok  = hrd_prob >= HRD_HIGH
    mmr_ok  = bool(mmr_sigs_found)
    hrr_ok  = bool(hrd_sigs_found)

    def pill(label, ok, val=''):
        col = '#22c55e' if ok else '#64748b'
        icon = '&#10003;' if ok else '&#10007;'
        v = (' ' + html_lib.escape(str(val))) if val else ''
        return ('<span style="display:inline-flex;align-items:center;gap:4px;'
                'padding:3px 10px;border-radius:16px;font-size:.75rem;font-weight:600;'
                'background:rgba(15,23,42,.7);border:1px solid ' + col + ';color:' + col + ';">'
                + icon + ' ' + label + v + '</span>')

    bio_bar = ('<div style="display:flex;flex-wrap:wrap;gap:8px;margin:12px 0 18px;">'
               + pill('TMB', tmb_ok, str(round(tmb, 1)) + ' mut/Mb')
               + pill('MSI', msi_ok, str(round(msi, 1)) + '%')
               + pill('HRD prob', hrd_ok, str(round(hrd_prob, 4)))
               + pill('MMR sigs', mmr_ok, ', '.join(mmr_sigs_found) if mmr_ok else
                      ('none' if sbs_calls['assessable'] else 'n/a (' + str(int(sbs_calls['total'])) + ' SBS)'))
               + pill('HRD sigs', hrr_ok, ', '.join(hrd_sigs_found) if hrr_ok else
                      ('none' if sbs_calls['assessable'] else 'n/a (' + str(int(sbs_calls['total'])) + ' SBS)'))
               + '</div>')

    # --- legend ---
    legend = ('<div style="font-size:.75rem;color:rgba(226,232,240,.7);margin-bottom:10px;">'
              '<strong>Interpretation codes:</strong>&nbsp;&nbsp;'
              + _badge('CONCORDANT', 'badge-P') + ' Both biallelic inactivation and pathway biomarkers confirmed&ensp;'
              + _badge('Biallelic / Pathway active', 'badge-LP') + ' Partial evidence&ensp;'
              + _badge('VUS', 'badge-VUS') + ' Uncertain significance&ensp;'
              + _badge('Germline only', 'badge-III') + ' No somatic confirmation'
              + '</div>')

    # --- main correlation table ---
    th_style = 'padding:.4rem .55rem;text-align:left;color:#fff;font-size:.75rem;'
    td_style = 'padding:.35rem .55rem;border-bottom:1px solid rgba(148,163,184,.2);font-size:.75rem;color:#e5e7eb;vertical-align:top;'

    thead = ('<thead style="background:rgba(249,115,22,.55);">'
             '<tr>'
             '<th style="' + th_style + '">Gene</th>'
             '<th style="' + th_style + '">Germline variant</th>'
             '<th style="' + th_style + '">ACMG</th>'
             '<th style="' + th_style + '">Pathway</th>'
             '<th style="' + th_style + '">Somatic 2nd hit</th>'
             '<th style="' + th_style + '">Pathway biomarkers</th>'
             '<th style="' + th_style + '">Interpretation</th>'
             '</tr></thead>')

    tbody = '<tbody>'
    for r in results:
        # Gene cell
        gene_cell = ('<td style="' + td_style + 'font-weight:700;">'
                     + html_lib.escape(r['gene']) + '</td>')

        # Variant cell — display differs for SNV vs CNV
        vtype = r.get('variant_type', 'SNV')
        if vtype == 'CNV':
            cnv_badge = ('<span style="font-size:.63rem;background:#7c3aed;color:#fff;'
                         'padding:1px 6px;border-radius:8px;margin-left:5px;">CNV</span>')
            var_cell = ('<td style="' + td_style + '">'
                        + '<div style="font-weight:600;">' + html_lib.escape(r['hgvs_c']) + cnv_badge + '</div>'
                        + '</td>')
        else:
            hgvs_p = r['hgvs_p'] or r['hgvs_c'] or '?'
            hgvs_c = r['hgvs_c'] if r['hgvs_p'] else ''
            var_cell = ('<td style="' + td_style + '">'
                        '<div style="font-weight:600;">' + html_lib.escape(hgvs_p) + '</div>'
                        + ('<div style="font-size:.68rem;color:#94a3b8;">' + html_lib.escape(hgvs_c) + '</div>' if hgvs_c else '')
                        + '</td>')

        # ACMG badge
        acmg = r['acmg']
        acmg_cell = '<td style="' + td_style + '">' + _badge(acmg, acmg) + '</td>'

        # Pathway cell
        pathway = r['pathway']
        pw_colors = {'MMR': '#6366f1', 'HRR': '#0ea5e9', 'Other TSG': '#8b5cf6'}
        pw_col = pw_colors.get(pathway, '#64748b')
        pw_cell = ('<td style="' + td_style + '">'
                   '<span style="color:' + pw_col + ';font-weight:600;">'
                   + html_lib.escape(pathway) + '</span>'
                   + '<div style="font-size:.65rem;color:#64748b;">'
                   + html_lib.escape(PATHWAY_GENES_DISPLAYED.get(pathway, '')) + '</div>'
                   + '</td>')

        # Second hit cell
        hits = r['second_hits']
        if hits:
            hit_lines = ''.join(
                '<div style="margin-bottom:3px;">'
                + '<span style="color:#38bdf8;">' + SECOND_HIT_ICON.get(h['type'], '&#9679;') + '</span> '
                + html_lib.escape(h['label'])
                + '</div>'
                for h in hits
            )
            sh_cell = '<td style="' + td_style + '">' + hit_lines + '</td>'
        else:
            sh_cell = ('<td style="' + td_style + 'color:#64748b;">'
                       '<em>No second hit detected</em></td>')

        # Biomarker cell
        bm_lines = ''
        for b in r['biomarkers']:
            icon = STATUS_ICON.get(b['status'], '')
            bm_lines += ('<div style="margin-bottom:2px;">'
                         + icon + ' <strong>' + html_lib.escape(b['key']) + ':</strong> '
                         + html_lib.escape(b['val'])
                         + ' <span style="font-size:.65rem;color:#64748b;">(thr: '
                         + html_lib.escape(b['threshold']) + ')</span>'
                         + '</div>')
        bm_cell = '<td style="' + td_style + '">' + bm_lines + '</td>'

        # Interpretation cell
        meaning = INTERP_MEANING.get(r['interp_key'], '')
        interp_cell = ('<td style="' + td_style + '">'
                       + _badge(r['interp_label'], r['interp_class'])
                       + ('<div style="font-size:.67rem;color:#94a3b8;margin-top:4px;">'
                          + html_lib.escape(meaning) + '</div>' if meaning else '')
                       + '</td>')

        tbody += ('<tr>' + gene_cell + var_cell + acmg_cell + pw_cell
                  + sh_cell + bm_cell + interp_cell + '</tr>\n')

    tbody += '</tbody>'

    table = ('<div style="overflow-x:auto;border-radius:10px;border:1px solid rgba(148,163,184,.2);margin-top:10px;">'
             '<table style="border-collapse:collapse;width:100%;font-size:.78rem;">'
             + thead + tbody + '</table></div>')

    # Methodology note
    method_note = ('<div style="font-size:.72rem;color:rgba(226,232,240,.55);'
                   'border-top:1px solid rgba(148,163,184,.2);margin-top:14px;padding-top:8px;">'
                   '<strong>Methodology:</strong> ACMG/AMP classification uses ClinVar annotation (primary), '
                   'population frequency (gnomAD), functional impact (SnpEff, MANE Select), and in silico '
                   'predictors (REVEL, SpliceAI). '
                   'Second hit detection uses somatic SNV (Mutect2, HIGH/MODERATE impact) and CNV/LOH '
                   '(ASCAT segments: cnLOH, HETDEL+LOH, HOMDEL). '
                   'Pathway biomarkers: MMR — MSI-high, or TMB &ge; 10 mut/Mb with an MMR signature '
                   '(SBS6/15/20/21/26/44 or ID1/ID2); HRR — HRProfiler probability &ge; 0.45 or SBS3/SBS8. '
                   'Signatures count only with &ge;10% of the assigned mutations, &ge;10 mutations and '
                   '&ge;100 SBS / 50 indels in the sample. '
                   'This analysis is intended as a decision-support tool and must be reviewed by a '
                   'certified clinical geneticist.</div>')

    # --- somatic biallelic events section ---
    somatic_bial_html = ''
    if biallelic_results:
        sb_th = 'padding:.4rem .55rem;text-align:left;color:#fff;font-size:.75rem;'
        sb_td = ('padding:.35rem .55rem;border-bottom:1px solid rgba(148,163,184,.2);'
                 'font-size:.75rem;color:#e5e7eb;vertical-align:top;')

        sb_thead = ('<thead style="background:rgba(139,92,246,.55);">'
                    '<tr>'
                    '<th style="' + sb_th + '">Gene</th>'
                    '<th style="' + sb_th + '">Pathway</th>'
                    '<th style="' + sb_th + '">Somatic SNV hit</th>'
                    '<th style="' + sb_th + '">Somatic CNV/LOH hit</th>'
                    '<th style="' + sb_th + '">Mechanism / Interpretation</th>'
                    '</tr></thead>')

        sb_tbody = '<tbody>'
        pw_colors = {'MMR': '#6366f1', 'HRR': '#0ea5e9', 'Other TSG': '#8b5cf6'}
        for b in biallelic_results:
            pw_col = pw_colors.get(b['pathway'], '#64748b')
            sb_tbody += (
                '<tr>'
                '<td style="' + sb_td + 'font-weight:700;">' + html_lib.escape(b['gene']) + '</td>'
                '<td style="' + sb_td + '">'
                + '<span style="color:' + pw_col + ';font-weight:600;">' + html_lib.escape(b['pathway']) + '</span>'
                + '</td>'
                '<td style="' + sb_td + '">'
                + '<span style="color:#38bdf8;">&#9632;</span> ' + html_lib.escape(b['snv_label'])
                + '</td>'
                '<td style="' + sb_td + '">'
                + '<span style="color:#f59e0b;">&#9651;</span> ' + html_lib.escape(b['cnv_label'])
                + '</td>'
                '<td style="' + sb_td + 'color:#cbd5e1;">' + html_lib.escape(b['mechanism']) + '</td>'
                '</tr>\n'
            )
        sb_tbody += '</tbody>'

        somatic_bial_html = (
            '<h4 style="color:#c4b5fd;font-size:.85rem;margin:22px 0 6px;font-weight:600;">'
            '&#9670; Somatic Biallelic Events (two somatic hits, no germline variant)</h4>'
            '<p style="font-size:.72rem;color:rgba(226,232,240,.6);margin:0 0 8px;">'
            'Genes with a somatic coding variant (HIGH/MODERATE) <em>and</em> a somatic '
            'copy-number loss (HETDEL, HOMDEL, or LOH) — independent of any germline variant. '
            'These represent candidate somatic biallelic inactivation events.</p>'
            '<div style="overflow-x:auto;border-radius:10px;border:1px solid rgba(148,163,184,.2);">'
            '<table style="border-collapse:collapse;width:100%;font-size:.78rem;">'
            + sb_thead + sb_tbody
            + '</table></div>'
        )

    return bio_bar + legend + table + method_note + somatic_bial_html


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(description='Germline-Somatic Correlation Analysis')
    ap.add_argument('--germline',  default='',  help='Germline SNV TSV (filtered, classified)')
    ap.add_argument('--somatic',   default='',  help='Somatic SNV TSV')
    ap.add_argument('--somatic-second-hit', default='', help='Somatic SNV TSV of the second-hit track (germline genes outside the core panel)')
    ap.add_argument('--cnv',           default='', help='Somatic CNV report TSV (parse_somatic_cnv.py output)')
    ap.add_argument('--cnv-germline',  default='', help='GATK gCNV germline annotated TSV (segments_annotated.tsv)')
    ap.add_argument('--sig',       default='',  help='SBS96 assignment activities TXT')
    ap.add_argument('--sig-id',    default='',  help='ID (indel) assignment activities TXT')
    ap.add_argument('--tumor',     default='',  help='Tumor sample name (for signature lookup)')
    ap.add_argument('--tmb',       default='0', help='TMB value (mut/Mb)')
    ap.add_argument('--msi',       default='0', help='MSI value (%)')
    ap.add_argument('--hrd',       default='0', help='HRD probability (0-1)')
    ap.add_argument('--html-out',  default='',  help='Output HTML fragment file')
    ap.add_argument('--germline-panel', default='', help='Germline gene panel file (one gene per line)')
    ap.add_argument('--json-out',  default='',  help='Output JSON file for PDF')
    return ap.parse_args()


def main():
    args = parse_args()

    def safe_float(s, default=0.0):
        try:
            return float(str(s).strip())
        except (ValueError, TypeError):
            return default

    tmb      = safe_float(args.tmb)
    msi      = safe_float(args.msi)
    hrd_prob = safe_float(args.hrd)

    results, biallelic_somatic, sbs_calls, id_calls = run_correlation(args)
    detected_sigs, detected_id_sigs = sbs_calls['detected'], id_calls['detected']

    if args.html_out:
        html_fragment = build_html(results, sbs_calls, tmb, msi, hrd_prob,
                                   biallelic_results=biallelic_somatic)
        with open(args.html_out, 'w') as fh:
            fh.write(html_fragment)
        print('[INFO] HTML correlation written: ' + args.html_out)

    if args.json_out:
        # Serialize (sets not JSON-serializable, convert to lists)
        out = {
            'results':           results,
            'biallelic_somatic': biallelic_somatic,
            'detected_sigs':     sorted(detected_sigs),
            'detected_id_sigs':  sorted(detected_id_sigs),
            # signature calls are only made with enough mutations (acmg_shared.signature_calls)
            'sig_sbs_assessable': sbs_calls['assessable'], 'sig_sbs_total': sbs_calls['total'],
            'sig_id_assessable':  id_calls['assessable'],  'sig_id_total':  id_calls['total'],
            'metrics':           {'tmb': tmb, 'msi': msi, 'hrd_prob': hrd_prob},
        }
        with open(args.json_out, 'w') as fh:
            json.dump(out, fh, indent=2)
        print('[INFO] JSON correlation written: ' + args.json_out)

    if not args.html_out and not args.json_out:
        # Print summary to stdout
        print('Germline-Somatic Correlation Summary')
        print('=' * 60)
        for r in results:
            print(r['gene'] + ' [' + r['acmg'] + '] ' + r['pathway']
                  + ' -> ' + r['interp_label'])
            for h in r['second_hits']:
                print('  2nd hit: ' + h['label'])
            for b in r['biomarkers']:
                print('  ' + b['key'] + ': ' + b['val'] + ' [' + b['status'] + ']')
        if biallelic_somatic:
            print('\nSomatic Biallelic Events')
            print('-' * 40)
            for b in biallelic_somatic:
                print(b['gene'] + ' [' + b['pathway'] + ']')
                print('  ' + b['snv_label'])
                print('  ' + b['cnv_label'])
                print('  => ' + b['mechanism'])
        print('Detected SBS signatures: ' + (', '.join(sorted(detected_sigs)) if detected_sigs else 'none'))
        print('Detected ID  signatures: ' + (', '.join(sorted(detected_id_sigs)) if detected_id_sigs else 'none'))


if __name__ == '__main__':
    main()

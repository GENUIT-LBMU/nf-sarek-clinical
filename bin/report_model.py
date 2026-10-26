#!/usr/bin/env python3
"""
report_model.py — the report model of one tumor/normal pair.

Reads every report input once and makes every interpretation the reports show:
ACMG classes, AMP tiers and CIViC evidence, biomarker statuses (MSI, HRD, low
coverage), ASCAT alerts, somatic CNV / LOH (with the events the clinical PDF
leaves out), germline CNV segments, signature activities, calls and aetiologies,
CPIC pharmacogenomics, per-gene coverage, clinical trials, analyses not run,
tool and resource versions. Writes {pair}.report_model.json.

The two reports only render this model (bin/render_html.py: laboratory report,
everything; bin/clinical_pdf_report.py: clinical subset), so they cannot
diverge. The JSON is published next to the reports as the record of what they
said.
"""

import argparse
import csv
import json
import os
import sys
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import amp_tier as _amp                                       # noqa: E402
from acmg_shared import (acmg_criteria_list, getter_for_dict,  # noqa: E402
                         msi_status, signature_calls)
from cnv_shared import (ascat_alerts, harmonize_cnv,           # noqa: E402
                        loh_relevant_genes, ploidy_note)
import report_shared as rs                                     # noqa: E402

_E = os.environ.get
SCHEMA_VERSION = 1

# Thresholds in force (nextflow.config env <- conf/filters.config), recorded in the model
THRESHOLDS = {
    'msi_high_pct':          float(_E('BIO_MSI_HIGH', '3.5')),
    'msi_borderline_pct':    float(_E('BIO_MSI_BORDERLINE', '1.5')),
    'tmb_high':              float(_E('BIO_TMB_HIGH', '10')),
    'tmb_borderline':        float(_E('BIO_TMB_BORDERLINE', '6.0')),
    'hrd_high':              float(_E('BIO_HRD_HIGH', '0.45')),
    'hrd_borderline':        float(_E('BIO_HRD_BORDERLINE', '0.30')),
    'low_tumor_coverage':    float(_E('QC_LOW_COV_THRESHOLD', '50')),
    'ascat_min_purity':      float(_E('QC_ASCAT_MIN_PURITY', '0.20')),
    'ascat_min_gof':         float(_E('QC_ASCAT_MIN_GOF', '90')),
    'ascat_ploidy_note':     float(_E('QC_ASCAT_PLOIDY_NOTE', '2.7')),
    'somatic_max_pop_af':    float(_E('FILT_SOMATIC_MAX_POP_AF', '0.01')),
    'germline_cnv_min_qa':   int(_E('FILT_GERMLINE_CNV_MIN_QA', '4')),
    'germline_cnv_lowmap_frac': float(_E('FILT_GERMLINE_CNV_LOWMAP_FRAC', '0.5')),
    'somatic_cnv_amp_cn':    5,
}

# ── COSMIC v3.4 aetiologies: (text, category). The category gives the colour. ─
SIG_AETIOLOGY = {
    'SBS1':  ('Deamination of 5-methylcytosine (ageing / clock-like)', 'clock'),
    'SBS2':  ('APOBEC cytidine deaminase activity', 'apobec'),
    'SBS3':  ('Defective homologous recombination DNA repair (HRD)', 'hrd'),
    'SBS4':  ('Tobacco smoking mutagen exposure', 'exposure'),
    'SBS5':  ('Clock-like (unknown mechanism, correlates with age)', 'clock'),
    'SBS6':  ('Mismatch repair deficiency (MMR-d)', 'mmr'),
    'SBS7a': ('UV light exposure (C>T at dipyrimidines)', 'exposure'),
    'SBS7b': ('UV light exposure (T>A at CC dinucleotides)', 'exposure'),
    'SBS7c': ('UV light exposure', 'exposure'),
    'SBS7d': ('UV light exposure', 'exposure'),
    'SBS8':  ('Defective homologous recombination / unknown', 'hrd'),
    'SBS9':  ('AID / somatic hypermutation (lymphoid)', 'other'),
    'SBS10a': ('POLE exonuclease domain mutation', 'pol'),
    'SBS10b': ('POLE exonuclease domain mutation', 'pol'),
    'SBS10c': ('POLE exonuclease domain mutation', 'pol'),
    'SBS10d': ('POLE exonuclease domain mutation', 'pol'),
    'SBS11': ('Temozolomide alkylating chemotherapy', 'treatment'),
    'SBS12': ('Unknown aetiology (found in liver cancer)', 'unknown'),
    'SBS13': ('APOBEC cytidine deaminase activity', 'apobec'),
    'SBS14': ('MMR deficiency + POLE exonuclease domain mutation', 'mmr'),
    'SBS15': ('MMR deficiency (associated with MLH1 silencing)', 'mmr'),
    'SBS16': ('Unknown aetiology (found in liver cancer)', 'unknown'),
    'SBS17a': ('Unknown (found in oesophagus/stomach cancer)', 'unknown'),
    'SBS17b': ('Unknown (found in oesophagus/stomach cancer)', 'unknown'),
    'SBS18': ('Reactive oxygen species (ROS) DNA damage', 'other'),
    'SBS19': ('Unknown aetiology', 'unknown'),
    'SBS20': ('MMR deficiency + POLD1 exonuclease domain mutation', 'mmr'),
    'SBS21': ('MMR deficiency (MSI-H tumours)', 'mmr'),
    'SBS22': ('Aristolochic acid exposure', 'exposure'),
    'SBS22a': ('Aristolochic acid exposure', 'exposure'),
    'SBS22b': ('Aristolochic acid exposure', 'exposure'),
    'SBS23': ('Unknown aetiology', 'unknown'),
    'SBS24': ('Aflatoxin exposure', 'exposure'),
    'SBS25': ('Chemotherapy (unknown agent)', 'treatment'),
    'SBS26': ('Mismatch repair deficiency', 'mmr'),
    'SBS28': ('Unknown aetiology', 'unknown'),
    'SBS29': ('Tobacco chewing (gutkha / betel quid)', 'exposure'),
    'SBS30': ('Defective base excision repair (NTHL1 biallelic)', 'other'),
    'SBS31': ('Platinum chemotherapy (cisplatin)', 'treatment'),
    'SBS32': ('Azathioprine immunosuppression', 'treatment'),
    'SBS33': ('Unknown aetiology', 'unknown'),
    'SBS34': ('Unknown aetiology', 'unknown'),
    'SBS35': ('Platinum chemotherapy (cisplatin)', 'treatment'),
    'SBS36': ('Defective base excision repair (MUTYH biallelic)', 'other'),
    'SBS37': ('Unknown aetiology', 'unknown'),
    'SBS38': ('UV light exposure (indirect)', 'exposure'),
    'SBS39': ('Unknown aetiology', 'unknown'),
    'SBS40': ('Unknown (clock-like, correlates with age)', 'clock'),
    'SBS41': ('Unknown aetiology', 'unknown'),
    'SBS42': ('Haloalkane exposure (occupational)', 'exposure'),
    'SBS44': ('Mismatch repair deficiency', 'mmr'),
    'SBS84': ('AID activity (lymphoid)', 'other'),
    'SBS85': ('AID activity (lymphoid)', 'other'),
    'SBS86': ('Unknown chemotherapy', 'treatment'),
    'SBS87': ('Thiopurine chemotherapy', 'treatment'),
    'SBS88': ('Colibactin (pks+ E. coli)', 'exposure'),
    'SBS89': ('Unknown aetiology', 'unknown'),
    'SBS90': ('Duocarmycin exposure', 'treatment'),
    'SBS92': ('Tobacco smoking', 'exposure'),
    'SBS93': ('Unknown aetiology', 'unknown'),
    'SBS94': ('Unknown aetiology', 'unknown'),
    'ID1':   ('Replication slippage at T homopolymers (clock-like)', 'mmr'),
    'ID2':   ('Replication slippage at T homopolymers — MMR-related', 'mmr'),
    'ID3':   ('Tobacco smoking', 'exposure'),
    'ID4':   ('Unknown aetiology', 'unknown'),
    'ID5':   ('Unknown aetiology', 'unknown'),
    'ID6':   ('Homologous recombination deficiency (HRD)', 'hrd'),
    'ID7':   ('Replication slippage at microsatellites — MMR-related', 'mmr'),
    'ID8':   ('Non-homologous end joining (NHEJ) repair', 'other'),
    'ID9':   ('Unknown aetiology', 'unknown'),
    'ID10':  ('Unknown aetiology', 'unknown'),
    'ID11':  ('Unknown aetiology', 'unknown'),
    'ID12':  ('Unknown aetiology', 'unknown'),
    'ID13':  ('UV light exposure', 'exposure'),
    'ID14':  ('Unknown aetiology', 'unknown'),
    'ID17':  ('Unknown aetiology', 'unknown'),
    'ID18':  ('Colibactin (pks+ E. coli)', 'exposure'),
}
_SIG_UNKNOWN = ('Unknown aetiology (see COSMIC v3.4)', 'unknown')

# ── CPIC: the rsIDs of assets/pharmacogenomics (DPYD, UGT1A1) ─────────────────
CPIC_TABLE = {
    'rs3918290':  {'gene': 'DPYD', 'allele': '*2A (c.1905+1G>A)', 'activity': '0',
                   'phenotype': 'No Function — Poor Metabolizer', 'drugs': '5-FU, Capecitabine, Tegafur',
                   'recommendation': 'Avoid fluoropyrimidines. If unavoidable, reduce starting dose by ≥50% and perform therapeutic drug monitoring (TDM). Risk of severe/fatal toxicity.',
                   'alert': 'HIGH'},
    'rs55886062': {'gene': 'DPYD', 'allele': '*13 (c.1679T>G, p.Ile560Ser)', 'activity': '0',
                   'phenotype': 'No Function — Poor Metabolizer', 'drugs': '5-FU, Capecitabine, Tegafur',
                   'recommendation': 'Avoid fluoropyrimidines. If unavoidable, reduce starting dose by ≥50% and perform TDM.',
                   'alert': 'HIGH'},
    'rs67376798': {'gene': 'DPYD', 'allele': 'c.2846A>T (p.Asp949Val)', 'activity': '0.5',
                   'phenotype': 'Decreased Function — Intermediate Met.', 'drugs': '5-FU, Capecitabine, Tegafur',
                   'recommendation': 'Reduce starting dose by 25–50%; titrate based on toxicity and pharmacokinetics. Monitor for hematologic and gastrointestinal toxicity.',
                   'alert': 'MODERATE'},
    'rs56038477': {'gene': 'DPYD', 'allele': 'HapB3 (c.1236G>A, E412E)', 'activity': '0.5',
                   'phenotype': 'Decreased Function — Intermediate Met.', 'drugs': '5-FU, Capecitabine, Tegafur',
                   'recommendation': 'Reduce starting dose by 25–50%; titrate based on toxicity and pharmacokinetics.',
                   'alert': 'MODERATE'},
    'rs4148323':  {'gene': 'UGT1A1', 'allele': '*6 (c.211G>A, p.Gly71Arg)', 'activity': 'Reduced',
                   'phenotype': 'Reduced UGT1A1 Activity', 'drugs': 'Irinotecan, Belinostat, Nilotinib',
                   'recommendation': 'Decreased-function allele, like *28 (more frequent in East Asian populations). Two decreased-function alleles (*6/*6, *6/*28): reduce irinotecan starting dose by ≥30%. One allele: standard dose acceptable; monitor for diarrhea and neutropenia.',
                   'alert': 'MODERATE'},
    'rs3064744':  {'gene': 'UGT1A1', 'allele': '*28 (A(TA)7TAA promoter repeat)', 'activity': 'Reduced',
                   'phenotype': 'Reduced UGT1A1 Activity', 'drugs': 'Irinotecan',
                   'recommendation': 'Consider reduced irinotecan starting dose; monitor for SN-38-related toxicity.',
                   'alert': 'MODERATE'},
}


# ── File helpers ──────────────────────────────────────────────────────────────

def present(path):
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def read_tsv(path):
    """(headers, rows as dicts); ([], []) when missing."""
    if not present(path):
        return [], []
    with open(path, encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def tsv_val(path, col):
    """Column `col` (case-insensitive) of the first data row, 'NA' when missing."""
    headers, rows = read_tsv(path)
    if not rows:
        return 'NA'
    for h in headers:
        if h.lower() == col.lower():
            v = (rows[0].get(h) or '').strip()
            return v if v else 'NA'
    return 'NA'


def msi_pct(path):
    """MSIsensor-pro: 3rd whitespace field of line 2, without '%'."""
    if not present(path):
        return 'NA'
    with open(path) as fh:
        lines = fh.read().splitlines()
    if len(lines) < 2 or len(lines[1].split()) < 3:
        return 'NA'
    return lines[1].split()[2].replace('%', '')


def cov_mean(path):
    """mosdepth summary: mean of the total_region row."""
    if not present(path):
        return 'NA'
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    if not rows:
        return 'NA'
    hdr = [h.lower() for h in rows[0]]
    if 'mean' not in hdr:
        return 'NA'
    for r in rows[1:]:
        if r and r[0] == 'total_region':
            return r[hdr.index('mean')]
    return 'NA'


def load_panel(path, upper=False):
    """Gene list (first column, '#' comments skipped), in file order."""
    genes = []
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                g = line.strip().split('\t')[0].split()[0] if line.strip() else ''
                if g and not g.startswith('#'):
                    genes.append(g.upper() if upper else g)
    return genes


def _num(v):
    try:
        f = float(str(v).strip())
        return None if f != f else f
    except (ValueError, TypeError):
        return None


# ── Metrics and statuses ──────────────────────────────────────────────────────

def build_metrics(a):
    purity_raw = tsv_val(a.ascat, 'purity')
    p = _num(purity_raw)
    purity_pct = 'NA' if p is None else f'{(p if p > 1.5 else p * 100):.1f}'
    m = {
        'tmb':    tsv_val(a.tmb, 'TMB_mut_per_mb'),
        'msi':    msi_pct(a.msi),
        'loh':    tsv_val(a.ascat, 'LOH'),
        'gi':     tsv_val(a.ascat, 'GI'),
        'wgd':    tsv_val(a.ascat, 'WGD'),
        'ploidy': tsv_val(a.ascat, 'ploidy'),
        'gof':    tsv_val(a.ascat, 'goodness_of_fit'),
        'purity': purity_pct,
        'cov_n':  cov_mean(a.cov_n),
        'cov_t':  cov_mean(a.cov_t),
    }
    m['msi_status'] = msi_status(m['msi'])
    # HRD: probability >= filters.biomarkers.hrd_high (the threshold also given to
    # HRProfiler); without a probability, the reason recorded by HRD_PROFILER
    prob = tsv_val(a.hrd, 'hrd.prob')
    m['hrd_prob'] = prob
    pf = _num(prob) if prob not in ('', 'NA', '.') else None
    if pf is None:
        why = tsv_val(a.hrd, 'status')
        m['hrd_status'] = 'NA' + (f' ({why})' if why not in ('', 'NA', 'ok') else '')
    else:
        m['hrd_status'] = 'HRD+' if pf >= THRESHOLDS['hrd_high'] else 'HRD-'
    m['hrd_prob_display'] = rs.fmt_hrd_prob(prob)
    tmb = _num(m['tmb'])
    m['tmb_status'] = 'NA' if tmb is None else ('TMB-H' if tmb >= THRESHOLDS['tmb_high'] else 'TMB-L')
    cov_t = _num(m['cov_t'])
    # Tumor coverage below QC_LOW_COV_THRESHOLD: TMB, MSI, HRD, signatures, CNV indeterminate
    m['low_tumor_coverage'] = cov_t is not None and cov_t < THRESHOLDS['low_tumor_coverage']
    return m


# ── Germline ──────────────────────────────────────────────────────────────────

def load_germline(path):
    """Every candidate row with its ACMG class and criteria (the PDF shows P/LP/VUS)."""
    _, rows = read_tsv(path)
    for row in rows:
        cls, crit = acmg_criteria_list(getter_for_dict(row))
        row['_ACMG'] = cls
        row['_ACMG_CRITERIA'] = crit
    return rows


def germline_cnv(path, panel_genes):
    """GATK gCNV events (ALT <DEL>/<DUP>, else CN != 2) grouped per segment, panel
    genes only, QA >= min; segments with >= lowmap_frac of their exon bp in GIAB
    low-mappability regions are listed apart."""
    _, rows = read_tsv(path)
    min_qa = THRESHOLDS['germline_cnv_min_qa']

    def is_event(row):
        if 'ALT' in row:
            return (row.get('ALT') or '').strip() in ('<DEL>', '<DUP>')
        try:
            return int((row.get('CN') or '2').strip()) != 2
        except ValueError:
            return False

    segments, order = {}, []
    for row in rows:
        if not is_event(row):
            continue
        key = (row.get('CHROM', ''), row.get('POS', ''), row.get('END', ''))
        if key not in segments:
            segments[key] = {'chrom': row.get('CHROM', ''), 'pos': row.get('POS', ''),
                             'end': row.get('END', ''), 'alt': row.get('ALT', '.'),
                             'cn': row.get('CN', '.'), 'np': row.get('NP', '.'),
                             'qa': row.get('QA', '0'), 'qs': row.get('QS', '.'),
                             'lowmap': _num(row.get('LOWMAP_SEG')), 'genes': []}
            order.append(key)
        gene = (row.get('gen') or '').strip()
        if gene and gene not in segments[key]['genes']:
            segments[key]['genes'].append(gene)

    kept = []
    for key in order:
        seg = segments[key]
        try:
            qa = int(float(seg['qa'] or ''))
        except (ValueError, TypeError, OverflowError):
            continue                    # non-numeric QA: below threshold
        hits = [g for g in seg['genes'] if g in panel_genes]
        if qa >= min_qa and hits:
            seg['genes'] = hits
            kept.append(seg)
    frac = THRESHOLDS['germline_cnv_lowmap_frac']
    lowmap = [s for s in kept if s['lowmap'] is not None and s['lowmap'] >= frac]
    return {'segments': [s for s in kept if s not in lowmap], 'lowmap': lowmap,
            'n_rows': len(rows)}


# ── Somatic ───────────────────────────────────────────────────────────────────

def load_gene_ctx(path):
    ctx = {}
    for row in read_tsv(path)[1]:
        gene = (row.get('GENE') or '').strip()
        if gene:
            ctx[gene] = {'role': (row.get('ROLE') or '').strip(),
                         'tier': (row.get('CGC_TIER') or '').strip(),
                         'pathway': (row.get('PATHWAY') or '').strip(),
                         'summary': (row.get('GENE_SUMMARY') or '').strip()}
    return ctx


def gene_consequence(role, impact, effect):
    """Clinical consequence sentence from gene role x impact x effect."""
    role, impact, effect = (role or '').upper(), (impact or '').upper(), (effect or '').lower()
    is_both = role == 'BOTH'
    is_tsg = is_both or 'TSG' in role or role == 'TUMOR SUPPRESSOR'
    is_onco = is_both or 'ONCOGENE' in role
    if impact == 'HIGH':
        if is_tsg or is_both:
            return ('Loss-of-function variant (frameshift, nonsense or splice). '
                    'Truncating mutations in a tumor suppressor eliminate functional protein; '
                    'biallelic inactivation through LOH or a second somatic hit drives tumorigenesis.')
        if is_onco:
            return ('High-impact variant in an oncogene — uncommon. '
                    'May represent loss of a regulatory domain or a rare neomorphic event. '
                    'Functional significance requires validation.')
        return 'High-impact variant (frameshift/nonsense/splice) of uncertain functional consequence.'
    if impact == 'MODERATE':
        if 'missense' in effect or 'inframe' in effect:
            if is_onco and not is_both:
                return ('Missense variant in an oncogene. '
                        'Many oncogenic missense mutations are activating, locking the protein '
                        'in a constitutively active state and driving downstream pathway signaling.')
            if is_tsg and not is_both:
                return ('Missense variant in a tumor suppressor. '
                        'May act as a dominant-negative allele (interfering with the wild-type copy) '
                        'or a hypomorphic allele with partial loss of function.')
            if is_both:
                return ('Missense variant in a gene with dual TSG/oncogenic roles. '
                        'Functional impact depends on variant position; activating and inactivating '
                        'effects are both described.')
            return 'Missense variant of uncertain functional consequence.'
        return 'Moderate-impact variant; functional significance depends on the specific effect and gene context.'
    return 'Variant of uncertain functional consequence in this gene.'


def load_somatic(path, civic_path, terms, track, gene_ctx):
    """Rows with AMP tier, tier-defining CIViC record, all CIViC records, the
    evidence summary, the track and the gene context."""
    _, rows = read_tsv(path)
    civic = _amp.load_civic(civic_path)
    for row in rows:
        records = civic.get((row.get('CHROM', ''), row.get('POS', ''),
                             row.get('REF', ''), row.get('ALT', '')), [])
        tier, rec = _amp.variant_tier(records, terms, row.get)
        row['_AMP_TIER'] = tier
        row['_CIVIC'] = rec or {}
        row['_CIVIC_ALL'] = records
        row['_EVIDENCE'] = _amp.evidence_summary(records, terms)
        row['_TRACK'] = track
        ctx = gene_ctx.get((row.get('ANN[*].GENE') or '').strip(), {})
        row['_CTX'] = ctx
        effect = row.get('ANN[*].EFFECT') or row.get('ANN[*].ANNOTATION') or ''
        row['_CONSEQUENCE'] = (gene_consequence(ctx.get('role'), row.get('ANN[*].IMPACT'), effect)
                               if ctx.get('role') else '')
    return rows


# ── Signatures ────────────────────────────────────────────────────────────────

def signatures(path, pair, kind):
    """This pair's assigned signatures (activity > 0, by count), with aetiology,
    and the MMR / HRD calls (acmg_shared.signature_calls)."""
    if not present(path):
        return None
    _, rows = read_tsv(path)
    row = next((r for r in rows if (r.get('Samples') or '') == pair), None)
    if row is None:
        print(f'[WARN] {pair} not found in {path}', file=sys.stderr)
        return {'path': path, 'assigned': [], 'calls': None}
    assigned = []
    for col, val in row.items():
        if not col or col.lower() in ('samples', 'sample'):
            continue
        v = _num(val)
        if v and v > 0:
            text, cat = SIG_AETIOLOGY.get(col, _SIG_UNKNOWN)
            assigned.append({'signature': col, 'count': int(round(v)), 'aetiology': text, 'category': cat})
    assigned.sort(key=lambda s: -s['count'])
    total = sum(s['count'] for s in assigned)
    for s in assigned:
        s['percent'] = round(100.0 * s['count'] / total, 1) if total else 0.0
    calls = None
    if kind in ('SBS', 'ID'):
        c = signature_calls(path, pair, kind)
        calls = {'detected': sorted(c['detected']), 'assessable': c['assessable'], 'total': c['total']}
    return {'path': path, 'assigned': assigned, 'total': total, 'calls': calls}


# ── Pharmacogenomics, coverage ────────────────────────────────────────────────

def pharmaco(path):
    out = []
    for row in read_tsv(path)[1]:
        for rsid in (row.get('ID') or '').split(';'):
            rsid = rsid.strip()
            if rsid in CPIC_TABLE:
                out.append(dict(CPIC_TABLE[rsid], rsid=rsid,
                                gt=(row.get('GEN[0].GT') or '').strip(),
                                dp=(row.get('GEN[0].DP') or row.get('DP') or '').strip(),
                                hgvs_c=(row.get('ANN[*].HGVS_C') or '').strip()))
    return out


def coverage(path):
    """{gene: {mean, min, max, n_regions}} from GENE_COVERAGE."""
    return {(r.get('gene') or '').strip(): {'mean': r.get('mean_coverage', 'NA'),
                                             'min': r.get('min_coverage', 'NA'),
                                             'max': r.get('max_coverage', 'NA'),
                                             'n_regions': r.get('n_regions', '')}
            for r in read_tsv(path)[1] if (r.get('gene') or '').strip()}


# ── Identification (extra samplesheet columns, all optional) ──────────────────
# Patient columns repeat on every row of the patient; sample columns per sample.
# sarek ignores them ("unidentified headers"); the pipeline reads the same file
# with --samplesheet.
PATIENT_FIELDS = ['patient_name', 'birth_date', 'sex', 'medical_record',
                  'ordering_physician', 'institution', 'diagnosis']
SAMPLE_FIELDS = ['accession_id', 'specimen_type', 'specimen_site', 'collection_date',
                 'received_date', 'tumor_content']


def identification(path, tumor, normal):
    """{'patient': {...}, 'tumor': {...}, 'normal': {...}, 'warnings': [...]} from the
    samplesheet rows of the two samples; a field missing everywhere is ''."""
    out = {'patient': {f: '' for f in PATIENT_FIELDS + ['patient']},
           'tumor': {f: '' for f in SAMPLE_FIELDS}, 'normal': {f: '' for f in SAMPLE_FIELDS},
           'warnings': [], 'source': os.path.basename(path) if present(path) else ''}
    if not present(path):
        return out
    with open(path, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        rows = [{k.strip(): (v or '').strip() for k, v in r.items() if k} for r in reader]
        header = {h.strip() for h in reader.fieldnames or []}
    # a plain sarek samplesheet (no identification column): reports stay compact
    if not header & (set(PATIENT_FIELDS + SAMPLE_FIELDS) - {'sex'}):
        out['source'] = ''
        return out
    rows = [r for r in rows if r.get('sample') in (tumor, normal)]
    if not rows:
        out['warnings'].append(f'{tumor} / {normal} not found in {out["source"]}')
        return out
    for f in PATIENT_FIELDS + ['patient']:
        vals = list(dict.fromkeys(r.get(f, '') for r in rows if r.get(f, '')))
        out['patient'][f] = vals[0] if vals else ''
        if len(vals) > 1:
            out['warnings'].append(f'{f}: different values in the samplesheet ({" / ".join(vals)})')
    for role, sample in (('tumor', tumor), ('normal', normal)):
        srows = [r for r in rows if r.get('sample') == sample]
        for f in SAMPLE_FIELDS:
            vals = list(dict.fromkeys(r.get(f, '') for r in srows if r.get(f, '')))
            out[role][f] = vals[0] if vals else ''
            if len(vals) > 1:
                out['warnings'].append(f'{sample} {f}: different values in the samplesheet ({" / ".join(vals)})')
    return out


# ── Model ─────────────────────────────────────────────────────────────────────

def build(a):
    terms, label = _amp.tumor_terms_for([a.patient, a.tumor, a.pair], a.tumor_type, a.tumor_types)
    metrics = build_metrics(a)
    gene_ctx = load_gene_ctx(a.gene_ctx)
    germline_panel = load_panel(a.germline_panel)
    somatic_panel = load_panel(a.somatic_panel)

    germline = load_germline(a.germline)
    core = load_somatic(a.somatic, a.civic, terms, 'core', gene_ctx)
    tso = load_somatic(a.somatic_tso500, a.civic_tso500, terms, 'TSO500', gene_ctx)
    # Main table: core panel + TSO500 Tier I/II; TSO500 Tier III apart; Tier IV
    # (common polymorphisms, probable germline) apart from the somatic findings
    somatic_main = [v for v in core + tso if v['_AMP_TIER'] != 'IV'
                    and not (v['_TRACK'] == 'TSO500' and v['_AMP_TIER'] == 'III')]
    tier_iv = [v for v in core + tso if v['_AMP_TIER'] == 'IV']
    tso500_vus = [v for v in tso if v['_AMP_TIER'] == 'III']
    second_hit = [v for v in read_tsv(a.somatic_second_hit)[1] if not _amp.population_common(v.get)]
    actionable = [v.get('ANN[*].GENE', '') for v in somatic_main if v['_AMP_TIER'] in rs.ACTIONABLE_TIERS]

    cnv = harmonize_cnv(a.cnv_somatic, a.cnv_gatk_somatic, a.somatic_cnv_panel,
                        a.germline_panel, THRESHOLDS['somatic_cnv_amp_cn'])
    relevant = loh_relevant_genes(a.germline or '', [a.somatic or '', a.somatic_second_hit or ''])
    cnv['loh_relevant'] = sorted(r['gene'] for r in cnv['loh'] if r['gene'] in relevant)
    cnv['loh_hidden'] = sorted(r['gene'] for r in cnv['loh'] if r['gene'] not in relevant)

    corr = {}
    if a.corr_json:
        with open(a.corr_json) as fh:
            corr = json.load(fh)
    corr_html = ''
    if a.corr_html:
        with open(a.corr_html) as fh:
            corr_html = fh.read()

    queries = rs.trial_queries(corr, actionable, metrics['tmb'], metrics['msi_status'],
                               metrics['hrd_status'], label)
    trials = rs.clinical_trials(queries, metrics['msi_status']) if queries else None
    if trials is not None:
        trials['queries'] = [q[0] for q in queries]

    qc = {'purity': metrics['purity'], 'gof': metrics['gof'],
          'ploidy': metrics['ploidy'], 'wgd': metrics['wgd']}
    keys_paths = {
        'germline_snv': [a.germline], 'germline_cnv': [a.cnv_germline_gatk], 'pharmaco': [a.pharmaco],
        'somatic_snv': [a.somatic], 'somatic_cnv': [a.cnv_somatic, a.cnv_gatk_somatic],
        'somatic_cnv_gatk': [a.cnv_gatk_somatic], 'tmb': [a.tmb],
        'signatures': [a.sig_sbs96, a.sig_id], 'hrd': [a.hrd],
        'coverage': [a.gene_cov_n, a.gene_cov_t],
    }
    info = rs.report_info()
    ident = identification(a.samplesheet, a.tumor, a.normal)
    return {
        'schema_version': SCHEMA_VERSION,
        'ids': {'pair': a.pair, 'tumor': a.tumor, 'normal': a.normal,
                'patient': ident['patient']['patient'] or a.patient,
                'run_name': a.run_name, 'date': date.today().isoformat()},
        'identification': ident,
        'versions': {'pipeline': rs.pipeline_label(), 'sarek': rs.sarek_label(),
                     'rows': rs.version_rows(), 'tools': rs.sarek_tool_versions(),
                     'info': info},
        'tumor_type': {'terms': sorted(terms), 'label': label},
        'thresholds': THRESHOLDS,
        'not_run': sorted(rs.not_run()),
        'unavailable': {k: rs.unavailable_text(k, *p) for k, p in keys_paths.items()},
        'metrics': metrics,
        'alerts': {'ascat': ascat_alerts(qc), 'ploidy_note': ploidy_note(qc),
                   'tert_promoter': _amp.tert_promoter_note(a.capture_bed)},
        'germline': germline,
        'germline_cnv': germline_cnv(a.cnv_germline_gatk, set(germline_panel)),
        'somatic': somatic_main,
        'somatic_tier_iv': tier_iv,
        'somatic_tso500_vus': tso500_vus,
        'somatic_second_hit': second_hit,
        'actionable_genes': list(dict.fromkeys(actionable)),
        'somatic_cnv': cnv,
        'signatures': {'SBS96': signatures(a.sig_sbs96, a.pair, 'SBS'),
                       'DBS': signatures(a.sig_dbs, a.pair, 'DBS'),
                       'ID': signatures(a.sig_id, a.pair, 'ID')},
        'correlation': corr,
        'correlation_html': corr_html,
        'pharmaco': pharmaco(a.pharmaco),
        'coverage': {'normal': coverage(a.gene_cov_n), 'tumor': coverage(a.gene_cov_t),
                     'germline_panel': germline_panel, 'somatic_panel': somatic_panel},
        'tables': {'tmb': read_tsv(a.tmb), 'hrd': read_tsv(a.hrd)},
        'trials': trials,
    }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('pair', 'run-name', 'patient', 'tumor', 'normal'):
        p.add_argument('--' + name, required=True)
    for name in ('germline', 'somatic', 'civic', 'somatic-second-hit', 'somatic-tso500', 'civic-tso500',
                 'capture-bed', 'tumor-type', 'tumor-types', 'tmb', 'msi', 'ascat', 'hrd', 'cov-n', 'cov-t',
                 'cnv-somatic', 'cnv-gatk-somatic', 'cnv-germline-gatk', 'somatic-cnv-panel',
                 'sig-sbs96', 'sig-dbs', 'sig-id', 'gene-ctx', 'gene-cov-n', 'gene-cov-t',
                 'germline-panel', 'somatic-panel', 'pharmaco', 'corr-json', 'corr-html', 'samplesheet'):
        p.add_argument('--' + name, default='')
    p.add_argument('--out', required=True)
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    model = build(a)
    with open(a.out, 'w') as fh:
        json.dump(model, fh, indent=1, ensure_ascii=False, default=sorted)
    print(f'[INFO] report model written: {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()

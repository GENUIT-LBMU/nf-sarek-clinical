"""
amp_tier.py — AMP/ASCO/CAP 2017 somatic tiering from CIViC evidence.

Single source for the HTML report, the PDF report and its findings card.

Only accepted evidence whose direction is SUPPORTS sets the tier. Submitted
evidence and DOES_NOT_SUPPORT records are listed, never tier-defining.

    CIViC level   same tumor type   other tumor type
    A             IA                IIC
    B             IB                IIC
    C             IIC               —
    D, E          IID               —

"Same tumor type" = the CIViC disease contains one of the patient's terms
(--tumor_type / --tumor_types; 'Breast Cancer' covers 'HER2 Positive Breast
Cancer') or is a pan-cancer term. Without a tumor type
every disease counts as "other", so the highest tier is IIC.
A CIViC assertion with an AMP category for the tumor type gives that tier.

Tier IV (benign / likely benign): a common population polymorphism, i.e. gnomAD
FAF95 (else AF) >= filters.somatic_snv.max_pop_af (1%), marked POP_COMMON=1 by
SNPSIFT_EXTRACT_SOMATIC. It wins over any CIViC evidence and the variant is kept
apart from the somatic findings (probable germline origin).
"""

import csv
import os
import re
from urllib.parse import unquote

TIER_ORDER = {'IA': 0, 'IB': 1, 'IIC': 2, 'IID': 3, 'III': 4, 'IV': 5}

PAN_CANCER = {'cancer', 'solid tumor', 'solid tumors', 'solid tumors, advanced',
              'advanced solid tumor', 'pan-cancer'}

_SAME_TUMOR  = {'A': 'IA', 'B': 'IB', 'C': 'IIC', 'D': 'IID', 'E': 'IID'}
_OTHER_TUMOR = {'A': 'IIC', 'B': 'IIC'}

_SIG_LABEL = {
    'SENSITIVITYRESPONSE': 'Sensitivity/Response',
    'RESISTANCE':          'Resistance',
    'REDUCED_SENSITIVITY': 'Reduced sensitivity',
    'ADVERSE_RESPONSE':    'Adverse response',
    'POOR_OUTCOME':        'Poor outcome',
    'BETTER_OUTCOME':      'Better outcome',
    'POSITIVE':            'Diagnostic (positive)',
    'NEGATIVE':            'Diagnostic (negative)',
    'ONCOGENIC':           'Oncogenic',
}


def clean(value):
    """Decode CIViC's URL-escaped text (e.g. 'Solid Tumors%2C Advanced')."""
    return unquote(str(value or '')).strip()


def _terms(spec):
    return {t.strip().lower() for t in clean(spec).split('|') if t.strip()}


def tumor_terms_for(ids, tumor_type='', tumor_types_tsv=''):
    """
    Disease terms for one patient. `ids` are the names the TSV may use
    (patient, tumor sample, pair). The TSV (id<TAB>term|term…) wins over the
    run-wide --tumor_type. Returns (set_of_terms, display_string).
    """
    if tumor_types_tsv and os.path.isfile(tumor_types_tsv):
        wanted = {str(i) for i in ids if i}
        with open(tumor_types_tsv) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[0].strip() in wanted:
                    return _terms(parts[1]), clean(parts[1])
    if tumor_type:
        return _terms(tumor_type), clean(tumor_type)
    return set(), ''


def _record(r):
    return {
        'gene':         clean(r.get('SYMBOL')),
        'profile':      clean(r.get('CIViC Molecular Profile Name') or r.get('CIViC Variant Name')),
        'entity_type':  clean(r.get('CIViC Entity Type')).lower(),
        'status':       clean(r.get('CIViC Entity Status')).lower(),
        'direction':    clean(r.get('CIViC Entity Direction')).upper(),
        'significance': clean(r.get('CIViC Entity Significance')).upper(),
        'disease':      re.sub(r'\s*\(DOID \d+\)$', '', clean(r.get('CIViC Entity Disease'))),
        'therapies':    re.sub(r'\s*\(NCIt ID [^)]*\)', '', clean(r.get('CIViC Entity Therapies'))).replace('&', ' + '),
        'level':        clean(r.get('CIViC Evidence Level')).upper(),
        'amp_category': clean(r.get('CIViC Assertion AMP Category')).upper(),
    }


def load_civic(path):
    """(CHROM, POS, REF, ALT) -> list of normalised CIViC records."""
    index = {}
    if not path or not os.path.isfile(path):
        return index
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            key = (r.get('CHROM', ''), r.get('POS', ''), r.get('REF', ''), r.get('ALT', ''))
            index.setdefault(key, []).append(_record(r))
    return index


def same_tumor(rec, terms):
    """Pan-cancer, or a tumor-type term contained in the disease name, so that
    'Breast Cancer' also covers 'HER2 Positive Breast Cancer'."""
    d = rec['disease'].lower()
    return d in PAN_CANCER or any(t in d for t in terms)


def _assertion_tier(category):
    m = re.search(r'TIER_(I{1,3}|IV)_LEVEL_([A-D])', category)
    if not m:
        return None
    tier = m.group(1) + m.group(2) if m.group(1) in ('I', 'II') else None
    return tier if tier in TIER_ORDER else None


def record_tier(rec, terms):
    """Tier this record supports, or None when it cannot define a tier."""
    if rec['status'] != 'accepted' or rec['direction'] != 'SUPPORTS':
        return None
    same = same_tumor(rec, terms)
    if rec['entity_type'] == 'assertion' and rec['amp_category']:
        tier = _assertion_tier(rec['amp_category'])
        if tier and same:
            return tier
        return 'IIC' if tier in ('IA', 'IB') else None
    return (_SAME_TUMOR if same else _OTHER_TUMOR).get(rec['level'])


def classify(records, terms):
    """Returns (tier, defining_record_or_None)."""
    best, best_rec = 'III', None
    level_rank = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}
    for rec in sorted(records or [], key=lambda r: level_rank.get(r['level'], 9)):
        tier = record_tier(rec, terms)
        if tier and TIER_ORDER[tier] < TIER_ORDER[best]:
            best, best_rec = tier, rec
    return best, best_rec


def population_common(get_field):
    """True for a common polymorphism (POP_COMMON=1): AMP Tier IV."""
    return (get_field('POP_COMMON') or '').strip() == '1'


def variant_tier(records, terms, get_field):
    """(tier, defining_record_or_None): IV for common polymorphisms, else classify()."""
    if population_common(get_field):
        return 'IV', None
    return classify(records, terms)


def significance_label(sig):
    return _SIG_LABEL.get(sig, sig.replace('_', ' ').capitalize() if sig else 'Other')


def therapy_cell(rec):
    """'Cetuximab (Resistance)' for the tier-defining record of a table row."""
    if not rec:
        return ''
    what = rec['therapies'] or rec['profile'] or ''
    return f"{what} ({significance_label(rec['significance'])})" if what else significance_label(rec['significance'])


def evidence_summary(records, terms):
    """
    Evidence grouped for display:
      {'groups': [(label, [item, …]), …], 'not_supporting': [item, …],
       'n_submitted': int}
    Each item: 'Therapy (B, Colorectal Cancer)' — accepted SUPPORTS only in
    groups; DOES_NOT_SUPPORT accepted in not_supporting; submitted counted.
    """
    groups, not_supporting, n_submitted = {}, [], 0
    for rec in records or []:
        if rec['status'] != 'accepted':
            n_submitted += 1
            continue
        what = rec['therapies'] or rec['profile'] or '—'
        flag = '' if same_tumor(rec, terms) else ', other tumor'
        item = f"{what} ({rec['level']}, {rec['disease'] or '?'}{flag})"
        if rec['direction'] != 'SUPPORTS':
            not_supporting.append(f"{significance_label(rec['significance'])}: {item}")
            continue
        items = groups.setdefault(significance_label(rec['significance']), [])
        if item not in items:
            items.append(item)
    order = sorted(groups.items(), key=lambda kv: kv[0])
    return {'groups': order, 'not_supporting': not_supporting, 'n_submitted': n_submitted}

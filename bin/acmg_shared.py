"""
ACMG germline classification — shared module.

Single source of truth for:
  * 9-tool in silico predictor panel (dbNSFP)
  * ACMG criteria scoring (ClinVar, BA1, PVS1, PM2, PP3/BP4, BS1)
  * Resolved REVEL vs in silico panel (no double counting per ClinGen 2022)
  * MMR / HRR gene sets used for germline pathway classification

Public API
----------
Field access is abstracted through a `get_field(name) -> str` callable, so the
same logic serves three different row representations in the codebase:
  1. (row: list, hdr: dict[str,int])      — final_report_acmg.py, snpsift_germline heredoc
  2. row: dict[str,str]                   — clinical_pdf_report.py, germline_somatic_correlate.py

Helpers `getter_for_dict()` and `getter_for_row_hdr()` build the callable from
each representation.

Functions
---------
count_insil(get_field)            -> (dam:int, ben:int)
gnomad_af(get_field)              -> float | None  (None = frequency unknown)
gnomad_supports_pm2(get_field)    -> bool  (False for sites failing gnomAD QC)
acmg_criteria_list(get_field)     -> (cls:str, criteria:list[str])
acmg_class(get_field)             -> str   (P / LP / VUS / LB / B)

Thresholds are read from environment variables (set by Nextflow from
conf/filters.config) with hardcoded ACMG defaults as fallback.
"""

import os as _os

_E = _os.environ.get


# =============================================================================
# Field-getter adapters
# =============================================================================
def getter_for_dict(row_dict):
    """Build a get_field callable from a row dict."""
    def _g(name):
        return str(row_dict.get(name, '') or '').strip()
    return _g


def getter_for_row_hdr(row, hdr):
    """Build a get_field callable from a (row, hdr) pair where row is a list
    and hdr maps column-name -> column-index."""
    def _g(name):
        i = hdr.get(name, -1)
        return row[i].strip() if 0 <= i < len(row) else '.'
    return _g


# =============================================================================
# In silico predictor panel (9 tools from dbNSFP)
# =============================================================================
INSIL_DAMAGING = {
    'dbNSFP_SIFT_pred':             ('D',),
    'dbNSFP_SIFT4G_pred':           ('D',),
    'dbNSFP_Polyphen2_HDIV_pred':   ('D', 'P'),  # P = possibly damaging
    'dbNSFP_Polyphen2_HVAR_pred':   ('D', 'P'),
    'dbNSFP_LRT_pred':              ('D',),
    'dbNSFP_MutationTaster_pred':   ('A', 'D'),  # A = disease causing automatic
    'dbNSFP_MutationAssessor_pred': ('H', 'M'),  # H = high, M = medium impact
    'dbNSFP_FATHMM_pred':           ('D',),
    'dbNSFP_PROVEAN_pred':          ('D',),
}

INSIL_BENIGN = {
    'dbNSFP_SIFT_pred':             ('T',),
    'dbNSFP_SIFT4G_pred':           ('T',),
    'dbNSFP_Polyphen2_HDIV_pred':   ('B',),
    'dbNSFP_Polyphen2_HVAR_pred':   ('B',),
    'dbNSFP_LRT_pred':              ('N',),
    'dbNSFP_MutationTaster_pred':   ('N', 'P'),  # N = neutral, P = polymorphism auto
    'dbNSFP_MutationAssessor_pred': ('L', 'N'),  # L = low, N = neutral
    'dbNSFP_FATHMM_pred':           ('T',),
    'dbNSFP_PROVEAN_pred':          ('N',),
}


def count_insil(get_field):
    """Count damaging vs benign votes across the 9-tool in silico panel.

    Only predictors with a clear D or B verdict are counted; missing values
    ('.', 'NA', '') abstain. dbNSFP stores multi-transcript predictions as
    comma-separated values; the first (primary transcript) is used.

    Returns (damaging_count, benign_count).
    """
    dam = ben = 0
    for field in INSIL_DAMAGING:
        v = get_field(field).split(',')[0].strip().upper()
        if v in ('', '.', 'NA'):
            continue
        if v in INSIL_DAMAGING[field]:
            dam += 1
        elif v in INSIL_BENIGN.get(field, ()):
            ben += 1
    return dam, ben


# =============================================================================
# ClinVar
# =============================================================================
def clinvar_sig(get_field):
    """ClinVar germline classification used by the ACMG engine and filters.

    Reads the canonical `ClinVar_sig` column written by SNPSIFT_EXTRACT_GERMLINE
    (local ClinVar release with --clinvar_tsv, else dbNSFP 4.1a). Falls back to
    `dbNSFP_clinvar_clnsig` only when the column is missing (older TSVs).
    """
    return (get_field('ClinVar_sig') or get_field('dbNSFP_clinvar_clnsig') or '').strip()


def clinvar_display(get_field):
    """'Benign/Likely benign (2*)' for reports; '-' when not in ClinVar."""
    sig = clinvar_sig(get_field)
    if sig in ('', '.'):
        return '-'
    st = (get_field('ClinVar_stars') or '').strip()
    txt = sig.replace('_', ' ')
    return f'{txt} ({st}*)' if st not in ('', '.') else txt


# =============================================================================
# Population frequency
# =============================================================================
def gnomad_af(get_field):
    """Population allele frequency used by BA1 / BS1 / PM2.

    Reads the canonical `gnomAD_AF` column written by SNPSIFT_EXTRACT_GERMLINE
    (local gnomAD v4.1 resource, or the legacy dbNSFP/web-API value) and falls
    back to `dbNSFP_gnomAD_exomes_AF` for tables produced before that column
    existed.

    Returns a float, or None when the frequency is UNKNOWN (local mode: site
    outside the resource regions with no dbNSFP value). Unknown must not be
    read as "absent": callers skip BA1, BS1 and PM2 in that case. In legacy
    mode a missing value keeps the historical behaviour (0.0, i.e. absent).
    See gnomad_supports_pm2() for sites failing gnomAD QC.
    """
    src = get_field('gnomAD_AF_source')
    if src == 'unknown':
        return None
    for name in ('gnomAD_AF', 'dbNSFP_gnomAD_exomes_AF'):
        vals = []
        for x in get_field(name).split(','):   # "0.53,0.53" = one value per transcript
            try:
                vals.append(float(x))
            except (TypeError, ValueError):
                pass
        if vals:
            return max(vals)
    return 0.0


def gnomad_supports_pm2(get_field):
    """False when the frequency comes from a site that fails gnomAD QC in both
    exomes and genomes (source *_qc_fail): low AF there reflects poor calling
    (paralogs, repeats), not rarity. Such AFs still count for BA1/BS1."""
    return not get_field('gnomAD_AF_source').endswith('_qc_fail')


# =============================================================================
# ACMG classification
# =============================================================================
def clinvar_stars(get_field):
    """Review stars of the ClinVar classification (0-4), or None when unknown.

    None (no `ClinVar_stars` column, older TSVs) keeps the pre-threshold
    behaviour: the ClinVar class is accepted.
    """
    st = (get_field('ClinVar_stars') or '').strip()
    try:
        return int(st)
    except ValueError:
        return None


def clinvar_class(clnsig):
    """ClinVar significance -> 'P' / 'LP' / 'VUS' / 'LB' / 'B', or None.

    None = not usable as a class (absent, conflicting, risk_factor, ...).
    """
    clnsig = (clnsig or '').lower()
    if not clnsig or clnsig in ('.', 'na', 'nan'):
        return None
    # "Conflicting_classifications_of_pathogenicity" must NOT be read as P
    if 'conflicting' in clnsig:
        return None
    if ('pathogenic' in clnsig
            and 'likely_pathogenic' not in clnsig
            and 'likely pathogenic' not in clnsig):
        return 'P' if 'likely' not in clnsig else None
    if 'likely_pathogenic' in clnsig or 'likely pathogenic' in clnsig:
        return 'LP'
    if 'uncertain' in clnsig:
        return 'VUS'
    if 'likely_benign' in clnsig or 'likely benign' in clnsig:
        return 'LB'
    if 'benign' in clnsig and 'pathogenic' not in clnsig:
        return 'B'
    return None


def clinvar_is_final(get_field):
    """True when the ClinVar review stars reach `clinvar_min_stars` (default 2)."""
    st = clinvar_stars(get_field)
    return st is None or st >= int(_E('ACMG_CLINVAR_MIN_STARS', '2'))


def acmg_criteria_list(get_field):
    """Classify a germline variant into P / LP / VUS / LB / B with ACMG criteria.

    Logic priority:
      1. ClinVar class, final only with >= clinvar_min_stars review stars:
           P / LP / LB / B  -> returned as is
           VUS              -> VUS, unless BA1 / BS1 bring it down to B / LB
         Below the threshold:
           P / LP  -> ACMG score, never below VUS (except BA1), criteria
                      carry 'ClinVar:LP(1*)' for manual review
           others  -> ACMG score
         Conflicting / other ClinVar terms -> ACMG score
      2. gnomAD BA1 rule (>5% AF -> Benign)
      3. Evidence scoring: PVS1 (HIGH impact) + PM2 (rare) + PP3/BP4 (REVEL OR
         9-tool panel, RESOLVED — only the strongest evidence is counted, no
         double counting per ClinGen 2022) + BS1 (>1% AF)

    Returns (cls, criteria_list) where criteria_list contains the activated
    ACMG evidence code strings (e.g. ['PVS1', 'PM2', 'PP3(7/9)']).
    """
    cv = clinvar_class(clinvar_sig(get_field))
    if cv is None:
        return _acmg_scored(get_field)
    if clinvar_is_final(get_field):
        if cv != 'VUS':
            return cv, ['ClinVar:' + cv]
        cls, criteria = _acmg_scored(get_field)
        if cls in ('LB', 'B') and ('BA1' in criteria or 'BS1' in criteria):
            return cls, ['ClinVar:VUS'] + criteria
        return 'VUS', ['ClinVar:VUS']
    cls, criteria = _acmg_scored(get_field)
    if cv in ('P', 'LP'):
        criteria = criteria + ['ClinVar:{}({}*)'.format(cv, clinvar_stars(get_field))]
        if cls in ('LB', 'B') and 'BA1' not in criteria:
            cls = 'VUS'
    return cls, criteria


def _acmg_scored(get_field):
    """BA1 + evidence scoring (no ClinVar)."""
    # ── 2. Population frequency — BA1 (common polymorphism) ────────────────
    gnomad = gnomad_af(get_field)          # None = unknown → no BA1/BS1/PM2
    if gnomad is not None and gnomad > float(_E('ACMG_BA1_AF', '0.05')):
        return 'B', ['BA1']

    # ── 3. Evidence scoring ────────────────────────────────────────────────
    score = 0
    criteria = []

    impact = get_field('ANN[*].IMPACT').upper()
    if impact == 'HIGH':
        score += int(_E('ACMG_PVS1_SCORE', '8'))
        criteria.append('PVS1')

    if (gnomad is not None and gnomad_supports_pm2(get_field)
            and gnomad < float(_E('ACMG_PM2_AF', '0.001'))):
        score += int(_E('ACMG_PM2_SCORE', '2'))
        criteria.append('PM2')

    # REVEL evidence
    revel_dir = 0
    revel_str = 0
    revel_lbl = ''
    try:
        revel = float(get_field('dbNSFP_REVEL_score'))
        if revel >= float(_E('ACMG_REVEL_PP3_STRONG', '0.932')):
            revel_dir, revel_str, revel_lbl = 1, 4, 'PP3_St'
        elif revel >= float(_E('ACMG_REVEL_PP3', '0.644')):
            revel_dir, revel_str, revel_lbl = 1, 2, 'PP3'
        elif revel >= float(_E('ACMG_REVEL_PP3_SUP', '0.5')):
            revel_dir, revel_str, revel_lbl = 1, 1, 'PP3_Sup'
        elif revel < float(_E('ACMG_REVEL_BP4_STRONG', '0.15')):
            revel_dir, revel_str, revel_lbl = -1, 4, 'BP4_St'
        elif revel < float(_E('ACMG_REVEL_BP4', '0.29')):
            revel_dir, revel_str, revel_lbl = -1, 2, 'BP4'
    except (ValueError, TypeError):
        pass

    # In silico panel evidence
    dam, ben = count_insil(get_field)
    total = dam + ben
    panel_dir = 0
    panel_str = 0
    panel_lbl = ''
    insil_maj = float(_E('ACMG_INSIL_MAJORITY', '0.75'))
    if total > 0:
        if ben > dam:
            panel_dir = -1
            ratio = ben / total
            if ratio >= insil_maj:
                panel_str, panel_lbl = 2, 'BP4({}/{})'.format(ben, total)
            else:
                panel_str, panel_lbl = 1, 'BP4_Sup({}/{})'.format(ben, total)
        elif dam > ben:
            panel_dir = 1
            ratio = dam / total
            if ratio >= insil_maj:
                panel_str, panel_lbl = 2, 'PP3({}/{})'.format(dam, total)
            else:
                panel_str, panel_lbl = 1, 'PP3_Sup({}/{})'.format(dam, total)

    # Resolve REVEL vs panel — never double-count (ClinGen 2022)
    if revel_dir != 0 and panel_dir != 0:
        if revel_dir == panel_dir:
            if revel_str >= panel_str:
                score += revel_dir * revel_str
                criteria.append(revel_lbl)
            else:
                score += panel_dir * panel_str
                criteria.append(panel_lbl)
        # opposite directions -> conflicting in silico evidence, both abstain
    elif revel_dir != 0:
        score += revel_dir * revel_str
        criteria.append(revel_lbl)
    elif panel_dir != 0:
        score += panel_dir * panel_str
        criteria.append(panel_lbl)

    # BS1 — frequency above benign threshold
    if gnomad is not None and gnomad > float(_E('ACMG_BS1_AF', '0.01')):
        score += int(_E('ACMG_BS1_SCORE', '-4'))
        criteria.append('BS1')

    sp = int(_E('ACMG_SCORE_P', '10'))
    slp = int(_E('ACMG_SCORE_LP', '6'))
    slb = int(_E('ACMG_SCORE_LB', '-2'))
    if score >= sp:
        return 'P', criteria
    if score >= slp:
        return 'LP', criteria
    if score <= slb:
        return 'LB', criteria
    return 'VUS', criteria


def acmg_class(get_field):
    """Thin wrapper — returns only the classification string (P/LP/VUS/LB/B)."""
    cls, _criteria = acmg_criteria_list(get_field)
    return cls


# =============================================================================
# Gene pathway sets (for germline-somatic pathway classification)
# =============================================================================
MMR_GENES = frozenset({
    'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM', 'MSH3', 'MLH3', 'PMS1',
})

HRR_GENES = frozenset({
    'BRCA1', 'BRCA2', 'PALB2', 'RAD51C', 'RAD51D', 'ATM', 'CHEK2',
    'BARD1', 'BRIP1', 'NBN', 'RAD50', 'MRE11', 'FANCA', 'FANCC',
    'FANCD2', 'FANCE', 'FANCF', 'FANCG', 'FANCL', 'FANCM',
    'CDK12', 'ABRAXAS1', 'ABRAXAS2',
    'SLX4', 'RBBP8', 'RTEL1',
})


# =============================================================================
# Mutational signature sets (env-overridable)
# =============================================================================
MMR_SBS_SIGNATURES = frozenset(_E('SIG_MMR_SBS', 'SBS6,SBS15,SBS20,SBS21,SBS26,SBS44').split(','))
HRD_SBS_SIGNATURES = frozenset(_E('SIG_HRD_SBS', 'SBS3,SBS8').split(','))
MMR_ID_SIGNATURES  = frozenset(_E('SIG_MMR_ID',  'ID1,ID2').split(','))

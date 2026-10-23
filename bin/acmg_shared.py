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
gnomad_faf95(get_field)           -> float | None  (max FAF95 across ancestry groups)
benign_freq(get_field)            -> float | None  (BA1/BS1 frequency: FAF95, else AF)
gnomad_display(get_field)         -> str   (report text)
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


def gnomad_faf95(get_field):
    """Max filtering allele frequency (95% CI) across gnomAD ancestry groups.

    Reads the canonical `gnomAD_FAF95` column (local gnomAD resource). None when
    it is missing / "." (legacy mode, outside the resource regions), or when
    frequency_metric is 'af': callers then use gnomad_af().
    """
    if _E('ACMG_FREQ_METRIC', 'faf95').strip().lower() == 'af':
        return None
    try:
        return float(get_field('gnomAD_FAF95'))
    except (ValueError, TypeError):
        return None


def benign_freq(get_field):
    """Frequency for BA1 / BS1: FAF95 when available, else AF (None = unknown)."""
    faf = gnomad_faf95(get_field)
    return faf if faf is not None else gnomad_af(get_field)


def gnomad_display(get_field):
    """'0.0048 / FAF 0.113' for reports (AF alone when there is no FAF95)."""
    af = (get_field('gnomAD_AF') or get_field('dbNSFP_gnomAD_exomes_AF') or '').strip()
    if af in ('', '.'):
        return '-'
    faf = (get_field('gnomAD_FAF95') or '').strip()
    try:
        af_txt = '{:.4g}'.format(float(af.split(',')[0]))
    except ValueError:
        af_txt = af
    try:
        return '{} / FAF {:.3g}'.format(af_txt, float(faf))
    except ValueError:
        return af_txt


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
      3. Evidence points (Tavtigian 2020): PVS1 by gene mechanism and variant
         type (pvs1_evidence) + PM2_Supporting (rare) + PP3/BP4 from REVEL
         (Pejaver 2022) and SpliceAI (Walker 2023) (pp3_bp4_evidence)
         + PS1 / PM5 from ClinVar at the same residue (ps1_pm5_evidence)
         + BS1 (FAF95 > 1%).
         P >= 10, LP 6..9, VUS 0..5, LB -1..-6, B <= -7

    Returns (cls, criteria_list) where criteria_list contains the activated
    ACMG evidence code strings (e.g. ['PVS1', 'PM2_Sup', 'PP3_Mod']).
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
    faf = gnomad_faf95(get_field)          # None → AF is used
    benign = benign_freq(get_field)        # BA1 / BS1: FAF95, else AF
    if benign is not None and benign > float(_E('ACMG_BA1_AF', '0.05')):
        return 'B', ['BA1']

    # ── 3. Evidence scoring (points: Tavtigian 2020) ─────────────────────
    score = 0
    criteria = []

    pts, code = pvs1_evidence(get_field)
    score += pts
    if code:
        criteria.append(code)

    # PM2_Supporting (ClinGen SVI 2020): rare overall AND in every ancestry group
    pm2_af = float(_E('ACMG_PM2_AF', '0.001'))
    if (gnomad is not None and gnomad_supports_pm2(get_field)
            and gnomad < pm2_af and (faf is None or faf < pm2_af)):
        score += int(_E('ACMG_PM2_SCORE', '1'))
        criteria.append('PM2_Sup')

    # PP3 / BP4: REVEL (Pejaver 2022) and SpliceAI (Walker 2023); not for null
    # variants, whose splicing / LOF effect is already weighed by PVS1
    if code is None:
        pts, code = pp3_bp4_evidence(get_field)
        score += pts
        if code:
            criteria.append(code)
        pts, code = ps1_pm5_evidence(get_field)
        score += pts
        if code:
            criteria.append(code)

    # BS1 — frequency above benign threshold
    if benign is not None and benign > float(_E('ACMG_BS1_AF', '0.01')):
        score += int(_E('ACMG_BS1_SCORE', '-4'))
        criteria.append('BS1')

    if score >= int(_E('ACMG_SCORE_P', '10')):
        return 'P', criteria
    if score >= int(_E('ACMG_SCORE_LP', '6')):
        return 'LP', criteria
    if score <= int(_E('ACMG_SCORE_B', '-7')):
        return 'B', criteria
    if score <= int(_E('ACMG_SCORE_LB', '-1')):
        return 'LB', criteria
    return 'VUS', criteria


# Genes where loss of function is not the cancer-predisposition mechanism
# (oncogenes / gain of function, or LOF causing a non-cancer phenotype).
# Overridden by filters.acmg.pvs1_not_applicable (env ACMG_PVS1_NA_GENES).
PVS1_NA_DEFAULT = (
    'ALK,AR,BRAF,CCND1,CCND2,CCND3,CCNE1,CDK4,CDK6,EGFR,EPCAM,ERBB2,FGFR1,FGFR2,'
    'FGFR3,GREM1,HOXB13,KIF1B,KIT,KRAS,MDM2,MDM4,MET,MITF,MYC,PDGFRA,PIK3CA,'
    'POLD1,POLE,RECQL,RET'
)


def pvs1_not_applicable_genes():
    raw = _E('ACMG_PVS1_NA_GENES', PVS1_NA_DEFAULT)
    return frozenset(g.strip().upper() for g in raw.split(',') if g.strip())


def _int_or_none(x):
    try:
        return int(str(x).split('/')[0])
    except (ValueError, TypeError):
        return None


def pvs1_evidence(get_field):
    """(points, code) for null variants — simplified ClinGen PVS1 tree
    (Abou Tayoun 2018). Only HIGH-impact variants are considered.

      gene in pvs1_not_applicable          -> 0, 'PVS1_NA'
      stop_gained / frameshift:
        not in the last exon               -> PVS1        (8)
        last exon, >= 10% of protein lost  -> PVS1_St     (4)
        last exon, <  10% of protein lost  -> PVS1_Mod    (2)
      splice donor / acceptor (+-1/2)      -> PVS1        (8)
      start_lost                           -> PVS1_Mod    (2)
      stop_lost                            -> PM4         (2)
      other HIGH                           -> PVS1        (8)

    The "last 50 nt of the penultimate exon" NMD rule is not applied (exon
    coordinates are not in the table): such variants keep full PVS1.
    """
    if get_field('ANN[*].IMPACT').upper() != 'HIGH':
        return 0, None
    if get_field('ANN[*].GENE').strip().upper() in pvs1_not_applicable_genes():
        return 0, 'PVS1_NA'
    full = int(_E('ACMG_PVS1_SCORE', '8'))
    strong = int(_E('ACMG_PVS1_STRONG_SCORE', '4'))
    moderate = int(_E('ACMG_PVS1_MODERATE_SCORE', '2'))
    effect = get_field('ANN[*].EFFECT').lower()
    if 'stop_gained' in effect or 'frameshift' in effect:
        rank = get_field('ANN[*].RANK').split('/')
        exon, n_exons = (_int_or_none(rank[0]), _int_or_none(rank[1])) if len(rank) == 2 else (None, None)
        if exon is not None and n_exons and exon == n_exons:
            aa_pos = _int_or_none(get_field('ANN[*].AA_POS'))
            aa_len = _int_or_none(get_field('ANN[*].AA_LEN'))
            min_lost = float(_E('ACMG_PVS1_LAST_EXON_MIN_LOST', '0.10'))
            if aa_pos and aa_len and (aa_len - aa_pos + 1) / aa_len < min_lost:
                return moderate, 'PVS1_Mod'
            return strong, 'PVS1_St'
        return full, 'PVS1'
    if 'splice_acceptor' in effect or 'splice_donor' in effect:
        return full, 'PVS1'
    if 'start_lost' in effect:
        return moderate, 'PVS1_Mod'
    if 'stop_lost' in effect:
        return int(_E('ACMG_PM4_SCORE', '2')), 'PM4'
    return full, 'PVS1'


def revel_score(get_field):
    """REVEL score: canonical `REVEL_score` (local REVEL v1.3, MANE Select
    transcript, or dbNSFP 4.1a), else `dbNSFP_REVEL_score` for older tables."""
    for name in ('REVEL_score', 'dbNSFP_REVEL_score'):
        v = (get_field(name) or '').strip()
        if v not in ('', '.'):
            try:
                return float(v.split(',')[0])
            except ValueError:
                pass
    return None


def spliceai_ds(get_field):
    """SpliceAI max delta score (None when SpliceAI did not run)."""
    try:
        return float(get_field('SpliceAI_DS'))
    except (ValueError, TypeError):
        return None


def revel_evidence(get_field):
    """(points, code) from REVEL with the ClinGen cut-offs (Pejaver 2022);
    (0, None) when REVEL is missing, (0, '') when it is indeterminate."""
    revel = revel_score(get_field)
    if revel is None:
        return 0, None
    for env, default, pts, code in (
            ('ACMG_REVEL_PP3_STRONG', '0.932', 4, 'PP3_St'),
            ('ACMG_REVEL_PP3_MOD',    '0.773', 2, 'PP3_Mod'),
            ('ACMG_REVEL_PP3_SUP',    '0.644', 1, 'PP3_Sup')):
        if revel >= float(_E(env, default)):
            return pts, code
    for env, default, pts, code in (
            ('ACMG_REVEL_BP4_STRONG', '0.016', -4, 'BP4_St'),
            ('ACMG_REVEL_BP4_MOD',    '0.183', -2, 'BP4_Mod'),
            ('ACMG_REVEL_BP4_SUP',    '0.290', -1, 'BP4_Sup')):
        if revel <= float(_E(env, default)):
            return pts, code
    return 0, ''


def pp3_bp4_evidence(get_field):
    """(points, code) for non-null variants — ClinGen SVI (Pejaver 2022 for REVEL,
    Walker 2023 for splicing):

      SpliceAI >= spliceai_pp3 (0.2)       -> PP3_Sup(splice x), unless REVEL gives a
                                              stronger PP3 (never added together)
      missense, REVEL PP3                   -> REVEL PP3 (Sup / Mod / St)
      missense, REVEL BP4                   -> BP4 only if SpliceAI <= spliceai_bp4 (0.1)
                                              or SpliceAI did not run
      no REVEL (in-frame indels, ...)       -> splicing PP3 only
    """
    r_pts, r_code = revel_evidence(get_field)
    ds = spliceai_ds(get_field)
    splice_pp3 = ds is not None and ds >= float(_E('ACMG_SPLICEAI_PP3', '0.2'))
    if splice_pp3 and r_pts < 1:
        return 1, 'PP3_Sup(splice {:.2f})'.format(ds)
    if r_pts < 0 and ds is not None and ds > float(_E('ACMG_SPLICEAI_BP4', '0.1')):
        return 0, ''
    return r_pts, (r_code or '')


def ps1_pm5_evidence(get_field):
    """(points, code) from ClinVar P/LP missense at the same residue
    (columns ClinVar_PS1 / ClinVar_PM5, SNPSIFT_EXTRACT_GERMLINE Step 1.9):

      same amino-acid change, other nucleotide  -> PS1 (+4), only if the variant
                                                   has no splicing prediction
                                                   (SpliceAI < spliceai_pp3)
      other amino-acid change, same residue     -> PM5 (+2), when PS1 does not apply
    """
    if 'missense' not in get_field('ANN[*].EFFECT').lower():
        return 0, None
    has = lambda c: (get_field(c) or '').strip() not in ('', '.')
    ds = spliceai_ds(get_field)
    if has('ClinVar_PS1') and (ds is None or ds < float(_E('ACMG_SPLICEAI_PP3', '0.2'))):
        return int(_E('ACMG_PS1_SCORE', '4')), 'PS1'
    if has('ClinVar_PM5'):
        return int(_E('ACMG_PM5_SCORE', '2')), 'PM5'
    return 0, None


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

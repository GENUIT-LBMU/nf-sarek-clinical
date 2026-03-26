"""
ACMG germline classification helpers for the nf-sarek-clinical final report.
Imported by modules/local/final_report/main.nf at run time via sys.path.
"""

_BENIGN_CRIT = {'BA1', 'BS1', 'BP4', 'BP4_St', 'BP4_Sup'}

# ---------------------------------------------------------------------------
# In silico predictor tables (9 tools from dbNSFP)
# ---------------------------------------------------------------------------
# Values that count as a "damaging" vote for each predictor
_INSIL_DAMAGING = {
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

# Values that count as a "benign/tolerated" vote for each predictor
_INSIL_BENIGN = {
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


def _count_insil(col_fn):
    """Count damaging vs benign votes across all in silico predictors.

    Returns (damaging_count, benign_count, total_with_vote).
    Only predictors that return a clear D or B verdict are counted;
    missing / unknown values ('.', 'NA', '') are abstentions.
    """
    dam = ben = 0
    for field in _INSIL_DAMAGING:
        # dbNSFP stores multi-transcript predictions as comma-separated values;
        # take the first one (primary transcript).
        v = col_fn(field).split(',')[0].strip().upper()
        if v in ('', '.', 'NA'):
            continue
        if v in _INSIL_DAMAGING[field]:
            dam += 1
        elif v in _INSIL_BENIGN.get(field, ()):
            ben += 1
    return dam, ben


def acmg_criteria_list(row, hdr):
    """Returns (cls, criteria_list) where cls is P/LP/VUS/LB/B and
    criteria_list contains the activated ACMG evidence code strings."""
    def col(name):
        i = hdr.get(name, -1)
        return row[i].strip() if i >= 0 and i < len(row) else '.'

    # ── ClinVar annotation (primary — highest weight) ────────────────────────
    clnsig = col('dbNSFP_clinvar_clnsig').lower()
    if clnsig and clnsig != '.':
        # Bug fix: "Conflicting_interpretations_of_pathogenicity" must NOT be
        # classified as P — it contains the substring "pathogenic" but is not P.
        if 'conflicting' in clnsig:
            # Fall through to score-based classification
            pass
        elif ('pathogenic' in clnsig
              and 'likely_pathogenic' not in clnsig
              and 'likely pathogenic' not in clnsig):
            if 'likely' not in clnsig:
                return 'P', ['ClinVar:P']
        elif 'likely_pathogenic' in clnsig or 'likely pathogenic' in clnsig:
            return 'LP', ['ClinVar:LP']
        elif 'uncertain' in clnsig:
            return 'VUS', ['ClinVar:VUS']
        elif 'likely_benign' in clnsig or 'likely benign' in clnsig:
            return 'LB', ['ClinVar:LB']
        elif 'benign' in clnsig and 'pathogenic' not in clnsig:
            return 'B', ['ClinVar:B']

    # ── Population frequency — BA1 (common polymorphism) ────────────────────
    try:
        gnomad = float(col('dbNSFP_gnomAD_exomes_AF'))
    except (ValueError, TypeError):
        gnomad = 0.0
    if gnomad > 0.05:
        return 'B', ['BA1']

    # ── Score-based classification ────────────────────────────────────────────
    score = 0
    criteria = []
    impact = col('ANN[*].IMPACT').upper()
    if impact == 'HIGH':
        score += 8
        criteria.append('PVS1')
    if gnomad < 0.001:
        score += 2
        criteria.append('PM2')

    # -- Computational evidence (PP3 / BP4) --
    # ClinGen SVI: REVEL and the 9-tool panel are evaluated together.
    # If they disagree, neither PP3 nor BP4 is applied.
    revel_dir = 0
    revel_str = 0
    revel_lbl = ''
    try:
        revel = float(col('dbNSFP_REVEL_score'))
        if revel >= 0.932:
            revel_dir = 1; revel_str = 4; revel_lbl = 'PP3_St'
        elif revel >= 0.644:
            revel_dir = 1; revel_str = 2; revel_lbl = 'PP3'
        elif revel >= 0.5:
            revel_dir = 1; revel_str = 1; revel_lbl = 'PP3_Sup'
        elif revel < 0.15:
            revel_dir = -1; revel_str = 4; revel_lbl = 'BP4_St'
        elif revel < 0.29:
            revel_dir = -1; revel_str = 2; revel_lbl = 'BP4'
    except (ValueError, TypeError):
        pass

    dam, ben = _count_insil(col)
    total = dam + ben
    panel_dir = 0
    panel_str = 0
    panel_lbl = ''
    if total > 0:
        if ben > dam:
            panel_dir = -1
            ratio = ben / total
            if ratio >= 0.75:
                panel_str = 2; panel_lbl = 'BP4({}/{})'.format(ben, total)
            else:
                panel_str = 1; panel_lbl = 'BP4_Sup({}/{})'.format(ben, total)
        elif dam > ben:
            panel_dir = 1
            ratio = dam / total
            if ratio >= 0.75:
                panel_str = 2; panel_lbl = 'PP3({}/{})'.format(dam, total)
            else:
                panel_str = 1; panel_lbl = 'PP3_Sup({}/{})'.format(dam, total)

    if revel_dir != 0 and panel_dir != 0:
        if revel_dir == panel_dir:
            if revel_str >= panel_str:
                score += revel_dir * revel_str; criteria.append(revel_lbl)
            else:
                score += panel_dir * panel_str; criteria.append(panel_lbl)
    elif revel_dir != 0:
        score += revel_dir * revel_str; criteria.append(revel_lbl)
    elif panel_dir != 0:
        score += panel_dir * panel_str; criteria.append(panel_lbl)

    if gnomad > 0.01:
        score -= 4
        criteria.append('BS1')

    if score >= 10:
        return 'P', criteria
    if score >= 6:
        return 'LP', criteria
    if score <= -2:
        return 'LB', criteria
    return 'VUS', criteria


def acmg_class(row, hdr):
    """Thin wrapper — returns only the classification string."""
    cls, _ = acmg_criteria_list(row, hdr)
    return cls

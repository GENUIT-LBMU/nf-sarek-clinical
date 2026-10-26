"""
ACMG germline classification — shared module.

Single source of truth for:
  * ACMG criteria scoring (ClinVar, BA1, PVS1, PM2, PP3/BP4, BS1)
  * PP3 / BP4 from REVEL (Pejaver 2022) and SpliceAI (Walker 2023)
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
gnomad_af(get_field)              -> float | None  (None = frequency unknown)
gnomad_faf95(get_field)           -> float | None  (max FAF95 across ancestry groups)
benign_freq(get_field)            -> float | None  (BA1/BS1 frequency: FAF95, else AF)
gnomad_display(get_field)         -> str   (report text)
gnomad_supports_pm2(get_field)    -> bool  (False for sites failing gnomAD QC)
vcep_frequency(get_field)         -> (spec, freqs) ClinGen VCEP gene rules or (None, None)
ba1_reason(get_field)             -> str   (BA1 text for the Step 2 filter, '' = not met)
pm1_evidence(get_field)           -> (points, code)  PM1 for VCEP genes
acmg_spec(get_field)              -> str   (rule set used: ClinGen GNxxx or general)
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
# ClinVar
# =============================================================================
def clinvar_sig(get_field):
    """ClinVar germline classification used by the ACMG engine and filters.

    Reads the canonical `ClinVar_sig` column written by SNPSIFT_EXTRACT_GERMLINE
    from the local ClinVar release (--clinvar_tsv).
    """
    return (get_field('ClinVar_sig') or '').strip()


def clinvar_display(get_field):
    """'Benign/Likely benign (2*)' for reports; '-' when not in ClinVar."""
    sig = clinvar_sig(get_field)
    if sig in ('', '.'):
        return '-'
    st = (get_field('ClinVar_stars') or '').strip()
    txt = sig.replace('_', ' ')
    return f'{txt} ({st}*)' if st not in ('', '.') else txt


def hgvs_display(get_field):
    """
    HGVS for reports: (hgvs_c, transcript, hgvs_p, plus). `plus` lists the MANE
    Plus Clinical HGVS of the same gene ('NM_001128425.2: p.Gly396Asp') when it
    differs from the one shown (SNPEFF_MANE, INFO/MANE_PLUS_CLINICAL).
    """
    def clean(x):
        x = (x or '').strip()
        return '' if x == '.' else x
    gene, tx = clean(get_field('ANN[*].GENE')), clean(get_field('ANN[*].FEATUREID'))
    hc, hp = clean(get_field('ANN[*].HGVS_C')), clean(get_field('ANN[*].HGVS_P'))
    plus = []
    for entry in clean(get_field('MANE_PLUS_CLINICAL')).split('&'):
        f = entry.split('|')
        if len(f) != 4 or f[0] != gene or f[1] == tx:
            continue
        shown, other = (hp, f[3]) if f[3] else (hc, f[2])
        if other and other != shown:
            plus.append(f'{f[1]}: {other}')
    return hc, tx, hp, plus


# =============================================================================
# Population frequency
# =============================================================================
def gnomad_af(get_field):
    """Population allele frequency used by BA1 / BS1 / PM2.

    Reads the canonical `gnomAD_AF` column written by SNPSIFT_EXTRACT_GERMLINE
    from the local gnomAD v4.1 resource (--gnomad_tsv).

    Returns a float, or None when the frequency is UNKNOWN (site outside the
    resource regions). Unknown must not be read as "absent": callers skip BA1,
    BS1 and PM2 in that case. See gnomad_supports_pm2() for sites failing
    gnomAD QC.
    """
    src = get_field('gnomAD_AF_source')
    if src == 'unknown':
        return None
    for name in ('gnomAD_AF',):
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
    af = (get_field('gnomAD_AF') or '').strip()
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
# ClinGen VCEP gene specifications (frequency rules + PM1)
# =============================================================================
# Values copied verbatim from the ClinGen Criteria Specification Registry
# (https://cspec.genome.network, API RuleSet/<id>?detail=high), Released
# specifications only, downloaded VCEP_SPECS_DATE. Frequencies refer to gnomAD
# v4.1 joint: 'faf' = grpmax FAF95 (gnomAD_FAF95; 0 when gnomAD computes none,
# e.g. AC = 1), 'grpaf' / 'grpac' / 'grpan' = AF / AC / AN of the group with the
# highest AF (gnomAD_AF_GRPMAX, ...). gnomAD v4 grpmax and FAF exclude the
# bottlenecked groups (ASJ, FIN, AMI, remaining), as the specifications ask.
# Not automated (documented in README): BRCA1/2 use the non-cancer gnomAD
# v2.1/v3.1 subsets (here: AC = 0 in v4.1); MMR BA1/BS1 exclude founder
# pathogenic variants (here: ClinVar P/LP >= clinvar_min_stars is never
# discarded by BA1 and keeps its ClinVar class); TP53 CHIP recalculation by VAF.
# Genes not listed keep the general rules (ba1_af / bs1_af / pm2_af).
VCEP_SPECS_DATE = '2026-10-23'

_MMR_PM2 = lambda f: f['af'] < 0.00002
_ABSENT = lambda f: f['ac'] == 0

VCEP_SPECS = {
    'APC':    dict(spec='GN089 v2.1', ba1=('faf', '>=', 0.001), bs1=('faf', '>=', 0.00001),
                   pm2=lambda f: (f['af'] <= 0.000003) if f['ac'] > 1 else (f['af'] < 0.00001)),
    'ATM':    dict(spec='GN020 v1.6', ba1=('faf', '>', 0.005), bs1=('faf', '>', 0.0005),
                   pm2=lambda f: f['grpaf'] <= 0.00001 or f['ac'] == 1),
    'BRCA1':  dict(spec='GN092 v1.2', ba1=('faf', '>', 0.001), bs1=('faf', '>', 0.0001),
                   bs1_sup=('faf', '>', 0.00002), pm2=_ABSENT),
    'BRCA2':  dict(spec='GN097 v1.2', ba1=('faf', '>', 0.001), bs1=('faf', '>', 0.0001),
                   bs1_sup=('faf', '>', 0.00002), pm2=_ABSENT),
    # CDH1: CSpec says "MAF cutoff"; grpmax FAF95 is used (ClinGen SVI)
    'CDH1':   dict(spec='GN007 v3.1', ba1=('faf', '>=', 0.002), bs1=('faf', '>=', 0.001),
                   pm2=lambda f: f['af'] <= 0.00001 and (f['grpac'] < 2 or f['grpaf'] <= 0.00002)),
    'DICER1': dict(spec='GN024 v1.4', ba1=('grpaf', '>', 0.003, 2001, 5),
                   bs1=('grpaf', '>', 0.0003, 2001, 5),
                   pm2=lambda f: f['af'] < 0.000005 and f['grpac'] <= 1),
    'ENG':    dict(spec='GN136 v1.1', ba1=('faf', '>=', 0.01), bs1=('faf', '>', 0.002),
                   bs1_sup=('faf', '>', 0.0008), bs1_sup_hom=2,
                   pm2=lambda f: f['ac'] < 6 or f['grpaf'] < 0.00004),
    'KRAS':   dict(spec='GN044 v2.3', ba1=('faf', '>=', 0.0005), bs1=('faf', '>=', 0.00025), pm2=_ABSENT),
    'BRAF':   dict(spec='GN049 v2.3', ba1=('faf', '>=', 0.0005), bs1=('faf', '>=', 0.00025), pm2=_ABSENT),
    'MLH1':   dict(spec='GN115 v2.0', ba1=('faf', '>=', 0.001), bs1=('faf', '>=', 0.0001), pm2=_MMR_PM2),
    'MSH2':   dict(spec='GN137 v2.0', ba1=('faf', '>=', 0.001), bs1=('faf', '>=', 0.0001), pm2=_MMR_PM2),
    'MSH6':   dict(spec='GN138 v2.0', ba1=('faf', '>=', 0.0022), bs1=('faf', '>=', 0.00022), pm2=_MMR_PM2),
    'PMS2':   dict(spec='GN139 v2.0', ba1=('faf', '>=', 0.0028), bs1=('faf', '>=', 0.00028), pm2=_MMR_PM2),
    'PALB2':  dict(spec='GN077 v1.3', ba1=('faf', '>', 0.001), bs1=('faf', '>', 0.0001),
                   pm2=lambda f: f['grpaf'] <= 1.0 / 300000 or f['ac'] == 1),
    'PTEN':   dict(spec='GN003 v3.2', ba1=('faf', '>', 0.00056), bs1=('faf', '>=', 0.000043),
                   bs1_sup=('faf', '>=', 0.0000043),
                   pm2=lambda f: f['af'] < 0.00001 and (f['grpac'] < 2 or f['grpaf'] < 0.00002)),
    'TP53':   dict(spec='GN009 v2.4', ba1=('faf', '>=', 0.001, 2000, 2), bs1=('faf', '>=', 0.0003, 2000, 2),
                   pm2=lambda f: f['af'] < 0.00003 and (f['grpac'] < 2 or f['grpaf'] < 0.00004)),
    'VHL':    dict(spec='GN078 v1.1', ba1=('faf', '>=', 0.000156), bs1=('faf', '>=', 0.0000156),
                   pm2=lambda f: f['faf'] <= 0.00000156),
}

# PM1 (missense only, MANE Select protein numbering). 'mod' / 'sup': residue
# ranges (start, end); 'exons_mod': exon numbers; 'hotspots': cancerhotspots.org
# counts for the same amino-acid change, (min count Moderate, min count
# Supporting). The VHL germline hotspot table of GN078 is not in the CSpec API
# and is not used.
VCEP_PM1 = {
    'TP53':   dict(mod=[(175, 175), (245, 245), (248, 249), (273, 273), (282, 282)], hotspots=(10, 2)),
    'PTEN':   dict(mod=[(90, 94), (123, 130), (166, 168)]),
    'DICER1': dict(mod=[(1344, 1344), (1705, 1705), (1709, 1709), (1713, 1713),
                        (1809, 1810), (1813, 1813)], sup=[(1682, 1846)]),
    'ENG':    dict(mod=[(278, 278), (282, 282), (207, 207), (363, 363), (382, 382), (412, 412),
                        (549, 549), (350, 350), (394, 394), (516, 516), (582, 582)]),
    'VHL':    dict(hotspots=(10, 1)),
    'KRAS':   dict(mod=[(10, 17), (25, 40), (57, 64), (145, 156)]),
    'BRAF':   dict(mod=[(459, 474), (594, 627)], exons_mod=(6, 11)),
}

# Generated by scripts/build_vcep_hotspots.py from cancerhotspots.org
# hotspots_v2.xls (md5 c12cea615260ed24bb8f3bbe7121524a), missense only.
VCEP_HOTSPOTS = {   # cancerhotspots.org v2: "p.Arg175His" -> somatic count
    'TP53': {
        'p.Gln104His': 1, 'p.Gly105Arg': 2, 'p.Gly105Asp': 4, 'p.Gly105Cys': 13,
        'p.Gly105Ser': 2, 'p.Gly105Val': 6, 'p.Tyr107Asp': 6, 'p.Phe109Cys': 5,
        'p.Phe109Ile': 1, 'p.Phe109Ser': 6, 'p.Phe109Val': 3, 'p.Arg110Cys': 6,
        'p.Arg110His': 1, 'p.Arg110Leu': 25, 'p.Arg110Pro': 11, 'p.Leu111Arg': 6,
        'p.Leu111Gln': 4, 'p.Leu111Met': 1, 'p.Leu111Pro': 7, 'p.Phe113Cys': 6,
        'p.Phe113Leu': 1, 'p.Phe113Pro': 1, 'p.Phe113Ser': 2, 'p.Phe113Val': 6,
        'p.Lys120Asn': 2, 'p.Lys120Glu': 5, 'p.Lys120Met': 1, 'p.Thr125Ala': 2,
        'p.Thr125Arg': 6, 'p.Thr125Lys': 3, 'p.Thr125Met': 16, 'p.Thr125Pro': 3,
        'p.Tyr126Asn': 2, 'p.Tyr126Asp': 8, 'p.Tyr126Cys': 11, 'p.Tyr126His': 2,
        'p.Tyr126Ser': 3, 'p.Ser127Cys': 1, 'p.Ser127Phe': 31, 'p.Ser127Pro': 4,
        'p.Ser127Thr': 1, 'p.Ser127Tyr': 10, 'p.Leu130Arg': 4, 'p.Leu130His': 1,
        'p.Leu130Phe': 10, 'p.Leu130Pro': 2, 'p.Leu130Val': 13, 'p.Asn131Ile': 8,
        'p.Asn131Lys': 1, 'p.Asn131Ser': 3, 'p.Lys132Arg': 27, 'p.Lys132Asn': 30,
        'p.Lys132Gln': 2, 'p.Lys132Glu': 14, 'p.Lys132Met': 5, 'p.Lys132Thr': 4,
        'p.Met133Arg': 3, 'p.Met133Ile': 1, 'p.Met133Lys': 5, 'p.Met133Thr': 1,
        'p.Phe134Cys': 4, 'p.Phe134Ile': 2, 'p.Phe134Leu': 12, 'p.Phe134Ser': 1,
        'p.Phe134Val': 5, 'p.Cys135Arg': 9, 'p.Cys135Gly': 2, 'p.Cys135Phe': 24,
        'p.Cys135Ser': 4, 'p.Cys135Trp': 4, 'p.Cys135Tyr': 27, 'p.Gln136Glu': 10,
        'p.Gln136His': 1, 'p.Gln136Leu': 1, 'p.Gln136Pro': 3, 'p.Cys141Arg': 5,
        'p.Cys141Gly': 4, 'p.Cys141Phe': 1, 'p.Cys141Ser': 1, 'p.Cys141Trp': 9,
        'p.Cys141Tyr': 25, 'p.Val143Ala': 7, 'p.Val143Glu': 2, 'p.Val143Gly': 5,
        'p.Val143Met': 12, 'p.Gln144His': 1, 'p.Gln144Leu': 2, 'p.Gln144Pro': 3,
        'p.Trp146Cys': 1, 'p.Trp146Gly': 1, 'p.Trp146Ser': 1, 'p.Pro151Ala': 11,
        'p.Pro151Arg': 9, 'p.Pro151His': 13, 'p.Pro151Leu': 1, 'p.Pro151Phe': 1,
        'p.Pro151Ser': 36, 'p.Pro151Thr': 8, 'p.Pro152Arg': 2, 'p.Pro152Gln': 1,
        'p.Pro152Leu': 28, 'p.Pro152Lys': 1, 'p.Pro152Ser': 2, 'p.Pro152Thr': 2,
        'p.Gly154Asp': 1, 'p.Gly154Cys': 1, 'p.Gly154Val': 23, 'p.Thr155Asn': 8,
        'p.Thr155Ile': 2, 'p.Thr155Pro': 7, 'p.Val157Ala': 1, 'p.Val157Asp': 1,
        'p.Val157Gly': 4, 'p.Val157Leu': 3, 'p.Val157Phe': 67, 'p.Arg158Cys': 5,
        'p.Arg158Gly': 9, 'p.Arg158His': 36, 'p.Arg158Leu': 42, 'p.Arg158Phe': 1,
        'p.Arg158Pro': 11, 'p.Arg158Ser': 2, 'p.Ala159Asp': 1, 'p.Ala159Pro': 23,
        'p.Ala159Val': 27, 'p.Ala161Asp': 7, 'p.Ala161Pro': 1, 'p.Ala161Ser': 4,
        'p.Ala161Thr': 12, 'p.Ala161Val': 2, 'p.Tyr163Asn': 9, 'p.Tyr163Asp': 3,
        'p.Tyr163Cys': 54, 'p.Tyr163His': 11, 'p.Lys164Glu': 9, 'p.Gln165His': 1,
        'p.His168Arg': 12, 'p.His168Leu': 4, 'p.His168Pro': 2, 'p.His168Tyr': 1,
        'p.Val172Asp': 5, 'p.Val172Gly': 1, 'p.Val172Phe': 12, 'p.Val173Ala': 4,
        'p.Val173Glu': 1, 'p.Val173Gly': 4, 'p.Val173Leu': 42, 'p.Val173Met': 27,
        'p.Arg175Cys': 8, 'p.Arg175Gly': 14, 'p.Arg175His': 386, 'p.Arg175Leu': 8,
        'p.Cys176Arg': 5, 'p.Cys176Gly': 5, 'p.Cys176Phe': 54, 'p.Cys176Ser': 9,
        'p.Cys176Trp': 7, 'p.Cys176Tyr': 30, 'p.His178Asn': 1, 'p.His178Asp': 5,
        'p.His178Gln': 2, 'p.His178Pro': 3, 'p.His179Arg': 63, 'p.His179Asn': 8,
        'p.His179Asp': 8, 'p.His179Gln': 12, 'p.His179Leu': 14, 'p.His179Pro': 3,
        'p.His179Tyr': 41, 'p.Arg181Cys': 9, 'p.Arg181His': 9, 'p.Arg181Pro': 8,
        'p.Pro190Arg': 1, 'p.Pro190Leu': 12, 'p.Pro190Ser': 1, 'p.Pro190Thr': 3,
        'p.Gln192His': 1, 'p.Gln192Leu': 1, 'p.His193Arg': 54, 'p.His193Asn': 5,
        'p.His193Asp': 5, 'p.His193Gln': 1, 'p.His193Leu': 19, 'p.His193Pro': 8,
        'p.His193Tyr': 31, 'p.Leu194Arg': 30, 'p.Leu194His': 4, 'p.Leu194Phe': 12,
        'p.Leu194Pro': 3, 'p.Ile195Asn': 10, 'p.Ile195Lys': 1, 'p.Ile195Met': 2,
        'p.Ile195Phe': 19, 'p.Ile195Ser': 4, 'p.Ile195Thr': 42, 'p.Arg196Gln': 1,
        'p.Arg196Leu': 1, 'p.Arg196Pro': 8, 'p.Val197Glu': 3, 'p.Val197Gly': 1,
        'p.Val197Leu': 5, 'p.Val197Met': 3, 'p.Tyr205Asn': 2, 'p.Tyr205Asp': 3,
        'p.Tyr205Cys': 35, 'p.Tyr205His': 8, 'p.Tyr205Phe': 6, 'p.Tyr205Ser': 6,
        'p.Asp208Asn': 2, 'p.Asp208Glu': 1, 'p.Asp208Gly': 2, 'p.Asp208Val': 7,
        'p.Arg213Gln': 25, 'p.Arg213Gly': 4, 'p.Arg213Leu': 17, 'p.Arg213Pro': 2,
        'p.His214Arg': 35, 'p.His214Gln': 2, 'p.His214Leu': 2, 'p.His214Pro': 1,
        'p.Ser215Arg': 13, 'p.Ser215Asn': 7, 'p.Ser215Cys': 1, 'p.Ser215Gly': 8,
        'p.Ser215Ile': 20, 'p.Ser215Thr': 1, 'p.Val216Glu': 4, 'p.Val216Gly': 5,
        'p.Val216Leu': 9, 'p.Val216Met': 24, 'p.Val218Glu': 1, 'p.Val218Gly': 5,
        'p.Val218Leu': 1, 'p.Val218Met': 1, 'p.Tyr220Asn': 3, 'p.Tyr220Asp': 4,
        'p.Tyr220Cys': 127, 'p.Tyr220His': 6, 'p.Tyr220Ser': 8, 'p.Ile232Asn': 3,
        'p.Ile232Phe': 3, 'p.Ile232Ser': 4, 'p.Ile232Thr': 3, 'p.Tyr234Asn': 6,
        'p.Tyr234Asp': 2, 'p.Tyr234Cys': 42, 'p.Tyr234His': 9, 'p.Tyr234Ser': 3,
        'p.Tyr236Asn': 5, 'p.Tyr236Asp': 5, 'p.Tyr236Cys': 40, 'p.Tyr236His': 2,
        'p.Tyr236Ser': 2, 'p.Met237Ile': 55, 'p.Met237Lys': 5, 'p.Met237Val': 4,
        'p.Cys238Arg': 5, 'p.Cys238Gly': 4, 'p.Cys238Phe': 23, 'p.Cys238Ser': 4,
        'p.Cys238Trp': 1, 'p.Cys238Tyr': 40, 'p.Asn239Asp': 18, 'p.Asn239His': 1,
        'p.Asn239Lys': 2, 'p.Asn239Ser': 19, 'p.Asn239Thr': 1, 'p.Ser241Ala': 2,
        'p.Ser241Cys': 19, 'p.Ser241Phe': 51, 'p.Ser241Pro': 6, 'p.Ser241Thr': 1,
        'p.Ser241Tyr': 15, 'p.Cys242Arg': 2, 'p.Cys242Gly': 3, 'p.Cys242Phe': 27,
        'p.Cys242Ser': 9, 'p.Cys242Trp': 4, 'p.Cys242Tyr': 12, 'p.Gly244Ala': 1,
        'p.Gly244Arg': 1, 'p.Gly244Asp': 21, 'p.Gly244Cys': 19, 'p.Gly244Phe': 1,
        'p.Gly244Ser': 15, 'p.Gly244Val': 6, 'p.Gly245Ala': 3, 'p.Gly245Arg': 7,
        'p.Gly245Asn': 1, 'p.Gly245Asp': 50, 'p.Gly245Cys': 22, 'p.Gly245Ser': 128,
        'p.Gly245Val': 31, 'p.Met246Arg': 2, 'p.Met246Ile': 10, 'p.Met246Leu': 3,
        'p.Met246Lys': 3, 'p.Met246Thr': 3, 'p.Met246Val': 16, 'p.Asn247Ile': 8,
        'p.Arg248Gln': 292, 'p.Arg248Gly': 8, 'p.Arg248Leu': 21, 'p.Arg248Pro': 11,
        'p.Arg248Trp': 222, 'p.Arg249Gly': 15, 'p.Arg249Ile': 1, 'p.Arg249Lys': 2,
        'p.Arg249Met': 17, 'p.Arg249Ser': 45, 'p.Arg249Thr': 11, 'p.Arg249Trp': 11,
        'p.Pro250Arg': 1, 'p.Pro250Leu': 22, 'p.Ile251Asn': 5, 'p.Ile251Met': 1,
        'p.Ile251Phe': 5, 'p.Ile251Ser': 3, 'p.Ile254Asn': 1, 'p.Ile254Leu': 1,
        'p.Ile254Met': 1, 'p.Ile254Phe': 1, 'p.Ile254Ser': 7, 'p.Ile254Val': 2,
        'p.Ile255Asn': 7, 'p.Ile255Phe': 12, 'p.Ile255Ser': 6, 'p.Ile255Thr': 6,
        'p.Leu257Arg': 4, 'p.Leu257Gln': 7, 'p.Leu257Pro': 4, 'p.Leu257Val': 1,
        'p.Glu258Ala': 4, 'p.Glu258Asp': 5, 'p.Glu258Gln': 4, 'p.Glu258Gly': 4,
        'p.Glu258Lys': 14, 'p.Glu258Val': 1, 'p.Asp259Asn': 2, 'p.Asp259His': 2,
        'p.Asp259Tyr': 19, 'p.Asp259Val': 7, 'p.Gly262Val': 17, 'p.Leu265Arg': 4,
        'p.Leu265Met': 1, 'p.Leu265Pro': 8, 'p.Gly266Arg': 40, 'p.Gly266Glu': 35,
        'p.Gly266Lys': 1, 'p.Gly266Thr': 1, 'p.Gly266Val': 26, 'p.Arg267Gln': 4,
        'p.Arg267Gly': 3, 'p.Arg267Leu': 5, 'p.Arg267Pro': 8, 'p.Arg267Trp': 15,
        'p.Phe270Cys': 4, 'p.Phe270Ile': 4, 'p.Phe270Leu': 8, 'p.Phe270Ser': 6,
        'p.Phe270Val': 8, 'p.Glu271Gln': 7, 'p.Glu271Lys': 27, 'p.Glu271Val': 5,
        'p.Val272Ala': 1, 'p.Val272Glu': 3, 'p.Val272Gly': 5, 'p.Val272Leu': 20,
        'p.Val272Lys': 2, 'p.Val272Met': 34, 'p.Arg273Cys': 283, 'p.Arg273Gly': 2,
        'p.Arg273His': 251, 'p.Arg273Leu': 56, 'p.Arg273Pro': 5, 'p.Arg273Ser': 12,
        'p.Val274Ala': 12, 'p.Val274Asp': 7, 'p.Val274Gly': 8, 'p.Val274Leu': 5,
        'p.Val274Phe': 14, 'p.Cys275Arg': 7, 'p.Cys275Gly': 2, 'p.Cys275Phe': 25,
        'p.Cys275Ser': 6, 'p.Cys275Trp': 3, 'p.Cys275Tyr': 24, 'p.Cys277Gly': 1,
        'p.Cys277Phe': 20, 'p.Cys277Trp': 1, 'p.Cys277Tyr': 3, 'p.Pro278Ala': 13,
        'p.Pro278Arg': 12, 'p.Pro278His': 6, 'p.Pro278Leu': 17, 'p.Pro278Phe': 3,
        'p.Pro278Ser': 37, 'p.Pro278Thr': 12, 'p.Arg280Gly': 17, 'p.Arg280Ile': 12,
        'p.Arg280Lys': 36, 'p.Arg280Ser': 8, 'p.Arg280Thr': 34, 'p.Asp281Ala': 1,
        'p.Asp281Asn': 14, 'p.Asp281Glu': 14, 'p.Asp281Gly': 4, 'p.Asp281His': 9,
        'p.Asp281Pro': 1, 'p.Asp281Tyr': 13, 'p.Asp281Val': 9, 'p.Arg282Gln': 3,
        'p.Arg282Gly': 9, 'p.Arg282Leu': 1, 'p.Arg282Pro': 4, 'p.Arg282Trp': 201,
        'p.Glu285Gln': 3, 'p.Glu285Lys': 68, 'p.Glu285Val': 9, 'p.Glu286Ala': 2,
        'p.Glu286Gln': 7, 'p.Glu286Gly': 6, 'p.Glu286Lys': 41, 'p.Glu286Val': 3,
        'p.Glu298Gln': 1, 'p.Arg306Gln': 1, 'p.Gln331His': 8, 'p.Arg337Cys': 36,
        'p.Arg337Gly': 1, 'p.Arg337His': 4, 'p.Arg337Leu': 20, 'p.Arg337Pro': 3,
        'p.Arg337Ser': 4, 'p.Arg342Pro': 8,
    },
    'VHL': {
        'p.Ser65Leu': 6, 'p.Ser65Pro': 2, 'p.Ser65Thr': 1, 'p.Ser65Trp': 1, 'p.Ser68Pro': 1,
        'p.Val74Asp': 5, 'p.Val74Gly': 1, 'p.Asn78Asp': 2, 'p.Asn78Lys': 1, 'p.Asn78Ser': 4,
        'p.Asn78Tyr': 2, 'p.Ser80Arg': 5, 'p.Ser80Asn': 4, 'p.Ser80Ile': 1, 'p.Pro86Arg': 1,
        'p.Pro86Leu': 2, 'p.Pro86Ser': 4, 'p.Pro86Thr': 2, 'p.Trp88Arg': 2, 'p.Trp88Leu': 2,
        'p.Leu89His': 9, 'p.Leu89Pro': 5, 'p.Ser111Arg': 2, 'p.Ser111Asn': 6,
        'p.Ser111Pro': 1, 'p.Tyr112Asp': 1, 'p.Tyr112Cys': 1, 'p.Gly114Ala': 4,
        'p.Gly114Arg': 2, 'p.Gly114Cys': 1, 'p.His115Asn': 4, 'p.His115Gln': 1,
        'p.His115Tyr': 2, 'p.Trp117Cys': 1, 'p.Trp117Leu': 3, 'p.Asp121Gly': 2,
        'p.Asp121Tyr': 6, 'p.Asp121Val': 1, 'p.Leu135Phe': 2, 'p.Ile151Asn': 2,
        'p.Ile151Phe': 1, 'p.Ile151Ser': 1, 'p.Ile151Thr': 2, 'p.Leu158Gln': 1,
        'p.Leu158Pro': 3, 'p.Leu158Val': 4, 'p.Arg161Gln': 2, 'p.Arg161Leu': 1,
        'p.Cys162Arg': 2, 'p.Cys162Phe': 1, 'p.Cys162Trp': 1, 'p.Cys162Tyr': 3,
        'p.Leu169Pro': 5,
    },
}


def vcep_enabled():
    return _E('ACMG_VCEP_SPECS', 'true').strip().lower() not in ('false', '0', 'no')


def vcep_spec(get_field):
    """VCEP_SPECS entry of the variant gene, or None (general rules)."""
    if not vcep_enabled():
        return None
    return VCEP_SPECS.get(get_field('ANN[*].GENE').strip().upper())


def _num0(v):
    vals = []
    for x in str(v or '').split(','):
        try:
            vals.append(float(x))
        except ValueError:
            pass
    return max(vals) if vals else None


def gnomad_vcep_freqs(get_field):
    """gnomAD values for the VCEP rules, or None when the grpmax columns are
    missing / "." (legacy mode, older resource, site outside the resource
    regions): the general rules apply then."""
    grpaf = _num0(get_field('gnomAD_AF_GRPMAX'))
    af = gnomad_af(get_field)
    if grpaf is None or af is None:
        return None
    return {
        'af': af,
        'ac': _num0(get_field('gnomAD_AC')) or 0,
        'faf': _num0(get_field('gnomAD_FAF95')) or 0.0,
        'grpaf': grpaf,
        'grpac': _num0(get_field('gnomAD_AC_GRPMAX')) or 0,
        'grpan': _num0(get_field('gnomAD_AN_GRPMAX')) or 0,
        'nhom': _num0(get_field('gnomAD_NHOMALT')) or 0,
    }


def _freq_rule(rule, f):
    """True when frequency rule (metric, op, threshold[, min_grpan, min_grpac]) is met."""
    metric, op, thr = rule[:3]
    if len(rule) > 3 and (f['grpan'] < rule[3] or f['grpac'] < rule[4]):
        return False
    v = f[metric]
    return v > thr if op == '>' else v >= thr


def vcep_frequency(get_field):
    """(spec, freqs) when the gene has a VCEP specification and the gnomAD
    grpmax data are available, else (None, None)."""
    spec = vcep_spec(get_field)
    if spec is None:
        return None, None
    f = gnomad_vcep_freqs(get_field)
    return (spec, f) if f is not None else (None, None)


def ba1_reason(get_field):
    """Text for the Step 2 BA1 filter ('' = BA1 not met), same rule as the
    ACMG engine: VCEP gene threshold when available, else ba1_af."""
    spec, f = vcep_frequency(get_field)
    if spec is not None:
        if _freq_rule(spec['ba1'], f):
            m, op, thr = spec['ba1'][:3]
            return 'gnomad_BA1 ({}={:.3g}; AF={:.3g}; ClinGen {} {} {:g})'.format(
                'FAF95' if m == 'faf' else 'AF_grpmax', f[m], f['af'], spec['spec'], op, thr)
        return ''
    benign = benign_freq(get_field)
    if benign is not None and benign > float(_E('ACMG_BA1_AF', '0.05')):
        faf = gnomad_faf95(get_field)
        if faf is not None:
            return 'gnomad_BA1 (FAF95={:.4f}; AF={:.4f})'.format(faf, gnomad_af(get_field) or 0)
        return 'gnomad_BA1 (AF={:.4f})'.format(benign)
    return ''


_AA3 = r'(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)'


def pm1_evidence(get_field):
    """(points, code) for PM1 from the VCEP_PM1 gene rules (missense only):
    Moderate +2 (PM1), Supporting +1 (PM1_Sup)."""
    import re
    gene = get_field('ANN[*].GENE').strip().upper()
    rule = VCEP_PM1.get(gene) if vcep_enabled() else None
    if rule is None or 'missense' not in get_field('ANN[*].EFFECT').lower():
        return 0, None
    m = re.match(r'^p\.(' + _AA3 + r')(\d+)(' + _AA3 + r')$', get_field('ANN[*].HGVS_P').strip())
    if not m:
        return 0, None
    pos = int(m.group(2))
    mod = int(_E('ACMG_PM1_SCORE', '2'))
    sup = int(_E('ACMG_PM1_SUP_SCORE', '1'))
    inside = lambda ranges: any(a <= pos <= b for a, b in ranges)
    if inside(rule.get('mod', ())):
        return mod, 'PM1'
    exon = _int_or_none(get_field('ANN[*].RANK'))
    if exon is not None and exon in rule.get('exons_mod', ()):
        return mod, 'PM1'
    if 'hotspots' in rule:
        n = VCEP_HOTSPOTS.get(gene, {}).get(m.group(0), 0)
        if n >= rule['hotspots'][0]:
            return mod, 'PM1(hotspot n={})'.format(n)
        if n >= rule['hotspots'][1]:
            return sup, 'PM1_Sup(hotspot n={})'.format(n)
    if inside(rule.get('sup', ())):
        return sup, 'PM1_Sup'
    return 0, None


def acmg_spec(get_field):
    """Rule set used for the variant (column ACMG_spec)."""
    spec = vcep_spec(get_field)
    if spec is None:
        return 'ACMG/SVI general'
    if gnomad_vcep_freqs(get_field) is None:
        return 'ClinGen {} (frequency: general rules, no gnomAD grpmax)'.format(spec['spec'])
    return 'ClinGen ' + spec['spec']


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
      2. gnomAD BA1 rule (FAF95 > 5%, or the ClinGen VCEP gene threshold)
      3. Evidence points (Tavtigian 2020): PVS1 by gene mechanism and variant
         type (pvs1_evidence) + PM2_Supporting (rare) + PP3/BP4 from REVEL
         (Pejaver 2022) and SpliceAI (Walker 2023) (pp3_bp4_evidence)
         + PM1 for VCEP genes (pm1_evidence; PP3 + PM1 <= 4 points)
         + PS1 / PM5 from ClinVar at the same residue (ps1_pm5_evidence)
         + BS1 (FAF95 > 1%) / BS1_Sup. BA1, BS1, BS1_Sup and PM2 follow the
         ClinGen VCEP gene specification when there is one (VCEP_SPECS).
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
    spec, freqs = vcep_frequency(get_field)   # ClinGen VCEP gene rules, or general
    if ba1_reason(get_field):
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
    if spec is not None:
        pm2 = gnomad_supports_pm2(get_field) and spec['pm2'](freqs)
    else:
        pm2 = (gnomad is not None and gnomad_supports_pm2(get_field)
               and gnomad < pm2_af and (faf is None or faf < pm2_af))
    if pm2:
        score += int(_E('ACMG_PM2_SCORE', '1'))
        criteria.append('PM2_Sup')

    # PP3 / BP4: REVEL (Pejaver 2022) and SpliceAI (Walker 2023); not for null
    # variants, whose splicing / LOF effect is already weighed by PVS1
    if code is None:
        pp3_pts, code = pp3_bp4_evidence(get_field)
        score += pp3_pts
        if code:
            criteria.append(code)
        # PM1 (VCEP genes); PP3 + PM1 together at most Strong (Pejaver 2022)
        pts, code = pm1_evidence(get_field)
        if code:
            if pp3_pts > 0 and pp3_pts + pts > 4:
                pts = max(0, 4 - pp3_pts)
                code += '(cap PP3+PM1={})'.format(pp3_pts + pts)
            score += pts
            criteria.append(code)
        pts, code = ps1_pm5_evidence(get_field)
        score += pts
        if code:
            criteria.append(code)

    # BS1 — frequency above benign threshold (VCEP gene rule, else bs1_af)
    if spec is not None:
        bs1 = _freq_rule(spec['bs1'], freqs)
        bs1_sup = not bs1 and 'bs1_sup' in spec and _freq_rule(spec['bs1_sup'], freqs)
        if bs1_sup and freqs['nhom'] >= spec.get('bs1_sup_hom', float('inf')):
            bs1, bs1_sup = True, False      # ENG: BS1_Sup + >= 2 homozygotes -> BS1
    else:
        bs1 = benign is not None and benign > float(_E('ACMG_BS1_AF', '0.01'))
        bs1_sup = False
    if bs1:
        score += int(_E('ACMG_BS1_SCORE', '-4'))
        criteria.append('BS1')
    elif bs1_sup:
        score += int(_E('ACMG_BS1_SUP_SCORE', '-1'))
        criteria.append('BS1_Sup')

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
    """REVEL score: canonical `REVEL_score` (local REVEL v1.3, MANE Select transcript)."""
    for name in ('REVEL_score',):
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
    (0, None) when REVEL is missing or the MANE effect is not missense (a
    synonymous / intronic MANE variant can carry the score of another
    transcript, REVEL_v1.3_other_tx), (0, '') when it is indeterminate."""
    revel = revel_score(get_field)
    if revel is None or 'missense' not in (get_field('ANN[*].EFFECT') or '').lower():
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
# Generated from assets/gene_panels/genes.tsv (scripts/build_gene_lists.py)
from gene_sets import MMR_GENES, HRR_GENES  # noqa: E402,F401


# =============================================================================
# Mutational signature sets (env-overridable)
# =============================================================================
MMR_SBS_SIGNATURES = frozenset(_E('SIG_MMR_SBS', 'SBS6,SBS15,SBS20,SBS21,SBS26,SBS44').split(','))
HRD_SBS_SIGNATURES = frozenset(_E('SIG_HRD_SBS', 'SBS3,SBS8').split(','))
MMR_ID_SIGNATURES  = frozenset(_E('SIG_MMR_ID',  'ID1,ID2').split(','))

# A signature is "detected" only with enough support (filters.signatures):
#   share of the sample's assigned mutations >= SIG_MIN_FRACTION,
#   >= SIG_MIN_COUNT mutations assigned to it, and a sample total of
#   >= SIG_MIN_TOTAL_SBS (SBS) / SIG_MIN_TOTAL_ID (ID) mutations; below that
#   total the sample is "not assessable" for that signature type.
SIG_MIN_FRACTION  = float(_E('SIG_MIN_FRACTION', '0.10'))
SIG_MIN_COUNT     = float(_E('SIG_MIN_COUNT', '10'))
SIG_MIN_TOTAL_SBS = float(_E('SIG_MIN_TOTAL_SBS', '100'))
SIG_MIN_TOTAL_ID  = float(_E('SIG_MIN_TOTAL_ID', '50'))


def signature_calls(activities_path, sample, kind='SBS'):
    """
    SigProfiler Assignment_Solution_Activities.txt -> dict:
      activities  {signature: mutations} (> 0) of the sample row
      total       mutations assigned in the sample
      assessable  total >= SIG_MIN_TOTAL_SBS / SIG_MIN_TOTAL_ID (kind 'SBS' / 'ID')
      detected    signatures passing SIG_MIN_FRACTION and SIG_MIN_COUNT (empty
                  when not assessable)
    `sample` matches the row exactly or as the tumor of '<tumor>_vs_<normal>'.
    """
    import csv
    out = {'activities': {}, 'total': 0.0, 'assessable': False, 'detected': set()}
    if not activities_path or not _os.path.isfile(activities_path):
        return out
    with open(activities_path, newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    if len(rows) < 2:
        return out
    header = rows[0]
    for row in rows[1:]:
        name = row[0] if row else ''
        if name != sample and not name.startswith(sample + '_vs_'):
            continue
        for i, sig in enumerate(header[1:], 1):
            try:
                v = float(row[i])
            except (IndexError, ValueError, TypeError):
                continue
            if v > 0:
                out['activities'][sig] = v
        break
    total = sum(out['activities'].values())
    out['total'] = total
    out['assessable'] = total >= (SIG_MIN_TOTAL_ID if kind == 'ID' else SIG_MIN_TOTAL_SBS)
    if out['assessable']:
        out['detected'] = {s for s, v in out['activities'].items()
                           if v >= SIG_MIN_COUNT and v / total >= SIG_MIN_FRACTION}
    return out


def variant_label(get_field):
    """HGVS.p, else HGVS.c, ignoring '.' placeholders ('' when neither)."""
    for k in ('ANN[*].HGVS_P', 'ANN[*].HGVS_C'):
        v = (get_field(k) or '').strip()
        if v and v != '.':
            return v
    return ''

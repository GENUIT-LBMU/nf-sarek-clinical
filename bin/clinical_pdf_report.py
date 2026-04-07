#!/usr/bin/env python3
"""
clinical_pdf_report.py
======================
Generates a clinical-grade PDF report for the GENUIT nf-sarek-clinical pipeline.

Sections:
  Page 1  – Cover / KPI Dashboard
  Page 2+ – Germline SNVs (ACMG P / LP / VUS)
  Page N+ – Somatic SNVs (AMP Tiers IA-IIC and III)

Dependencies: reportlab (Platypus layout engine)
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import date

# ---------------------------------------------------------------------------
# ReportLab imports
# ---------------------------------------------------------------------------
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    HRFlowable,
    KeepTogether,
    PageBreak,
)
from reportlab.platypus.flowables import Flowable

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
C_DARK_HEADER  = colors.HexColor('#1e293b')
C_ACCENT       = colors.HexColor('#f97316')
C_SLATE        = colors.HexColor('#334155')
C_ROW_ALT      = colors.HexColor('#f8fafc')
C_WHITE        = colors.white
C_TEXT         = colors.HexColor('#1e293b')
C_LIGHT_BORDER = colors.HexColor('#e2e8f0')

# ACMG badge colours
C_ACMG_P   = colors.HexColor('#ef4444')   # Pathogenic – red
C_ACMG_LP  = colors.HexColor('#f97316')   # Likely Pathogenic – orange
C_ACMG_VUS = colors.HexColor('#eab308')   # VUS – yellow

# AMP tier colours
C_AMP_IA_IB  = colors.HexColor('#dc2626')  # Tier IA / IB – dark red
C_AMP_IIA_IIB = colors.HexColor('#ea580c') # Tier IIA / IIB – orange
C_AMP_IIC    = colors.HexColor('#d97706')  # Tier IIC – amber
C_AMP_III    = colors.HexColor('#6b7280')  # Tier III – gray

# HRD status colours
C_HRD_POS = colors.HexColor('#ef4444')
C_HRD_NEG = colors.HexColor('#22c55e')

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
PAGE_W, PAGE_H = A4          # 595.28 x 841.89 pt
MARGIN_L = 1.8 * cm
MARGIN_R = 1.8 * cm
MARGIN_T = 2.0 * cm
MARGIN_B = 2.2 * cm
USABLE_W = PAGE_W - MARGIN_L - MARGIN_R   # ~553 pt

# ---------------------------------------------------------------------------
# Internationalisation (i18n) — English / Spanish
# ---------------------------------------------------------------------------
_LANG = 'en'

_TR = {
    'report_title':      {'en': 'GENUIT Molecular Report',          'es': 'Informe Molecular GENUIT'},
    'meta_run':          {'en': 'Run',      'es': 'Corrida'},
    'meta_patient':      {'en': 'Patient',  'es': 'Paciente'},
    'meta_tumor':        {'en': 'Tumor',    'es': 'Tumor'},
    'meta_normal':       {'en': 'Normal',   'es': 'Normal'},
    'meta_date':         {'en': 'Date',     'es': 'Fecha'},
    'meta_pipeline':     {'en': 'Pipeline', 'es': 'Pipeline'},
    'footer_conf':       {'en': 'CONFIDENTIAL', 'es': 'CONFIDENCIAL'},
    # Section headers
    'sect_global_bm':    {'en': 'Global Biomarkers',        'es': 'Biomarcadores Globales'},
    'sect_major_finds':  {'en': 'Major Clinical Findings',  'es': 'Hallazgos Cl\u00ednicos Principales'},
    'sect_secondary_bm': {'en': 'Secondary Biomarkers',     'es': 'Biomarcadores Secundarios'},
    'sect_quality':      {'en': 'Genomic Quality Metrics',  'es': 'M\u00e9tricas de Calidad Gen\u00f3mica'},
    'sect_germline':     {'en': 'Germline Variants (ACMG Classification)',
                          'es': 'Variantes Germinales (Clasificaci\u00f3n ACMG)'},
    'sect_somatic':      {'en': 'Somatic Variants (AMP/ASCO/CAP Tier Classification)',
                          'es': 'Variantes Som\u00e1ticas (Clasificaci\u00f3n AMP/ASCO/CAP)'},
    'sect_germ_cnv':     {'en': 'Germline CNVs (GATK gCNV) \u2014 CN \u2260 2',
                          'es': 'CNVs Germinales (GATK gCNV) \u2014 CN \u2260 2'},
    'sect_som_cnv':      {'en': 'Somatic CNVs \u2014 Reportable Events',
                          'es': 'CNVs Som\u00e1ticos \u2014 Eventos Reportables'},
    'sect_seqq':         {'en': 'Sequencing Quality \u2014 Per-Gene Coverage',
                          'es': 'Calidad de Secuenciaci\u00f3n \u2014 Cobertura por Gen'},
    'sect_correlation':  {'en': 'Germline\u2013Somatic Correlation Analysis',
                          'es': 'An\u00e1lisis de Correlaci\u00f3n Germinal\u2013Som\u00e1tica'},
    'sect_mutsig':       {'en': 'Mutational Signatures',
                          'es': 'Firmas Mutacionales'},
    # KPI labels
    'kpi_tmb':      {'en': 'TMB (mut/Mb)',        'es': 'TMB (mut/Mb)'},
    'kpi_msi':      {'en': 'MSI',                 'es': 'MSI'},
    'kpi_hrd_prob': {'en': 'HRD Probability',     'es': 'Probabilidad HRD'},
    'kpi_hrd_stat': {'en': 'HRD Status',          'es': 'Estado HRD'},
    'kpi_wgd':      {'en': 'WGD',                 'es': 'WGD'},
    'kpi_gi':       {'en': 'Genomic Instability', 'es': 'Inestabilidad Gen\u00f3mica'},
    'kpi_loh':      {'en': 'LOH',                 'es': 'LOH'},
    'kpi_ploidy':   {'en': 'Ploidy',              'es': 'Ploid\u00eda'},
    'kpi_purity':   {'en': 'Tumor Purity (%)',    'es': 'Pureza Tumoral (%)'},
    'kpi_cov_n':    {'en': 'Coverage Normal (X)', 'es': 'Cobertura Normal (X)'},
    'kpi_cov_t':    {'en': 'Coverage Tumor (X)',  'es': 'Cobertura Tumor (X)'},
    # MSI status strings
    'msi_detected': {'en': 'MSI detected',           'es': 'MSI detectada'},
    'msi_not_det':  {'en': 'MSI not detected (MSS)', 'es': 'MSI no detectada (MSS)'},
    # Major findings — card titles
    'find_hrd':     {'en': 'HRD Pathway',               'es': 'V\u00eda HRD'},
    'find_mmr':     {'en': 'MMR Deficiency',            'es': 'Deficiencia MMR'},
    'find_dsh':     {'en': 'Double Somatic Hits',       'es': 'Doble Impacto Som\u00e1tico'},
    'find_cnv':     {'en': 'High-Confidence CNV Events','es': 'Eventos CNV de Alta Confianza'},
    'find_civic':   {'en': 'CIViC Actionable Variants', 'es': 'Variantes Accionables CIViC'},
    'civic_resist': {'en': 'Resistance',                'es': 'Resistencia'},
    'civic_sensit': {'en': 'Sensitivity/Response',      'es': 'Sensibilidad/Respuesta'},
    # Low-coverage warning
    'low_cov_title': {'en': 'Low Tumor Sequencing Depth',
                      'es': 'Baja Profundidad de Secuenciación Tumoral'},
    'low_cov_body':  {
        'en': ('Tumor mean coverage is <b>{cov}X</b> (recommended ≥100X for WES). '
               'The following results may have <b>reduced sensitivity</b> and should be '
               'interpreted with caution:'),
        'es': ('La cobertura media del tumor es <b>{cov}X</b> (recomendado ≥100X para WES). '
               'Los siguientes resultados pueden tener <b>sensibilidad reducida</b> y deben '
               'interpretarse con precaución:'),
    },
    'low_cov_tmb':    {'en': 'TMB — underestimated due to fewer detected variants',
                       'es': 'TMB — subestimado por menor detección de variantes'},
    'low_cov_msi':    {'en': 'MSI — microsatellite instability detection may be unreliable',
                       'es': 'MSI — la detección de inestabilidad microsatelital puede ser poco fiable'},
    'low_cov_mutsig': {'en': 'Mutational Signatures — reduced mutation count limits signature decomposition accuracy',
                       'es': 'Firmas Mutacionales — el bajo número de mutaciones limita la precisión de la descomposición'},
    'low_cov_somatic':{'en': 'Somatic Variants — low-VAF variants may be missed',
                       'es': 'Variantes Somáticas — pueden no detectarse variantes de baja VAF'},
    'low_cov_cnv':    {'en': 'Somatic CNVs — copy-number segmentation has reduced resolution',
                       'es': 'CNVs Somáticos — la segmentación de número de copias tiene resolución reducida'},
    'low_cov_hrd':    {'en': 'HRD — probability derived from mutational signatures is unreliable with low mutation counts',
                       'es': 'HRD — la probabilidad derivada de firmas mutacionales no es fiable con bajo conteo de mutaciones'},
    'indeterminate':  {'en': 'Indeterminate',    'es': 'Indeterminado'},
    'low_cov_indet':  {
        'en': ('Result indeterminate due to low tumor sequencing depth. '
               'Insufficient coverage compromises the quality of this analysis.'),
        'es': ('Resultado indeterminado debido a la baja profundidad de secuenciación tumoral. '
               'La cobertura insuficiente compromete la calidad de este análisis.'),
    },
    # Major findings — content labels
    'hrd_status_lbl': {'en': 'HRD Status:',         'es': 'Estado HRD:'},
    'prob_lbl':       {'en': 'probability',          'es': 'probabilidad'},
    'germ_hit_lbl':   {'en': 'Germline hit:',        'es': 'Variante germinal:'},
    '2nd_hits_lbl':   {'en': '2nd hit(s):',          'es': '2do impacto:'},
    '2nd_hit_lbl':    {'en': '2nd hit:',             'es': '2do impacto:'},
    'hrd_sigs_lbl':   {'en': 'HRD signatures:',      'es': 'Firmas HRD:'},
    'sbs3_not_det':   {'en': 'SBS3/SBS8 not detected','es': 'SBS3/SBS8 no detectadas'},
    'loh_det':        {'en': 'LOH detected (somatic 2nd hit)',  'es': 'LOH detectado (2do impacto)'},
    'homdel_hit':     {'en': 'HOMDEL (CN=0, somatic 2nd hit)', 'es': 'HOMDEL (CN=0, 2do impacto)'},
    'somatic_pfx':    {'en': 'Somatic:',             'es': 'Som\u00e1tico:'},
    'not_detected':   {'en': 'not detected',         'es': 'no detectado'},
    'sbs_mmr_lbl':    {'en': 'SBS MMR signatures:',  'es': 'Firmas SBS MMR:'},
    'id_mmr_lbl':     {'en': 'ID MMR signatures:',   'es': 'Firmas ID MMR:'},
    'not_det_f':      {'en': 'not detected',         'es': 'no detectadas'},
    'biallelic_lbl':  {'en': 'Biallelic somatic:',   'es': 'Bial\u00e9lico som\u00e1tico:'},
    'somatic_var_lbl':{'en': 'Somatic variant:',     'es': 'Variante som\u00e1tica:'},
    # "No findings" message
    'no_findings': {
        'en': ('No major clinical findings meeting reporting criteria were identified. '
               'Detailed variant classifications are provided in the sections below.'),
        'es': ('No se identificaron hallazgos cl\u00ednicos principales que cumplan los criterios '
               'de reporte. Las clasificaciones detalladas se proporcionan en las secciones siguientes.'),
    },
    # Disclaimer
    'disclaimer': {
        'en': ('<b>DISCLAIMER:</b> This report is generated for research and clinical '
               'decision-support purposes only. Findings should be interpreted by a qualified '
               'clinical geneticist or molecular pathologist in the context of the patient\'s '
               'clinical history. Variant classifications follow ACMG/AMP 2015 guidelines '
               'and are subject to revision as knowledge evolves.'),
        'es': ('<b>AVISO LEGAL:</b> Este informe ha sido generado con fines de apoyo a la '
               'investigaci\u00f3n y a la toma de decisiones cl\u00ednicas. Los hallazgos deben '
               'ser interpretados por un genetista cl\u00ednico o pat\u00f3logo molecular '
               'calificado en el contexto de la historia cl\u00ednica del paciente. Las '
               'clasificaciones de variantes siguen las gu\u00edas ACMG/AMP 2015 y est\u00e1n '
               'sujetas a revisi\u00f3n conforme avanza el conocimiento.'),
    },
    # Germline no-variants message
    'germ_no_vars': {
        'en': ('No pathogenic, likely pathogenic, or uncertain significance germline '
               'variants were identified in this sample.'),
        'es': ('No se identificaron variantes germinales patog\u00e9nicas, probablemente '
               'patog\u00e9nicas o de significado incierto en esta muestra.'),
    },
    # ACMG legend labels
    'acmg_p':   {'en': 'Pathogenic',                       'es': 'Patog\u00e9nica'},
    'acmg_lp':  {'en': 'Likely Pathogenic',                'es': 'Prob. Patog\u00e9nica'},
    'acmg_vus': {'en': 'Variant of Uncertain Significance', 'es': 'Variante de Signif. Incierto'},
}


def T(key):
    """Return translated string for the current language, falling back to English."""
    entry = _TR.get(key, {})
    return entry.get(_LANG) or entry.get('en') or key


def set_lang(lang):
    """Set the report language globally ('en' or 'es')."""
    global _LANG
    _LANG = (lang or 'en').lower()


def _fmt_hrd_prob(val):
    """Format HRD probability: '<0.01' if very small, else 2 decimal places."""
    try:
        f = float(str(val).strip())
        if f != f:
            return 'NA'
        if f < 0.01:
            return '<0.01'
        return f'{f:.2f}'
    except (ValueError, TypeError):
        return str(val) if val not in ('NA', None, '') else 'NA'


def _make_kpi_cell(label, value):
    """Create a KPI Paragraph: small grey label above large coloured value.
    HRD+ gets red, HRD- gets green; everything else gets orange.
    """
    hrd_val = str(value).upper() if value else ''
    if hrd_val.startswith('HRD+'):
        val_colour = C_HRD_POS.hexval()
    elif hrd_val.startswith('HRD-'):
        val_colour = C_HRD_NEG.hexval()
    else:
        val_colour = '#f97316'
    return Paragraph(
        f'<font size="7" color="#64748b"><b>{label}</b></font><br/>'
        f'<font size="13" color="{val_colour}"><b>{value}</b></font>',
        styles['Normal'],
    )


def _build_kpi_table(kpi_list, ncols=3):
    """Build a uniform KPI grid from [(label, value), ...] with ncols columns."""
    cells = [_make_kpi_cell(lbl, val) for lbl, val in kpi_list]
    while len(cells) % ncols != 0:
        cells.append(Paragraph('', styles['Normal']))
    rows = [cells[i:i + ncols] for i in range(0, len(cells), ncols)]
    col_w = USABLE_W / ncols
    tbl = Table(rows, colWidths=[col_w] * ncols)
    tbl.setStyle(TableStyle([
        ('BACKGROUND',     (0, 0), (-1, -1), colors.HexColor('#f8fafc')),
        ('ROWBACKGROUNDS', (0, 0), (-1, -1),
         [colors.HexColor('#ffffff'), colors.HexColor('#f1f5f9')]),
        ('BOX',            (0, 0), (-1, -1), 0.5, C_LIGHT_BORDER),
        ('GRID',           (0, 0), (-1, -1), 0.3, C_LIGHT_BORDER),
        ('TOPPADDING',     (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING',  (0, 0), (-1, -1), 10),
        ('LEFTPADDING',    (0, 0), (-1, -1), 12),
        ('RIGHTPADDING',   (0, 0), (-1, -1), 8),
        ('VALIGN',         (0, 0), (-1, -1), 'MIDDLE'),
    ]))
    return tbl


# ---------------------------------------------------------------------------
# Gene context helpers
# ---------------------------------------------------------------------------
def load_gene_ctx(path):
    """Load somatic_gene_context.tsv → {gene: {role, tier, pathway, summary}}."""
    ctx = {}
    if not path or not os.path.isfile(path):
        return ctx
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                gene = (row.get('GENE') or '').strip()
                if gene:
                    ctx[gene] = {
                        'role':    (row.get('ROLE') or '').strip(),
                        'tier':    (row.get('CGC_TIER') or '').strip(),
                        'pathway': (row.get('PATHWAY') or '').strip(),
                        'summary': (row.get('GENE_SUMMARY') or '').strip(),
                    }
    except Exception as exc:
        print(f'[WARN] Could not load gene context TSV {path}: {exc}', file=sys.stderr)
    return ctx


def gene_consequence(role, impact, effect):
    """Return a clinical consequence sentence based on gene role × impact × effect."""
    role   = (role   or '').upper()
    impact = (impact or '').upper()
    effect = (effect or '').lower()
    is_both = role == 'BOTH'
    is_tsg  = is_both or 'TSG' in role or role == 'TUMOR SUPPRESSOR'
    is_onco = is_both or 'ONCOGENE' in role

    if impact == 'HIGH':
        if is_tsg or is_both:
            return ('Loss-of-function variant (frameshift, nonsense or splice). '
                    'Truncating mutations in a tumor suppressor eliminate functional protein; '
                    'biallelic inactivation through LOH or a second somatic hit drives tumorigenesis.')
        elif is_onco:
            return ('High-impact variant in an oncogene — uncommon. '
                    'May represent loss of a regulatory domain or a rare neomorphic event. '
                    'Functional significance requires validation.')
        else:
            return 'High-impact variant (frameshift/nonsense/splice) of uncertain functional consequence.'

    elif impact == 'MODERATE':
        if 'missense' in effect or 'inframe' in effect:
            if is_onco and not is_both:
                return ('Missense variant in an oncogene. '
                        'Many oncogenic missense mutations are activating, locking the protein '
                        'in a constitutively active state and driving downstream pathway signaling.')
            elif is_tsg and not is_both:
                return ('Missense variant in a tumor suppressor. '
                        'May act as a dominant-negative allele (interfering with the wild-type copy) '
                        'or a hypomorphic allele with partial loss of function.')
            elif is_both:
                return ('Missense variant in a gene with dual TSG/oncogenic roles. '
                        'Functional impact depends on variant position; activating and inactivating '
                        'effects are both described.')
            else:
                return 'Missense variant of uncertain functional consequence.'
        else:
            return 'Moderate-impact variant; functional significance depends on the specific effect and gene context.'

    return 'Variant of uncertain functional consequence in this gene.'


def trunc(text, n=40):
    """Truncate text to n characters, appending '...' if needed."""
    s = str(text).strip() if text else ''
    if s in ('.', '', 'NA', 'N/A', 'nan'):
        return '-'
    return s if len(s) <= n else s[:n - 3] + '...'


def safe_float(val, default=None):
    """Return float or default if conversion fails."""
    try:
        f = float(str(val).strip())
        if f != f:   # nan check
            return default
        return f
    except (ValueError, TypeError):
        return default


def read_tsv(path):
    """
    Read a TSV file.  Returns (headers: list[str], rows: list[dict]) or ([], [])
    on failure.
    """
    if not path or not os.path.isfile(path):
        return [], []
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            rows = list(reader)
            headers = reader.fieldnames or []
        return list(headers), rows
    except Exception as exc:
        print(f'[WARN] Could not read {path}: {exc}', file=sys.stderr)
        return [], []


# ---------------------------------------------------------------------------
# ACMG classification — 9-tool in silico panel (mirrors final_report_acmg.py)
# ---------------------------------------------------------------------------
_INSIL_DAMAGING_PDF = {
    'dbNSFP_SIFT_pred':             ('D',),
    'dbNSFP_SIFT4G_pred':           ('D',),
    'dbNSFP_Polyphen2_HDIV_pred':   ('D', 'P'),
    'dbNSFP_Polyphen2_HVAR_pred':   ('D', 'P'),
    'dbNSFP_LRT_pred':              ('D',),
    'dbNSFP_MutationTaster_pred':   ('A', 'D'),
    'dbNSFP_MutationAssessor_pred': ('H', 'M'),
    'dbNSFP_FATHMM_pred':           ('D',),
    'dbNSFP_PROVEAN_pred':          ('D',),
}
_INSIL_BENIGN_PDF = {
    'dbNSFP_SIFT_pred':             ('T',),
    'dbNSFP_SIFT4G_pred':           ('T',),
    'dbNSFP_Polyphen2_HDIV_pred':   ('B',),
    'dbNSFP_Polyphen2_HVAR_pred':   ('B',),
    'dbNSFP_LRT_pred':              ('N',),
    'dbNSFP_MutationTaster_pred':   ('N', 'P'),
    'dbNSFP_MutationAssessor_pred': ('L', 'N'),
    'dbNSFP_FATHMM_pred':           ('T',),
    'dbNSFP_PROVEAN_pred':          ('N',),
}


def _count_insil_pdf(row):
    """Count damaging vs benign votes from 9-tool panel for a row dict.
    dbNSFP stores multi-transcript predictions comma-separated; we take the first.
    """
    dam = ben = 0
    for field in _INSIL_DAMAGING_PDF:
        v = str(row.get(field, '') or '').split(',')[0].strip().upper()
        if v in ('', '.', 'NA'):
            continue
        if v in _INSIL_DAMAGING_PDF[field]:
            dam += 1
        elif v in _INSIL_BENIGN_PDF.get(field, ()):
            ben += 1
    return dam, ben


def classify_acmg_with_criteria(row):
    """
    Classify a germline variant row into: P, LP, VUS, LB, B
    Returns (cls, criteria_list) where criteria_list contains the activated
    ACMG evidence code strings (e.g. ['PVS1', 'PM2', 'PP3']).
    Logic priority:
      1. ClinVar annotation
      2. gnomAD BA1 rule
      3. Evidence scoring (9-tool in silico panel)
    """
    clnsig = str(row.get('dbNSFP_clinvar_clnsig', '') or '').lower()
    gnomad_af_raw = row.get('dbNSFP_gnomAD_exomes_AF', None)
    gnomad_af = safe_float(gnomad_af_raw)
    impact = str(row.get('ANN[*].IMPACT', '') or '').upper()
    revel = safe_float(row.get('dbNSFP_REVEL_score', None))

    # 1. ClinVar rule (conflicting = fall through to score-based)
    if clnsig and clnsig not in ('.', '', 'na', 'nan'):
        if 'conflicting' in clnsig:
            pass  # fall through to evidence scoring
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

    # 2. BA1 – common variant
    if gnomad_af is not None and gnomad_af > 0.05:
        return 'B', ['BA1']

    # 3. Evidence scoring
    score = 0
    criteria = []

    # Pathogenic evidence
    if impact == 'HIGH':
        score += 8
        criteria.append('PVS1')
    if gnomad_af is None or gnomad_af < 0.001:
        score += 2
        criteria.append('PM2')
    # -- Computational evidence (PP3 / BP4) --
    # ClinGen SVI: REVEL and the 9-tool panel are evaluated together.
    # If they disagree, neither PP3 nor BP4 is applied.
    revel_dir = 0
    revel_str = 0
    revel_lbl = ''
    if revel is not None:
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

    dam, ben = _count_insil_pdf(row)
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

    # Benign evidence
    if gnomad_af is not None and gnomad_af > 0.01:
        score -= 4
        criteria.append('BS1')

    if score >= 10:
        cls = 'P'
    elif score >= 6:
        cls = 'LP'
    elif score <= -2:
        cls = 'LB'
    else:
        cls = 'VUS'
    return cls, criteria


def classify_acmg(row):
    """Thin wrapper — returns only the classification string (P/LP/VUS/LB/B)."""
    cls, _ = classify_acmg_with_criteria(row)
    return cls


def load_germline_variants(path):
    """
    Read germline TSV, run ACMG classification, and return list of
    dicts for variants classified as P, LP, or VUS only.
    Each dict gains '_ACMG' (classification) and '_ACMG_CRITERIA' (list of codes).
    """
    _, rows = read_tsv(path)
    results = []
    for row in rows:
        cls, criteria = classify_acmg_with_criteria(row)
        if cls in ('P', 'LP', 'VUS'):
            row['_ACMG'] = cls
            row['_ACMG_CRITERIA'] = criteria
            results.append(row)
    return results


# ---------------------------------------------------------------------------
# AMP tier classification
# ---------------------------------------------------------------------------
def load_civic_map(civic_path):
    """
    Build lookup: (CHROM, POS, REF, ALT) -> list[dict] of CIViC records.
    """
    _, rows = read_tsv(civic_path)
    civic_map = {}
    for row in rows:
        key = (
            str(row.get('CHROM', '')),
            str(row.get('POS', '')),
            str(row.get('REF', '')),
            str(row.get('ALT', '')),
        )
        civic_map.setdefault(key, []).append(row)
    return civic_map


def classify_amp(row, civic_map):
    """
    Classify a somatic variant into AMP tier using CIViC evidence.
    Returns (tier_str, civic_record_or_None).
    """
    key = (
        str(row.get('CHROM', '')),
        str(row.get('POS', '')),
        str(row.get('REF', '')),
        str(row.get('ALT', '')),
    )
    impact = str(row.get('ANN[*].IMPACT', '') or '').upper()
    civic_records = civic_map.get(key, [])

    if civic_records:
        # Pick best record by CIViC evidence level (A > B > C > D)
        level_order = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}
        def record_sort_key(r):
            lvl = str(r.get('CIViC Evidence Level', 'Z') or 'Z').upper()
            return level_order.get(lvl, 99)
        civic_records_sorted = sorted(civic_records, key=record_sort_key)
        best = civic_records_sorted[0]
        lvl = str(best.get('CIViC Evidence Level', '') or '').upper()
        sig = str(best.get('CIViC Entity Significance', '') or '').upper()

        if lvl == 'A':
            if 'SENSITIVITYRESPONSE' in sig or sig == 'SENSITIVITY/RESPONSE':
                return 'IA', best
            return 'IIA', best
        if lvl == 'B':
            if 'SENSITIVITYRESPONSE' in sig or sig == 'SENSITIVITY/RESPONSE' or 'RESISTANCE' in sig:
                return 'IB', best
            return 'IIB', best
        if lvl in ('C', 'D'):
            return 'IIC', best

    # No CIViC match – use impact
    if impact in ('HIGH', 'MODERATE'):
        return 'III', None
    return 'III', None


def load_somatic_variants(somatic_path, civic_path):
    """
    Read somatic TSV, classify each variant with AMP tiers.
    Returns list of variant dicts with '_AMP_TIER' and '_CIVIC' keys.
    """
    _, rows = read_tsv(somatic_path)
    civic_map = load_civic_map(civic_path) if (civic_path and os.path.isfile(civic_path or '')) else {}
    results = []
    for row in rows:
        tier, civic_rec = classify_amp(row, civic_map)
        row['_AMP_TIER'] = tier
        row['_CIVIC'] = civic_rec or {}
        results.append(row)
    return results


# ---------------------------------------------------------------------------
# Colour helpers for badges
# ---------------------------------------------------------------------------
def acmg_colour(cls):
    return {
        'P':   C_ACMG_P,
        'LP':  C_ACMG_LP,
        'VUS': C_ACMG_VUS,
    }.get(cls, C_AMP_III)


def amp_colour(tier):
    if tier in ('IA', 'IB'):
        return C_AMP_IA_IB
    if tier in ('IIA', 'IIB'):
        return C_AMP_IIA_IIB
    if tier == 'IIC':
        return C_AMP_IIC
    return C_AMP_III


# ---------------------------------------------------------------------------
# Document template with page numbers
# ---------------------------------------------------------------------------
class FooterCanvas(BaseDocTemplate):
    """
    Custom BaseDocTemplate that draws a persistent footer on every page.
    """
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.page_count = 0

    def handle_documentEnd(self):
        # Two-pass: first pass counts pages, second pass draws.
        pass

    def afterPage(self):
        self.page_count = self.page

    def build(self, flowables, onFirstPage=None, onLaterPages=None):
        self._calc()
        super().build(flowables, onFirstPage=onFirstPage, onLaterPages=onLaterPages)


def build_footer(canvas, doc):
    """Draw footer on every page."""
    canvas.saveState()
    footer_text = 'GENUIT  \u00b7  nf-sarek-clinical  \u00b7  ' + T('footer_conf')
    page_text = f'Page {doc.page}'
    style = ParagraphStyle(
        'footer',
        fontName='Helvetica',
        fontSize=7,
        textColor=colors.HexColor('#64748b'),
    )
    canvas.setFont('Helvetica', 7)
    canvas.setFillColor(colors.HexColor('#64748b'))
    y_pos = MARGIN_B - 0.8 * cm
    canvas.drawString(MARGIN_L, y_pos, footer_text)
    canvas.drawRightString(PAGE_W - MARGIN_R, y_pos, page_text)
    # thin separator line
    canvas.setStrokeColor(colors.HexColor('#cbd5e1'))
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN_L, y_pos + 0.4 * cm, PAGE_W - MARGIN_R, y_pos + 0.4 * cm)
    canvas.restoreState()


# ---------------------------------------------------------------------------
# Style definitions
# ---------------------------------------------------------------------------
styles = getSampleStyleSheet()

STYLE_H1 = ParagraphStyle(
    'h1_report',
    fontName='Helvetica-Bold',
    fontSize=18,
    textColor=C_WHITE,
    spaceAfter=4,
)
STYLE_H2 = ParagraphStyle(
    'h2_report',
    fontName='Helvetica-Bold',
    fontSize=12,
    textColor=C_DARK_HEADER,
    spaceBefore=12,
    spaceAfter=4,
)
STYLE_H3 = ParagraphStyle(
    'h3_report',
    fontName='Helvetica-Bold',
    fontSize=10,
    textColor=C_DARK_HEADER,
    spaceBefore=8,
    spaceAfter=3,
)
STYLE_SUBTITLE = ParagraphStyle(
    'subtitle',
    fontName='Helvetica',
    fontSize=9,
    textColor=colors.HexColor('#94a3b8'),
    spaceAfter=2,
)
STYLE_BODY = ParagraphStyle(
    'body_report',
    fontName='Helvetica',
    fontSize=9,
    textColor=C_TEXT,
    spaceAfter=4,
)
STYLE_CELL = ParagraphStyle(
    'cell',
    fontName='Helvetica',
    fontSize=7,
    textColor=C_TEXT,
    leading=9,
)
STYLE_CELL_BOLD = ParagraphStyle(
    'cell_bold',
    fontName='Helvetica-Bold',
    fontSize=7,
    textColor=C_TEXT,
    leading=9,
)
STYLE_BADGE = ParagraphStyle(
    'badge',
    fontName='Helvetica-Bold',
    fontSize=7,
    textColor=C_WHITE,
    alignment=TA_CENTER,
    leading=9,
)
STYLE_DISCLAIMER = ParagraphStyle(
    'disclaimer',
    fontName='Helvetica-Oblique',
    fontSize=7,
    textColor=colors.HexColor('#64748b'),
    spaceAfter=4,
    leading=10,
)
STYLE_WARN = ParagraphStyle(
    'warn',
    fontName='Helvetica-Oblique',
    fontSize=8,
    textColor=colors.HexColor('#92400e'),
    spaceAfter=6,
)


# ---------------------------------------------------------------------------
# Coloured rectangle badge flowable
# ---------------------------------------------------------------------------
class Badge(Flowable):
    """A small coloured badge with white text, used for ACMG/AMP classification."""

    def __init__(self, label, bg_colour, width=32, height=14):
        super().__init__()
        self.label = label
        self.bg_colour = bg_colour
        self.width = width
        self.height = height

    def wrap(self, available_w, available_h):
        return self.width, self.height

    def draw(self):
        self.canv.setFillColor(self.bg_colour)
        self.canv.roundRect(0, 0, self.width, self.height, 3, fill=1, stroke=0)
        self.canv.setFillColor(C_WHITE)
        self.canv.setFont('Helvetica-Bold', 7)
        self.canv.drawCentredString(self.width / 2, 3, self.label)


# ---------------------------------------------------------------------------
# Cover page builder
# ---------------------------------------------------------------------------
def build_cover_page(args, metrics, story):
    """Append cover-page flowables: header + metadata + Global Biomarkers KPI.
    Does NOT append a PageBreak — caller continues with Major Findings.
    """
    # ── Header banner ────────────────────────────────────────────────────────
    header_data = [[
        Paragraph(
            f'<font color="white" size="16"><b>{T("report_title")}</b></font>',
            styles['Normal'],
        ),
        Paragraph(
            f'<font color="#f97316" size="9"><b>{args.pair}</b></font>',
            ParagraphStyle('hdr_right', fontName='Helvetica', fontSize=9,
                           textColor=C_ACCENT, alignment=TA_RIGHT),
        ),
    ]]
    header_table = Table(header_data, colWidths=[USABLE_W * 0.7, USABLE_W * 0.3])
    header_table.setStyle(TableStyle([
        ('BACKGROUND',  (0, 0), (-1, -1), C_DARK_HEADER),
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), [C_DARK_HEADER]),
        ('TOPPADDING',  (0, 0), (-1, -1), 12),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
        ('LEFTPADDING',  (0, 0), (-1, -1), 14),
        ('RIGHTPADDING', (0, 0), (-1, -1), 14),
        ('VALIGN',       (0, 0), (-1, -1), 'MIDDLE'),
        ('ROUNDEDCORNERS', [6]),
    ]))
    story.append(header_table)
    story.append(Spacer(1, 0.4 * cm))

    # ── Metadata table ────────────────────────────────────────────────────────
    meta_rows = [
        [T('meta_run'),   args.run_name or 'N/A', T('meta_patient'), args.patient or 'N/A'],
        [T('meta_tumor'), args.tumor or 'N/A',    T('meta_normal'),  args.normal or 'N/A'],
        [T('meta_date'),  date.today().strftime('%Y-%m-%d'), T('meta_pipeline'), 'nf-sarek-clinical'],
    ]
    col_w = USABLE_W / 4
    meta_table = Table(meta_rows, colWidths=[col_w] * 4)
    meta_table.setStyle(TableStyle([
        ('BACKGROUND',    (0, 0), (-1, -1), colors.HexColor('#f1f5f9')),
        ('TEXTCOLOR',     (0, 0), (0, -1), C_SLATE),
        ('TEXTCOLOR',     (2, 0), (2, -1), C_SLATE),
        ('FONTNAME',      (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTNAME',      (2, 0), (2, -1), 'Helvetica-Bold'),
        ('FONTNAME',      (1, 0), (1, -1), 'Helvetica'),
        ('FONTNAME',      (3, 0), (3, -1), 'Helvetica'),
        ('FONTSIZE',      (0, 0), (-1, -1), 8),
        ('TOPPADDING',    (0, 0), (-1, -1), 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
        ('LEFTPADDING',   (0, 0), (-1, -1), 10),
        ('RIGHTPADDING',  (0, 0), (-1, -1), 8),
        ('GRID',          (0, 0), (-1, -1), 0.3, C_LIGHT_BORDER),
        ('ROWBACKGROUNDS', (0, 0), (-1, -1),
         [colors.HexColor('#f8fafc'), colors.HexColor('#f1f5f9')]),
    ]))
    story.append(meta_table)
    story.append(Spacer(1, 0.5 * cm))

    # ── Determine low-coverage flag ────────────────────────────────────────────
    _LOW_COV_THRESHOLD = 50   # mean X — below this, flag quality concern
    try:
        _cov_t_f = float(metrics.get('cov_t', 'NA'))
    except (ValueError, TypeError):
        _cov_t_f = None
    is_low_cov = _cov_t_f is not None and _cov_t_f < _LOW_COV_THRESHOLD
    metrics['is_low_cov'] = is_low_cov   # propagate to other sections

    # ── Global Biomarkers (4-column) ──────────────────────────────────────────
    story.append(Paragraph(T('sect_global_bm'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=6))

    _indet = f'<font color="#b91c1c"><b>{T("indeterminate")}</b></font>'

    _msi_raw_pdf = metrics.get('msi', 'NA')
    if is_low_cov:
        _msi_disp_pdf = _indet
    else:
        try:
            _msi_pf = float(_msi_raw_pdf)
            _msi_st = T('msi_detected') if _msi_pf >= 5.0 else T('msi_not_det')
            _msi_disp_pdf = _msi_st + '<br/>' + _msi_raw_pdf + '%'
        except (ValueError, TypeError):
            _msi_disp_pdf = _msi_raw_pdf

    _hrd_prob_disp = _indet if is_low_cov else _fmt_hrd_prob(metrics.get('hrd_prob', 'NA'))

    global_kpis = [
        (T('kpi_tmb'),      _indet if is_low_cov else metrics.get('tmb', 'N/A')),
        (T('kpi_msi'),      _msi_disp_pdf),
        (T('kpi_hrd_prob'), _hrd_prob_disp),
        (T('kpi_hrd_stat'), _indet if is_low_cov else metrics.get('hrd_status', 'N/A')),
    ]
    story.append(_build_kpi_table(global_kpis, ncols=4))
    story.append(Spacer(1, 0.4 * cm))

    # ── Low tumor-coverage warning ───────────────────────────────────────────
    if is_low_cov:
        _cov_str = f'{_cov_t_f:.0f}'
        _warn_items = ''.join(
            f'<br/>&nbsp;&nbsp;• {T(k)}'
            for k in ('low_cov_tmb', 'low_cov_msi', 'low_cov_hrd',
                      'low_cov_mutsig', 'low_cov_somatic', 'low_cov_cnv')
        )
        _warn_body = T('low_cov_body').format(cov=_cov_str) + _warn_items
        _warn_s = ParagraphStyle(
            'LowCovWarn', fontName='Helvetica', fontSize=7.5,
            textColor=colors.HexColor('#7c2d12'), leading=11)
        _warn_data = [[
            Paragraph(
                f'<font size="9"><b>⚠ {T("low_cov_title")}</b></font><br/>'
                f'{_warn_body}',
                _warn_s),
        ]]
        _warn_tbl = Table(_warn_data, colWidths=[USABLE_W])
        _warn_tbl.setStyle(TableStyle([
            ('BACKGROUND',   (0, 0), (-1, -1), colors.HexColor('#fef3c7')),
            ('TOPPADDING',   (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING',(0, 0), (-1, -1), 8),
            ('LEFTPADDING',  (0, 0), (-1, -1), 12),
            ('RIGHTPADDING', (0, 0), (-1, -1), 12),
            ('ROUNDEDCORNERS', [4]),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#f59e0b')),
        ]))
        story.append(_warn_tbl)
        story.append(Spacer(1, 0.3 * cm))

    # Story continues immediately with Major Clinical Findings (no PageBreak).


# ---------------------------------------------------------------------------
# Secondary Biomarkers + Disclaimer + PageBreak  (end of page 1)
# ---------------------------------------------------------------------------
def build_secondary_biomarkers(metrics, story):
    """Append Secondary Biomarkers section and Disclaimer."""
    story.append(PageBreak())
    story.append(Paragraph(T('sect_secondary_bm'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=6))

    secondary_kpis = [
        (T('kpi_wgd'),    metrics.get('wgd',    'N/A')),
        (T('kpi_gi'),     metrics.get('gi',     'N/A')),
        (T('kpi_loh'),    metrics.get('loh',    'N/A')),
        (T('kpi_ploidy'), metrics.get('ploidy', 'N/A')),
    ]
    story.append(_build_kpi_table(secondary_kpis, ncols=4))
    story.append(Spacer(1, 0.5 * cm))
    story.append(Paragraph(T('disclaimer'), STYLE_DISCLAIMER))


# ---------------------------------------------------------------------------
# Genomic Quality Metrics page
# ---------------------------------------------------------------------------
def build_quality_metrics_page(metrics, story):
    """Append Genomic Quality Metrics section (Purity + Coverage)."""
    story.append(PageBreak())
    story.append(Paragraph(T('sect_quality'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=6))

    quality_kpis = [
        (T('kpi_purity'), metrics.get('purity', 'N/A')),
        (T('kpi_cov_n'),  metrics.get('cov_n',  'N/A')),
        (T('kpi_cov_t'),  metrics.get('cov_t',  'N/A')),
    ]
    story.append(_build_kpi_table(quality_kpis, ncols=3))
    story.append(Spacer(1, 0.6 * cm))




# ---------------------------------------------------------------------------
# Table helper with alternating row colours and header
# ---------------------------------------------------------------------------
def build_variant_table(col_headers, col_widths, data_rows, colour_col_idx=None, colour_fn=None):
    """
    Build a ReportLab Table with styled header row and alternating body rows.

    colour_col_idx – index of the column that gets a coloured badge background
    colour_fn      – callable(cell_value) -> Color
    """
    # Header row
    header_cells = [
        Paragraph(h, ParagraphStyle(
            'th', fontName='Helvetica-Bold', fontSize=7,
            textColor=C_WHITE, leading=9
        ))
        for h in col_headers
    ]
    table_data = [header_cells]

    for i, row_vals in enumerate(data_rows):
        row_cells = []
        for j, val in enumerate(row_vals):
            if j == colour_col_idx and colour_fn is not None:
                colour = colour_fn(str(val))
                cell = Paragraph(
                    f'<font color="white"><b>{val}</b></font>',
                    ParagraphStyle('badge_cell', fontName='Helvetica-Bold',
                                   fontSize=7, textColor=C_WHITE, alignment=TA_CENTER, leading=9),
                )
                row_cells.append(cell)
            else:
                # If val is already a ReportLab Paragraph/Flowable, use it directly.
                # Calling Paragraph(str(val), …) on a Paragraph object would produce
                # the garbage text "<Paragraph object at 0x…>" in the PDF.
                if isinstance(val, Paragraph):
                    row_cells.append(val)
                else:
                    row_cells.append(Paragraph(str(val), STYLE_CELL))
        table_data.append(row_cells)

    tbl = Table(table_data, colWidths=col_widths, repeatRows=1)

    style_cmds = [
        ('BACKGROUND',    (0, 0), (-1, 0),  C_SLATE),
        ('TEXTCOLOR',     (0, 0), (-1, 0),  C_WHITE),
        ('FONTNAME',      (0, 0), (-1, 0),  'Helvetica-Bold'),
        ('FONTSIZE',      (0, 0), (-1, -1), 7),
        ('TOPPADDING',    (0, 0), (-1, -1), 4),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('LEFTPADDING',   (0, 0), (-1, -1), 5),
        ('RIGHTPADDING',  (0, 0), (-1, -1), 5),
        ('GRID',          (0, 0), (-1, -1), 0.3, C_LIGHT_BORDER),
        ('VALIGN',        (0, 0), (-1, -1), 'TOP'),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [C_WHITE, C_ROW_ALT]),
    ]

    # Colour background for badge column in data rows
    if colour_col_idx is not None and colour_fn is not None:
        for i, row_vals in enumerate(data_rows):
            val = str(row_vals[colour_col_idx])
            bg = colour_fn(val)
            r = i + 1  # offset for header
            style_cmds.append(('BACKGROUND', (colour_col_idx, r), (colour_col_idx, r), bg))

    tbl.setStyle(TableStyle(style_cmds))
    return tbl


# ---------------------------------------------------------------------------
# Germline section
# ---------------------------------------------------------------------------
def build_germline_section(germline_variants, story):
    """Append germline SNV section flowables to story."""
    story.append(PageBreak())
    story.append(Paragraph(T('sect_germline'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))

    # Legend
    legend_items = [
        ('P',   C_ACMG_P,   T('acmg_p')),
        ('LP',  C_ACMG_LP,  T('acmg_lp')),
        ('VUS', C_ACMG_VUS, T('acmg_vus')),
    ]
    legend_parts = '&nbsp;&nbsp;&nbsp;'.join(
        f'<font color="{c.hexval()}"><b>[{lbl}]</b></font> {desc}'
        for lbl, c, desc in legend_items
    )
    story.append(Paragraph(f'<font size="7">{legend_parts}</font>', styles['Normal']))
    story.append(Spacer(1, 0.2 * cm))

    if not germline_variants:
        story.append(Paragraph(T('germ_no_vars'), STYLE_WARN))
        return

    # Table columns (widths sum to USABLE_W ~553 pt)
    # Keep every header short enough to fit its column at font-size 7 without wrapping.
    # 'gnomAD AF' (9 chars) overflowed 40 pt → shortened to 'gnomAD'.
    # 'ACMG' bumped from 28→32; 'Criteria' from 86→90 to give criteria badges room.
    col_headers = [
        'Gene', 'HGVS_c', 'HGVS_p', 'Impact', 'DP', 'VAF',
        'REVEL', 'gnomAD', 'ClinVar', 'ACMG', 'Criteria'
    ]
    col_widths = [44, 84, 72, 36, 22, 26, 28, 38, 79, 32, 92]
    # sum = 44+84+72+36+22+26+28+38+79+32+92 = 553 ✓

    _crit_style = ParagraphStyle(
        'CritCell',
        parent=styles['Normal'],
        fontSize=6,
        textColor=colors.HexColor('#475569'),
        leading=8,
    )

    data_rows = []
    for v in germline_variants:
        gene   = trunc(v.get('ANN[*].GENE', '-'), 20)
        hgvsc  = trunc(v.get('ANN[*].HGVS_C', '-'), 30)
        hgvsp  = trunc(v.get('ANN[*].HGVS_P', '-'), 26)
        impact = trunc(v.get('ANN[*].IMPACT', '-'), 12)
        # DP and VAF
        dp_raw = v.get('GEN[0].DP', '')
        ad_raw = v.get('GEN[0].AD', '')
        try:
            dp_int = int(dp_raw)
            alt_reads = int(ad_raw.split(',')[1])
            vaf_str = '{:.2f}'.format(alt_reads / dp_int) if dp_int > 0 else '.'
        except Exception:
            dp_int  = dp_raw or '-'
            vaf_str = '-'
        revel  = trunc(v.get('dbNSFP_REVEL_score', '-'), 8)
        gnomad = trunc(v.get('dbNSFP_gnomAD_exomes_AF', '-'), 12)
        clnsig = trunc(v.get('dbNSFP_clinvar_clnsig', '-'), 42)
        cls    = v.get('_ACMG', 'VUS')
        criteria_codes = v.get('_ACMG_CRITERIA', [])
        # Colour-code each criterion: pathogenic (red) vs benign (blue)
        _BENIGN_CRIT = {'BA1', 'BS1', 'BP4', 'BP4_St'}
        parts = []
        for code in criteria_codes:
            if code.startswith('ClinVar'):
                colour = '#7c3aed'
            elif code in _BENIGN_CRIT:
                colour = '#3b82f6'
            else:
                colour = '#dc2626'
            parts.append(f'<font color="{colour}"><b>{code}</b></font>')
        crit_para = Paragraph(' · '.join(parts) if parts else '—', _crit_style)
        data_rows.append([gene, hgvsc, hgvsp, impact, str(dp_int), vaf_str,
                          revel, gnomad, clnsig, cls, crit_para])

    def _acmg_bg(val):
        return acmg_colour(val.strip() if isinstance(val, str) else 'VUS')

    tbl = build_variant_table(
        col_headers, col_widths, data_rows,
        colour_col_idx=9, colour_fn=_acmg_bg,
    )
    story.append(tbl)
    story.append(Spacer(1, 0.3 * cm))
    story.append(Paragraph(
        f'<font size="7" color="#64748b">'
        f'{len(germline_variants)} variant(s) shown (P/LP/VUS only). '
        f'Benign and likely benign variants are excluded from this report.</font>',
        styles['Normal'],
    ))


# ---------------------------------------------------------------------------
# Somatic section
# ---------------------------------------------------------------------------
def _somatic_row(v, include_civic=False):
    """Extract display fields from a somatic variant dict."""
    civic = v.get('_CIVIC') or {}
    gene   = trunc(v.get('ANN[*].GENE', '-'), 20)
    hgvsc  = trunc(v.get('ANN[*].HGVS_C', '-'), 32)
    hgvsp  = trunc(v.get('ANN[*].HGVS_P', '-'), 28)
    effect = trunc(v.get('ANN[*].EFFECT', '-'), 22)
    af_raw = v.get('GEN[1].AF', v.get('GEN[1].AD', '-'))
    try:
        af_val = float(str(af_raw).split(',')[0])
        af_str = f'{af_val:.3f}'
        low_vaf = af_val < 0.05
    except (ValueError, TypeError):
        af_str = trunc(af_raw, 8)
        low_vaf = False
    # VAF cell with low-VAF warning
    if low_vaf:
        vaf_cell = Paragraph(
            f'<font color="#b91c1c"><b>{af_str}</b></font>'
            f'<br/><font size="5" color="#b91c1c">\u26a0 low</font>',
            ParagraphStyle('vaf_warn', fontSize=7, fontName='Helvetica',
                           textColor=C_TEXT, leading=8, alignment=TA_CENTER))
    else:
        vaf_cell = af_str
    callers_raw = v.get('Callers', 'Mutect2')
    if 'Strelka2' in callers_raw:
        callers_cell = Paragraph(
            '<font size="5.5" color="#166534">M+S</font>',
            ParagraphStyle('cal2', fontSize=5.5, alignment=TA_CENTER, leading=7))
    else:
        callers_cell = Paragraph(
            '<font size="5.5" color="#92400e">M</font>',
            ParagraphStyle('cal1', fontSize=5.5, alignment=TA_CENTER, leading=7))
    tier = v.get('_AMP_TIER', 'III')
    if include_civic:
        disease = trunc(civic.get('CIViC Entity Disease', '-'), 30)
        therapy = trunc(civic.get('CIViC Entity Therapies', '-'), 30)
        evid    = trunc(civic.get('CIViC Evidence Level', '-'), 3)
        return [gene, hgvsc, hgvsp, effect, vaf_cell, callers_cell, evid, disease, therapy, tier]
    return [gene, hgvsc, hgvsp, effect, vaf_cell, callers_cell, tier]


def build_somatic_section(somatic_variants, story, gene_ctx=None):
    """Append somatic SNV section: CIViC actionable table + other variants table."""
    story.append(PageBreak())
    story.append(Paragraph(T('sect_somatic'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))

    if not somatic_variants:
        story.append(Paragraph(
            'No somatic variants were identified in this sample.',
            STYLE_WARN,
        ))
        return

    # Split into CIViC actionable (Tier I/II) vs other (Tier III)
    _TIER_ORDER = {'IA': 0, 'IB': 1, 'IIA': 2, 'IIB': 3, 'IIC': 4, 'III': 5}
    actionable = [v for v in somatic_variants if v.get('_AMP_TIER', 'III') not in ('III',)]
    other      = [v for v in somatic_variants if v.get('_AMP_TIER', 'III') == 'III']
    actionable.sort(key=lambda v: _TIER_ORDER.get(v.get('_AMP_TIER', 'III'), 5))

    # Legend
    _legend_style = ParagraphStyle('som_legend', fontSize=6, leading=8,
                                    textColor=colors.HexColor('#64748b'))
    story.append(Paragraph(
        '<b>Callers:</b> <font color="#166534">M+S</font> = Mutect2 + Strelka2 (consensus) '
        '&nbsp; <font color="#92400e">M</font> = Mutect2 only '
        '&nbsp;&nbsp; <font color="#b91c1c">\u26a0 low</font> = VAF &lt; 5%',
        _legend_style))
    story.append(Spacer(1, 0.2 * cm))

    # ── Table 1: CIViC Actionable Variants (Tier I/II) ─────────────────────
    if actionable:
        story.append(Paragraph(
            '<b>Clinically Actionable Variants</b> (CIViC Level A\u2013D, AMP Tier I\u2013II)',
            ParagraphStyle('sub_h', fontSize=8, fontName='Helvetica-Bold',
                           textColor=C_TEXT, spaceBefore=2, spaceAfter=4)))

        act_headers = ['Gene', 'HGVS_c', 'HGVS_p', 'Effect', 'VAF',
                        'Call', 'Evid', 'Disease', 'Therapy', 'Tier']
        act_widths  = [45, 85, 75, 65, 32, 25, 25, 80, 88, 30]

        act_rows = [_somatic_row(v, include_civic=True) for v in actionable]

        tbl_act = build_variant_table(
            act_headers, act_widths, act_rows,
            colour_col_idx=9, colour_fn=lambda v: amp_colour(v.strip()),
        )
        story.append(tbl_act)
        story.append(Spacer(1, 0.15 * cm))
        story.append(Paragraph(
            f'<font size="6" color="#64748b">'
            f'{len(actionable)} actionable variant(s). '
            f'CIViC annotations based on database snapshot at analysis date.</font>',
            styles['Normal']))
    else:
        story.append(Paragraph(
            '<b>Clinically Actionable Variants:</b> None identified (no CIViC Level A\u2013D matches).',
            ParagraphStyle('no_act', fontSize=7.5, textColor=colors.HexColor('#64748b'),
                           spaceBefore=2, spaceAfter=4)))

    # ── Table 2: Other Somatic Variants (Tier III) ─────────────────────────
    if other:
        story.append(Spacer(1, 0.3 * cm))
        story.append(Paragraph(
            '<b>Other Somatic Variants</b> (AMP Tier III \u2014 uncertain clinical significance)',
            ParagraphStyle('sub_h2', fontSize=8, fontName='Helvetica-Bold',
                           textColor=C_TEXT, spaceBefore=2, spaceAfter=4)))

        oth_headers = ['Gene', 'HGVS_c', 'HGVS_p', 'Effect', 'VAF', 'Call', 'Tier']
        oth_widths  = [50, 120, 105, 95, 38, 30, 30]

        oth_rows = [_somatic_row(v, include_civic=False) for v in other]

        tbl_oth = build_variant_table(
            oth_headers, oth_widths, oth_rows,
            colour_col_idx=6, colour_fn=lambda v: amp_colour(v.strip()),
        )
        story.append(tbl_oth)
        story.append(Spacer(1, 0.15 * cm))

    # Summary
    story.append(Paragraph(
        f'<font size="6.5" color="#64748b">'
        f'{len(somatic_variants)} somatic variant(s) total: '
        f'{len(actionable)} actionable (Tier I/II), {len(other)} Tier III.</font>',
        styles['Normal'],
    ))

    # ── Gene Annotation Notes ────────────────────────────────────────────────
    if gene_ctx:
        # Collect unique genes in order of appearance
        seen_genes = []
        for v in somatic_variants:
            gene = (v.get('ANN[*].GENE') or '').strip()
            if gene and gene not in seen_genes and gene in gene_ctx:
                seen_genes.append(gene)

        if seen_genes:
            story.append(Spacer(1, 0.4 * cm))
            story.append(Paragraph('Gene Annotation Notes', STYLE_H2))
            story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))

            # Style for gene note entries
            _style_note_header = ParagraphStyle(
                'note_header', fontName='Helvetica-Bold', fontSize=7.5,
                textColor=C_TEXT, spaceAfter=2, spaceBefore=0,
            )
            _style_note_body = ParagraphStyle(
                'note_body', fontName='Helvetica', fontSize=7,
                textColor=colors.HexColor('#334155'), spaceAfter=3, spaceBefore=0,
                leading=9.5,
            )

            for gene in seen_genes:
                ctx = gene_ctx[gene]
                role    = ctx.get('role', '')
                tier    = ctx.get('tier', '')
                pathway = ctx.get('pathway', '')
                summary = ctx.get('summary', '')

                # Find the worst impact observed for this gene across variants
                impacts_order = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
                best_impact = ''
                best_effect = ''
                for v in somatic_variants:
                    if (v.get('ANN[*].GENE') or '').strip() == gene:
                        imp = (v.get('ANN[*].IMPACT') or '').upper()
                        eff = (v.get('ANN[*].EFFECT') or v.get('ANN[*].ANNOTATION') or '').lower()
                        if impacts_order.get(imp, -1) > impacts_order.get(best_impact, -1):
                            best_impact = imp
                            best_effect = eff

                conseq = gene_consequence(role, best_impact, best_effect)

                tier_str = f' · CGC Tier {tier}' if tier else ''
                path_str = f' · {pathway}' if pathway else ''
                header_text = f'{gene}  ({role}{path_str}{tier_str})'

                # Wrap in a light-background block
                note_data = [[
                    Paragraph(header_text, _style_note_header),
                ]]
                note_body_parts = []
                if summary:
                    note_body_parts.append(Paragraph(
                        f'<b>Function:</b> {summary}', _style_note_body))
                if conseq:
                    note_body_parts.append(Paragraph(
                        f'<b>Tumor consequence:</b> {conseq}', _style_note_body))

                block_content = [Paragraph(header_text, _style_note_header)] + note_body_parts
                block_tbl = Table(
                    [[block_content]],
                    colWidths=[USABLE_W - 12],
                )
                block_tbl.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#f8fafc')),
                    ('BOX',        (0, 0), (-1, -1), 0.5, colors.HexColor('#e2e8f0')),
                    ('LEFTPADDING',  (0, 0), (-1, -1), 8),
                    ('RIGHTPADDING', (0, 0), (-1, -1), 8),
                    ('TOPPADDING',   (0, 0), (-1, -1), 5),
                    ('BOTTOMPADDING',(0, 0), (-1, -1), 5),
                    ('VALIGN',      (0, 0), (-1, -1), 'TOP'),
                ]))
                story.append(KeepTogether([block_tbl, Spacer(1, 0.15 * cm)]))


# ---------------------------------------------------------------------------
# Germline-Somatic Correlation section (reads JSON from germline_somatic_correlate.py)
# ---------------------------------------------------------------------------
INTERP_COLORS = {
    'CONCORDANT':             colors.HexColor('#ef4444'),
    'VUS_UPGRADED':           colors.HexColor('#f97316'),
    'BIOMARKER_NO_2ND_HIT':   colors.HexColor('#f97316'),
    'BIALLELIC_NO_BIOMARKER': colors.HexColor('#f97316'),
    'BIALLELIC_OTHER':        colors.HexColor('#f97316'),
    'VUS_UNCERTAIN':          colors.HexColor('#eab308'),
    'GERMLINE_ONLY':          colors.HexColor('#6b7280'),
}

PATHWAY_COLORS = {
    'MMR':       colors.HexColor('#6366f1'),
    'HRR':       colors.HexColor('#0ea5e9'),
    'Other TSG': colors.HexColor('#8b5cf6'),
}


def _run_correlation_script(args, tmb_val, msi_val, hrd_prob_val):
    corr_json = getattr(args, 'corr_json', '') or ''
    if corr_json and os.path.isfile(corr_json):
        try:
            with open(corr_json) as fh:
                return json.load(fh)
        except Exception:
            pass

    # Try to run the correlation script inline
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'germline_somatic_correlate.py')
    if not os.path.isfile(script):
        return None

    tmp_json = args.out + '.corr.tmp.json'
    cmd = [
        sys.executable, script,
        '--germline',  getattr(args, 'germline', '') or '',
        '--somatic',   getattr(args, 'somatic', '') or '',
        '--cnv',       getattr(args, 'cnv_somatic', '') or '',
        '--sig',       getattr(args, 'sig_sbs96', '') or '',
        '--tumor',     getattr(args, 'tumor', '') or '',
        '--tmb',       str(tmb_val),
        '--msi',       str(msi_val),
        '--hrd',       str(hrd_prob_val),
        '--json-out',  tmp_json,
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        if os.path.isfile(tmp_json):
            with open(tmp_json) as fh:
                data = json.load(fh)
            os.unlink(tmp_json)
            return data
    except Exception as exc:
        print('[WARN] Correlation script failed: ' + str(exc), file=sys.stderr)
    return None


def build_signature_plots_pdf(act_path, sig_type, story):
    """Embed signature plots into PDF story.

    Priority order per plot type:
      1. Custom PNG files generated by generate_signature_plots.py
         (matrix profile + etiology activity plot)
      2. SigProfiler PDF outputs (fallback via fitz)
    """
    if not act_path or not os.path.isfile(act_path):
        return

    import io as _io
    from reportlab.platypus import Image as _RLImage

    act_dir = os.path.dirname(act_path)
    caption_style = ParagraphStyle('SigCaption', parent=styles['Normal'],
                                   fontSize=7, textColor=colors.Color(0.6, 0.7, 0.8),
                                   spaceAfter=2)

    def _png_size_px(png_path):
        """Read PNG width/height from the IHDR chunk — no PIL needed."""
        import struct
        with open(png_path, 'rb') as f:
            f.read(8)   # PNG signature
            f.read(8)   # IHDR length + 'IHDR'
            w = struct.unpack('>I', f.read(4))[0]
            h = struct.unpack('>I', f.read(4))[0]
        return w, h

    def _embed_png(png_path, caption, max_w):
        """Embed a PNG file directly (no fitz / PIL needed)."""
        if not os.path.isfile(png_path):
            return False
        try:
            img_w_px, img_h_px = _png_size_px(png_path)
            dpi = 150.0
            img_w_pt = img_w_px / dpi * 72
            img_h_pt = img_h_px / dpi * 72
            scale = min(1.0, max_w / img_w_pt)
            story.append(Paragraph(caption, caption_style))
            story.append(_RLImage(png_path, width=img_w_pt * scale, height=img_h_pt * scale))
            return True
        except Exception as e:
            print('[WARN] Could not embed PNG ' + png_path + ': ' + str(e), file=sys.stderr)
            return False

    def _embed_pdf_page(pdf_path, caption, max_w):
        """Embed first page of a PDF as a rasterised image (requires fitz)."""
        if not os.path.isfile(pdf_path):
            return False
        try:
            import fitz as _fitz
            doc = _fitz.open(pdf_path)
            if doc.page_count == 0:
                return False
            mat = _fitz.Matrix(2.0, 2.0)
            pix = doc[0].get_pixmap(matrix=mat)
            png_bytes = pix.tobytes('png')
            doc.close()
            img_w = pix.width / 2.0
            img_h = pix.height / 2.0
            scale = min(1.0, max_w / img_w)
            story.append(Paragraph(caption, caption_style))
            story.append(_RLImage(_io.BytesIO(png_bytes), width=img_w * scale, height=img_h * scale))
            return True
        except Exception as e:
            print('[WARN] Could not embed PDF plot ' + pdf_path + ': ' + str(e), file=sys.stderr)
            return False

    # ── Derive the matrix_profiles/output/plots/ directory ──────────────────
    # act_dir: .../somatic/signatures/assignment/<proj>_SBS96/Assignment_Solution/Activities/
    # sig_root: .../somatic/signatures/
    sig_root = os.path.normpath(os.path.join(act_dir, '..', '..', '..', '..'))
    mg_plots_dir = os.path.join(sig_root, 'matrix_profiles', 'output', 'plots')

    # ── 1. Canonical SigProfilerMatrixGenerator spectrum plot ───────────────
    # Pattern: SBS_96_plots_*.pdf  /  ID_83_plots_*.pdf  (1 page = observed spectrum)
    mg_pdf = None
    if os.path.isdir(mg_plots_dir):
        prefix = 'SBS_96_plots_' if sig_type == 'SBS96' else 'ID_83_plots_'
        for fn in sorted(os.listdir(mg_plots_dir)):
            if fn.startswith(prefix) and fn.endswith('.pdf'):
                mg_pdf = os.path.join(mg_plots_dir, fn)
                break

    caption_obs = ('SBS96 Mutational Spectrum — Observed' if sig_type == 'SBS96'
                   else 'ID83 Indel Spectrum — Observed')
    if mg_pdf and not _embed_pdf_page(mg_pdf, caption_obs, USABLE_W):
        mg_pdf = None  # fall through to custom PNG

    if not mg_pdf:
        # Fallback: custom PNG from generate_signature_plots.py
        profile_png = 'SBS96_matrix_profile.png' if sig_type == 'SBS96' else 'ID83_matrix_profile.png'
        _embed_png(os.path.join(act_dir, profile_png), caption_obs, USABLE_W)

    story.append(Spacer(1, 0.2 * cm))

    # ── 2. Etiology-annotated activity bar chart ────────────────────────────
    if not _embed_png(os.path.join(act_dir, 'etiology_activity_plot.png'),
                      'COSMIC Signature Assignment \u2014 Aetiology (' + sig_type + ')', USABLE_W):
        _embed_pdf_page(
            os.path.join(act_dir, 'Assignment_Solution_Activity_Plots.pdf'),
            'COSMIC Signature Assignment \u2014 ' + sig_type, USABLE_W * 0.55)


# ---------------------------------------------------------------------------
# Mutational Signatures section
# ---------------------------------------------------------------------------
_SIG_AETIOLOGY = {
    'SBS1':  'Deamination of 5-methylcytosine (ageing)',
    'SBS2':  'APOBEC cytidine deaminase',
    'SBS3':  'Defective homologous recombination (HRD)',
    'SBS4':  'Tobacco mutagen exposure',
    'SBS5':  'Clock-like replication errors',
    'SBS6':  'Mismatch repair deficiency (MMR)',
    'SBS7a': 'UV light exposure',
    'SBS7b': 'UV light exposure',
    'SBS8':  'Defective homologous recombination',
    'SBS13': 'APOBEC cytidine deaminase',
    'SBS15': 'MMR deficiency (MLH1 hypermethylation)',
    'SBS18': 'Reactive oxygen species damage',
    'SBS20': 'MMR deficiency + POLD1 mutation',
    'SBS21': 'MMR deficiency (MSI-H)',
    'SBS26': 'MMR deficiency',
    'SBS30': 'Defective base excision repair (NTHL1)',
    'SBS44': 'MMR deficiency',
    'ID1':   'Replication slippage at homopolymers',
    'ID2':   'Replication slippage — MMR-related',
    'ID3':   'Tobacco smoking',
    'ID4':   'Unknown aetiology',
    'ID5':   'Unknown aetiology',
    'ID6':   'Homologous recombination deficiency',
    'ID7':   'Replication slippage at microsatellites',
    'ID8':   'Non-homologous end joining',
    'ID13':  'UV light exposure',
}


def build_mutational_signatures_section(args, story):
    """Mutational Signatures section: detected-signature summary table + spectrum plots."""
    story.append(PageBreak())
    story.append(Paragraph(T('sect_mutsig'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))

    sig_sbs96 = getattr(args, 'sig_sbs96', None) or ''
    sig_id    = getattr(args, 'sig_id', None) or ''

    def _read_activities(act_path):
        """Return list of (sig_name, activity_int) for detected sigs (activity > 0)."""
        if not act_path or not os.path.isfile(act_path):
            return []
        _, rows = read_tsv(act_path)
        if not rows:
            return []
        row = rows[0]  # single sample per run
        detected = []
        for col, val in row.items():
            if not col or col.lower() in ('samples', 'sample', ''):
                continue
            try:
                v = float(val)
            except (ValueError, TypeError):
                continue
            if v > 0:
                detected.append((col, int(round(v))))
        return sorted(detected, key=lambda x: -x[1])

    sbs_activities = _read_activities(sig_sbs96)
    id_activities  = _read_activities(sig_id)

    if sbs_activities or id_activities:
        col_headers = ['Signature', 'Type', 'Mutations Attributed', 'COSMIC Aetiology']
        col_widths  = [USABLE_W * w for w in [0.13, 0.07, 0.20, 0.60]]
        rows_data   = []
        for sig, act in sbs_activities:
            rows_data.append([sig, 'SBS', str(act),
                               _SIG_AETIOLOGY.get(sig, 'See COSMIC mutational signatures database')])
        for sig, act in id_activities:
            rows_data.append([sig, 'ID', str(act),
                               _SIG_AETIOLOGY.get(sig, 'See COSMIC mutational signatures database')])
        tbl = build_variant_table(col_headers, col_widths, rows_data)
        story.append(tbl)
        story.append(Spacer(1, 0.3 * cm))
        n_total = sum(a for _, a in sbs_activities) + sum(a for _, a in id_activities)
        story.append(Paragraph(
            f'<font size="7" color="#64748b">'
            f'{len(sbs_activities) + len(id_activities)} signature(s) detected '
            f'({len(sbs_activities)} SBS, {len(id_activities)} ID) — '
            f'{n_total} total mutations attributed. '
            f'COSMIC v3.4 reference signatures.</font>',
            styles['Normal'],
        ))
    else:
        story.append(Paragraph(
            'No mutational signatures detected above threshold.', STYLE_WARN))

    story.append(Spacer(1, 0.4 * cm))
    if sig_sbs96 and os.path.isfile(sig_sbs96):
        build_signature_plots_pdf(sig_sbs96, 'SBS96', story)
    if sig_id and os.path.isfile(sig_id):
        build_signature_plots_pdf(sig_id, 'ID', story)


def _load_cnv_panel(path):
    genes = set()
    if not path or not os.path.isfile(path):
        return genes
    with open(path) as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith('#'):
                genes.add(g.upper())
    return genes


# ---------------------------------------------------------------------------
# LOH filter helpers
# ---------------------------------------------------------------------------
def _germline_plp_genes(path):
    """Return set of gene names (uppercase) with ACMG P or LP classification."""
    genes = set()
    if not path or not os.path.isfile(path):
        return genes
    _, rows = read_tsv(path)
    for row in rows:
        cls = classify_acmg(row)
        if cls in ('P', 'LP'):
            gene = (row.get('ANN[*].GENE') or '').strip().upper()
            if gene:
                genes.add(gene)
    return genes


def _somatic_variant_genes(path):
    """Return set of gene names (uppercase) present in somatic variant TSV."""
    genes = set()
    if not path or not os.path.isfile(path):
        return genes
    _, rows = read_tsv(path)
    for row in rows:
        gene = (row.get('ANN[*].GENE') or '').strip().upper()
        if gene:
            genes.add(gene)
    return genes


def _harmonize_cnv_pdf(ascat_cnvkit_path, gatk_somatic_path,
                       somatic_panel_path, germline_panel_path, amp_threshold=5):
    # Returns dict: {'homdel': [...], 'amp': [...], 'loh': [...]}
    import csv as _csv
    somatic_panel  = _load_cnv_panel(somatic_panel_path)
    germline_panel = _load_cnv_panel(germline_panel_path)

    a_homdel, a_amp, a_loh = {}, {}, {}
    if ascat_cnvkit_path and os.path.isfile(ascat_cnvkit_path):
        with open(ascat_cnvkit_path) as fh:
            for row in _csv.DictReader(fh, delimiter='\t'):
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
                    a_homdel[gene] = {'cn': 0, 'chrom': chrom, 'log2': log2}
                elif etype == 'AMP' and gene in somatic_panel and cn >= amp_threshold:
                    if gene not in a_amp or cn > a_amp[gene]['cn']:
                        a_amp[gene] = {'cn': cn, 'chrom': chrom, 'log2': log2}
                elif etype == '2nd_hit_LOH' and gene in germline_panel:
                    a_loh[gene] = {'cn': cn, 'loh_type': loh_type, 'chrom': chrom}

    g_homdel, g_amp, g_loh = {}, {}, {}
    if gatk_somatic_path and os.path.isfile(gatk_somatic_path):
        with open(gatk_somatic_path) as fh:
            for row in _csv.DictReader(fh, delimiter='\t'):
                gene  = row.get('GENE', '').strip().upper()
                etype = row.get('EVENT_TYPE', '').strip()
                try:
                    cn = int(float(row.get('CN_ABSOLUTE', -1)))
                except (ValueError, TypeError):
                    cn = -1
                chrom = row.get('CONTIG', '')
                if etype == 'DELETION' and cn == 0 and gene in somatic_panel:
                    g_homdel[gene] = {'cn': 0, 'chrom': chrom}
                elif etype == 'AMPLIFICATION' and cn >= amp_threshold and gene in somatic_panel:
                    if gene not in g_amp or cn > g_amp[gene]['cn']:
                        g_amp[gene] = {'cn': cn, 'chrom': chrom}
                elif etype in ('CN-LOH', 'DEL+LOH') and gene in germline_panel:
                    g_loh[gene] = {
                        'cn': cn,
                        'loh_type': 'cnLOH' if etype == 'CN-LOH' else 'HETDEL+LOH',
                        'chrom': chrom,
                    }

    homdel_list, amp_list, loh_list = [], [], []
    for gene in sorted(set(a_homdel) | set(g_homdel)):
        a, g = a_homdel.get(gene), g_homdel.get(gene)
        tools = '+'.join(t for t in ['CNVkit+ASCAT' if a else '', 'GATK' if g else ''] if t)
        homdel_list.append({'gene': gene, 'cn_p': 0,
                            'cn_s': g['cn'] if g else '',
                            'chrom': (a or g)['chrom'],
                            'tools': tools, 'conf': 'HIGH' if (a and g) else 'MODERATE'})
    for gene in sorted(set(a_amp) | set(g_amp)):
        a, g = a_amp.get(gene), g_amp.get(gene)
        tools = '+'.join(t for t in ['CNVkit+ASCAT' if a else '', 'GATK' if g else ''] if t)
        amp_list.append({'gene': gene,
                         'cn_p': a['cn'] if a else (g['cn'] if g else ''),
                         'cn_s': g['cn'] if g else '',
                         'chrom': (a or g)['chrom'],
                         'tools': tools, 'conf': 'HIGH' if (a and g) else 'MODERATE'})
    for gene in sorted(set(a_loh) | set(g_loh)):
        a, g = a_loh.get(gene), g_loh.get(gene)
        tools = '+'.join(t for t in ['ASCAT' if a else '', 'GATK' if g else ''] if t)
        loh_type = a['loh_type'] if a else (g['loh_type'] if g else 'LOH')
        loh_list.append({'gene': gene, 'loh_type': loh_type,
                         'cn_p': a['cn'] if a else '',
                         'cn_s': g['cn'] if g else '',
                         'chrom': (a or g)['chrom'],
                         'tools': tools, 'conf': 'HIGH' if (a and g) else 'MODERATE'})
    _loh_order = {'HETDEL+LOH': 0, 'cnLOH': 1}
    loh_list.sort(key=lambda r: (_loh_order.get(r['loh_type'], 2), r['gene']))

    return {'homdel': homdel_list, 'amp': amp_list, 'loh': loh_list}


def build_somatic_cnv_section(ascat_cnvkit_path, gatk_somatic_path,
                               somatic_panel_path, germline_panel_path, story,
                               amp_threshold=5,
                               germline_snv_path='', somatic_snv_path='',
                               is_low_cov=False):
    # Unified somatic CNV section: HOMDEL, AMP, LOH
    story.append(PageBreak())
    story.append(Paragraph(T('sect_som_cnv'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))

    if is_low_cov:
        _indet_s = ParagraphStyle('CnvIndet', parent=styles['Normal'],
                                   fontSize=8, textColor=colors.HexColor('#7c2d12'), leading=12)
        _indet_data = [[Paragraph(
            f'<b>⚠ {T("indeterminate")}</b><br/>{T("low_cov_indet")}',
            _indet_s)]]
        _indet_tbl = Table(_indet_data, colWidths=[USABLE_W])
        _indet_tbl.setStyle(TableStyle([
            ('BACKGROUND',   (0, 0), (-1, -1), colors.HexColor('#fef3c7')),
            ('TOPPADDING',   (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING',(0, 0), (-1, -1), 8),
            ('LEFTPADDING',  (0, 0), (-1, -1), 12),
            ('RIGHTPADDING', (0, 0), (-1, -1), 12),
            ('ROUNDEDCORNERS', [4]),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#f59e0b')),
        ]))
        story.append(_indet_tbl)
        return

    data = _harmonize_cnv_pdf(ascat_cnvkit_path, gatk_somatic_path,
                               somatic_panel_path, germline_panel_path, amp_threshold)

    # ── Filter LOH: only show if gene has germline P/LP or somatic variant ────
    if germline_snv_path or somatic_snv_path:
        plp_genes = _germline_plp_genes(germline_snv_path or '')
        som_genes = _somatic_variant_genes(somatic_snv_path or '')
        loh_before = len(data['loh'])
        data['loh'] = [r for r in data['loh']
                       if r['gene'] in plp_genes or r['gene'] in som_genes]
        if loh_before != len(data['loh']):
            print(f'[INFO] LOH filter: {loh_before} → {len(data["loh"])} events '
                  f'(kept genes with germline P/LP or somatic variant)', file=sys.stderr)

    note_s = ParagraphStyle('CnvNote', parent=styles['Normal'],
                             fontSize=6, textColor=colors.Color(0.6, 0.7, 0.8))

    C_HOMDEL = colors.HexColor('#ef4444')
    C_AMP    = colors.HexColor('#f97316')
    C_LOHH   = colors.HexColor('#dc2626')   # HETDEL+LOH
    C_LOHC   = colors.HexColor('#ca8a04')   # cnLOH
    C_LOHO   = colors.HexColor('#7c3aed')   # other LOH
    C_HIGH   = colors.HexColor('#16a34a')
    C_MOD    = colors.HexColor('#6b7280')

    def conf_color(c):
        return C_HIGH if c == 'HIGH' else C_MOD

    # ── Subsection: HOMDEL ────────────────────────────────────────────────────
    story.append(Paragraph('Homozygous Deletions (CN = 0)', STYLE_H3))
    if data['homdel']:
        hdrs = ['Gene', 'Chr', 'CN (ASCAT)', 'CN (GATK)', 'Tools', 'Confidence']
        wids = [USABLE_W * w for w in [0.13, 0.10, 0.12, 0.12, 0.36, 0.17]]
        rows_d = []
        for r in data['homdel']:
            cn_s_str = str(r['cn_s']) if r['cn_s'] != '' else '\u2014'
            conf_c   = '#16a34a' if r['conf'] == 'HIGH' else '#6b7280'
            rows_d.append([
                f"<b>{r['gene']}</b>",
                r['chrom'],
                '<font color="#ef4444"><b>0</b></font>',
                cn_s_str,
                r['tools'],
                f'<font color="{conf_c}">{r["conf"]}</font>',
            ])
        story.append(build_variant_table(hdrs, wids, rows_d,
                                          colour_col_idx=5,
                                          colour_fn=lambda v: C_HIGH if 'HIGH' in str(v) else C_MOD))
    else:
        story.append(Paragraph(
            'No homozygous deletions detected in the somatic panel.', note_s))

    story.append(Spacer(1, 8))

    # ── Subsection: AMP ───────────────────────────────────────────────────────
    story.append(Paragraph(f'Amplifications (CN \u2265 {amp_threshold})', STYLE_H3))
    if data['amp']:
        hdrs = ['Gene', 'Chr', 'CN (ASCAT)', 'CN (GATK)', 'Tools', 'Confidence']
        wids = [USABLE_W * w for w in [0.13, 0.10, 0.12, 0.12, 0.36, 0.17]]
        amp_sorted = sorted(data['amp'],
                            key=lambda r: -(r['cn_p'] if isinstance(r['cn_p'], int) else 0))
        rows_a = []
        for r in amp_sorted:
            cn_s_str = str(r['cn_s']) if r['cn_s'] != '' else '\u2014'
            conf_c   = '#16a34a' if r['conf'] == 'HIGH' else '#6b7280'
            rows_a.append([
                f"<b>{r['gene']}</b>",
                r['chrom'],
                f'<font color="#f97316"><b>{r["cn_p"]}</b></font>',
                cn_s_str,
                r['tools'],
                f'<font color="{conf_c}">{r["conf"]}</font>',
            ])
        story.append(build_variant_table(hdrs, wids, rows_a,
                                          colour_col_idx=5,
                                          colour_fn=lambda v: C_HIGH if 'HIGH' in str(v) else C_MOD))
    else:
        story.append(Paragraph(
            f'No amplifications (CN \u2265 {amp_threshold}) detected in the somatic panel.', note_s))

    story.append(Spacer(1, 8))

    # ── Subsection: LOH ───────────────────────────────────────────────────────
    story.append(Paragraph('Loss of Heterozygosity (LOH) \u2014 Germline Panel', STYLE_H3))
    if data['loh']:
        hdrs = ['Gene', 'Chr', 'LOH Type', 'CN (ASCAT)', 'CN (GATK)', 'Tools', 'Confidence']
        wids = [USABLE_W * w for w in [0.11, 0.08, 0.14, 0.09, 0.09, 0.30, 0.19]]
        rows_l = []
        loh_colors = {'HETDEL+LOH': '#dc2626', 'cnLOH': '#ca8a04'}
        for r in data['loh']:
            cn_p_str = str(r['cn_p']) if r['cn_p'] != '' else '\u2014'
            cn_s_str = str(r['cn_s']) if r['cn_s'] != '' else '\u2014'
            loh_c    = loh_colors.get(r['loh_type'], '#7c3aed')
            conf_c   = '#16a34a' if r['conf'] == 'HIGH' else '#6b7280'
            rows_l.append([
                f"<b>{r['gene']}</b>",
                r['chrom'],
                f'<font color="{loh_c}"><b>{r["loh_type"]}</b></font>',
                cn_p_str,
                cn_s_str,
                r['tools'],
                f'<font color="{conf_c}">{r["conf"]}</font>',
            ])
        story.append(build_variant_table(hdrs, wids, rows_l,
                                          colour_col_idx=6,
                                          colour_fn=lambda v: C_HIGH if 'HIGH' in str(v) else C_MOD))
    else:
        story.append(Paragraph(
            'No LOH events detected in the germline panel.', note_s))

    n_total = len(data['homdel']) + len(data['amp']) + len(data['loh'])
    story.append(Paragraph(
        f'{n_total} reportable CNV event(s): '
        f'{len(data["homdel"])} HOMDEL \u2022 '
        f'{len(data["amp"])} AMP \u2022 '
        f'{len(data["loh"])} LOH',
        note_s))


def build_germline_cnv_section(path, story, germline_panel_path=''):
    """Append GATK germline CNV section: show segments with CN != 2.

    Only panel genes (103-gene cancer-predisposition panel) are reported.
    Filter: QA >= 4 (no NP filter — small CNVs of 1-3 exons are clinically critical).
    """
    story.append(PageBreak())
    story.append(Paragraph(T('sect_germ_cnv'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))
    MIN_QA = 4

    # Load panel genes
    panel_genes = set()
    if germline_panel_path and os.path.isfile(germline_panel_path):
        with open(germline_panel_path) as _pf:
            for _line in _pf:
                g = _line.strip().split('\t')[0].strip()
                if g and not g.startswith('#'):
                    panel_genes.add(g)

    if not path or not os.path.isfile(path):
        story.append(Paragraph('Germline CNV (GATK) file not found.', STYLE_WARN))
        return

    import csv as _csv
    with open(path) as fh:
        rows = list(_csv.DictReader(fh, delimiter='\t'))

    # Group by segment (CHROM, POS, END), filter CN != 2
    segments = {}
    seg_order = []
    for row in rows:
        try:
            cn_int = int(row.get('CN', '2').strip())
        except ValueError:
            cn_int = 2
        if cn_int == 2:
            continue
        key = (row.get('CHROM', ''), row.get('POS', ''), row.get('END', ''))
        if key not in segments:
            segments[key] = {
                'chrom': row.get('CHROM', ''),
                'pos':   row.get('POS', ''),
                'end':   row.get('END', ''),
                'alt':   row.get('ALT', '.'),
                'cn':    row.get('CN', '.'),
                'np':    row.get('NP', '.'),
                'qa':    row.get('QA', '0'),
                'qs':    row.get('QS', '.'),
                'genes': [],
            }
            seg_order.append(key)
        gene = row.get('gen', '').strip()
        if gene and gene not in segments[key]['genes']:
            segments[key]['genes'].append(gene)

    # Filter: only panel genes with QA >= MIN_QA
    filtered = []
    for key in seg_order:
        seg = segments[key]
        try:
            qa_val = int(seg.get('qa', '0') or '0')
        except (ValueError, TypeError):
            continue
        if qa_val < MIN_QA:
            continue
        # Keep only panel genes in the segment
        panel_hits = [g for g in seg['genes'] if g in panel_genes]
        if panel_hits:
            seg['genes'] = panel_hits
            filtered.append(seg)

    if not filtered:
        # Styled callout box
        callout_style = ParagraphStyle(
            'CnvCallout',
            parent=styles['Normal'],
            fontSize=8,
            textColor=colors.HexColor('#1e3a5f'),
            leading=12,
        )
        callout_data = [[
            Paragraph(
                '<b>Germline CNV Analysis \u2014 No Reportable Events</b><br/>'
                'The GATK germline copy-number variant analysis did not identify '
                'any copy-number alterations (CN \u2260 2) meeting the quality '
                'threshold (QA \u2265 4) in the 103-gene germline '
                'cancer predisposition panel. The raw annotated segments are '
                'available for manual review in:<br/>'
                '\u00a0\u00a0<i>clinical_results/germline/cnv_gatk/annotated/</i>',
                callout_style,
            )
        ]]
        callout_tbl = Table(callout_data, colWidths=[USABLE_W])
        callout_tbl.setStyle(TableStyle([
            ('BACKGROUND',    (0, 0), (-1, -1), colors.HexColor('#dbeafe')),
            ('LEFTPADDING',   (0, 0), (-1, -1), 10),
            ('RIGHTPADDING',  (0, 0), (-1, -1), 10),
            ('TOPPADDING',    (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('LINEBEFORE',    (0, 0), (0, -1), 4, colors.HexColor('#3b82f6')),
            ('ROUNDEDCORNERS', [4]),
        ]))
        story.append(callout_tbl)
        return

    def cn_type(alt, cn_str):
        try:
            cn = int(cn_str)
        except (ValueError, TypeError):
            cn = -1
        if cn == 0:   return 'HOMDEL'
        elif cn == 1: return 'DEL (CN=1)'
        elif cn == 3: return 'DUP (CN=3)'
        elif cn >= 4: return f'AMP (CN={cn})'
        return alt or 'CNV'

    col_headers = ['Locus', 'Gene(s)', 'Type', 'CN', 'NP', 'QA', 'QS', 'Quality']
    col_widths  = [USABLE_W * w for w in [0.20, 0.16, 0.12, 0.05, 0.06, 0.06, 0.06, 0.09]]

    data_rows = []
    for seg in filtered:
        locus = seg['chrom'] + ':' + seg['pos'] + '-' + seg['end']
        genes = ', '.join(seg['genes']) if seg['genes'] else '.'
        ctype = cn_type(seg['alt'], seg['cn'])
        try:
            qs_f = float(seg['qs'])
            qs_label = 'High' if qs_f >= 30 else ('Medium' if qs_f >= 10 else 'Low')
        except (ValueError, TypeError):
            qs_label = 'Unknown'
        data_rows.append([locus, genes, ctype, seg['cn'], seg['np'], seg['qa'], seg['qs'], qs_label])

    tbl = build_variant_table(
        col_headers, col_widths, data_rows,
        colour_col_idx=7,
        colour_fn=lambda val: (C_ACMG_P if val == 'High' else
                               (C_ACMG_LP if val == 'Medium' else
                                (C_ACMG_VUS if val == 'Low' else colors.grey))),
    )
    note_style = ParagraphStyle('CnvNote', parent=styles['Normal'],
                                fontSize=6, textColor=colors.Color(0.6, 0.7, 0.8))
    story.append(Paragraph(
        str(len(filtered)) + ' segment(s) shown (103-gene panel only, QA \u2265 4). Quality: '
        'High (QS \u226530), Medium (QS 10\u201329), Low (QS <10 \u2014 manual review recommended).',
        note_style))
    story.append(tbl)


def build_sequencing_quality_section(path_n, path_t, normal_name, tumor_name, story,
                                     germline_panel_path='', somatic_panel_path=''):
    # Per-gene mean coverage as justified grid: gene name + colored value, two subsections.
    story.append(PageBreak())
    story.append(Paragraph(T('sect_seqq'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=4))

    def load_panel(path):
        genes = []
        if not path or not os.path.isfile(path):
            return genes
        with open(path) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith('#'):
                    genes.append(g.split()[0])
        return genes

    def load_cov(path):
        data = {}
        if not path or not os.path.isfile(path):
            return data
        _, rows = read_tsv(path)
        for row in rows:
            gene = (row.get('gene') or '').strip()
            if gene:
                data[gene] = row
        return data

    def cov_color_rl(val_str):
        try:
            v = float(val_str)
        except (ValueError, TypeError):
            return colors.HexColor('#6b7280')
        if v >= 100:
            return colors.HexColor('#16a34a')
        elif v >= 50:
            return colors.HexColor('#ca8a04')
        elif v >= 20:
            return colors.HexColor('#ea580c')
        else:
            return colors.HexColor('#dc2626')

    def fmt_cov(v):
        try:
            return f'{float(v):.0f}'
        except (ValueError, TypeError):
            return '-'

    # Grid renderer: takes a list of (gene, mean_cov_str) tuples
    NCOLS = 9
    CELL_W = USABLE_W / NCOLS  # ~57.3pt

    gene_name_style = ParagraphStyle(
        'gcov_gene', fontName='Helvetica-Bold', fontSize=6.5,
        textColor=C_TEXT, alignment=TA_CENTER, leading=8,
    )
    cov_val_base = ParagraphStyle(
        'gcov_val', fontName='Helvetica-Bold', fontSize=9,
        alignment=TA_CENTER, leading=11, spaceAfter=0,
    )

    def render_grid(gene_cov_list, cov_file_name):
        if not gene_cov_list:
            story.append(Paragraph(
                f'<font size="7" color="#64748b">No coverage data available ({cov_file_name}).</font>',
                styles['Normal'],
            ))
            return
        # Pad to full rows
        padded = list(gene_cov_list)
        while len(padded) % NCOLS != 0:
            padded.append(('', ''))

        rows = []
        for r in range(len(padded) // NCOLS):
            row_cells = []
            for c in range(NCOLS):
                gene, cov_str = padded[r * NCOLS + c]
                if not gene:
                    row_cells.append('')
                    continue
                # ReportLab doesn't support ParagraphStyle per-cell with color, so build
                # two paragraphs stacked inside the cell using a nested Table.
                cov_c = cov_color_rl(cov_str)
                gene_p = Paragraph(gene, gene_name_style)
                val_p  = Paragraph(
                    f'<font color="{cov_c.hexval()}">{fmt_cov(cov_str)}</font>',
                    cov_val_base,
                )
                inner = Table(
                    [[gene_p], [val_p]],
                    colWidths=[CELL_W - 4],
                )
                inner.setStyle(TableStyle([
                    ('TOPPADDING',    (0, 0), (-1, -1), 1),
                    ('BOTTOMPADDING',(0, 0), (-1, -1), 1),
                    ('LEFTPADDING',  (0, 0), (-1, -1), 0),
                    ('RIGHTPADDING', (0, 0), (-1, -1), 0),
                    ('ALIGN',        (0, 0), (-1, -1), 'CENTER'),
                ]))
                row_cells.append(inner)
            rows.append(row_cells)

        tbl = Table(rows, colWidths=[CELL_W] * NCOLS)
        tbl.setStyle(TableStyle([
            ('GRID',          (0, 0), (-1, -1), 0.3, C_LIGHT_BORDER),
            ('ROWBACKGROUNDS', (0, 0), (-1, -1), [C_WHITE, C_ROW_ALT]),
            ('TOPPADDING',    (0, 0), (-1, -1), 2),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('LEFTPADDING',   (0, 0), (-1, -1), 2),
            ('RIGHTPADDING',  (0, 0), (-1, -1), 2),
            ('VALIGN',        (0, 0), (-1, -1), 'MIDDLE'),
        ]))
        story.append(tbl)

    legend_style = ParagraphStyle(
        'cov_leg', fontName='Helvetica', fontSize=7, textColor=C_TEXT, spaceAfter=6,
    )
    story.append(Paragraph(
        '<font color="#16a34a"><b>\u25a0</b></font> \u2265100X &nbsp;&nbsp;'
        '<font color="#ca8a04"><b>\u25a0</b></font> 50\u2013100X &nbsp;&nbsp;'
        '<font color="#ea580c"><b>\u25a0</b></font> 20\u201350X &nbsp;&nbsp;'
        '<font color="#dc2626"><b>\u25a0</b></font> &lt;20X',
        legend_style,
    ))

    # ── Germline panel ──────────────────────────────────────────────────────
    cov_n = load_cov(path_n)
    germline_genes = load_panel(germline_panel_path)
    if not germline_genes:
        germline_genes = sorted(cov_n.keys())

    subhdr = ParagraphStyle(
        'cov_subhdr', fontName='Helvetica-Bold', fontSize=8.5,
        textColor=C_ACCENT, spaceBefore=6, spaceAfter=4,
    )
    story.append(Paragraph(
        f'Germline Panel \u2014 {normal_name} ({len([g for g in germline_genes if g in cov_n])} genes)',
        subhdr,
    ))
    germ_pairs = [(g, cov_n[g].get('mean_coverage', '-')) for g in germline_genes if g in cov_n]
    render_grid(germ_pairs, path_n or 'normal')

    story.append(Spacer(1, 0.3 * cm))

    # ── Somatic panel ───────────────────────────────────────────────────────
    cov_t = load_cov(path_t)
    somatic_genes = load_panel(somatic_panel_path)
    if not somatic_genes:
        somatic_genes = sorted(cov_t.keys())

    story.append(Paragraph(
        f'Somatic Panel \u2014 {tumor_name} ({len([g for g in somatic_genes if g in cov_t])} genes)',
        subhdr,
    ))
    som_pairs = [(g, cov_t[g].get('mean_coverage', '-')) for g in somatic_genes if g in cov_t]
    render_grid(som_pairs, path_t or 'tumor')

    story.append(Spacer(1, 0.2 * cm))
    story.append(Paragraph(
        f'<font size="6.5" color="#64748b">Mean coverage (X) per gene over panel capture regions. '
        f'Thresholds: \u226520X adequate \u2022 \u226550X good \u2022 \u2265100X excellent.</font>',
        styles['Normal'],
    ))


def build_correlation_section(corr_data, story):
    story.append(PageBreak())
    story.append(Paragraph(T('sect_correlation'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=6))

    results           = (corr_data or {}).get('results', [])
    biallelic_somatic = (corr_data or {}).get('biallelic_somatic', [])

    if not corr_data or (not results and not biallelic_somatic):
        story.append(Paragraph(
            'No germline P/LP/VUS variants or somatic biallelic events available '
            'for correlation analysis.',
            STYLE_WARN))
        return

    metrics = corr_data.get('metrics', {})
    detected_sigs = sorted(corr_data.get('detected_sigs', []))

    # Biomarker summary line
    tmb_s  = str(round(float(metrics.get('tmb', 0) or 0), 1))
    msi_s  = str(round(float(metrics.get('msi', 0) or 0), 1))
    hrd_s  = str(round(float(metrics.get('hrd_prob', 0) or 0), 4))
    mmr_found = [s for s in detected_sigs if s in ('SBS6','SBS15','SBS20','SBS21','SBS26','SBS44')]
    hrd_found = [s for s in detected_sigs if s in ('SBS3','SBS8')]
    try:
        _msi_corr_status = 'detected' if float(msi_s or 0) >= 5.0 else 'MSS'
    except (ValueError, TypeError):
        _msi_corr_status = ''
    bio_line = (
        'Global biomarkers:  TMB=' + tmb_s + ' mut/Mb'
        + '  |  MSI=' + msi_s + '% (' + _msi_corr_status + ')'
        + '  |  HRD prob=' + hrd_s
        + ('  |  MMR signatures: ' + ', '.join(mmr_found) if mmr_found else '  |  MMR signatures: none')
        + ('  |  HRD signatures: ' + ', '.join(hrd_found) if hrd_found else '  |  HRD signatures: none')
    )
    story.append(Paragraph('<font size="7" color="#94a3b8">' + bio_line + '</font>',
                           styles['Normal']))
    story.append(Spacer(1, 0.25 * cm))

    # Column headers
    col_headers = ['Gene', 'Germline variant', 'ACMG', 'Pathway',
                   'Somatic 2nd hit', 'Biomarkers', 'Interpretation']
    col_widths   = [1.3*cm, 3.4*cm, 1.0*cm, 1.5*cm, 3.8*cm, 3.8*cm, 3.2*cm]

    header_style = ParagraphStyle('ch', fontSize=6.5, fontName='Helvetica-Bold',
                                  textColor=colors.white, leading=8)
    cell_style   = ParagraphStyle('cc', fontSize=6.5, fontName='Helvetica',
                                  textColor=colors.HexColor('#1e293b'), leading=8)
    small_style  = ParagraphStyle('cs', fontSize=5.5, fontName='Helvetica',
                                  textColor=colors.HexColor('#64748b'), leading=7)

    def trunc(s, n=38):
        s = str(s or '')
        return s if len(s) <= n else s[:n-1] + '\u2026'

    table_data = [[Paragraph(h, header_style) for h in col_headers]]

    for r in results:
        # Gene
        gene_p = Paragraph('<b>' + r['gene'] + '</b>', cell_style)

        # Variant — different display for CNV vs SNV
        vtype = r.get('variant_type', 'SNV')
        if vtype == 'CNV':
            var_p = Paragraph(
                trunc(r.get('hgvs_c', '?'), 35)
                + '<br/><font size="5" color="#7c3aed"><b>[CNV germline hit]</b></font>',
                cell_style)
        else:
            hp = r.get('hgvs_p') or r.get('hgvs_c') or '?'
            hc = r.get('hgvs_c', '') if r.get('hgvs_p') else ''
            var_p = Paragraph(trunc(hp, 28)
                              + ('<br/><font size="5" color="#94a3b8">' + trunc(hc, 28) + '</font>'
                                 if hc else ''), cell_style)

        # ACMG badge
        acmg_colors_map = {'P': '#ef4444', 'LP': '#f97316', 'VUS': '#eab308',
                           'LB': '#3b82f6', 'B': '#22c55e'}
        ac = r.get('acmg', 'VUS')
        ac_col = acmg_colors_map.get(ac, '#6b7280')
        acmg_p = Paragraph('<font color="' + ac_col + '"><b>' + ac + '</b></font>', cell_style)

        # Pathway
        pw = r.get('pathway', 'Other TSG')
        pw_col_hex = {
            'MMR': '#6366f1', 'HRR': '#0ea5e9', 'Other TSG': '#8b5cf6'
        }.get(pw, '#8b5cf6')
        pw_p = Paragraph('<font color="' + pw_col_hex + '"><b>' + pw + '</b></font>', cell_style)

        # Second hits
        hits = r.get('second_hits', [])
        if hits:
            hit_lines = '<br/>'.join(
                '\u25a0 ' + trunc(h.get('label', ''), 36) for h in hits)
        else:
            hit_lines = '<font color="#94a3b8"><i>Not detected</i></font>'
        hits_p = Paragraph(hit_lines, cell_style)

        # Biomarkers
        bm_items = r.get('biomarkers', [])
        status_sym = {'HIGH': '\u2713', 'borderline': '\u25cf', 'low': '\u2717',
                      'absent': '\u2717', 'detected': '\u2713', 'n/a': '\u2014'}
        bm_lines = '<br/>'.join(
            status_sym.get(b.get('status', ''), '\u25cf') + ' '
            + trunc(b.get('key', '') + ': ' + b.get('val', ''), 34)
            for b in bm_items
        )
        bm_p = Paragraph(bm_lines, cell_style)

        # Interpretation
        ik  = r.get('interp_key', 'GERMLINE_ONLY')
        il  = r.get('interp_label', 'Germline only')
        ic  = INTERP_COLORS.get(ik, colors.HexColor('#6b7280'))
        ic_hex = '#{:02x}{:02x}{:02x}'.format(
            int(ic.red * 255), int(ic.green * 255), int(ic.blue * 255))
        interp_p = Paragraph('<font color="' + ic_hex + '"><b>' + trunc(il, 22) + '</b></font>',
                             cell_style)

        table_data.append([gene_p, var_p, acmg_p, pw_p, hits_p, bm_p, interp_p])

    tbl = Table(table_data, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(TableStyle([
        ('BACKGROUND',   (0, 0), (-1, 0),  C_SLATE),
        ('TEXTCOLOR',    (0, 0), (-1, 0),  colors.white),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, C_ROW_ALT]),
        ('GRID',         (0, 0), (-1, -1), 0.4, colors.HexColor('#e2e8f0')),
        ('VALIGN',       (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING',  (0, 0), (-1, -1), 4),
        ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING',   (0, 0), (-1, -1), 3),
        ('BOTTOMPADDING',(0, 0), (-1, -1), 3),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 0.3 * cm))
    story.append(Paragraph(
        '<font size="6" color="#94a3b8">Methodology: ACMG/AMP 2015 classification using ClinVar, '
        'gnomAD, REVEL, SIFT, PolyPhen2, and SnpEff impact. Second hit detection from somatic SNV '
        '(Mutect2, HIGH/MODERATE) and ASCAT LOH events (cnLOH, HETDEL+LOH, HOMDEL). '
        'Pathway biomarkers: MMR=TMB\u226510 mut/Mb or MSI\u22653.5% + SBS6/15/20/21/26/44; '
        'HRR=HRD prob\u22650.5 + SBS3/SBS8. For decision support only &#8212; '
        'must be reviewed by a certified clinical geneticist.</font>',
        styles['Normal']))

    # --- Somatic Biallelic Events sub-section ---
    if biallelic_somatic:
        story.append(Spacer(1, 0.4 * cm))
        story.append(Paragraph(
            '\u25c6 Somatic Biallelic Events (two somatic hits, no germline variant)',
            ParagraphStyle('sbh', fontSize=8, fontName='Helvetica-Bold',
                           textColor=colors.HexColor('#7c3aed'), leading=10)))
        story.append(Paragraph(
            '<font size="6" color="#94a3b8">'
            'Genes with a somatic coding variant (HIGH/MODERATE) and a somatic copy-number loss '
            '(HETDEL, HOMDEL, or LOH) — independent of any germline variant. '
            'These are candidate somatic biallelic inactivation events.</font>',
            styles['Normal']))
        story.append(Spacer(1, 0.15 * cm))

        sb_header_style = ParagraphStyle('sbhdr', fontSize=6.5, fontName='Helvetica-Bold',
                                         textColor=colors.white, leading=8)
        sb_cell_style   = ParagraphStyle('sbcell', fontSize=6.5, fontName='Helvetica',
                                         textColor=colors.HexColor('#1e293b'), leading=8)

        sb_col_headers = ['Gene', 'Pathway', 'Somatic SNV', 'Somatic CNV/LOH', 'Mechanism']
        sb_col_widths  = [1.2*cm, 1.4*cm, 3.2*cm, 3.2*cm, 9.0*cm]

        sb_data = [[Paragraph(h, sb_header_style) for h in sb_col_headers]]
        pw_colors_pdf = {'MMR': '#6366f1', 'HRR': '#0ea5e9', 'Other TSG': '#8b5cf6'}
        def _trunc_sb(s, n=55):
            s = str(s or '')
            return s if len(s) <= n else s[:n-1] + '\u2026'

        for b in biallelic_somatic:
            pw_hex = pw_colors_pdf.get(b['pathway'], '#8b5cf6')
            sb_data.append([
                Paragraph('<b>' + b['gene'] + '</b>', sb_cell_style),
                Paragraph('<font color="' + pw_hex + '"><b>' + b['pathway'] + '</b></font>',
                          sb_cell_style),
                Paragraph('\u25a0 ' + _trunc_sb(b['snv_label'], 38), sb_cell_style),
                Paragraph('\u25b3 ' + _trunc_sb(b['cnv_label'], 38), sb_cell_style),
                Paragraph(_trunc_sb(b['mechanism'], 110), sb_cell_style),
            ])

        sb_tbl = Table(sb_data, colWidths=sb_col_widths, repeatRows=1)
        sb_tbl.setStyle(TableStyle([
            ('BACKGROUND',    (0, 0), (-1, 0),  colors.HexColor('#5b21b6')),
            ('TEXTCOLOR',     (0, 0), (-1, 0),  colors.white),
            ('ROWBACKGROUNDS',(0, 1), (-1, -1), [colors.HexColor('#faf5ff'),
                                                  colors.HexColor('#f3e8ff')]),
            ('GRID',          (0, 0), (-1, -1), 0.4, colors.HexColor('#e9d5ff')),
            ('VALIGN',        (0, 0), (-1, -1), 'TOP'),
            ('LEFTPADDING',   (0, 0), (-1, -1), 4),
            ('RIGHTPADDING',  (0, 0), (-1, -1), 4),
            ('TOPPADDING',    (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ]))
        story.append(sb_tbl)


# ---------------------------------------------------------------------------
# Major Findings section
# ---------------------------------------------------------------------------
_HRR_GENES = {
    'BRCA1', 'BRCA2', 'RAD51C', 'RAD51D', 'PALB2', 'BRIP1',
    'ATM', 'CHEK2', 'NBN', 'RAD50', 'MRE11', 'BARD1', 'CDK12',
    'FANCD2', 'FANCA', 'FANCM', 'SLX4', 'RBBP8', 'RTEL1',
}
_MMR_GENES = {
    'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM', 'MSH3', 'PMS1', 'MLH3',
}


def build_major_findings_section(args, metrics, germline_variants, somatic_variants,
                                  corr_data, story):
    """Append Major Clinical Findings section."""
    story.append(PageBreak())
    story.append(Paragraph(T('sect_major_finds'), STYLE_H2))
    story.append(HRFlowable(width=USABLE_W, thickness=1, color=C_ACCENT, spaceAfter=6))

    def _sf(v, d=None):
        try:
            return float(str(v).strip())
        except (ValueError, TypeError):
            return d

    # Detected signatures from correlation JSON
    detected_sigs    = set(corr_data.get('detected_sigs', [])) if corr_data else set()
    detected_id_sigs = set(corr_data.get('detected_id_sigs', [])) if corr_data else set()
    mmr_sbs_sigs = [s for s in ('SBS6', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44')
                    if s in detected_sigs]
    mmr_id_sigs  = [s for s in ('ID1', 'ID2') if s in detected_id_sigs]
    mmr_sigs     = mmr_sbs_sigs  # keep for legacy compatibility
    hrd_sigs     = [s for s in ('SBS3', 'SBS8') if s in detected_sigs]
    biallelic_somatic = (corr_data or {}).get('biallelic_somatic', [])

    # Germline P/LP lookup by gene
    germ_plp = {}
    for v in germline_variants:
        if v.get('_ACMG') in ('P', 'LP'):
            g = (v.get('ANN[*].GENE') or '').strip().upper()
            if g:
                germ_plp.setdefault(g, []).append(v)

    # Somatic variants by gene
    som_by_gene = {}
    for v in somatic_variants:
        g = (v.get('ANN[*].GENE') or '').strip().upper()
        if g:
            som_by_gene.setdefault(g, []).append(v)

    # CNV data for second hits
    cnv_data = _harmonize_cnv_pdf(
        getattr(args, 'cnv_somatic', '') or '',
        getattr(args, 'cnv_gatk_somatic', '') or '',
        getattr(args, 'somatic_cnv_panel', '') or '',
        getattr(args, 'germline_panel', '') or '',
    )
    loh_genes    = {r['gene'] for r in cnv_data['loh']}
    homdel_genes = {r['gene'] for r in cnv_data['homdel']}

    # Metric values
    hrd_prob   = _sf(metrics.get('hrd_prob', 'NA'))
    hrd_status = str(metrics.get('hrd_status', 'NA'))
    tmb_val    = _sf(metrics.get('tmb', 'NA'))
    msi_val    = _sf(metrics.get('msi', 'NA'))

    findings = []  # list of (title, colour, [line_markup_str, ...])

    # ── Finding 1: HRD / HRR ────────────────────────────────────────────────
    hrr_germ = [g for g in germ_plp if g in _HRR_GENES]
    hrr_som  = [g for g in som_by_gene if g in _HRR_GENES]
    hrd_pos  = hrd_prob is not None and hrd_prob >= 0.45
    if hrd_pos or hrr_germ:
        lines = []
        prob_str  = f'{hrd_prob:.3f}' if hrd_prob is not None else 'NA'
        stat_col  = '#ef4444' if hrd_pos else '#22c55e'
        lines.append(f'<b>{T("hrd_status_lbl")}</b> <font color="{stat_col}"><b>{hrd_status}</b></font>'
                     f'&nbsp;&nbsp;({T("prob_lbl")}&nbsp;=&nbsp;{prob_str})')
        for g in hrr_germ:
            for v in germ_plp[g]:
                hgvsp = v.get('ANN[*].HGVS_P', '') or v.get('ANN[*].HGVS_C', '') or ''
                ac    = v.get('_ACMG', '')
                ac_c  = '#ef4444' if ac == 'P' else '#f97316'
                lines.append(f'<b>{T("germ_hit_lbl")}</b> {g} '
                              f'<font color="{ac_c}"><b>[{ac}]</b></font> '
                              f'{hgvsp[:45] if hgvsp else ""}')
                hits = []
                if g in loh_genes:
                    hits.append(T('loh_det'))
                if g in homdel_genes:
                    hits.append(T('homdel_hit'))
                if g in som_by_gene:
                    for sv in som_by_gene[g][:2]:
                        hp2 = sv.get('ANN[*].HGVS_P', '') or sv.get('ANN[*].HGVS_C', '') or 'variant'
                        hits.append(f'{T("somatic_pfx")} {hp2[:30]}')
                if hits:
                    lines.append(f'<b>{T("2nd_hits_lbl")}</b> ' + ' | '.join(hits[:3]))
                else:
                    lines.append(f'<b>{T("2nd_hit_lbl")}</b> <font color="#6b7280">{T("not_detected")}</font>')
        if hrd_sigs:
            lines.append(f'<b>{T("hrd_sigs_lbl")}</b> {", ".join(hrd_sigs)}')
        else:
            lines.append(f'<b>{T("hrd_sigs_lbl")}</b> <font color="#6b7280">{T("sbs3_not_det")}</font>')
        findings.append((T('find_hrd'), colors.HexColor('#0ea5e9'), lines))

    # ── Finding 2: MMR Deficiency ────────────────────────────────────────────
    mmr_germ = [g for g in germ_plp if g in _MMR_GENES]
    mmr_som  = [g for g in som_by_gene if g in _MMR_GENES]
    msi_pos  = msi_val is not None and msi_val >= 3.5
    # MMR biallelic somatic events (e.g. MSH2 somatic SNV + cnLOH)
    mmr_biallelic_som = [b for b in biallelic_somatic
                         if b.get('pathway') == 'MMR']
    if mmr_germ or mmr_sbs_sigs or mmr_id_sigs or msi_pos or mmr_biallelic_som:
        lines = []
        msi_s  = f'{msi_val:.1f}%' if msi_val is not None else 'NA'
        tmb_s  = f'{tmb_val:.1f} mut/Mb' if tmb_val is not None else 'NA'
        msi_st = 'MSI-H' if (msi_val is not None and msi_val >= 5) else 'MSI-L/MSS'
        msi_c  = '#ef4444' if (msi_val is not None and msi_val >= 5) else '#6b7280'
        tmb_c  = '#ef4444' if (tmb_val is not None and tmb_val >= 10) else '#6b7280'
        lines.append(
            f'<b>MSI:</b> <font color="{msi_c}">{msi_s} ({msi_st})</font>'
            f'&nbsp;&nbsp;<b>TMB:</b> <font color="{tmb_c}">{tmb_s}</font>')
        for g in mmr_germ:
            for v in germ_plp[g]:
                hgvsp = v.get('ANN[*].HGVS_P', '') or v.get('ANN[*].HGVS_C', '') or ''
                ac    = v.get('_ACMG', '')
                ac_c  = '#ef4444' if ac == 'P' else '#f97316'
                lines.append(f'<b>{T("germ_hit_lbl")}</b> {g} '
                              f'<font color="{ac_c}"><b>[{ac}]</b></font> '
                              f'{hgvsp[:45] if hgvsp else ""}')
                hits = []
                if g in loh_genes:
                    hits.append(T('loh_det'))
                if g in som_by_gene:
                    for sv in som_by_gene[g][:2]:
                        hp2 = sv.get('ANN[*].HGVS_P', '') or sv.get('ANN[*].HGVS_C', '') or 'variant'
                        hits.append(f'{T("somatic_pfx")} {hp2[:30]}')
                if hits:
                    lines.append(f'<b>{T("2nd_hits_lbl")}</b> ' + ' | '.join(hits[:3]))
                else:
                    lines.append(f'<b>{T("2nd_hit_lbl")}</b> <font color="#6b7280">{T("not_detected")}</font>')
        # Biallelic somatic MMR events (no germline anchor — e.g. MSH2 cnLOH + splice)
        for b in mmr_biallelic_som:
            g = b.get('gene', '')
            snv_l = b.get('snv_label', '')
            cnv_l = b.get('cnv_label', '')
            # Extract HGVS from all_snvs if available
            snvs = b.get('all_snvs', [])
            hgvs_c = snvs[0].get('hgvs_c', '') if snvs else ''
            vaf_s  = (' VAF=' + snvs[0]['vaf']) if snvs and snvs[0].get('vaf') not in ('', '.', None) else ''
            lines.append(
                f'<font color="#ef4444"><b>{T("biallelic_lbl")}</b></font> '
                f'<b>{g}</b> {hgvs_c}{vaf_s}')
            lines.append(f'&nbsp;&nbsp;&nbsp;{T("2nd_hit_lbl")} {cnv_l or snv_l}')
        if mmr_som and not mmr_germ and not mmr_biallelic_som:
            for g in mmr_som:
                for sv in som_by_gene[g][:2]:
                    hp2  = sv.get('ANN[*].HGVS_P', '') or sv.get('ANN[*].HGVS_C', '') or ''
                    tier = sv.get('_AMP_TIER', 'III')
                    lines.append(f'<b>{T("somatic_var_lbl")}</b> {g} {hp2[:35]} [Tier {tier}]')
        # Separate SBS vs ID signature lines
        if mmr_sbs_sigs:
            lines.append(f'<b>{T("sbs_mmr_lbl")}</b> <font color="#ef4444">{", ".join(mmr_sbs_sigs)}</font>')
        else:
            lines.append(f'<b>{T("sbs_mmr_lbl")}</b> <font color="#6b7280">{T("not_det_f")}</font>')
        if mmr_id_sigs:
            lines.append(f'<b>{T("id_mmr_lbl")}</b> <font color="#ef4444">{", ".join(mmr_id_sigs)}</font>')
        else:
            lines.append(f'<b>{T("id_mmr_lbl")}</b> <font color="#6b7280">{T("not_det_f")}</font>')
        findings.append((T('find_mmr'), colors.HexColor('#6366f1'), lines))

    # ── Finding 3: Double Somatic Hits (only genes with biomarker support) ───
    # Only include genes already covered in correlation results or biallelic_somatic,
    # i.e., genes for which there is somatic + pathway biomarker evidence.
    corr_results     = (corr_data or {}).get('results', [])
    corr_genes       = {r['gene'] for r in corr_results if r.get('biomarker_match')}
    biallelic_genes  = {b['gene'] for b in biallelic_somatic}
    supported_genes  = corr_genes | biallelic_genes
    double_hits = {
        g: vl for g, vl in som_by_gene.items()
        if len(vl) >= 2 and g in supported_genes
    }
    if double_hits:
        lines = []
        for g, vlist in sorted(double_hits.items()):
            descs = []
            for sv in vlist[:3]:
                hp    = sv.get('ANN[*].HGVS_P', '') or sv.get('ANN[*].HGVS_C', '') or '?'
                tier  = sv.get('_AMP_TIER', 'III')
                imp   = sv.get('ANN[*].IMPACT', '')
                descs.append(f'{hp[:25]} (Tier {tier}, {imp})')
            lines.append(f'<b>{g}:</b> ' + ' &nbsp;/&nbsp; '.join(descs))
        findings.append((T('find_dsh'), colors.HexColor('#8b5cf6'), lines))

    # ── Finding 4: High-Confidence Actionable CNVs ───────────────────────────
    hc_homdel = [r for r in cnv_data['homdel'] if r['conf'] == 'HIGH']
    hc_amp    = [r for r in cnv_data['amp']    if r['conf'] == 'HIGH']
    if hc_homdel or hc_amp:
        lines = []
        for r in hc_homdel:
            lines.append(
                f'<b>HOMDEL (CN=0):</b> <font color="#ef4444"><b>{r["gene"]}</b></font>'
                f' &nbsp; chr{r["chrom"]} — confirmed by {r["tools"]}')
        for r in hc_amp:
            lines.append(
                f'<b>AMP (CN={r["cn_p"]}):</b> <font color="#f97316"><b>{r["gene"]}</b></font>'
                f' &nbsp; chr{r["chrom"]} — confirmed by {r["tools"]}')
        findings.append((T('find_cnv'), colors.HexColor('#ea580c'), lines))

    # ── Finding 5: CIViC Actionable Variants (Level A or B, SUPPORTS) ────────
    civic_path = getattr(args, 'civic', '') or ''
    if civic_path and os.path.isfile(civic_path):
        import csv as _csv_civ
        civic_findings = {}  # (gene, var_name) -> {'level': X, 'resist': set, 'sensit': set}
        with open(civic_path) as _cf:
            # Fix header: first column may contain literal \t (bug in CIVIC_EXTRACT)
            raw_header = _cf.readline().rstrip('\n')
            # Expand literal \t in the first column so DictReader keys are correct
            real_header = raw_header.replace('\\t', '\t')
            cols = real_header.split('\t')
            reader = _csv_civ.DictReader(_cf, fieldnames=cols, delimiter='\t')
            for cr in reader:
                lvl = str(cr.get('CIViC Evidence Level', '') or '').strip().upper()
                if lvl not in ('A', 'B'):
                    continue
                direction = str(cr.get('CIViC Entity Direction', '') or '').strip().upper()
                if direction != 'SUPPORTS':
                    continue
                sig = str(cr.get('CIViC Entity Significance', '') or '').strip().upper()
                if 'SENSITIVITY' not in sig and 'RESISTANCE' not in sig:
                    continue
                gene     = str(cr.get('SYMBOL', '') or '').strip()
                var_name = str(cr.get('CIViC Variant Name', '') or '').strip()
                if not gene or not var_name:
                    continue
                therapy  = str(cr.get('CIViC Entity Therapies', '') or '—').strip()
                key = (gene, var_name)
                if key not in civic_findings:
                    civic_findings[key] = {'level': lvl, 'resist': set(), 'sensit': set()}
                # Keep highest level (A > B)
                if lvl < civic_findings[key]['level']:
                    civic_findings[key]['level'] = lvl
                if 'RESISTANCE' in sig:
                    civic_findings[key]['resist'].add(therapy[:40])
                else:
                    civic_findings[key]['sensit'].add(therapy[:40])

        civ_lines = []
        for (gene, var_name), data in sorted(civic_findings.items()):
            lvl = data['level']
            if data['resist']:
                ts = ', '.join(sorted(data['resist'])[:3])
                civ_lines.append(
                    f'<b>{gene} {var_name}</b> [Level {lvl}] '
                    f'— <font color="#ef4444"><b>{T("civic_resist")}:</b></font> {ts}')
            if data['sensit']:
                ts = ', '.join(sorted(data['sensit'])[:3])
                civ_lines.append(
                    f'<b>{gene} {var_name}</b> [Level {lvl}] '
                    f'— <font color="#16a34a"><b>{T("civic_sensit")}:</b></font> {ts}')
        if civ_lines:
            findings.append((T('find_civic'), colors.HexColor('#0f766e'), civ_lines))

    # ── Render findings as styled cards ──────────────────────────────────────
    _card_title_s = ParagraphStyle(
        'card_title', fontName='Helvetica-Bold', fontSize=9,
        textColor=C_WHITE, leading=12, spaceAfter=0, spaceBefore=0,
    )
    _card_body_s = ParagraphStyle(
        'card_body', fontName='Helvetica', fontSize=7.5,
        textColor=C_TEXT, leading=10, spaceAfter=2, spaceBefore=0,
    )

    if not findings:
        story.append(Paragraph(T('no_findings'), STYLE_BODY))
        return

    for title, hdr_color, lines in findings:
        body_content = [Paragraph(ln, _card_body_s) for ln in lines]
        card_data = [
            [Paragraph(title, _card_title_s)],
            [body_content],
        ]
        card_tbl = Table(card_data, colWidths=[USABLE_W])
        card_tbl.setStyle(TableStyle([
            ('BACKGROUND',    (0, 0), (0, 0), hdr_color),
            ('BACKGROUND',    (0, 1), (0, 1), colors.HexColor('#f8fafc')),
            ('BOX',           (0, 0), (-1, -1), 0.8, hdr_color),
            ('LEFTPADDING',   (0, 0), (-1, -1), 10),
            ('RIGHTPADDING',  (0, 0), (-1, -1), 10),
            ('TOPPADDING',    (0, 0), (0, 0), 6),
            ('BOTTOMPADDING', (0, 0), (0, 0), 6),
            ('TOPPADDING',    (0, 1), (0, 1), 6),
            ('BOTTOMPADDING', (0, 1), (0, 1), 8),
            ('VALIGN',        (0, 0), (-1, -1), 'TOP'),
        ]))
        story.append(KeepTogether([card_tbl, Spacer(1, 0.2 * cm)]))
    story.append(Spacer(1, 0.4 * cm))
    # No PageBreak — story continues with Secondary Biomarkers.


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
# CPIC Pharmacogenomics table
# ---------------------------------------------------------------------------
_CPIC_TABLE = {
    'rs3918290':   {'gene':'DPYD',    'allele':'*2A (c.1905+1G>A)',        'activity':'0',       'phenotype':'No Function — Poor Metabolizer',         'drugs':'5-FU, Capecitabine, Tegafur',                         'recommendation':'Avoid fluoropyrimidines. If unavoidable, reduce starting dose by ≥50% and perform therapeutic drug monitoring (TDM). Risk of severe/fatal toxicity.',                     'alert':'HIGH'},
    'rs55886062':  {'gene':'DPYD',    'allele':'*13 (c.1679T>G, p.Ile560Ser)', 'activity':'0',   'phenotype':'No Function — Poor Metabolizer',         'drugs':'5-FU, Capecitabine, Tegafur',                         'recommendation':'Avoid fluoropyrimidines. If unavoidable, reduce starting dose by ≥50% and perform TDM.',                                                                                'alert':'HIGH'},
    'rs67376798':  {'gene':'DPYD',    'allele':'c.2846A>T (p.Asp949Val)',   'activity':'0.5',    'phenotype':'Decreased Function — Intermediate Met.', 'drugs':'5-FU, Capecitabine, Tegafur',                         'recommendation':'Reduce starting dose by 25–50%; titrate based on toxicity and pharmacokinetics. Monitor for hematologic and gastrointestinal toxicity.',                              'alert':'MODERATE'},
    'rs56038477':  {'gene':'DPYD',    'allele':'HapB3 (c.1236G>A, E412E)', 'activity':'0.5',    'phenotype':'Decreased Function — Intermediate Met.', 'drugs':'5-FU, Capecitabine, Tegafur',                         'recommendation':'Reduce starting dose by 25–50%; titrate based on toxicity and pharmacokinetics.',                                                                                     'alert':'MODERATE'},
    'rs4148323':   {'gene':'UGT1A1',  'allele':'*28 ([TA]7)',               'activity':'Reduced','phenotype':'Reduced UGT1A1 Activity',                'drugs':'Irinotecan, Belinostat, Nilotinib',                    'recommendation':'*28/*28 homozygous: reduce irinotecan starting dose by ≥30%. *1/*28 heterozygous: standard dose acceptable; monitor for diarrhea and neutropenia.',                  'alert':'MODERATE'},
    'rs3064744':   {'gene':'UGT1A1',  'allele':'*28 TA7 repeat',           'activity':'Reduced','phenotype':'Reduced UGT1A1 Activity',                'drugs':'Irinotecan',                                          'recommendation':'Consider reduced irinotecan starting dose; monitor for SN-38-related toxicity.',                                                                                    'alert':'MODERATE'},
}

C_PGX_HIGH     = colors.HexColor('#dc2626')
C_PGX_MODERATE = colors.HexColor('#d97706')
C_PGX_LOW      = colors.HexColor('#22c55e')
C_PGX_HEADER   = colors.HexColor('#1e3a8a')
C_PGX_ROW_H    = colors.HexColor('#eff6ff')
C_PGX_ROW_M    = colors.HexColor('#fffbeb')
C_PGX_ROW_N    = colors.HexColor('#f0fdf4')


def build_pharmaco_section(pharmaco_path, story):
    """Pharmacogenomics section: CPIC rsID → phenotype → dose-modification table."""
    story.append(PageBreak())
    story.append(Paragraph('Pharmacogenomics (PGx) — CPIC Dose Recommendations',
                            ParagraphStyle('SH', parent=getSampleStyleSheet()['Normal'],
                                           fontSize=13, fontName='Helvetica-Bold',
                                           textColor=C_DARK_HEADER, spaceAfter=6)))
    story.append(HRFlowable(width='100%', thickness=1, color=C_ACCENT, spaceAfter=8))

    # ── Read pharmaco TSV and look up each rsID in CPIC table ─────────────────
    cpic_rows = []
    if pharmaco_path and os.path.isfile(pharmaco_path):
        import csv as _csv
        with open(pharmaco_path) as fh:
            reader = _csv.DictReader(fh, delimiter='\t')
            for row in reader:
                rsid = (row.get('ID') or '').strip()
                if not rsid or rsid == '.':
                    continue
                # rsID may be semicolon-separated (SnpSift multi-annotation)
                for r in rsid.split(';'):
                    r = r.strip()
                    cpic = _CPIC_TABLE.get(r)
                    if cpic:
                        gt  = (row.get('GEN[0].GT') or '').strip()
                        dp  = (row.get('GEN[0].DP') or row.get('DP') or '').strip()
                        hgvs_c = (row.get('ANN[*].HGVS_C') or '').strip()
                        cpic_rows.append({**cpic, 'rsid': r, 'gt': gt, 'dp': dp, 'hgvs_c': hgvs_c})

    STYLE_H = ParagraphStyle('PGH', parent=getSampleStyleSheet()['Normal'],
                              fontSize=7, fontName='Helvetica-Bold', textColor=C_WHITE)
    STYLE_C = ParagraphStyle('PGC', parent=getSampleStyleSheet()['Normal'],
                              fontSize=7, fontName='Helvetica', textColor=C_TEXT,
                              leading=9, wordWrap='CJK')
    STYLE_SM = ParagraphStyle('PGSM', parent=getSampleStyleSheet()['Normal'],
                               fontSize=6.5, fontName='Helvetica', textColor=C_TEXT,
                               leading=8, wordWrap='CJK')

    if cpic_rows:
        col_headers = ['Gene', 'rsID', 'Allele / Variant', 'Genotype', 'Phenotype', 'Activity\nScore', 'Drugs Affected', 'CPIC Recommendation']
        col_widths  = [38, 55, 80, 42, 70, 35, 75, 158]  # sum=553

        tbl_data = [[Paragraph(h, STYLE_H) for h in col_headers]]
        tbl_style = [
            ('BACKGROUND', (0, 0), (-1, 0), C_PGX_HEADER),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [C_WHITE, C_ROW_ALT]),
            ('GRID', (0, 0), (-1, -1), 0.4, C_LIGHT_BORDER),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTSIZE', (0, 1), (-1, -1), 7),
        ]

        for i, r in enumerate(cpic_rows):
            alert = r.get('alert', 'LOW')
            row_bg = C_PGX_ROW_H if alert == 'HIGH' else (C_PGX_ROW_M if alert == 'MODERATE' else C_PGX_ROW_N)
            tbl_style.append(('BACKGROUND', (0, i+1), (-1, i+1), row_bg))
            alert_color = C_PGX_HIGH if alert == 'HIGH' else (C_PGX_MODERATE if alert == 'MODERATE' else C_PGX_LOW)
            gene_para = Paragraph(
                f'<font color="{alert_color.hexval()}" size="8"><b>{r["gene"]}</b></font>',
                STYLE_C)
            recom_para = Paragraph(r['recommendation'], STYLE_SM)
            tbl_data.append([
                gene_para,
                Paragraph(r['rsid'], STYLE_C),
                Paragraph(r['allele'], STYLE_SM),
                Paragraph(r['gt'] or '—', STYLE_C),
                Paragraph(r['phenotype'], STYLE_SM),
                Paragraph(str(r['activity']), STYLE_C),
                Paragraph(r['drugs'], STYLE_SM),
                recom_para,
            ])

        t = Table(tbl_data, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle(tbl_style))
        story.append(t)
    else:
        # No actionable PGx variants detected
        callout_data = [[
            Paragraph(
                '<b>No Actionable Pharmacogenomic Variants Detected</b><br/>'
                'The pharmacogenomics panel (DPYD, UGT1A1) '
                'did not identify any germline variants with CPIC actionable dosing recommendations. '
                'Standard dosing for fluoropyrimidines and irinotecan '
                'is indicated based on pharmacogenomic genotyping.',
                ParagraphStyle('Callout', parent=getSampleStyleSheet()['Normal'],
                               fontSize=8, fontName='Helvetica', textColor=colors.HexColor('#1e3a8a'),
                               leading=11))
        ]]
        callout_t = Table(callout_data, colWidths=[USABLE_W])
        callout_t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#eff6ff')),
            ('LEFTPADDING', (0, 0), (0, 0), 12),
            ('RIGHTPADDING', (0, 0), (0, 0), 12),
            ('TOPPADDING', (0, 0), (0, 0), 10),
            ('BOTTOMPADDING', (0, 0), (0, 0), 10),
            ('LINEBEFORE', (0, 0), (0, -1), 4, colors.HexColor('#3b82f6')),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#bfdbfe')),
            ('ROUNDEDCORNERS', [6, 6, 6, 6]),
        ]))
        story.append(callout_t)

    story.append(Spacer(1, 6))
    story.append(Paragraph(
        '<i>CPIC = Clinical Pharmacogenomics Implementation Consortium. '
        'Recommendations based on CPIC guidelines for DPYD (2017) and UGT1A1 (2014). '
        'rsID lookup against PharmGKB / CPIC variant annotations.</i>',
        ParagraphStyle('Note', parent=getSampleStyleSheet()['Normal'],
                       fontSize=6.5, fontName='Helvetica-Oblique',
                       textColor=colors.HexColor('#64748b'), leading=8)))


# ---------------------------------------------------------------------------
# Clinical Trials section
# ---------------------------------------------------------------------------
def _query_clinical_trials(gene, max_results=3):
    """Query ClinicalTrials.gov API v2. Returns list of trial dicts or []."""
    import urllib.request
    import urllib.parse
    import json as _json
    try:
        params = urllib.parse.urlencode({
            'query.term': f'{gene} cancer',
            'filter.overallStatus': 'RECRUITING',
            'fields': 'NCTId,BriefTitle,Phase,OverallStatus,Condition,StartDate',
            'pageSize': max_results,
            'format': 'json',
        })
        url = f'https://clinicaltrials.gov/api/v2/studies?{params}'
        req = urllib.request.Request(url, headers={'User-Agent': 'nf-sarek-clinical/1.0'})
        with urllib.request.urlopen(req, timeout=8) as resp:
            data = _json.loads(resp.read().decode('utf-8'))
        studies = data.get('studies', [])
        results = []
        for s in studies:
            proto = s.get('protocolSection', {})
            id_mod = proto.get('identificationModule', {})
            stat_mod = proto.get('statusModule', {})
            design_mod = proto.get('designModule', {})
            cond_mod = proto.get('conditionsModule', {})
            nct  = id_mod.get('nctId', '')
            title = id_mod.get('briefTitle', '')
            phase = ', '.join(design_mod.get('phases', [])) or 'N/A'
            status = stat_mod.get('overallStatus', '')
            conds  = ', '.join(cond_mod.get('conditions', [])[:2])
            start  = (stat_mod.get('startDateStruct', {}) or {}).get('date', '')
            results.append({'nct': nct, 'title': title, 'phase': phase,
                            'status': status, 'conditions': conds, 'start': start})
        return results
    except Exception:
        return None  # None = API unavailable; [] = no results


def build_clinical_trials_section(corr_data, tmb_val, msi_val, hrd_status, story):
    """Clinical Trials section: query ClinicalTrials.gov for relevant recruiting trials."""
    story.append(PageBreak())
    story.append(Paragraph('Clinical Trials — Relevant Recruiting Studies',
                            ParagraphStyle('SH', parent=getSampleStyleSheet()['Normal'],
                                           fontSize=13, fontName='Helvetica-Bold',
                                           textColor=C_DARK_HEADER, spaceAfter=6)))
    story.append(HRFlowable(width='100%', thickness=1, color=C_ACCENT, spaceAfter=8))

    STYLE_H  = ParagraphStyle('CTH',  parent=getSampleStyleSheet()['Normal'],
                               fontSize=7, fontName='Helvetica-Bold', textColor=C_WHITE)
    STYLE_C  = ParagraphStyle('CTC',  parent=getSampleStyleSheet()['Normal'],
                               fontSize=7, fontName='Helvetica', textColor=C_TEXT, leading=9)
    STYLE_SM = ParagraphStyle('CTSM', parent=getSampleStyleSheet()['Normal'],
                               fontSize=6.5, fontName='Helvetica', textColor=C_TEXT,
                               leading=8, wordWrap='CJK')

    # Determine genes/biomarkers to query
    query_genes = []
    if corr_data:
        for r in corr_data.get('results', []):
            if r.get('acmg') in ('P', 'LP'):
                g = r.get('gene', '')
                if g and g not in query_genes:
                    query_genes.append(g)
        for r in corr_data.get('biallelic_somatic', []):
            g = r.get('gene', '')
            if g and g not in query_genes:
                query_genes.append(g)
    # Biomarker queries
    try:
        if float(str(tmb_val).replace('NA','0') or '0') >= 10:
            query_genes.append('TMB-H pembrolizumab')
    except (ValueError, TypeError):
        pass
    try:
        if float(str(msi_val).replace('NA','0') or '0') >= 5.0:
            query_genes.append('microsatellite instability MSI-H')
    except (ValueError, TypeError):
        pass
    if str(hrd_status) == 'HRD+':
        query_genes.append('BRCA homologous recombination deficiency PARP inhibitor')

    if not query_genes:
        story.append(Paragraph(
            'No reportable germline/somatic variants or elevated biomarkers identified '
            'to drive clinical trial queries.',
            ParagraphStyle('Warn', parent=getSampleStyleSheet()['Normal'],
                           fontSize=8, textColor=C_SLATE)))
        return

    all_trials = []
    api_ok = True
    for g in query_genes[:4]:   # limit to 4 queries to avoid latency
        res = _query_clinical_trials(g, max_results=3)
        if res is None:
            api_ok = False
            break
        for t in res:
            if t['nct'] not in {x['nct'] for x in all_trials}:
                t['query'] = g
                all_trials.append(t)

    if not api_ok:
        callout_data = [[Paragraph(
            '<b>ClinicalTrials.gov API Unavailable</b><br/>'
            'Could not connect to ClinicalTrials.gov at report generation time (network timeout). '
            'Please search manually at clinicaltrials.gov for recruiting studies involving: '
            + ', '.join(query_genes[:4]),
            ParagraphStyle('CWarn', parent=getSampleStyleSheet()['Normal'],
                           fontSize=8, textColor=colors.HexColor('#7c2d12'), leading=11))]]
        ct = Table(callout_data, colWidths=[USABLE_W])
        ct.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#fff7ed')),
            ('LINEBEFORE', (0, 0), (0, -1), 4, colors.HexColor('#f97316')),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#fed7aa')),
            ('LEFTPADDING', (0, 0), (0, 0), 12),
            ('TOPPADDING', (0, 0), (0, 0), 10),
            ('BOTTOMPADDING', (0, 0), (0, 0), 10),
        ]))
        story.append(ct)
        return

    if not all_trials:
        story.append(Paragraph(
            'No recruiting clinical trials found for the identified variants/biomarkers '
            'in the ClinicalTrials.gov registry at report generation time.',
            ParagraphStyle('Warn', parent=getSampleStyleSheet()['Normal'],
                           fontSize=8, textColor=C_SLATE)))
        return

    col_headers = ['NCT ID', 'Trial Title', 'Phase', 'Status', 'Condition(s)', 'Start Date', 'Query Gene/Biomarker']
    col_widths  = [60, 180, 45, 60, 80, 50, 78]  # sum=553

    tbl_data = [[Paragraph(h, STYLE_H) for h in col_headers]]
    for t in all_trials:
        tbl_data.append([
            Paragraph(t['nct'], STYLE_C),
            Paragraph(t['title'], STYLE_SM),
            Paragraph(t['phase'], STYLE_C),
            Paragraph(t['status'], STYLE_C),
            Paragraph(t['conditions'], STYLE_SM),
            Paragraph(t['start'], STYLE_C),
            Paragraph(t['query'], STYLE_SM),
        ])

    tbl = Table(tbl_data, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1e3a8a')),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [C_WHITE, C_ROW_ALT]),
        ('GRID', (0, 0), (-1, -1), 0.4, C_LIGHT_BORDER),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 6))
    story.append(Paragraph(
        '<i>Trials sourced from ClinicalTrials.gov at report generation date. '
        'Eligibility criteria must be verified at clinicaltrials.gov before referral. '
        'Only RECRUITING status trials shown. Queries limited to primary actionable findings.</i>',
        ParagraphStyle('Note', parent=getSampleStyleSheet()['Normal'],
                       fontSize=6.5, fontName='Helvetica-Oblique',
                       textColor=colors.HexColor('#64748b'), leading=8)))


# ---------------------------------------------------------------------------
# Methodology section
# ---------------------------------------------------------------------------
def build_methodology_section(story):
    """Final methodology section documenting all callers, tools and calculations."""
    story.append(PageBreak())
    story.append(Paragraph('Test Methodology',
                            ParagraphStyle('SH', parent=getSampleStyleSheet()['Normal'],
                                           fontSize=13, fontName='Helvetica-Bold',
                                           textColor=C_DARK_HEADER, spaceAfter=6)))
    story.append(HRFlowable(width='100%', thickness=1, color=C_ACCENT, spaceAfter=10))

    STYLE_BODY = ParagraphStyle('MB', parent=getSampleStyleSheet()['Normal'],
                                fontSize=8, fontName='Helvetica', textColor=C_TEXT,
                                leading=11, spaceAfter=6)
    STYLE_SUB  = ParagraphStyle('MSub', parent=getSampleStyleSheet()['Normal'],
                                fontSize=9, fontName='Helvetica-Bold',
                                textColor=C_SLATE, spaceAfter=3, spaceBefore=8)
    STYLE_TH   = ParagraphStyle('MTH', parent=getSampleStyleSheet()['Normal'],
                                fontSize=7, fontName='Helvetica-Bold', textColor=C_WHITE)
    STYLE_TC   = ParagraphStyle('MTC', parent=getSampleStyleSheet()['Normal'],
                                fontSize=7, fontName='Helvetica', textColor=C_TEXT,
                                leading=9, wordWrap='CJK')

    # ── Overview table ────────────────────────────────────────────────────────
    story.append(Paragraph('Bioinformatics Pipeline Overview', STYLE_SUB))
    overview_headers = ['Report Section', 'Variant Caller / Tool', 'Genome Build', 'Key Parameters / Thresholds']
    overview_data = [
        ['Germline SNVs',      'GATK HaplotypeCaller v4 (gVCF mode)',              'GRCh38', 'GATK Best Practices; SnpSift extract + SnpEff annotation; ACMG/AMP 2015 classification'],
        ['Somatic SNVs',       'GATK Mutect2 v4 (tumor–normal mode)',              'GRCh38', 'FilterMutectCalls PASS; SnpSift + CIViC annotation; AMP/ASCO/CAP Tier IA–III'],
        ['Somatic CNVs',       'ASCAT v3 + CNVkit v0.9',                           'GRCh38', 'Consensus HOMDEL (CN=0), AMP (CN>=5), 2nd-hit LOH; reported for 38-gene somatic CNV panel'],
        ['Germline CNVs',      'GATK gCNV (case mode)',                            'GRCh38', 'QA>=4 Phred quality; CN!=2 reported for 103-gene germline panel only'],
        ['MSI',                'MSIsensorPro v1.x (tumor–normal)',                 'N/A',    'Score = % unstable microsatellite loci; MSI-H threshold >=5%'],
        ['Somatic Consensus',  'bcftools isec (Mutect2 \u2229 Strelka2)',            'GRCh38', 'Intersection of PASS variants from both callers; Mutect2 FORMAT fields preserved'],
        ['HRD',                'HRProfiler (logistic regression model)',            'N/A',    'Consensus VCF input; SBS + ID + CN48 features; HRD+ threshold: probability >=0.45'],
        ['Mut. Signatures',    'SigProfilerMatrixGenerator + SigProfilerAssignment','N/A',   'Consensus VCF input; COSMIC v3.4 SBS96, DBS and ID signatures'],
        ['TMB',                'Consensus VCF / effective capture Mb',             'GRCh38', 'Non-synonymous PASS variants from consensus; TMB-H threshold >=10 mut/Mb'],
        ['Pharmacogenomics',   'SnpSift extract from HaplotypeCaller gVCF',        'GRCh38', 'DPYD, UGT1A1 rsIDs; CPIC annotation'],
        ['Germline–Somatic',   'Custom integration script (germline_somatic_correlate.py)', 'N/A', 'P/LP germline variants cross-referenced with somatic 2nd hits and HRD/MMR biomarkers'],
    ]

    col_w = [100, 120, 58, 275]
    tbl_data = [[Paragraph(h, STYLE_TH) for h in overview_headers]]
    for row in overview_data:
        tbl_data.append([Paragraph(c, STYLE_TC) for c in row])

    t = Table(tbl_data, colWidths=col_w, repeatRows=1)
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), C_DARK_HEADER),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [C_WHITE, C_ROW_ALT]),
        ('GRID', (0, 0), (-1, -1), 0.4, C_LIGHT_BORDER),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
    ]))
    story.append(t)

    # ── Detailed subsections ──────────────────────────────────────────────────
    story.append(Paragraph('Gene Panels', STYLE_SUB))
    story.append(Paragraph(
        '<b>Germline SNV + CNV panel (103 genes):</b> Used for germline SNV filtering '
        '(HIGH/MODERATE impact variants) and germline CNV reporting (CN\u22602, QA\u22654). '
        'Sources: ACMG SF v3.2, ESMO guidelines, published WES cancer panels.<br/>'
        '<font size="5.5">AIP, ALK, APC, AR, ATM, ATR, AXIN2, BAP1, BARD1, BLM, BMPR1A, BRAF, BRCA1, '
        'BRCA2, BRIP1, BUB1B, CASR, CCND1, CCND2, CCND3, CCNE1, CDC73, CDH1, CDK4, CDK6, CDKN1B, '
        'CDKN2A, CFTR, CHEK2, CTNNA1, CTRC, DICER1, EGFR, ENG, EPCAM, ERBB2, FANCC, FANCD2, FANCM, '
        'FGFR1, FGFR2, FGFR3, FH, FLCN, GREM1, HOXB13, KIF1B, KIT, KRAS, LZTR1, MAX, MBD4, MDM2, '
        'MDM4, MEN1, MET, MITF, MLH1, MLH3, MSH2, MSH3, MSH6, MUTYH, MYC, NBN, NF1, NF2, NTHL1, '
        'PALB2, PDGFRA, PIK3CA, PMS2, POLD1, POLE, POT1, PRKAR1A, PTCH1, PTEN, RAD51C, RB1, RECQL, '
        'RECQL4, RET, RNF43, RPS20, SDHA, SDHAF2, SDHB, SDHC, SDHD, SMAD4, SMARCA4, SMARCB1, '
        'SMARCE1, SPINK1, STK11, SUFU, TMEM127, TP53, TSC1, TSC2, VHL, WRN</font>',
        STYLE_BODY))
    story.append(Paragraph(
        '<b>Somatic SNV panel (60 genes):</b> Used for Mutect2 somatic SNV/indel filtering '
        '(HIGH/MODERATE impact). Covers actionable oncogenes, tumor suppressors, and DNA repair genes.<br/>'
        '<font size="5.5">APC, AR, ARAF, ATM, ATR, B2M, BLM, BRAF, BRCA1, BRCA2, BRIP1, CDH1, CDKN1B, '
        'CDKN2B, CDKN2C, CHEK2, CTNNB1, EGFR, ERBB2, FBXW7, FGFR2, GNAS, IDH1, IDH2, KIT, KRAS, '
        'MAP2K1, MEN1, MLH1, MLH3, MSH2, MSH6, MTOR, NF1, NOTCH3, NRAS, NTHL1, NTRK1, NTRK3, PALB2, '
        'PDGFRA, PDGFRB, PIK3CA, PIK3R1, PMS2, POLD1, POLE, PTEN, RECQL4, RET, RNF43, ROS1, SDHA, '
        'SDHD, SETD2, SMAD4, SMARCA4, STK11, TP53, WRN</font>',
        STYLE_BODY))
    story.append(Paragraph(
        '<b>Somatic CNV panel (38 genes):</b> Used for ASCAT/CNVkit somatic copy-number event '
        'reporting: HOMDEL (CN=0), AMP (CN\u22655), and 2nd-hit LOH.<br/>'
        '<font size="5.5">APC, B2M, BLM, BRAF, BRCA1, BRCA2, CTNNB1, EGFR, EPCAM, ERBB2, FBXW7, '
        'KRAS, MAP2K1, MEN1, MET, MLH1, MLH3, MSH2, MSH3, MSH6, MTOR, MYC, NOTCH3, NRAS, PALB2, '
        'PIK3CA, PIK3R1, PMS2, POLD1, POLE, PTEN, RECQL4, RNF43, SDHA, SDHD, SMAD4, TP53, WRN</font>',
        STYLE_BODY))

    story.append(Paragraph('TMB Calculation', STYLE_SUB))
    story.append(Paragraph(
        'Tumor Mutational Burden (TMB) is calculated as the number of somatic PASS variants '
        'with non-synonymous functional annotation (missense, nonsense, frameshift, splice) '
        'divided by the effective captured exome size (in megabases). '
        '<b>TMB-High (TMB-H) threshold: >=10 mut/Mb</b>, consistent with FDA-approved '
        'pembrolizumab (KEYTRUDA) indication for solid tumors (KEYNOTE-158). '
        'Variants in known germline polymorphism databases (gnomAD AF > 1%) are excluded.',
        STYLE_BODY))

    story.append(Paragraph('MSI — Microsatellite Instability', STYLE_SUB))
    story.append(Paragraph(
        'MSIsensorPro compares microsatellite loci repeat-length distributions between '
        'tumor and matched normal samples using a T/N somatic approach. The MSI score '
        'represents the percentage of unstable microsatellite sites. '
        '<b>MSI-High (MSI-H) threshold: score >=5%</b> (consistent with FDA guidelines). '
        'MSI-H tumors may benefit from immune checkpoint inhibitors (anti-PD-1/PD-L1).',
        STYLE_BODY))

    story.append(Paragraph('HRD — Homologous Recombination Deficiency', STYLE_SUB))
    story.append(Paragraph(
        'HRProfiler applies a logistic regression model trained on TCGA/PCAWG cohorts. '
        'Input features include: COSMIC SBS3, SBS8 (HRD-associated substitution signatures), '
        'ID6 (NHEJ deletion signature), and copy-number segment features. '
        '<b>HRD+ threshold: probability >=0.45.</b> '
        'HRD+ tumors may benefit from PARP inhibitors (olaparib, rucaparib, niraparib, talazoparib). '
        'Supporting biomarkers: SBS3/SBS8 detection and HRD probability are reported jointly.',
        STYLE_BODY))

    story.append(Paragraph('Mutational Signatures', STYLE_SUB))
    story.append(Paragraph(
        'Mutational signatures are derived from somatic variant patterns using the '
        'SigProfiler framework (COSMIC v3.4). '
        '<b>SigProfilerMatrixGenerator</b> constructs SBS96, DBS and ID mutation matrices. '
        '<b>SigProfilerAssignment</b> performs constrained attribution of COSMIC reference '
        'signatures to the tumour spectrum. '
        'Reported signatures: all signatures with activity >0.05 mutations/Mb. '
        'Aetiology annotations sourced from the COSMIC Mutational Signatures website '
        '(cancer.sanger.ac.uk/signatures).',
        STYLE_BODY))

    story.append(Paragraph('ACMG/AMP Variant Classification', STYLE_SUB))
    story.append(Paragraph(
        'Germline variants are classified per ACMG/AMP 2015 standards (Richards et al., '
        'Genetics in Medicine 2015) with updated evidence weights (Tavtigian et al. 2020, '
        'Biesecker et al. 2018). Evidence criteria include: PVS1 (+8), PM2 (+2), PP3 variants '
        '(Strong/Moderate/Supporting: +4/+2/+1), BS1 (-4), BP4 variants (-4/-2). '
        'ClinVar pathogenic/benign classifications override the evidence score when available. '
        'gnomAD AF > 5% triggers BA1 (stand-alone benign). '
        'Score thresholds: P >=10, LP 6-9, VUS -5 to 5, LB -9 to -6, B <=-10.',
        STYLE_BODY))

    story.append(Paragraph('Somatic Variant Interpretation (AMP/ASCO/CAP Tiers)', STYLE_SUB))
    story.append(Paragraph(
        'Somatic variants are tiered per AMP/ASCO/CAP 2017 guidelines using CIViC database '
        'evidence levels (A–E). Tier IA/IB: FDA-approved or guideline-supported therapeutic '
        'implications. Tier IIA/IIB: Clinical trials or emerging evidence. '
        'Tier IIC: Biologically plausible evidence. Tier III: variants of unknown significance '
        '(HIGH/MODERATE impact without CIViC annotation).',
        STYLE_BODY))

    story.append(Paragraph('Germline–Somatic Integration', STYLE_SUB))
    story.append(Paragraph(
        'The germline–somatic correlation module identifies P/LP/VUS germline variants in cancer '
        'predisposition genes and cross-references them with: (1) somatic second-hit variants '
        '(SNV/indel) in the same gene, (2) allele-specific copy-number loss/LOH events from ASCAT '
        'and GATK gCNV, and (3) pathway-level biomarkers (HRD for HRR genes; TMB/MSI for MMR genes). '
        'Genes without a somatic second hit and without pathway biomarker support are excluded '
        'from the correlation table to focus on actionable findings.',
        STYLE_BODY))

    story.append(Spacer(1, 8))
    story.append(Paragraph(
        '<i>This report was generated by the GENUIT nf-sarek-clinical pipeline based on '
        'nf-core/sarek v3.x. All variant interpretations are for research/clinical '
        'investigational use. Results should be confirmed in a certified clinical laboratory '
        'before clinical decision-making. References: GATK Best Practices (Van der Auwera 2013), '
        'ACMG/AMP 2015 (Richards et al.), CPIC guidelines (cpicpgx.org), '
        'COSMIC v3.4 (cancer.sanger.ac.uk/signatures), ClinicalTrials.gov.</i>',
        ParagraphStyle('Disc', parent=getSampleStyleSheet()['Normal'],
                       fontSize=6.5, fontName='Helvetica-Oblique',
                       textColor=colors.HexColor('#64748b'), leading=8)))


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description='Generate GENUIT clinical PDF report for nf-sarek-clinical pipeline.'
    )
    p.add_argument('--pair',       required=True,  help='TUMOR_vs_NORMAL string')
    p.add_argument('--run-name',   required=True,  help='Run/batch identifier')
    p.add_argument('--patient',    required=True,  help='Patient ID')
    p.add_argument('--tumor',      required=True,  help='Tumor sample name')
    p.add_argument('--normal',     required=True,  help='Normal sample name')
    p.add_argument('--germline',   default='',     help='Germline SNV TSV')
    p.add_argument('--somatic',    default='',     help='Somatic SNV TSV')
    p.add_argument('--civic',      default='',     help='CIViC exploded TSV')
    p.add_argument('--tmb',        default='',     help='TMB TSV')
    p.add_argument('--ascat',      default='',     help='ASCAT metrics TSV')
    p.add_argument('--hrd',        default='',     help='HRD predictions TSV')
    p.add_argument('--cnv-somatic',        default='', help='Somatic CNV report TSV (parse_somatic_cnv output)')
    p.add_argument('--cnv-gatk-somatic',  default='', help='GATK somatic CNV annotated TSV (MANE_ANNOTATE_SOMATIC_CNV output)')
    p.add_argument('--cnv-germline-gatk', default='', help='Germline CNV GATK annotated TSV (segments_annotated.tsv)')
    p.add_argument('--somatic-cnv-panel', default='', help='Somatic 38-gene CNV panel (somatic_38_oncogenes.txt)')
    p.add_argument('--sig-sbs96',  default='', help='SBS96 assignment activities TXT')
    p.add_argument('--sig-id',     default='', help='ID assignment activities TXT')
    p.add_argument('--corr-json',  default='',     help='Pre-computed correlation JSON (from germline_somatic_correlate.py)')
    p.add_argument('--gene-ctx',   default='',     help='Gene context TSV (somatic_gene_context.tsv)')
    p.add_argument('--gene-cov-n',       default='', help='Gene coverage TSV for normal sample')
    p.add_argument('--gene-cov-t',       default='', help='Gene coverage TSV for tumor sample')
    p.add_argument('--germline-panel',   default='', help='Germline gene panel list (one gene per line)')
    p.add_argument('--somatic-panel',    default='', help='Somatic gene panel list (one gene per line)')
    p.add_argument('--pharmaco',   default='',     help='Pharmacogenomics TSV (DPYD/UGT1A1 variants)')
    p.add_argument('--out',        required=True,  help='Output PDF path')
    p.add_argument('--language',   default='en',   help='Report language: en (default) or es (Spanish)')
    # Pre-computed metric strings (passed from bash awk extractions)
    p.add_argument('--tmb-val',    default='NA',   help='TMB value string')
    p.add_argument('--msi-val',    default='NA',   help='MSI value string')
    p.add_argument('--loh-val',    default='NA',   help='LOH value string')
    p.add_argument('--gi-val',     default='NA',   help='GI value string')
    p.add_argument('--wgd-val',    default='NA',   help='WGD value string')
    p.add_argument('--ploidy-val', default='NA',   help='Ploidy value string')
    p.add_argument('--purity-val', default='NA',   help='Tumor purity value string')
    p.add_argument('--cov-n-val',  default='NA',   help='Coverage Normal value string')
    p.add_argument('--cov-t-val',  default='NA',   help='Coverage Tumor value string')
    return p.parse_args()


def main():
    args = parse_args()

    # ── Derive HRD metrics from the HRD file if not provided on CLI ──────────
    hrd_prob_val   = 'NA'
    hrd_status_val = 'NA'
    if args.hrd and os.path.isfile(args.hrd):
        _, hrd_rows = read_tsv(args.hrd)
        if hrd_rows:
            row0 = hrd_rows[0]
            hrd_prob_val = str(row0.get('hrd.prob', 'NA') or 'NA').strip()
            raw_pred     = str(row0.get('prediction', '') or '').strip()
            if raw_pred == '1':
                hrd_status_val = 'HRD+'
            elif raw_pred == '0':
                hrd_status_val = 'HRD-'
            else:
                hrd_status_val = raw_pred if raw_pred else 'NA'
            # Override with custom probability threshold (>= 0.45 = HRD+)
            try:
                hrd_status_val = 'HRD+' if float(hrd_prob_val) >= 0.45 else 'HRD-'
            except (ValueError, TypeError):
                pass

    metrics = {
        'tmb':        args.tmb_val    or 'NA',
        'msi':        args.msi_val    or 'NA',
        'loh':        args.loh_val    or 'NA',
        'gi':         args.gi_val     or 'NA',
        'wgd':        args.wgd_val    or 'NA',
        'ploidy':     args.ploidy_val or 'NA',
        'purity':     args.purity_val or 'NA',
        'cov_n':      args.cov_n_val  or 'NA',
        'cov_t':      args.cov_t_val  or 'NA',
        'hrd_prob':   hrd_prob_val,
        'hrd_status': hrd_status_val,
    }

    # ── Variant loading ───────────────────────────────────────────────────────
    germline_variants = []
    if args.germline and os.path.isfile(args.germline):
        germline_variants = load_germline_variants(args.germline)
        print(f'[INFO] Germline: {len(germline_variants)} P/LP/VUS variant(s) after ACMG filtering',
              file=sys.stderr)
    else:
        print(f'[WARN] Germline TSV not found or not provided: {args.germline}', file=sys.stderr)

    somatic_variants = []
    if args.somatic and os.path.isfile(args.somatic):
        somatic_variants = load_somatic_variants(args.somatic, args.civic)
        print(f'[INFO] Somatic: {len(somatic_variants)} variant(s) loaded', file=sys.stderr)
    else:
        print(f'[WARN] Somatic TSV not found or not provided: {args.somatic}', file=sys.stderr)

    gene_ctx_data = {}
    gene_ctx_path = getattr(args, 'gene_ctx', '') or ''
    if gene_ctx_path and os.path.isfile(gene_ctx_path):
        gene_ctx_data = load_gene_ctx(gene_ctx_path)
        print(f'[INFO] Gene context: {len(gene_ctx_data)} gene(s) loaded from {gene_ctx_path}', file=sys.stderr)
    else:
        print(f'[WARN] Gene context TSV not found: {gene_ctx_path}', file=sys.stderr)

    # ── Document setup ────────────────────────────────────────────────────────
    doc = BaseDocTemplate(
        args.out,
        pagesize=A4,
        leftMargin=MARGIN_L,
        rightMargin=MARGIN_R,
        topMargin=MARGIN_T,
        bottomMargin=MARGIN_B,
        title=f'GENUIT Clinical Report – {args.pair}',
        author='GENUIT nf-sarek-clinical',
    )

    frame = Frame(
        MARGIN_L, MARGIN_B,
        USABLE_W, PAGE_H - MARGIN_T - MARGIN_B,
        id='main_frame',
    )
    pt = PageTemplate(
        id='standard',
        frames=[frame],
        onPage=build_footer,
    )
    doc.addPageTemplates([pt])

    # ── Germline-somatic correlation ──────────────────────────────────────────
    def safe_float_m(s, default=0.0):
        try:
            return float(str(s).strip())
        except (ValueError, TypeError):
            return default

    tmb_f      = safe_float_m(args.tmb_val)
    msi_f      = safe_float_m(args.msi_val)
    hrd_prob_f = safe_float_m(hrd_prob_val)
    corr_data  = _run_correlation_script(args, tmb_f, msi_f, hrd_prob_f)

    # ── Language ──────────────────────────────────────────────────────────────
    set_lang(getattr(args, 'language', 'en') or 'en')

    # ── Story ─────────────────────────────────────────────────────────────────
    story = []
    # Each section starts on a new page (leading PageBreak inside each build_* function).
    # Cover page is always page 1 (no leading PageBreak).
    build_cover_page(args, metrics, story)
    build_major_findings_section(args, metrics, germline_variants, somatic_variants,
                                  corr_data, story)
    build_secondary_biomarkers(metrics, story)
    build_quality_metrics_page(metrics, story)
    build_germline_section(germline_variants, story)
    build_germline_cnv_section(args.cnv_germline_gatk, story,
                               germline_panel_path=getattr(args, 'germline_panel', '') or '')
    build_somatic_section(somatic_variants, story, gene_ctx=gene_ctx_data)
    build_mutational_signatures_section(args, story)
    build_somatic_cnv_section(
        getattr(args, 'cnv_somatic', '') or '',
        getattr(args, 'cnv_gatk_somatic', '') or '',
        getattr(args, 'somatic_cnv_panel', '') or '',
        getattr(args, 'germline_panel', '') or '',
        story,
        germline_snv_path=args.germline or '',
        somatic_snv_path=args.somatic or '',
        is_low_cov=metrics.get('is_low_cov', False),
    )
    build_correlation_section(corr_data, story)
    # Pharmacogenomics with CPIC dose recommendations
    build_pharmaco_section(getattr(args, 'pharmaco', '') or '', story)
    # Clinical Trials — recruiting studies for identified variants/biomarkers
    build_clinical_trials_section(corr_data, metrics.get('tmb', 'NA'),
                                   metrics.get('msi', 'NA'), metrics.get('hrd_status', 'NA'), story)
    # Sequencing Quality
    build_sequencing_quality_section(
        getattr(args, 'gene_cov_n', '') or '',
        getattr(args, 'gene_cov_t', '') or '',
        args.normal, args.tumor, story,
        germline_panel_path=getattr(args, 'germline_panel', '') or '',
        somatic_panel_path=getattr(args, 'somatic_panel', '') or '',
    )
    # Methodology section — always last
    build_methodology_section(story)

    # Build PDF
    doc.build(story)
    print(f'[INFO] PDF written: {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()

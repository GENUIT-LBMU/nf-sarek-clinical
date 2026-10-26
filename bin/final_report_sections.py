#!/usr/bin/env python3
"""
final_report_sections.py
========================
Helper functions for the GENUIT nf-sarek-clinical HTML report:
  - Pharmacogenomics rendering
  - ClinicalTrials.gov API query + rendering
  - Test Methodology section rendering

These functions are imported by the FINAL_REPORT Nextflow module
to keep the heredoc small enough for the Groovy parser.
"""

import html
import os
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from acmg_shared import MSI_HIGH_PCT, MSI_BORDER_PCT  # noqa: E402

_GCNV_MIN_QA = str(int(os.environ.get('FILT_GERMLINE_CNV_MIN_QA', '4')))
_GCNV_LOWMAP_PCT = str(int(round(float(os.environ.get('FILT_GERMLINE_CNV_LOWMAP_FRAC', '0.5')) * 100)))
_HRD_HIGH = format(float(os.environ.get('BIO_HRD_HIGH', '0.45')), 'g')
_SOM_MAX_POP_AF = str(float(os.environ.get('FILT_SOMATIC_MAX_POP_AF', '0.01')))

# ---------------------------------------------------------------------------
# CPIC Pharmacogenomics lookup table
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# CPIC Pharmacogenomics lookup table
# ---------------------------------------------------------------------------
def _alert_badge(alert):
    colors_map = {'HIGH': '#dc2626', 'MODERATE': '#d97706', 'LOW': '#16a34a'}
    c = colors_map.get(alert, '#6b7280')
    return "<span style='background:" + c + ";color:#fff;padding:2px 7px;border-radius:10px;font-size:.62rem;font-weight:700;'>" + html.escape(alert) + "</span>"

def tsv_to_html_pharmaco(cpic_rows, not_assessable=(), assessed=()):
    # Report model 'pharmaco': genotype-based phenotype and recommendation per gene
    cpic_rows = [dict(r, gt=html.escape(r['gt']) or '&mdash;') for r in cpic_rows]
    na_html = ("<p class='warn'><strong>Not assessable (outside the capture regions of this assay):</strong> "
               + html.escape(', '.join(not_assessable)) + ". The phenotype of the gene is based on its other "
               "variants and may be underestimated (e.g. UGT1A1 *28 requires a dedicated assay).</p>"
               if not_assessable else "")
    if not cpic_rows:
        return (
            "<div class='callout-info'>"
            "<strong>No Actionable Pharmacogenomic Variants Detected</strong><br/>"
            "None of the assessed variants (" + html.escape(', '.join(assessed)) + ") was found: "
            "normal metabolizer for these variants, standard dosing from the pharmacogenomic "
            "point of view. Other variants of these genes are not tested."
            "</div>" + na_html
        )

    alert_row_bg = {'HIGH': 'rgba(220,38,38,0.08)', 'MODERATE': 'rgba(217,119,6,0.08)', 'LOW': 'rgba(22,163,74,0.06)'}
    th = (
        "<th>Gene</th><th>rsID</th><th>Allele / Variant</th>"
        "<th>Genotype</th><th>Phenotype</th><th>Activity Score</th>"
        "<th>Drugs Affected</th><th>Alert</th><th>Recommendation</th>"
    )
    tbody = ""
    for r in cpic_rows:
        bg = alert_row_bg.get(r['alert'], '')
        style = ("style='background:" + bg + ";'" if bg else "")
        tbody += (
            "<tr " + style + ">"
            "<td><strong>" + html.escape(r['gene']) + "</strong></td>"
            "<td style='font-family:monospace;font-size:.72rem;'>" + html.escape(r['rsid']) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(r['allele']) + "</td>"
            "<td style='font-family:monospace;'>" + r['gt'] + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(r['phenotype']) + "</td>"
            "<td style='text-align:center;'>" + html.escape(str(r['activity'])) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(r['drugs']) + "</td>"
            "<td>" + _alert_badge(r['alert']) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(r['recommendation']) + "<br><span class='tx'>" + html.escape(r['guideline']) + "</span></td>"
            "</tr>\n"
        )
    note = (
        "<p style='font-size:.72rem;color:#94a3b8;margin-top:6px;'>"
        "Phenotype from both alleles of each gene: DPYD activity score (no-function allele -1, decreased-function allele -0.5; two variants of unknown phase taken as on different alleles), CPIC 2017 for fluoropyrimidines; UGT1A1 decreased-function alleles (*28, *6), DPWG for irinotecan. Only the listed variants are genotyped (DPYD *2A, *13, c.2846A>T, HapB3; UGT1A1 *28, *6); a recommendation must be confirmed by a clinical pharmacologist.</p>"
    )
    return (
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody>" + tbody + "</tbody></table></div>" + na_html + note
    )

# ---------------------------------------------------------------------------
# Clinical Trials section (ClinicalTrials.gov API v2)
# ---------------------------------------------------------------------------
def clinical_trials_html(res):
    """Recruiting trials queried once in the report model (report_shared.trial_queries)."""
    if not res:
        return "<p class='warn'>No reportable findings to drive clinical trial queries.</p>"
    all_trials = res['trials']
    queries = [(q,) for q in res['queries']]

    if res['error'] and not all_trials:
        return (
            "<div class='callout-info' style='background:#fff7ed;border-color:#f97316;color:#7c2d12;'>"
            "<strong>ClinicalTrials.gov API Unavailable</strong><br/>"
            "Could not connect to ClinicalTrials.gov on " + res['date'] + " (no network or API down). "
            "Please search manually at <a href='https://clinicaltrials.gov' style='color:#f97316;'>clinicaltrials.gov</a> "
            "for recruiting studies involving: " + html.escape('; '.join(q[0] for q in queries)) +
            "</div>"
        )

    if not all_trials:
        return (
            "<div class='callout-info'>"
            "<strong>No Recruiting Trials Found</strong><br/>"
            "No recruiting clinical trials were found for the identified variants/biomarkers "
            "in the ClinicalTrials.gov registry on " + res['date'] + "."
            "</div>"
        )

    th = ("<th>NCT ID</th><th>Trial Title</th><th>Phase</th>"
          "<th>Status</th><th>Condition(s)</th><th>Start Date</th><th>Query (finding &middot; condition)</th>")
    tbody = ""
    for t in all_trials:
        tbody += (
            "<tr>"
            "<td style='font-family:monospace;font-size:.72rem;white-space:nowrap;'>"
            + html.escape(t['nct']) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(t['title']) + "</td>"
            "<td style='font-size:.72rem;white-space:nowrap;'>" + html.escape(t['phase']) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(t['status']) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(t['conds']) + "</td>"
            "<td style='font-size:.72rem;white-space:nowrap;'>" + html.escape(t['start']) + "</td>"
            "<td style='font-size:.72rem;'>" + html.escape(t['query']) + "</td>"
            "</tr>\n"
        )
    note = (
        "<p style='font-size:.72rem;color:#94a3b8;margin-top:6px;'>"
        "ClinicalTrials.gov queried on " + res['date'] + " (results change over time): actionable somatic "
        "genes (AMP Tier I/II), genes of the germline&ndash;somatic findings and biomarkers, with the tumor "
        "type as condition; MSI-H and TMB-H also in any solid tumor; trials restricted to the opposite MSI "
        "status are not shown. Only RECRUITING trials, up to 3 per query. Eligibility must be verified "
        "before referral.</p>"
    )
    return (
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody>" + tbody + "</tbody></table></div>" + note
    )

# ---------------------------------------------------------------------------
# Test Methodology section
# ---------------------------------------------------------------------------
def build_methodology_html(m):
    # Versions from the report model
    v = m['versions']
    pipeline_label = lambda: v['pipeline']
    sarek_label = lambda: v['sarek']
    version_rows = lambda: v['rows']

    def tool_label(name, key):
        ver = v['tools'].get(key, '')
        return f'{name} {ver}' if ver else name

    method_css = (
        "<style>"
        ".meth-table { border-collapse:collapse; width:100%; font-size:.78rem; margin-top:8px; }"
        ".meth-table th { background:#1e293b; color:#f8fafc; padding:.4rem .6rem; text-align:left; white-space:nowrap; }"
        ".meth-table td { padding:.35rem .6rem; border-bottom:1px solid rgba(148,163,184,0.15); color:#e5e7eb; vertical-align:top; }"
        ".meth-table tr:nth-child(even) td { background:rgba(15,23,42,0.4); }"
        ".meth-h3 { font-size:.9rem; font-weight:700; color:#f97316; margin:18px 0 6px; border-left:3px solid #f97316; padding-left:8px; }"
        ".meth-p { font-size:.8rem; color:#cbd5e1; line-height:1.5; margin:4px 0 10px; }"
        "</style>"
    )
    overview_rows = [
        ['Germline SNVs',     tool_label('GATK HaplotypeCaller', 'gatk4') + ' (gVCF mode)',              'GRCh38', 'GATK Best Practices; SnpEff HGVS on the transcript shown (MANE Select with --snpeff_mane_data; MANE Plus Clinical HGVS added when different); ACMG/AMP 2015 classification'],
        ['Somatic SNVs',      tool_label('GATK Mutect2', 'gatk4') + ' (tumor&ndash;normal mode)',        'GRCh38', 'FilterMutectCalls PASS; SnpEff HGVS on the transcript shown (MANE Select with --snpeff_mane_data; MANE Plus Clinical HGVS added when different); CIViC annotation; AMP/ASCO/CAP Tier IA&ndash;III'],
        ['Somatic CNVs',      tool_label('ASCAT', 'ascat') + ' + ' + tool_label('CNVkit', 'cnvkit') + ' (consensus)',               'GRCh38', 'HOMDEL (CN=0), AMP (CN&ge;5), 2nd-hit LOH; 38-gene somatic CNV panel'],
        ['Germline CNVs',     'GATK gCNV (case mode)',                            'GRCh38', 'QA&ge;' + _GCNV_MIN_QA + ' Phred; deletions / duplications vs contig ploidy; germline panel genes, somatic oncogenes excluded; GIAB low-mappability segments shown apart'],
        ['MSI',               tool_label('MSIsensor-pro', 'msisensor-pro') + ' (tumor&ndash;normal)',           'N/A',    'Score = % unstable microsatellite loci; MSI-H &ge;' + format(MSI_HIGH_PCT, 'g') + '%, borderline &ge;' + format(MSI_BORDER_PCT, 'g') + '%'],
        ['Somatic Consensus', 'bcftools isec (Mutect2 &cap; ' + tool_label('Strelka2', 'strelka') + ')',            'GRCh38', 'Intersection of PASS variants from both callers; Mutect2 FORMAT preserved'],
        ['HRD',               'HRProfiler (logistic regression on SBS+ID+CN)',    'N/A',    'Consensus VCF + ASCAT; HRD+ probability &ge;' + _HRD_HIGH + ' (models validated in breast / ovarian only); SBS3/SBS8'],
        ['Mut. Signatures',   'SigProfilerMatrixGenerator + SigProfilerAssignment','N/A',  'Consensus VCF input; COSMIC v3.4 SBS96, DBS, ID signatures; MMR / HRD calls: &ge;10% of the assigned mutations, &ge;10 mutations, &ge;100 SBS / 50 ID in the sample'],
        ['TMB',               'Consensus VCF / coding territory Mb',              'GRCh38', 'Non-synonymous consensus variants in capture &cap; MANE CDS; TMB-H &ge;10 mut/Mb'],
        ['Pharmacogenomics',  'SnpSift extract from HaplotypeCaller gVCF',        'GRCh38', 'DPYD, UGT1A1 genotypes; phenotype from both alleles (CPIC DPYD, DPWG UGT1A1-irinotecan)'],
        ['Germline&ndash;Somatic', 'germline_somatic_correlate.py',               'N/A',    'P/LP germline x somatic 2nd hits; HRD/MMR pathway biomarkers'],
    ]
    rows_html = ""
    for r in overview_rows:
        rows_html += "<tr>" + "".join("<td>" + c + "</td>" for c in r) + "</tr>\n"
    overview_table = (
        "<div class='tbl-wrap'>"
        "<table class='meth-table'><thead><tr>"
        "<th>Report Section</th><th>Caller / Tool</th><th>Genome Build</th><th>Key Parameters &amp; Thresholds</th>"
        "</tr></thead><tbody>" + rows_html + "</tbody></table></div>"
    )

    detail_sections = [
        ('Gene Panels',
         '<b>Germline SNV + CNV panel (105 genes)</b> &mdash; Used for germline SNV filtering '
         '(HIGH/MODERATE impact; splice-region, synonymous and intronic with SpliceAI &ge;0.2) and germline CNV reporting (deletions / duplications, somatic oncogenes excluded, QA&ge;' + _GCNV_MIN_QA + '; segments with &ge;' + _GCNV_LOWMAP_PCT + '% of exon bases in GIAB low-mappability regions are listed apart and need orthogonal confirmation). '
         'Sources: ACMG SF v3.2, ESMO guidelines, published WES cancer panels.<br/>'
         '<span style="font-size:0.75em; color:#94a3b8;">'
         'AIP, ALK, APC, AR, ATM, ATR, AXIN2, BAP1, BARD1, BLM, BMPR1A, BRAF, BRCA1, '
         'BRCA2, BRIP1, BUB1B, CASR, CCND1, CCND2, CCND3, CCNE1, CDC73, CDH1, CDK4, CDK6, CDKN1B, '
         'CDKN2A, CFTR, CHEK2, CTNNA1, CTRC, DICER1, EGFR, ENG, EPCAM, ERBB2, FANCC, FANCD2, FANCM, '
         'FGFR1, FGFR2, FGFR3, FH, FLCN, GREM1, HOXB13, KIF1B, KIT, KRAS, LZTR1, MAX, MBD4, MDM2, '
         'MDM4, MEN1, MET, MITF, MLH1, MLH3, MRE11, MSH2, MSH3, MSH6, MUTYH, MYC, NBN, NF1, NF2, NTHL1, '
         'PALB2, PDGFRA, PIK3CA, PMS2, POLD1, POLE, POT1, PRKAR1A, PTCH1, PTEN, RAD51C, RAD51D, RB1, RECQL, '
         'RECQL4, RET, RNF43, RPS20, SDHA, SDHAF2, SDHB, SDHC, SDHD, SMAD4, SMARCA4, SMARCB1, '
         'SMARCE1, SPINK1, STK11, SUFU, TMEM127, TP53, TSC1, TSC2, VHL, WRN'
         '</span><br/><br/>'
         '<b>Somatic SNV panel (60 genes)</b> &mdash; Used for Mutect2 somatic SNV/indel filtering '
         '(HIGH/MODERATE impact). Covers actionable oncogenes, tumor suppressors, and DNA repair genes.<br/>'
         '<span style="font-size:0.75em; color:#94a3b8;">'
         'APC, AR, ARAF, ATM, ATR, B2M, BLM, BRAF, BRCA1, BRCA2, BRIP1, CDH1, CDKN1B, '
         'CDKN2B, CDKN2C, CHEK2, CTNNB1, EGFR, ERBB2, FBXW7, FGFR2, GNAS, IDH1, IDH2, KIT, KRAS, '
         'MAP2K1, MEN1, MLH1, MLH3, MSH2, MSH6, MTOR, NF1, NOTCH3, NRAS, NTHL1, NTRK1, NTRK3, PALB2, '
         'PDGFRA, PDGFRB, PIK3CA, PIK3R1, PMS2, POLD1, POLE, PTEN, RECQL4, RET, RNF43, ROS1, SDHA, '
         'SDHD, SETD2, SMAD4, SMARCA4, STK11, TP53, WRN'
         '</span><br/><br/>'
         '<b>Somatic CNV panel (38 genes)</b> &mdash; Used for ASCAT/CNVkit somatic copy-number event '
         'reporting: HOMDEL (CN=0), AMP (CN&ge;5), and 2nd-hit LOH.<br/>'
         '<span style="font-size:0.75em; color:#94a3b8;">'
         'APC, B2M, BLM, BRAF, BRCA1, BRCA2, CTNNB1, EGFR, EPCAM, ERBB2, FBXW7, '
         'KRAS, MAP2K1, MEN1, MET, MLH1, MLH3, MSH2, MSH3, MSH6, MTOR, MYC, NOTCH3, NRAS, PALB2, '
         'PIK3CA, PIK3R1, PMS2, POLD1, POLE, PTEN, RECQL4, RNF43, SDHA, SDHD, SMAD4, TP53, WRN'
         '</span>'),
        ('TMB Calculation',
         'Tumor Mutational Burden (TMB) = somatic non-synonymous variants (missense, nonsense, frameshift, '
         'in-frame indel, splice site, start / stop loss; synonymous variants are not counted) passing in Mutect2 '
         'and Strelka2 (consensus), tumor DP &ge; 10, VAF &ge; 5% and &ge; 3 ALT reads, inside the coding territory '
         '(capture BED &cap; MANE Select CDS, computed per run) / that territory in Mb (WES non-synonymous TMB). '
         'Germline variants are removed by the matched normal. '
         '<b>TMB-High (TMB-H) threshold: &ge;10 mut/Mb</b> (FDA pembrolizumab indication, KEYNOTE-158, '
         'FoundationOne CDx); this WES definition is not calibrated against FoundationOne CDx.'),
        ('MSI &mdash; Microsatellite Instability',
         'MSIsensorPro compares microsatellite repeat-length distributions between tumor and matched normal '
         'using a somatic T/N model. MSI score = % unstable microsatellite sites. '
         '<b>MSI-High (MSI-H) threshold: &ge;' + format(MSI_HIGH_PCT, 'g') + '%</b> (MSIsensor tumor-normal exome cut-off, '
         'Niu et al. 2014); ' + format(MSI_BORDER_PCT, 'g') + '% to ' + format(MSI_HIGH_PCT, 'g') + '%: borderline. '
         'MSI-H tumors may benefit from immune checkpoint inhibitors (anti-PD-1/PD-L1).'),
        ('HRD &mdash; Homologous Recombination Deficiency',
         'HRProfiler (Chopra et al., 2024) is a logistic regression model on six features of the consensus '
         'somatic variants and ASCAT segments: the fractions of N[C&gt;T]G and N[C&gt;G]T substitutions, deletions '
         '&ge;5 bp with microhomology, LOH segments of 1&ndash;40 Mb and heterozygous segments. '
         '<b>Its models were trained and validated in breast and ovarian cancer only</b>; for other tumor types the '
         'breast (or ovarian) model is applied without specific validation, and the probability should be read as '
         'supporting evidence. '
         '<b>HRD+ threshold: probability &ge;' + _HRD_HIGH + '.</b> '
         'HRD+ tumors may benefit from PARP inhibitors (olaparib, rucaparib, niraparib, talazoparib). '
         'Supporting biomarkers (SBS3/SBS8 detection, HRD probability) are reported jointly.'),
        ('Mutational Signatures',
         'Signatures are derived using the SigProfiler framework (COSMIC v3.4). '
         '<b>SigProfilerMatrixGenerator</b> constructs SBS96, DBS, and ID mutation matrices from somatic VCFs. '
         '<b>SigProfilerAssignment</b> performs constrained attribution of COSMIC reference signatures. '
         'All assigned signatures are listed; for the MMR / HRD findings a signature counts as detected with '
         '&ge;10% of the assigned mutations, &ge;10 mutations and &ge;100 SBS / 50 indels in the sample (fewer: not assessable). '
         'Aetiology annotations from COSMIC Mutational Signatures (cancer.sanger.ac.uk/signatures).'),
        ('ACMG/AMP Variant Classification (Germline)',
         'Germline variants are classified per ACMG/AMP 2015 standards (Richards et al., '
         '<i>Genetics in Medicine</i> 2015) with updated evidence weights (Tavtigian et al. 2020). '
         'Evidence points: PVS1 (+8; +4 / +2 in the last exon, +2 for start-lost; not applied in genes where '
         'loss of function is not the cancer mechanism), PM4 stop-lost (+2), PM2_Supporting (+1), '
         'PP3 / BP4 from REVEL with the ClinGen calibrated cut-offs (Pejaver et al. 2022: '
         '&plusmn;1 / &plusmn;2 / &plusmn;4) and SpliceAI (Walker et al. 2023; splice-region, synonymous and intronic variants '
         'are reported only with SpliceAI &ge;0.2 or ClinVar P/LP), PS1 (+4) / PM5 (+2) from ClinVar '
         'P/LP missense at the same residue, PM1 (+2 / +1; PP3 + PM1 &le; 4), BS1 (&minus;4), BS1_Supporting (&minus;1). '
         'ClinVar classifications with &ge;2 review stars are final. gnomAD frequency&nbsp;&gt;&nbsp;5% triggers BA1 (stand-alone benign); BA1, BS1 and PM2 use the maximum filtering allele frequency (FAF95) across ancestry groups. '
         'For APC, ATM, BRAF, BRCA1, BRCA2, CDH1, DICER1, ENG, KRAS, MLH1, MSH2, MSH6, PALB2, PMS2, PTEN, TP53 and VHL, '
         'BA1 / BS1 / PM2 thresholds and PM1 hotspots / domains follow the ClinGen Variant Curation Expert Panel '
         'specifications (CSpec registry; PM1 hotspot counts from cancerhotspots.org v2); ClinVar P/LP with &ge;2 stars is never removed by BA1. '
         'Point thresholds: P &ge;10, LP 6&ndash;9, VUS 0&ndash;5, LB &minus;1 to &minus;6, B &le;&minus;7.'),
        ('Somatic Variant Interpretation (AMP/ASCO/CAP Tiers)',
         'Somatic variants are tiered per AMP/ASCO/CAP 2017 from CIViC evidence. Only accepted evidence whose '
         'direction supports the claim sets the tier; submitted evidence and evidence against are listed but never raise it. '
         'In the patient&#39;s tumor type (--tumor_type / --tumor_types; pan-cancer evidence included): CIViC level '
         '<b>A &rarr; IA</b>, <b>B &rarr; IB</b>, <b>C &rarr; IIC</b>, <b>D/E &rarr; IID</b>; level A/B from another '
         'tumor type &rarr; <b>IIC</b>. Without a tumor type the highest tier is IIC. '
         '<b>Tier III</b>: HIGH/MODERATE impact variants without tier-defining evidence. '
         '<b>Tier IV</b>: population polymorphisms (local gnomAD v4.1 FAF95, else AF, &ge; ' + _SOM_MAX_POP_AF + '), '
         'listed apart as probable germline variants and not used as somatic findings.'),
        ('Germline&ndash;Somatic Integration',
         'The correlation module cross-references P/LP/VUS germline variants in cancer-predisposition genes '
         'with (1) somatic second-hit SNVs/indels in the same gene, '
         '(2) allele-specific copy-number loss/LOH (from ASCAT+CNVkit and GATK gCNV), and '
         '(3) pathway-level biomarkers: HRD probability &ge;' + _HRD_HIGH + ' or SBS3/SBS8 for HRR genes; '
         'TMB &ge;10 mut/Mb, MSI &ge;3.5%, or SBS6/15/20/21/26/44 for MMR genes. '
         'Genes without a somatic second hit <i>and</i> without pathway biomarker support are excluded.'),
    ]
    details_html = ""
    for title, body in detail_sections:
        details_html += (
            "<div class='meth-h3'>" + title + "</div>"
            "<p class='meth-p'>" + body + "</p>"
        )

    disclaimer = (
        "<p style='margin-top:20px;font-size:.72rem;color:#64748b;border-top:1px solid rgba(148,163,184,0.2);padding-top:10px;'>"
        "<i>This report was generated by " + html.escape(pipeline_label()) + " on " + html.escape(sarek_label()) + ". "
        "All variant interpretations are for investigational use. Results must be confirmed in a "
        "certified clinical laboratory before clinical decision-making. "
        "References: GATK Best Practices; ACMG/AMP 2015 (Richards et al.); "
        "CPIC guidelines (cpicpgx.org); DPWG (UGT1A1-irinotecan); COSMIC v3.4 (cancer.sanger.ac.uk); ClinicalTrials.gov.</i>"
        "</p>"
    )

    versions_html = (
        "<div class='meth-h3'>Versions and Resources</div>"
        "<table class='meth-table'><thead><tr><th>Component</th><th>Version / file</th></tr></thead><tbody>"
        + ''.join("<tr><td>" + html.escape(k) + "</td><td>" + html.escape(str(v)) + "</td></tr>"
                  for k, v in version_rows())
        + "</tbody></table>"
    )
    return method_css + overview_table + versions_html + details_html + disclaimer


# ----
# Germline CNV GATK table: ALT <DEL>/<DUP> (CN vs contig ploidy), deduplicated by segment, quality-badged
# ----
def germline_cnv_is_event(row):
    # gCNV sets ALT to <DEL>/<DUP> against the contig ploidy (male chrX: 1 copy is
    # reference). Old tables without ALT fall back to CN != 2.
    if 'ALT' in row:
        return (row.get('ALT') or '').strip() in ('<DEL>', '<DUP>')
    try:
        return int((row.get('CN') or '2').strip()) != 2
    except ValueError:
        return False

def germline_cnv_lowmap(row):
    # LOWMAP_SEG = fraction of the segment's panel exon bp in low-mappability
    # regions (pseudogenes / paralogs); None for old tables without the column.
    try:
        return float((row.get('LOWMAP_SEG') or '').strip())
    except ValueError:
        return None

def germline_cnv_html(gcnv, thresholds):
    # Report model 'germline_cnv': segments grouped per (CHROM, POS, END), panel genes,
    # QA >= min; low-mappability segments in their own table, never in the main one.
    MIN_QA = thresholds['germline_cnv_min_qa']
    LOWMAP_FRAC = thresholds['germline_cnv_lowmap_frac']
    if not gcnv.get('n_rows'):
        return "<p class='warn'>Germline CNV (GATK): empty file</p>"
    filtered = gcnv['segments']
    lowmap_segs = gcnv['lowmap']
    main_html = ""
    if not filtered:
        main_html = (
            "<div class='callout-info'>"
            "<strong>No Reportable Germline CNV Events</strong><br>"
            "The GATK germline copy-number variant analysis did not identify any deletion or duplication "
            "meeting the quality threshold (QA &ge; " + str(MIN_QA) + ") in the germline panel genes "
            "(somatic oncogenes excluded)"
            + (" outside low-mappability regions (see below). " if lowmap_segs else ". ") +
            "Raw annotated segments are available for review in "
            "<code>germline/cnv_gatk/annotated/</code>."
            "</div>"
        )
    def qs_badge(qs_str):
        try:
            qs = float(qs_str)
        except (ValueError, TypeError):
            return "<span class='badge badge-III'>Unknown</span>"
        if qs >= 30:
            return "<span class='badge badge-P'>High</span>"
        elif qs >= 10:
            return "<span class='badge badge-LP'>Medium</span>"
        else:
            return "<span class='badge badge-VUS'>Low</span>"
    def cn_type(alt, cn_str):
        try:
            cn = int(cn_str)
        except (ValueError, TypeError):
            cn = -1
        if alt == '<DEL>':
            return "HOMDEL" if cn == 0 else "DEL (CN=" + str(cn) + ")"
        if alt == '<DUP>':
            return ("AMP (CN=" + str(cn) + ")") if cn >= 4 else ("DUP (CN=" + str(cn) + ")")
        if cn == 0: return "HOMDEL"
        elif cn == 1: return "DEL (CN=1)"
        elif cn == 3: return "DUP (CN=3)"
        elif cn >= 4: return "AMP (CN=" + str(cn) + ")"
        return alt or "CNV"
    def seg_rows(segs, lowmap):
        out = ""
        for seg in segs:
            locus = seg['chrom'] + ":" + seg['pos'] + "-" + seg['end']
            genes = ", ".join(seg['genes']) if seg['genes'] else "."
            ctype = cn_type(seg['alt'], seg['cn'])
            last = ("<td>" + str(int(round(seg['lowmap'] * 100))) + "%</td>") if lowmap \
                else ("<td>" + qs_badge(seg['qs']) + "</td>")
            out += (
                "<tr>"
                "<td>" + html.escape(locus) + "</td>"
                "<td><strong>" + html.escape(genes) + "</strong></td>"
                "<td>" + html.escape(ctype) + "</td>"
                "<td>" + html.escape(seg['cn']) + "</td>"
                "<td>" + html.escape(seg['np']) + "</td>"
                "<td>" + html.escape(seg['qa']) + "</td>"
                "<td>" + html.escape(seg['qs']) + "</td>"
                + last +
                "</tr>\n"
            )
        return out
    th = ("<th>Locus</th><th>Gene(s)</th><th>Type</th><th>CN</th>"
          "<th>NP (intervals)</th><th>QA</th><th>QS</th>")
    lowmap_html = ""
    if lowmap_segs:
        lowmap_html = (
            "<h4 style='margin:14px 0 4px;'>Low-mappability regions (pseudogenes / paralogs) &mdash; "
            "not reportable without orthogonal confirmation</h4>"
            "<div style='font-size:.75rem;color:rgba(226,232,240,.6);margin:0 0 4px;'>"
            + str(len(lowmap_segs)) + " segment(s) with &ge; " + str(int(round(LOWMAP_FRAC * 100)))
            + "% of their panel exon bases in GIAB low-mappability regions. Read counts there are "
            "unreliable (reads map to the pseudogene with low MAPQ); confirm by MLPA or another "
            "orthogonal method before clinical use. Not used in the germline&ndash;somatic correlation."
            "</div>"
            "<div class='tbl-wrap'><table><thead><tr>" + th + "<th>Low-map exon bases</th></tr></thead>"
            "<tbody>" + seg_rows(lowmap_segs, True) + "</tbody></table></div>"
        )
    if not filtered:
        return main_html + lowmap_html
    tbody = seg_rows(filtered, False)
    note = ("<div style='font-size:.75rem;color:rgba(226,232,240,.6);margin:6px 0 4px;'>"
            + str(len(filtered)) + " deletion / duplication segment(s) (germline panel genes, somatic oncogenes excluded; QA &ge; " + str(MIN_QA) + ")."
            " Quality: <span class='badge badge-P'>High</span> QS &ge; 30 &nbsp;"
            "<span class='badge badge-LP'>Medium</span> QS 10&ndash;29 &nbsp;"
            "<span class='badge badge-VUS'>Low</span> QS &lt; 10 (manual review recommended)."
            "</div>")
    return (
        note +
        "<div class='tbl-wrap'><table><thead><tr>" + th + "<th>Quality</th></tr></thead>"
        "<tbody>" + tbody + "</tbody></table></div>"
        + lowmap_html
    )

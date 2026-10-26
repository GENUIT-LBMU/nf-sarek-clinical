#!/usr/bin/env python3
"""
final_report_sections.py
========================
Helper functions for the GENUIT nf-sarek-clinical HTML report:
  - CPIC Pharmacogenomics lookup table + rendering
  - ClinicalTrials.gov API query + rendering
  - Test Methodology section rendering

These functions are imported by the FINAL_REPORT Nextflow module
to keep the heredoc small enough for the Groovy parser.
"""

import html
import os
import csv

_GCNV_MIN_QA = str(int(os.environ.get('FILT_GERMLINE_CNV_MIN_QA', '4')))
_GCNV_LOWMAP_PCT = str(int(round(float(os.environ.get('FILT_GERMLINE_CNV_LOWMAP_FRAC', '0.5')) * 100)))
_SOM_MAX_POP_AF = str(float(os.environ.get('FILT_SOMATIC_MAX_POP_AF', '0.01')))

# ---------------------------------------------------------------------------
# CPIC Pharmacogenomics lookup table
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# CPIC Pharmacogenomics lookup table
# ---------------------------------------------------------------------------
_CPIC_TABLE = {
    'rs3918290':   {'gene':'DPYD',    'allele':'*2A (c.1905+1G>A)',         'activity':'0',       'phenotype':'No Function -- Poor Metabolizer',          'drugs':'5-FU, Capecitabine, Tegafur',                       'recommendation':'Avoid fluoropyrimidines. If unavoidable, reduce starting dose by >=50% and perform therapeutic drug monitoring. Risk of severe/fatal toxicity.',           'alert':'HIGH'},
    'rs55886062':  {'gene':'DPYD',    'allele':'*13 (c.1679T>G)',            'activity':'0',       'phenotype':'No Function -- Poor Metabolizer',          'drugs':'5-FU, Capecitabine, Tegafur',                       'recommendation':'Avoid fluoropyrimidines. If unavoidable, reduce starting dose by >=50% and perform TDM.',                                                                'alert':'HIGH'},
    'rs67376798':  {'gene':'DPYD',    'allele':'c.2846A>T (p.Asp949Val)',    'activity':'0.5',     'phenotype':'Decreased Function -- Intermediate Met.',   'drugs':'5-FU, Capecitabine, Tegafur',                       'recommendation':'Reduce starting dose by 25-50%; titrate based on toxicity and pharmacokinetics.',                                                                        'alert':'MODERATE'},
    'rs56038477':  {'gene':'DPYD',    'allele':'HapB3 (c.1236G>A, E412E)',  'activity':'0.5',     'phenotype':'Decreased Function -- Intermediate Met.',   'drugs':'5-FU, Capecitabine, Tegafur',                       'recommendation':'Reduce starting dose by 25-50%; titrate based on toxicity and pharmacokinetics.',                                                                        'alert':'MODERATE'},
    'rs4148323':   {'gene':'UGT1A1',  'allele':'*28 ([TA]7)',                'activity':'Reduced', 'phenotype':'Reduced UGT1A1 Activity',                  'drugs':'Irinotecan, Belinostat, Nilotinib',                 'recommendation':'*28/*28: reduce irinotecan starting dose by >=30%. *1/*28: standard dose; monitor for diarrhea/neutropenia.',                                            'alert':'MODERATE'},
    'rs3064744':   {'gene':'UGT1A1',  'allele':'*28 TA7 repeat',            'activity':'Reduced', 'phenotype':'Reduced UGT1A1 Activity',                  'drugs':'Irinotecan',                                        'recommendation':'Consider reduced irinotecan starting dose; monitor for SN-38-related toxicity.',                                                                          'alert':'MODERATE'},
    'rs116855232': {'gene':'NUDT15',  'allele':'*2 (c.415C>T, p.Arg139Cys)','activity':'0',       'phenotype':'No Function -- Poor Metabolizer',          'drugs':'Mercaptopurine, Azathioprine, Thioguanine',         'recommendation':'Start thiopurine at 10% of standard dose; adjust based on hematological toxicity and TGN nucleotide levels.',                                           'alert':'HIGH'},
    'rs186364861': {'gene':'NUDT15',  'allele':'*3',                        'activity':'0',       'phenotype':'No Function -- Poor Metabolizer',          'drugs':'Mercaptopurine, Azathioprine, Thioguanine',         'recommendation':'Start thiopurine at 10% of standard dose; adjust based on toxicity monitoring.',                                                                        'alert':'HIGH'},
    'rs1800462':   {'gene':'TPMT',    'allele':'*2 (c.238G>C)',             'activity':'Decreased','phenotype':'Decreased Function',                       'drugs':'Mercaptopurine, Azathioprine, Thioguanine',         'recommendation':'Intermediate metabolizer: reduce thiopurine starting dose by 30-70%; titrate based on hematological toxicity.',                                         'alert':'MODERATE'},
    'rs1800460':   {'gene':'TPMT',    'allele':'*3B (c.460G>A)',            'activity':'Decreased','phenotype':'Decreased Function',                       'drugs':'Mercaptopurine, Azathioprine, Thioguanine',         'recommendation':'Intermediate metabolizer: reduce thiopurine starting dose by 30-70%.',                                                                                   'alert':'MODERATE'},
    'rs1142345':   {'gene':'TPMT',    'allele':'*3C (c.719A>G)',            'activity':'Decreased','phenotype':'Decreased Function',                       'drugs':'Mercaptopurine, Azathioprine, Thioguanine',         'recommendation':'Intermediate metabolizer: reduce thiopurine starting dose by 30-70%.',                                                                                   'alert':'MODERATE'},
    'rs3892097':   {'gene':'CYP2D6',  'allele':'*4 (c.506-1G>A)',           'activity':'0',       'phenotype':'No Function -- Poor Metabolizer',          'drugs':'Tamoxifen',                                         'recommendation':'Consider aromatase inhibitor (AI) rather than tamoxifen for postmenopausal ER+ patients. CYP2D6 PM have impaired conversion to active endoxifen.',       'alert':'HIGH'},
    'rs1065852':   {'gene':'CYP2D6',  'allele':'*10 (c.100C>T)',            'activity':'0.5',     'phenotype':'Decreased Function -- Intermediate Met.',   'drugs':'Tamoxifen',                                         'recommendation':'Consider alternative hormonal therapy (AI) if clinically feasible; reduced endoxifen levels may compromise tamoxifen efficacy.',                           'alert':'MODERATE'},
    'rs4244285':   {'gene':'CYP2C19', 'allele':'*2 (c.681G>A)',             'activity':'0',       'phenotype':'No Function -- Poor Metabolizer',          'drugs':'Clopidogrel, Voriconazole',                         'recommendation':'For clopidogrel: use alternative antiplatelet (prasugrel, ticagrelor). For voriconazole: dose reduction; PMs have substantially increased exposure.',      'alert':'MODERATE'},
}

def _alert_badge(alert):
    colors_map = {'HIGH': '#dc2626', 'MODERATE': '#d97706', 'LOW': '#16a34a'}
    c = colors_map.get(alert, '#6b7280')
    return "<span style='background:" + c + ";color:#fff;padding:2px 7px;border-radius:10px;font-size:.62rem;font-weight:700;'>" + html.escape(alert) + "</span>"

def tsv_to_html_pharmaco(path):
    cpic_rows = []
    if path and os.path.isfile(path):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                rsid_raw = (row.get('ID') or '').strip()
                if not rsid_raw or rsid_raw == '.':
                    continue
                for rsid in rsid_raw.split(';'):
                    rsid = rsid.strip()
                    cpic = _CPIC_TABLE.get(rsid)
                    if cpic:
                        cpic_rows.append({
                            'rsid':   rsid,
                            'gene':   cpic['gene'],
                            'allele': cpic['allele'],
                            'gt':     (row.get('GEN[0].GT') or '').strip() or '&mdash;',
                            'phenotype': cpic['phenotype'],
                            'activity':  cpic['activity'],
                            'drugs':     cpic['drugs'],
                            'recommendation': cpic['recommendation'],
                            'alert':         cpic['alert'],
                        })

    if not cpic_rows:
        return (
            "<div class='callout-info'>"
            "<strong>No Actionable Pharmacogenomic Variants Detected</strong><br/>"
            "The pharmacogenomics panel (DPYD, UGT1A1, NUDT15, TPMT, CYP2D6, CYP2C19) "
            "did not identify any germline variants with CPIC actionable dosing recommendations. "
            "Standard dosing for fluoropyrimidines, irinotecan, thiopurines, and tamoxifen "
            "is indicated based on pharmacogenomic genotyping."
            "</div>"
        )

    alert_row_bg = {'HIGH': 'rgba(220,38,38,0.08)', 'MODERATE': 'rgba(217,119,6,0.08)', 'LOW': 'rgba(22,163,74,0.06)'}
    th = (
        "<th>Gene</th><th>rsID</th><th>Allele / Variant</th>"
        "<th>Genotype</th><th>Phenotype</th><th>Activity Score</th>"
        "<th>Drugs Affected</th><th>Alert</th><th>CPIC Recommendation</th>"
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
            "<td style='font-size:.72rem;'>" + html.escape(r['recommendation']) + "</td>"
            "</tr>\n"
        )
    note = (
        "<p style='font-size:.72rem;color:#94a3b8;margin-top:6px;'>"
        "CPIC = Clinical Pharmacogenomics Implementation Consortium. "
        "Recommendations based on CPIC guidelines for DPYD (2017), UGT1A1 (2014), "
        "NUDT15/TPMT (2018), CYP2D6/Tamoxifen (2018), CYP2C19 (2013).</p>"
    )
    return (
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody>" + tbody + "</tbody></table></div>" + note
    )

# ---------------------------------------------------------------------------
# Clinical Trials section (ClinicalTrials.gov API v2)
# ---------------------------------------------------------------------------
def fetch_clinical_trials_html(corr_json_path, tmb_val, msi_val, hrd_status):
    import urllib.request, urllib.parse, json as _json

    # Determine query genes from correlation JSON
    query_terms = []
    if corr_json_path and os.path.isfile(corr_json_path):
        try:
            with open(corr_json_path) as fh:
                cdata = _json.load(fh)
            for r in cdata.get('results', []):
                if r.get('acmg') in ('P', 'LP'):
                    g = r.get('gene', '')
                    if g and g not in query_terms:
                        query_terms.append(g)
            for r in cdata.get('biallelic_somatic', []):
                g = r.get('gene', '')
                if g and g not in query_terms:
                    query_terms.append(g)
        except Exception:
            pass
    try:
        if float(str(tmb_val).replace('NA','0') or '0') >= 10:
            query_terms.append('TMB-H pembrolizumab')
    except (ValueError, TypeError):
        pass
    try:
        if float(str(msi_val).replace('NA','0') or '0') >= 5.0:
            query_terms.append('microsatellite instability MSI-H')
    except (ValueError, TypeError):
        pass
    if str(hrd_status) == 'HRD+':
        query_terms.append('BRCA homologous recombination deficiency PARP inhibitor')

    if not query_terms:
        return "<p class='warn'>No reportable findings to drive clinical trial queries.</p>"

    all_trials = []
    api_error = None
    for term in query_terms[:4]:
        try:
            params = urllib.parse.urlencode({
                'query.term': term + ' cancer',
                'filter.overallStatus': 'RECRUITING',
                'fields': 'NCTId,BriefTitle,Phase,OverallStatus,Condition,StartDate',
                'pageSize': '3',
                'format': 'json',
            })
            url = 'https://clinicaltrials.gov/api/v2/studies?' + params
            req = urllib.request.Request(url, headers={'User-Agent': 'nf-sarek-clinical/1.0'})
            with urllib.request.urlopen(req, timeout=8) as resp:
                data = _json.loads(resp.read().decode('utf-8'))
            existing_ncts = {t['nct'] for t in all_trials}
            for s in data.get('studies', []):
                proto  = s.get('protocolSection', {})
                id_mod = proto.get('identificationModule', {})
                stat   = proto.get('statusModule', {})
                des    = proto.get('designModule', {})
                cond   = proto.get('conditionsModule', {})
                nct    = id_mod.get('nctId', '')
                if nct and nct not in existing_ncts:
                    all_trials.append({
                        'nct':    nct,
                        'title':  id_mod.get('briefTitle', ''),
                        'phase':  ', '.join(des.get('phases', [])) or 'N/A',
                        'status': stat.get('overallStatus', ''),
                        'conds':  ', '.join(cond.get('conditions', [])[:2]),
                        'start':  (stat.get('startDateStruct') or {}).get('date', ''),
                        'query':  term,
                    })
        except Exception as e:
            api_error = str(e)
            break

    if api_error and not all_trials:
        return (
            "<div class='callout-info' style='background:#fff7ed;border-color:#f97316;color:#7c2d12;'>"
            "<strong>ClinicalTrials.gov API Unavailable</strong><br/>"
            "Could not connect to ClinicalTrials.gov at report generation time (timeout/network). "
            "Please search manually at <a href='https://clinicaltrials.gov' style='color:#f97316;'>clinicaltrials.gov</a> "
            "for recruiting studies involving: " + html.escape(', '.join(query_terms[:4])) +
            "</div>"
        )

    if not all_trials:
        return (
            "<div class='callout-info'>"
            "<strong>No Recruiting Trials Found</strong><br/>"
            "No recruiting clinical trials were found for the identified variants/biomarkers "
            "in the ClinicalTrials.gov registry at report generation time."
            "</div>"
        )

    th = ("<th>NCT ID</th><th>Trial Title</th><th>Phase</th>"
          "<th>Status</th><th>Condition(s)</th><th>Start Date</th><th>Query Gene/Biomarker</th>")
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
        "Sourced from ClinicalTrials.gov at report generation date. "
        "Only RECRUITING status trials shown. Eligibility must be verified before referral.</p>"
    )
    return (
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody>" + tbody + "</tbody></table></div>" + note
    )

# ---------------------------------------------------------------------------
# Test Methodology section
# ---------------------------------------------------------------------------
def build_methodology_html():
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
        ['Germline SNVs',     'GATK HaplotypeCaller v4 (gVCF mode)',              'GRCh38', 'GATK Best Practices; SnpEff HGVS on the transcript shown (MANE Select with --snpeff_mane_data; MANE Plus Clinical HGVS added when different); ACMG/AMP 2015 classification'],
        ['Somatic SNVs',      'GATK Mutect2 v4 (tumor&ndash;normal mode)',        'GRCh38', 'FilterMutectCalls PASS; SnpEff HGVS on the transcript shown (MANE Select with --snpeff_mane_data; MANE Plus Clinical HGVS added when different); CIViC annotation; AMP/ASCO/CAP Tier IA&ndash;III'],
        ['Somatic CNVs',      'ASCAT v3 + CNVkit v0.9 (consensus)',               'GRCh38', 'HOMDEL (CN=0), AMP (CN&ge;5), 2nd-hit LOH; 38-gene somatic CNV panel'],
        ['Germline CNVs',     'GATK gCNV (case mode)',                            'GRCh38', 'QA&ge;' + _GCNV_MIN_QA + ' Phred; deletions / duplications vs contig ploidy; germline panel genes, somatic oncogenes excluded; GIAB low-mappability segments shown apart'],
        ['MSI',               'MSIsensorPro v1.x (tumor&ndash;normal)',           'N/A',    'Score = % unstable microsatellite loci; MSI-H threshold &ge;5%'],
        ['Somatic Consensus', 'bcftools isec (Mutect2 &cap; Strelka2)',            'GRCh38', 'Intersection of PASS variants from both callers; Mutect2 FORMAT preserved'],
        ['HRD',               'HRProfiler (logistic regression on SBS+ID+CN)',    'N/A',    'Consensus VCF input; HRD+ threshold: probability &ge;0.45; SBS3/SBS8'],
        ['Mut. Signatures',   'SigProfilerMatrixGenerator + SigProfilerAssignment','N/A',  'Consensus VCF input; COSMIC v3.4 SBS96, DBS, ID signatures; &gt;0.05 mut/Mb'],
        ['TMB',               'Consensus VCF / coding territory Mb',              'GRCh38', 'Non-synonymous consensus variants in capture &cap; MANE CDS; TMB-H &ge;10 mut/Mb'],
        ['Pharmacogenomics',  'SnpSift extract from HaplotypeCaller gVCF',        'GRCh38', 'DPYD, UGT1A1, NUDT15, TPMT, CYP2D6, CYP2C19; CPIC guideline annotation'],
        ['Germline&ndash;Somatic', 'germline_somatic_correlate.py',               'N/A',    'P/LP/VUS x somatic 2nd hits; HRD/MMR pathway biomarkers'],
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
         '<b>MSI-High (MSI-H) threshold: &ge;5%.</b> '
         'MSI-H tumors may benefit from immune checkpoint inhibitors (anti-PD-1/PD-L1).'),
        ('HRD &mdash; Homologous Recombination Deficiency',
         'HRProfiler applies a logistic regression model trained on TCGA/PCAWG cohorts using '
         'COSMIC SBS3, SBS8, ID6 signatures and copy-number features as inputs. '
         '<b>HRD+ threshold: probability &ge;0.45.</b> '
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
         '(3) pathway-level biomarkers: HRD probability &ge;0.45 or SBS3/SBS8 for HRR genes; '
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
        "<i>This report was generated by the GENUIT nf-sarek-clinical pipeline (nf-core/sarek v3.x). "
        "All variant interpretations are for investigational use. Results must be confirmed in a "
        "certified clinical laboratory before clinical decision-making. "
        "References: GATK Best Practices; ACMG/AMP 2015 (Richards et al.); "
        "CPIC guidelines (cpicpgx.org); COSMIC v3.4 (cancer.sanger.ac.uk); ClinicalTrials.gov.</i>"
        "</p>"
    )

    return method_css + overview_table + details_html + disclaimer


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

def tsv_to_html_germline_cnv(path, _panel_genes):
    # One row per exon -> group by (CHROM, POS, END). Only panel genes, QA >= MIN_QA.
    # Low-mappability segments go to a separate table, never to the main one.
    MIN_QA = int(os.environ.get('FILT_GERMLINE_CNV_MIN_QA','4'))
    LOWMAP_FRAC = float(os.environ.get('FILT_GERMLINE_CNV_LOWMAP_FRAC','0.5'))
    if not path or not os.path.isfile(path):
        return "<p class='warn'>Germline CNV (GATK) file not found: " + html.escape(path or "N/A") + "</p>"
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = list(reader)
    if not rows:
        return "<p class='warn'>Germline CNV (GATK): empty file</p>"
    # Group rows by segment
    segments = {}
    seg_order = []
    for row in rows:
        cn_raw = row.get('CN', '2').strip()
        if not germline_cnv_is_event(row):
            continue
        key = (row.get('CHROM',''), row.get('POS',''), row.get('END',''))
        if key not in segments:
            segments[key] = {
                'chrom': row.get('CHROM',''),
                'pos':   row.get('POS',''),
                'end':   row.get('END',''),
                'alt':   row.get('ALT','.'),
                'qual':  row.get('QUAL','.'),
                'filter':row.get('FILTER','.'),
                'gt':    row.get('GT','.'),
                'cn':    cn_raw,
                'np':    row.get('NP','.'),
                'qa':    row.get('QA','0'),
                'qs':    row.get('QS','.'),
                'lowmap': germline_cnv_lowmap(row),
                'genes': [],
            }
            seg_order.append(key)
        gene = row.get('gen','').strip()
        if gene and gene not in segments[key]['genes']:
            segments[key]['genes'].append(gene)
    # Panel genes only, QA >= MIN_QA
    filtered = []
    for key in seg_order:
        seg = segments[key]
        try:
            qa_val = int(float(seg.get('qa') or ''))
        except (ValueError, TypeError, OverflowError):
            continue  # non-numeric QA: treated as below threshold
        if qa_val < MIN_QA:
            continue
        panel_hits = [g for g in seg['genes'] if g in _panel_genes]
        if panel_hits:
            seg['genes'] = panel_hits
            filtered.append(seg)
    lowmap_segs = [s for s in filtered if s['lowmap'] is not None and s['lowmap'] >= LOWMAP_FRAC]
    filtered = [s for s in filtered if s not in lowmap_segs]
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

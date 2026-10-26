#!/usr/bin/env python3
"""
render_html.py — laboratory HTML report of one pair, rendered from the report
model (bin/report_model.py). Shows the whole model: every germline class, every
AMP tier, Tier IV and TSO500 apart, all LOH events (those the clinical PDF lists
marked), DBS signatures, raw TMB / HRD tables. No classification or threshold
lives here.

    render_html.py --model PAIR.report_model.json --out PAIR_molecular_report.html
"""

import argparse
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from acmg_shared import clinvar_display, gnomad_display, hgvs_display  # noqa: E402
from report_shared import fmt_num as _fmt_num                           # noqa: E402
import amp_tier as _amp                                                 # noqa: E402
from final_report_sections import (tsv_to_html_pharmaco, clinical_trials_html,  # noqa: E402
                                   build_methodology_html, germline_cnv_html)


def _conf_badge(c):
    # CNV call confidence: HIGH blue, MODERATE grey (the same colours in the PDF)
    bg = '#1d4ed8' if c == 'HIGH' else '#6b7280'
    return ("<span class='badge' style='font-size:10px;background:" + bg + ";color:#fff'>"
            + ('HIGH' if c == 'HIGH' else 'MODERATE') + "</span>")

def _tools_html(t):
    out = ''
    for tool in t.split('+'):
        out += "<span class='badge badge-VUS' style='font-size:10px'>" + html.escape(tool.strip()) + "</span> "
    return out

def _loh_badge(lt):
    colors = {'HETDEL+LOH': 'badge-P', 'cnLOH': 'badge-LP'}
    cls = colors.get(lt, 'badge-III')
    return "<span class='badge " + cls + "' style='font-size:10px'>" + html.escape(lt) + "</span>"

def build_unified_cnv_html(data, amp_threshold):
    # data: report model 'somatic_cnv' (cnv_shared.harmonize_cnv); loh_relevant: genes
    # with a germline P/LP or somatic variant, whose LOH the clinical PDF lists
    loh_relevant = set(data['loh_relevant'])
    out  = "<div class='cnv-unified'>"

    # HOMDEL
    out += ("<h3 style='color:#ef4444;margin-top:16px'>"
            "Homozygous Deletions (CN = 0)</h3>")
    if data['homdel']:
        out += ("<table><thead><tr>"
                "<th>Gene</th><th>Chr</th>"
                "<th>CN ASCAT/CNVkit</th><th>CN GATK</th>"
                "<th>Tools</th><th>Confidence</th></tr></thead><tbody>")
        for r in data['homdel']:
            out += ("<tr>"
                    "<td><b>" + html.escape(r['gene']) + "</b></td>"
                    "<td>" + html.escape(r['chrom']) + "</td>"
                    "<td style='color:#ef4444;font-weight:700'>0</td>"
                    "<td>" + (str(r['cn_s']) if r['cn_s'] != '' else '&mdash;') + "</td>"
                    "<td>" + _tools_html(r['tools']) + "</td>"
                    "<td>" + _conf_badge(r['conf']) + "</td>"
                    "</tr>")
        out += "</tbody></table>"
    else:
        out += "<p class='warn'>No homozygous deletions detected in the somatic panel.</p>"

    # AMP
    out += ("<h3 style='color:#f97316;margin-top:16px'>"
            "Amplifications (CN &ge; " + str(amp_threshold) + ")</h3>")
    if data['amp']:
        amp_sorted = sorted(data['amp'], key=lambda r: -(r['cn_p'] if isinstance(r['cn_p'], int) else 0))
        out += ("<table><thead><tr>"
                "<th>Gene</th><th>Chr</th>"
                "<th>CN ASCAT/CNVkit</th><th>CN GATK</th>"
                "<th>Tools</th><th>Confidence</th></tr></thead><tbody>")
        for r in amp_sorted:
            out += ("<tr>"
                    "<td><b>" + html.escape(r['gene']) + "</b></td>"
                    "<td>" + html.escape(r['chrom']) + "</td>"
                    "<td style='color:#f97316;font-weight:700'>" + str(r['cn_p']) + "</td>"
                    "<td>" + (str(r['cn_s']) if r['cn_s'] != '' else '&mdash;') + "</td>"
                    "<td>" + _tools_html(r['tools']) + "</td>"
                    "<td>" + _conf_badge(r['conf']) + "</td>"
                    "</tr>")
        out += "</tbody></table>"
    else:
        out += ("<p class='warn'>No amplifications (CN &ge; "
                + str(amp_threshold) + ") detected in the somatic panel.</p>")

    # LOH
    out += ("<h3 style='color:#eab308;margin-top:16px'>"
            "Loss of Heterozygosity (LOH) &mdash; Germline Panel</h3>")
    if data['loh']:
        out += ("<table><thead><tr>"
                "<th>Gene</th><th>Chr</th><th>LOH Type</th>"
                "<th>CN ASCAT</th><th>CN GATK</th>"
                "<th>Tools</th><th>Confidence</th><th>Clinical</th></tr></thead><tbody>")
        for r in data['loh']:
            _rel = loh_relevant is not None and r['gene'] in loh_relevant
            out += ("<tr>"
                    "<td><b>" + html.escape(r['gene']) + "</b></td>"
                    "<td>" + html.escape(r['chrom']) + "</td>"
                    "<td>" + _loh_badge(r['loh_type']) + "</td>"
                    "<td>" + (str(r['cn_p']) if r['cn_p'] != '' else '&mdash;') + "</td>"
                    "<td>" + (str(r['cn_s']) if r['cn_s'] != '' else '&mdash;') + "</td>"
                    "<td>" + _tools_html(r['tools']) + "</td>"
                    "<td>" + _conf_badge(r['conf']) + "</td>"
                    "<td>" + ("<span class='badge badge-IB' title='Germline P/LP or somatic variant in the gene: listed in the PDF'>2nd hit</span>"
                              if _rel else "&mdash;") + "</td>"
                    "</tr>")
        out += "</tbody></table>"
        out += ("<p style='font-size:.75rem;color:#94a3b8;'>Clinical = gene with a germline P/LP or a somatic "
                "variant (core panel or second-hit track): the LOH listed in the clinical PDF.</p>")
    else:
        out += "<p class='warn'>No LOH events detected in the germline panel.</p>"

    out += "</div>"
    return out


def hgvs_cells(get_field):
    hc, tx, hp, plus = hgvs_display(get_field)
    sub = lambda t: "<br><span class='tx'>" + html.escape(t) + "</span>"
    return ("<td>" + html.escape(hc) + (sub(tx) if tx else "") + "</td>",
            "<td>" + html.escape(hp) + "".join(sub(p) for p in plus) + "</td>")


def tsv_to_html(table):
    # table: [headers, rows as dicts] from the report model
    header, rows = table
    if not header:
        return "<p class='warn'>Empty file</p>"
    th = "".join("<th>" + html.escape(c if c not in ('', '0') else '#') + "</th>" for c in header)
    tbody = ""
    for row in rows:
        cells = [row.get(c) or '' for c in header]
        td = "".join("<td>" + html.escape(_fmt_num(v)) + "</td>" for v in cells)
        tbody += "<tr>" + td + "</tr>\n"
    return ("<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
            "<tbody>" + tbody + "</tbody></table></div>")

# ----
# Germline SNV table with ACMG classification
# ----
def tsv_to_html_germline(rows):
    # rows: report model 'germline' (every candidate, with its ACMG class and criteria)

    # Display column map: display_name -> raw_column_name (or None for computed)
    disp_cols = [
        ('Locus',       None),
        ('Gene',        'ANN[*].GENE'),
        ('HGVS_c',      'ANN[*].HGVS_C'),
        ('HGVS_p',      'ANN[*].HGVS_P'),
        ('Impact',      'ANN[*].IMPACT'),
        ('DP',          None),
        ('VAF',         None),
        ('Cig.',        'ZYGOSITY'),
        ('REVEL',       None),
        ('SpliceAI',    'SpliceAI_DS'),
        ('gnomAD',      None),
        ('ClinVar',     None),
        ('ACMG',        None),
        ('Criteria',    None),
    ]

    acmg_colors = {
        'P':   'badge-P',
        'LP':  'badge-LP',
        'VUS': 'badge-VUS',
        'LB':  'badge-LB',
        'B':   'badge-B',
    }

    def crit_badge(code):
        if code.startswith('ClinVar'):
            cls = 'badge-crit-cv'
        elif code.endswith('_NA'):
            cls = 'badge-crit-na'
        elif code.startswith('BP') or code.startswith('BA') or code.startswith('BS'):
            cls = 'badge-crit-b'
        else:
            cls = 'badge-crit-p'
        return "<span class='badge-crit " + cls + "'>" + html.escape(code) + "</span>"

    def get_col(row, name):
        return (row.get(name) or '').strip()

    th = "".join("<th>" + html.escape(d) + "</th>" for d, _ in disp_cols)

    filter_bar = (
        "<div class='filter-bar'>"
        "<span style='font-size:.78rem;color:rgba(226,232,240,0.7);align-self:center;'>Germline filter:</span>"
        "<button class='filter-btn germ-filter-btn active' data-filter='all' onclick='filterGermline(this.dataset.filter)'>All</button>"
        "<button class='filter-btn germ-filter-btn' data-filter='plpvus' onclick='filterGermline(this.dataset.filter)'>P + LP + VUS</button>"
        "<button class='filter-btn germ-filter-btn' data-filter='plp' onclick='filterGermline(this.dataset.filter)'>P + LP only</button>"
        "</div>"
    )

    tbody = ""
    for row in rows:
        acmg, criteria = row['_ACMG'], row['_ACMG_CRITERIA']
        badge_cls = acmg_colors.get(acmg, 'badge-VUS')

        chrom = get_col(row, 'CHROM')
        pos   = get_col(row, 'POS')
        locus = chrom + ':' + pos if chrom and pos else (chrom or pos or '.')

        dp_raw  = get_col(row, 'GEN[0].DP')
        ad_raw  = get_col(row, 'GEN[0].AD')
        try:
            dp_val = int(dp_raw)
            alt_reads = int(ad_raw.split(',')[1])
            vaf_val = "{:.2f}".format(alt_reads / dp_val) if dp_val > 0 else '.'
        except Exception:
            dp_val  = dp_raw or '.'
            vaf_val = '.'

        tds = ""
        for disp_name, raw_col in disp_cols:
            if disp_name == 'Locus':
                tds += "<td>" + html.escape(locus) + "</td>"
            elif disp_name == 'ACMG':
                tds += "<td><span class='badge " + badge_cls + "'>" + html.escape(acmg) + "</span></td>"
            elif disp_name == 'Criteria':
                tds += "<td>" + "".join(crit_badge(c) for c in criteria) + "</td>"
            elif disp_name == 'DP':
                tds += "<td>" + html.escape(str(dp_val)) + "</td>"
            elif disp_name == 'HGVS_c':
                tds += hgvs_cells(lambda k: get_col(row, k))[0]
            elif disp_name == 'HGVS_p':
                tds += hgvs_cells(lambda k: get_col(row, k))[1]
            elif disp_name == 'ClinVar':
                tds += "<td>" + html.escape(clinvar_display(lambda k: get_col(row, k))) + "</td>"
            elif disp_name == 'REVEL':
                rv = get_col(row, 'REVEL_score')
                tds += "<td>" + html.escape(rv) + "</td>"
            elif disp_name == 'gnomAD':
                tds += "<td>" + html.escape(gnomad_display(lambda k: get_col(row, k))) + "</td>"
            elif disp_name == 'VAF':
                # QC_FLAG (Step 2 of SNPSIFT_EXTRACT_GERMLINE) shown inside the VAF cell
                qc = get_col(row, 'QC_FLAG')
                vaf_txt = str(vaf_val) + (" (" + qc + ")" if qc not in ('', '.') else "")
                tds += "<td>" + html.escape(vaf_txt) + "</td>"
            else:
                val = get_col(row, raw_col) if raw_col else ''
                tds += "<td>" + html.escape(val) + "</td>"

        tbody += "<tr data-acmg='" + acmg + "'>" + tds + "</tr>\n"

    return (
        filter_bar +
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody id='tbl-germline'>" + tbody + "</tbody></table></div>"
    )

# ----
# Somatic SNV table with AMP tier classification
# ----
def tsv_to_html_somatic(m):
    # Report model rows: main table (core panel + TSO500 Tier I/II), then Tier IV
    # (population polymorphisms) and the TSO500 Tier III apart, as classified there
    all_rows = ([(v['_TRACK'], v) for v in m['somatic']] +
                [(v['_TRACK'], v) for v in m['somatic_tier_iv']] +
                [(v['_TRACK'], v) for v in m['somatic_tso500_vus']])
    disp_cols = [
        ('Locus',        None),
        ('Gene',         'ANN[*].GENE'),
        ('HGVS_c',       'ANN[*].HGVS_C'),
        ('HGVS_p',       'ANN[*].HGVS_P'),
        ('Effect',       'ANN[*].EFFECT'),
        ('VAF',          None),
        ('gnomAD',       None),
        ('Call',         None),
        ('CIViC Level',  None),
        ('Disease',      None),
        ('Therapy',      None),
        ('AMP Tier',     None),
    ]

    tier_colors = {
        'IA':  'badge-IA',
        'IB':  'badge-IB',
        'IIC': 'badge-IIC',
        'IID': 'badge-IIC',
        'III': 'badge-III',
        'IV':  'badge-IV',
    }

    def get_col(row, name):
        return (row.get(name) or '').strip()

    th = "".join("<th>" + html.escape(d) + "</th>" for d, _ in disp_cols)

    filter_bar = (
        "<div class='filter-bar'>"
        "<span style='font-size:.78rem;color:rgba(226,232,240,0.7);'>Somatic:</span>"
        "<button class='filter-btn som-filter-btn active' data-filter='all' onclick='filterSomatic(this.dataset.filter)'>All</button>"
        "<button class='filter-btn som-filter-btn' data-filter='ti-tii' onclick='filterSomatic(this.dataset.filter)'>Tier I+II</button>"
        "<button class='filter-btn som-filter-btn' data-filter='ti' onclick='filterSomatic(this.dataset.filter)'>Tier I</button>"
        "</div>"
    )

    tbody = ""
    tbody_pop = ""
    tbody_tso = ""
    for track, row in all_rows:
        tier, civic = row['_AMP_TIER'], row['_CIVIC'] or {}
        track = 'tso500' if track == 'TSO500' else track
        row_start = len(tbody)
        badge_cls = tier_colors.get(tier, 'badge-III')

        chrom = get_col(row, 'CHROM')
        pos   = get_col(row, 'POS')
        locus = chrom + ':' + pos if chrom and pos else (chrom or pos or '.')

        civic_level   = civic.get('level', '')
        civic_disease = civic.get('disease', '')
        civic_therapy = _amp.therapy_cell(civic)
        summ = row['_EVIDENCE']
        civic_p = ''.join("<p class='ctx-p'><strong>" + html.escape(lbl) + ":</strong> " + html.escape('; '.join(it)) + "</p>"
                          for lbl, it in summ['groups'])
        if summ['not_supporting']:
            civic_p += "<p class='ctx-p'><strong>Evidence against:</strong> " + html.escape('; '.join(summ['not_supporting'])) + "</p>"
        if summ['n_submitted']:
            civic_p += "<p class='ctx-p' style='color:#888;'>" + str(summ['n_submitted']) + " submitted (unreviewed) CIViC record(s) not used for the tier.</p>"

        ctx = row['_CTX']
        ctx_role    = ctx.get('role', '')
        ctx_tier    = ctx.get('tier', '')
        ctx_pathway = ctx.get('pathway', '')
        ctx_summary = ctx.get('summary', '')
        conseq_text = row['_CONSEQUENCE']

        vaf_raw = get_col(row, 'GEN[1].AF')
        try:
            vaf_val = float(vaf_raw.split(',')[0])
            vaf_str = '{:.3f}'.format(vaf_val)
            low_vaf = vaf_val < 0.05
        except (ValueError, IndexError):
            vaf_str = vaf_raw or '.'
            low_vaf = False
        callers_raw = get_col(row, 'Callers') or 'Mutect2'
        callers_html = ("<span style='color:#166534;font-weight:bold;'>M+S</span>"
            if 'Strelka2' in callers_raw else "<span style='color:#92400e;font-weight:bold;'>M</span>")

        tds = ""
        for disp_name, raw_col in disp_cols:
            if disp_name == 'Locus':
                toggle_btn = ("<button class='ctx-toggle' onclick='toggleCtx(this)' "
                              "aria-label='Toggle gene context'>&#9654;</button> ")
                tds += "<td style='white-space:nowrap;'>" + toggle_btn + html.escape(locus) + "</td>"
            elif disp_name == 'VAF':
                if low_vaf:
                    tds += "<td style='color:#b91c1c;font-weight:bold;'>" + html.escape(vaf_str) + " &#9888;</td>"
                else:
                    tds += "<td>" + html.escape(vaf_str) + "</td>"
            elif disp_name == 'HGVS_c':
                tds += hgvs_cells(lambda k: get_col(row, k))[0]
            elif disp_name == 'HGVS_p':
                tds += hgvs_cells(lambda k: get_col(row, k))[1]
            elif disp_name == 'Gene' and track == 'tso500':
                tds += ("<td>" + html.escape(get_col(row, raw_col)) +
                        "<br><span class='tx'>TSO500</span></td>")
            elif disp_name == 'gnomAD':
                tds += "<td>" + html.escape(gnomad_display(lambda k: get_col(row, k))) + "</td>"
            elif disp_name == 'Call':
                tds += "<td style='text-align:center;'>" + callers_html + "</td>"
            elif disp_name == 'CIViC Level':
                tds += "<td>" + html.escape(civic_level) + "</td>"
            elif disp_name == 'Disease':
                tds += "<td>" + html.escape(civic_disease) + "</td>"
            elif disp_name == 'Therapy':
                tds += "<td>" + html.escape(civic_therapy) + "</td>"
            elif disp_name == 'AMP Tier':
                tds += "<td><span class='badge " + badge_cls + "'>" + html.escape(tier) + "</span></td>"
            else:
                val = get_col(row, raw_col) if raw_col else ''
                tds += "<td>" + html.escape(val) + "</td>"

        tbody += "<tr data-tier='" + tier + "'>" + tds + "</tr>\n"

        # Hidden context row
        if ctx_role or ctx_summary or civic_p:
            role_badge = "<span class='badge badge-role'>" + html.escape(ctx_role or 'Unknown') + "</span>"
            path_badge = "<span class='badge badge-path'>" + html.escape(ctx_pathway) + "</span>" if ctx_pathway else ""
            tier_badge = ("<span class='badge badge-cgc'>CGC Tier " + html.escape(ctx_tier) + "</span>") if ctx_tier else ""
            summary_p  = ("<p class='ctx-p'><strong>Gene function:</strong> " + html.escape(ctx_summary) + "</p>") if ctx_summary else ""
            conseq_p   = ("<p class='ctx-p'><strong>Tumor consequence:</strong> " + html.escape(conseq_text) + "</p>") if conseq_text else ""
            ctx_inner  = ("<div class='ctx-box'>"
                          "<div class='ctx-badges'>" + role_badge + path_badge + tier_badge + "</div>"
                          + summary_p + conseq_p + civic_p
                          + "</div>")
            tbody += ("<tr class='ctx-row' data-tier='" + tier + "' style='display:none;'>"
                      "<td colspan='12'>" + ctx_inner + "</td></tr>\n")
        else:
            tbody += ("<tr class='ctx-row' data-tier='" + tier + "' style='display:none;'>"
                      "<td colspan='12'><em style='color:#888;'>No annotation.</em></td></tr>\n")
        # AMP Tier IV (common population polymorphism): own table, not a somatic finding
        if tier == 'IV':
            tbody_pop += tbody[row_start:]
            tbody = tbody[:row_start]
        elif track == 'tso500' and tier == 'III':
            tbody_tso += tbody[row_start:]
            tbody = tbody[:row_start]

    TUMOR_LABEL = m['tumor_type']['label']
    tumor_note = ("<p style='font-size:.8rem;'><b>Tumor type for tiering:</b> " + html.escape(TUMOR_LABEL) + "</p>" if TUMOR_LABEL else
                  "<p class='warn'>Tumor type not provided (--tumor_type): CIViC evidence is treated as from another tumor type, so the highest tier is IIC.</p>")
    pop_html = ""
    if tbody_pop:
        pop_html = ("<h3 style='font-size:.95rem;margin-top:18px;'>Tier IV &mdash; population polymorphisms "
                    "(gnomAD &ge; " + str(m['thresholds']['somatic_max_pop_af']) + "; probable germline origin, not interpreted as somatic)</h3>"
                    "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
                    "<tbody id='tbl-somatic-pop'>" + tbody_pop + "</tbody></table></div>")
    tert_note = m['alerts']['tert_promoter']
    tert_html = ("<p class='warn' style='margin-top:14px;'>" + html.escape(tert_note) + "</p>") if tert_note else ""
    tso_html = ""
    if tbody_tso:
        n_tso = tbody_tso.count("<tr data-tier=")
        tso_html = ("<details style='margin-top:18px;'><summary style='cursor:pointer;font-size:.95rem;font-weight:600;'>"
                    "Pan-tumor (TSO500) &mdash; " + str(n_tso) + " variant(s) of uncertain significance "
                    "(AMP Tier III, TruSight Oncology 500 genes outside the core panel)</summary>"
                    "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
                    "<tbody id='tbl-somatic-tso500'>" + tbody_tso + "</tbody></table></div></details>")
    return (
        tumor_note + filter_bar +
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody id='tbl-somatic'>" + tbody + "</tbody></table></div>" + pop_html + tso_html + tert_html
    )

# ----
# Signature table helper
# ----

_SIG_COLOURS = {'clock': '#4e79a7', 'apobec': '#f28e2b', 'hrd': '#e15759', 'mmr': '#b07aa1',
                'pol': '#ff9da7', 'exposure': '#9c755f', 'treatment': '#59a14f', 'other': '#76b7b2',
                'unknown': '#888888'}


def sig_to_html(sig):
    # sig: report model signatures[SBS96|DBS|ID] (assigned signatures, COSMIC aetiology)
    if not sig or not sig['assigned']:
        return "<p class='warn'>No signatures assigned (all zero)</p>"
    th = "<th>Signature</th><th>Mutations Assigned</th><th>Proportion (%)</th><th>COSMIC Aetiology</th>"
    tbody = ""
    for s in sig['assigned']:
        badge = ("<span style='background:" + _SIG_COLOURS.get(s['category'], '#888888') + ";color:#fff;padding:2px 8px;"
                 "border-radius:4px;font-size:.72rem;'>" + html.escape(s['aetiology']) + "</span>")
        tbody += ("<tr><td><strong>" + html.escape(s['signature']) + "</strong></td>"
                  "<td>" + str(s['count']) + "</td><td>" + "{:.1f}".format(s['percent']) + "%</td>"
                  "<td>" + badge + "</td></tr>\n")
    return ("<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
            "<tbody>" + tbody + "</tbody></table></div>")

def hrd_status_color(v):
    if v.upper().startswith("HRD+"):  return "#ef4444"
    if v.upper().startswith("HRD-"):  return "#22c55e"
    return "#f97316"

# ----
# Signature plots: prefer new custom PNGs; fall back to SigProfiler PDFs
# ----
def sig_plots_to_html(act_path, sig_type):
    # act_path: path to Assignment_Solution_Activities.txt
    # sig_type: 'SBS96', 'DBS', or 'ID'
    if not act_path or not os.path.isfile(act_path):
        return ""
    import base64 as _b64m
    act_dir = os.path.dirname(act_path)

    _img_style = "max-width:100%;border-radius:8px;border:1px solid rgba(148,163,184,.2);margin-top:4px;"
    _lbl_style = "font-size:.78rem;color:rgba(226,232,240,.65);margin-top:10px;margin-bottom:2px;"
    parts = []

    def png_b64(png_path):
        if not png_path or not os.path.isfile(png_path):
            return None
        try:
            with open(png_path, "rb") as _f:
                return _b64m.b64encode(_f.read()).decode()
        except Exception:
            return None

    def pdf_b64(pdf_path, dpi=150):
        if not pdf_path or not os.path.isfile(pdf_path):
            return None
        try:
            import fitz as _fitz
            doc = _fitz.open(pdf_path)
            if doc.page_count == 0:
                return None
            mat = _fitz.Matrix(dpi / 72.0, dpi / 72.0)
            pix = doc[0].get_pixmap(matrix=mat)
            b64 = _b64m.b64encode(pix.tobytes("png")).decode()
            doc.close()
            return b64
        except Exception:
            return None

    # ── 1. Matrix profile plot (SBS96 or ID83) ────────────────────────────
    # New custom PNG generated by generate_signature_plots.py
    profile_png_name = "SBS96_matrix_profile.png" if sig_type == "SBS96" else "ID83_matrix_profile.png"
    b64_profile = png_b64(os.path.join(act_dir, profile_png_name))
    if b64_profile:
        label = "SBS96 Mutational Spectrum" if sig_type == "SBS96" else "ID83 Indel Spectrum"
        parts.append("<div><div style='" + _lbl_style + "'>" + label + " (reconstructed)</div>"
                     + "<img src='data:image/png;base64," + b64_profile + "' style='" + _img_style + "'/></div>")
    elif sig_type == "SBS96":
        # Fall back to SigProfiler reconstruction PDF
        b64_recon = pdf_b64(os.path.join(act_dir, "SampleReconstruction", "Reconstructed_Sample_Plots_96.pdf"))
        if b64_recon:
            parts.append("<div><div style='" + _lbl_style + "'>SBS96 Mutational Spectrum &amp; Reconstruction</div>"
                         + "<img src='data:image/png;base64," + b64_recon + "' style='" + _img_style + "'/></div>")

    # ── 2. Etiology-annotated activity bar chart ──────────────────────────
    # New custom PNG generated by generate_signature_plots.py
    b64_etiol = png_b64(os.path.join(act_dir, "etiology_activity_plot.png"))
    if b64_etiol:
        parts.append("<div><div style='" + _lbl_style + "'>COSMIC Signature Assignment — Aetiology (" + sig_type + ")</div>"
                     + "<img src='data:image/png;base64," + b64_etiol + "' style='" + _img_style + "max-width:800px;'/></div>")
    else:
        # Fall back to SigProfiler activity PDF
        b64_act = pdf_b64(os.path.join(act_dir, "Assignment_Solution_Activity_Plots.pdf"))
        if b64_act:
            parts.append("<div style='max-width:460px;'><div style='" + _lbl_style + "'>COSMIC Signature Assignment (" + sig_type + ")</div>"
                         + "<img src='data:image/png;base64," + b64_act + "' style='" + _img_style + "'/></div>")

    if not parts:
        return ""
    return "<div style='display:flex;flex-direction:column;gap:6px;margin-top:8px;'>" + "".join(parts) + "</div>"

# ----
# Precompute KPI cards
# ----

def build_coverage_quality_html(cov, normal_name, tumor_name):
    # Return HTML for the per-gene coverage quality section.

    def cov_badge(val_str):
        try:
            v = float(val_str)
        except (ValueError, TypeError):
            return "<span style='color:#888;'>NA</span>"
        if v >= 100:
            color = "#22c55e"
        elif v >= 50:
            color = "#eab308"
        elif v >= 20:
            color = "#f97316"
        else:
            color = "#ef4444"
        return "<span style='color:" + color + ";font-weight:700;'>" + html.escape(f'{v:.1f}') + "</span>"

    cov_n = {g: {'mean_coverage': r['mean'], 'min_coverage': r['min']} for g, r in cov['normal'].items()}
    cov_t = {g: {'mean_coverage': r['mean'], 'min_coverage': r['min']} for g, r in cov['tumor'].items()}
    all_genes = sorted(set(list(cov_n.keys()) + list(cov_t.keys())))

    if not all_genes:
        return "<p class='warn'>Gene coverage files not found.</p>"

    th = (
        "<th>Gene</th>"
        "<th title='Normal mean coverage'>Normal (X)</th>"
        "<th title='Normal min coverage'>N min</th>"
        "<th title='Tumor mean coverage'>Tumor (X)</th>"
        "<th title='Tumor min coverage'>T min</th>"
    )

    tbody = ""
    for gene in all_genes:
        rn = cov_n.get(gene, {})
        rt = cov_t.get(gene, {})
        mn  = rn.get('mean_coverage', 'NA')
        minn = rn.get('min_coverage', 'NA')
        mt  = rt.get('mean_coverage', 'NA')
        mint = rt.get('min_coverage', 'NA')
        tbody += (
            "<tr>"
            "<td style='font-weight:600;font-family:monospace;'>" + html.escape(gene) + "</td>"
            "<td style='text-align:right;'>" + cov_badge(mn) + "</td>"
            "<td style='text-align:right;color:#94a3b8;font-size:.75rem;'>" + html.escape(str(minn)) + "</td>"
            "<td style='text-align:right;'>" + cov_badge(mt) + "</td>"
            "<td style='text-align:right;color:#94a3b8;font-size:.75rem;'>" + html.escape(str(mint)) + "</td>"
            "</tr>\n"
        )

    legend = (
        "<div style='display:flex;gap:16px;font-size:.75rem;margin-bottom:8px;align-items:center;'>"
        "<span style='color:#22c55e;font-weight:700;'>&#9632;</span> &ge;100X &nbsp;"
        "<span style='color:#eab308;font-weight:700;'>&#9632;</span> 50&#8211;100X &nbsp;"
        "<span style='color:#f97316;font-weight:700;'>&#9632;</span> 20&#8211;50X &nbsp;"
        "<span style='color:#ef4444;font-weight:700;'>&#9632;</span> &lt;20X"
        "</div>"
    )
    sample_hdr = (
        "<div style='font-size:.8rem;margin-bottom:6px;'>"
        "<strong>Normal:</strong> " + html.escape(normal_name) + " &nbsp;|&nbsp; "
        "<strong>Tumor:</strong> " + html.escape(tumor_name)
        + "</div>"
    )

    return (
        sample_hdr + legend
        + "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        + "<tbody>" + tbody + "</tbody></table></div>"
    )


def kpis(m):
    """Dashboard values (report model metrics; 'Not run' / 'Indeterminate' as in the PDF)."""
    mt, nr = m['metrics'], set(m['not_run'])
    indet = 'Indeterminate (tumor coverage < ' + format(m['thresholds']['low_tumor_coverage'], 'g') + 'X)'
    low = mt['low_tumor_coverage']
    st = mt['msi_status']
    msi = 'NA' if st == 'NA' else {'MSI-H': 'MSI detected', 'MSI borderline': 'MSI borderline'}.get(
        st, 'MSI not detected (MSS)') + " \u2014 " + mt['msi'] + "%"
    nrun = lambda k, v: 'Not run' if k in nr else v
    return {
        "TMB (mut/Mb)":        indet if low else nrun('tmb', mt['tmb']),
        "MSI":                 indet if low else msi,
        "LOH":                 mt['loh'],
        "Genomic Instability": mt['gi'],
        "WGD":                 mt['wgd'],
        "Ploidy":              mt['ploidy'] + (" (WGD)" if mt['wgd'].strip() in ("1", "True", "TRUE") else ""),
        "Tumor Purity (%)":    mt['purity'],
        "Coverage Normal (X)": mt['cov_n'],
        "Coverage Tumor (X)":  mt['cov_t'],
        "HRD Probability":     indet if low else nrun('hrd', mt['hrd_prob_display']),
        "HRD Status":          indet if low else nrun('hrd', mt['hrd_status']),
    }


def unavail(m, key):
    """'<p class=warn>why there is no data</p>', or '' when the analysis has a result."""
    t = m['unavailable'].get(key, '')
    return ("<p class='warn'>" + html.escape(t) + "</p>") if t else ''


def render(m):
    ids = m['ids']
    pair, tumor, normal, patient_id, run_name = (ids['pair'], ids['tumor'], ids['normal'],
                                                 ids['patient'], ids['run_name'])
    sig = m['signatures']
    kpi_cards = ""
    for k, v in kpis(m).items():
        color = hrd_status_color(v) if k == "HRD Status" else "var(--accent)"
        kpi_cards += ("<div class='kpi'><div class='kpi-label'>" + html.escape(k) + "</div>"
                      "<div class='kpi-value' style='color:" + color + "'>" + html.escape(v) + "</div></div>")

    sections = ""
    sections += "<section><h2>Germline SNVs</h2>" + (unavail(m, 'germline_snv') or tsv_to_html_germline(m['germline'])) + "</section>\n"
    sections += "<section><h2>Somatic SNVs</h2>" + (unavail(m, 'somatic_snv') or tsv_to_html_somatic(m)) + "</section>\n"
    sections += "<section><h2>TMB</h2>" + (unavail(m, 'tmb') or tsv_to_html(m['tables']['tmb'])) + "</section>\n"
    sections += "<section><h2>HRD Predictions</h2>" + (unavail(m, 'hrd') or tsv_to_html(m['tables']['hrd'])) + "</section>\n"
    sections += ("<section><h2>Pharmacogenomics (PGx) &#8212; CPIC Dose Recommendations</h2>"
                 + (unavail(m, 'pharmaco') or tsv_to_html_pharmaco(m['pharmaco'])) + "</section>\n")
    sections += ("<section><h2>Sequencing Quality &#8212; Per-Gene Coverage</h2>"
                 + (unavail(m, 'coverage') or build_coverage_quality_html(m['coverage'], normal, tumor))
                 + "</section>\n")
    sbs = sig.get('SBS96')
    sections += ("<section><h2>Mutational Signatures &#8212; SBS96</h2>"
                 + (unavail(m, 'signatures') or (sig_to_html(sbs) + sig_plots_to_html((sbs or {}).get('path', ''), 'SBS96')))
                 + "</section>\n")
    for kind in ('DBS', 'ID'):
        if sig.get(kind):
            sections += ("<section><h2>Mutational Signatures &#8212; " + kind + "</h2>"
                         + sig_to_html(sig[kind]) + sig_plots_to_html(sig[kind]['path'], kind) + "</section>\n")
    sections += '<section><h2>Germline&#8211;Somatic Correlation</h2>' + m['correlation_html'] + '</section>\n'
    sections += ("<section><h2>Germline CNVs (GATK) &#8212; deletions / duplications</h2>"
                 + (unavail(m, 'germline_cnv') or germline_cnv_html(m['germline_cnv'], m['thresholds'])) + "</section>\n")
    _ploidy_txt = m['alerts']['ploidy_note']
    sections += ("<section><h2>Somatic CNVs &#8212; Reportable Events</h2>"
                 + (("<p class='warn'>" + html.escape(_ploidy_txt) + "</p>") if _ploidy_txt else "")
                 + (unavail(m, 'somatic_cnv') or build_unified_cnv_html(m['somatic_cnv'], m['thresholds']['somatic_cnv_amp_cn']))
                 + "</section>\n")
    sections += ("<section><h2>Clinical Trials &#8212; Relevant Recruiting Studies</h2>"
                 + clinical_trials_html(m['trials']) + "</section>\n")
    sections += "<section><h2>Test Methodology</h2>" + build_methodology_html(m) + "</section>\n"
    _ascat_alerts = m['alerts']['ascat']
    # Patient / specimen identification (samplesheet extra columns, report model)
    idn = m.get('identification') or {}
    ident_html = ''
    if idn.get('source'):
        pt, tu, no = idn['patient'], idn['tumor'], idn['normal']
        nv = lambda v: html.escape(v) if v else "<span style='opacity:.5'>Not provided</span>"
        spec = lambda d: ' · '.join(x for x in (d.get('specimen_type'), d.get('specimen_site'), d.get('accession_id')) if x)
        when = lambda d: ' / '.join(x or '—' for x in (d.get('collection_date'), d.get('received_date')))
        items = [('Patient name', pt.get('patient_name')), ('Birth date', pt.get('birth_date')),
                 ('Sex', pt.get('sex')), ('Medical record', pt.get('medical_record')),
                 ('Diagnosis', pt.get('diagnosis')), ('Ordering physician', pt.get('ordering_physician')),
                 ('Institution', pt.get('institution')),
                 ('Tumor specimen', spec(tu)), ('Tumor collected / received', when(tu)),
                 ('Tumor content (pathology)', (tu.get('tumor_content') + '%') if tu.get('tumor_content') else ''),
                 ('Normal specimen', spec(no)), ('Normal collected / received', when(no))]
        ident_html = ''.join("    <div><strong>" + k + ":</strong> " + nv(v) + "</div>\n" for k, v in items)
    ident_html += ''.join("    <div style='color:#f59e0b'><strong>Samplesheet:</strong> " + html.escape(w) + "</div>\n"
                          for w in idn.get('warnings', []))
    _pipeline_label = lambda: m['versions']['pipeline']
    _sarek_label = lambda: m['versions']['sarek']
    # ----
    # CSS
    # ----
    _css = (
        ':root { --bg:#020617; --fg:#f8fafc; --accent:#f97316; --line:rgba(148,163,184,0.2); }'
        'body { background:var(--bg); color:var(--fg); font-family:system-ui,sans-serif; margin:0; padding:0; }'
        '.container { max-width:1180px; margin:0 auto; padding:28px 18px 80px; }'
        'h1 { font-size:28px; } h2 { font-size:18px; margin-top:28px; border-bottom:1px solid var(--line); padding-bottom:6px; }'
        '.meta { display:grid; grid-template-columns:1fr 1fr; gap:8px 20px; background:rgba(15,23,42,0.6);'
        '        border:1px solid rgba(249,115,22,0.2); border-radius:12px; padding:14px 18px; margin:14px 0; }'
        '.kpi-grid { display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin:16px 0; }'
        '.kpi { background:rgba(15,23,42,0.7); border:1px solid var(--line); border-radius:12px; padding:14px 16px; }'
        '.kpi-label { font-size:11px; color:rgba(226,232,240,0.7); text-transform:uppercase; letter-spacing:.08em; }'
        '.kpi-value { font-size:22px; font-weight:700; color:var(--accent); margin-top:4px; }'
        'section { margin-top:20px; }'
        '.tbl-wrap { overflow-x:auto; border-radius:10px; border:1px solid var(--line); margin-top:10px; }'
        'table { border-collapse:collapse; width:100%; font-size:.78rem; }'
        'thead { background:rgba(249,115,22,0.6); } th { padding:.4rem .6rem; text-align:left; color:#fff; white-space:nowrap; }'
        'td { padding:.35rem .6rem; border-bottom:1px solid var(--line); color:#e5e7eb; }'
        'tr:hover td { background:rgba(249,115,22,0.15); }'
        '.warn { color:#fed7aa; background:rgba(124,45,18,0.3); border-radius:8px; padding:8px 12px; font-size:.85rem; }'
        '.callout-info { background:#eff6ff; border-left:4px solid #3b82f6; border-radius:6px; padding:10px 14px; font-size:.85rem; color:#1e3a5f; margin:8px 0; }'
        '.badge-crit { display:inline-block; padding:1px 5px; border-radius:10px; font-size:.62rem; font-weight:600; margin:1px; white-space:nowrap; }'
        '.badge-crit-p { background:rgba(239,68,68,0.15); color:#f87171; border:1px solid rgba(239,68,68,0.35); }'
        '.badge-crit-b { background:rgba(59,130,246,0.15); color:#60a5fa; border:1px solid rgba(59,130,246,0.35); }'
        '.badge-crit-cv { background:rgba(139,92,246,0.15); color:#a78bfa; border:1px solid rgba(139,92,246,0.35); }'
        '.badge-crit-na { background:rgba(107,114,128,0.15); color:#9ca3af; border:1px solid rgba(107,114,128,0.35); }'
        '.badge { display:inline-block; padding:2px 8px; border-radius:12px; font-size:.7rem; font-weight:700; }'
        '.badge-P { background:#ef4444; color:#fff; }'
        '.badge-LP { background:#f97316; color:#fff; }'
        '.badge-VUS { background:#eab308; color:#111; }'
        '.badge-LB { background:#3b82f6; color:#fff; }'
        '.badge-B { background:#22c55e; color:#fff; }'
        '.tx { font-size:.68rem; color:#94a3b8; white-space:nowrap; }'
        '.badge-IA, .badge-IB { background:#dc2626; color:#fff; }'
        '.badge-IIC { background:#d97706; color:#fff; }'
        '.badge-III { background:#6b7280; color:#fff; }'
        '.badge-IV { background:#16a34a; color:#fff; }'
        '.filter-bar { display:flex; gap:8px; margin:10px 0; flex-wrap:wrap; }'
        '.filter-btn { padding:4px 14px; border-radius:20px; border:1px solid var(--line); background:rgba(15,23,42,0.7); color:#e5e7eb; cursor:pointer; font-size:.78rem; }'
        '.filter-btn.active { background:var(--accent); color:#fff; border-color:var(--accent); }'
        '.ctx-toggle { background:none; border:none; color:rgba(249,115,22,0.8); cursor:pointer; font-size:.7rem; padding:0 4px 0 0; line-height:1; transition:transform .15s; }'
        '.ctx-toggle.open { transform:rotate(90deg); }'
        '.ctx-row td { background:rgba(15,23,42,0.9); padding:.5rem .8rem; }'
        '.ctx-box { border-left:3px solid rgba(249,115,22,0.5); padding:6px 12px; border-radius:4px; }'
        '.ctx-badges { display:flex; gap:6px; flex-wrap:wrap; margin-bottom:6px; }'
        '.badge-role { background:#7c3aed; color:#fff; }'
        '.badge-path { background:#0f766e; color:#fff; font-size:.65rem; }'
        '.badge-cgc  { background:#1e40af; color:#fff; font-size:.65rem; }'
        '.ctx-p { margin:4px 0; font-size:.76rem; color:#cbd5e1; line-height:1.45; }'
        '@media print {'
        '  .filter-bar { display:none !important; }'
        '  body { background:#fff !important; color:#111 !important; font-size:9pt; }'
        '  .kpi-grid { grid-template-columns:repeat(4,1fr) !important; }'
        '  .kpi { background:#f8fafc !important; border:1px solid #e2e8f0 !important; }'
        '  .kpi-value { color:#1e293b !important; font-size:16px !important; }'
        '  table { font-size:7pt !important; }'
        '  thead { background:#334155 !important; -webkit-print-color-adjust:exact; print-color-adjust:exact; }'
        '  tr[data-acmg]:not([data-acmg="P"]):not([data-acmg="LP"]):not([data-acmg="VUS"]) { display:none !important; }'
        '  tr[data-tier="IV"] { display:none !important; }'
        '  .ctx-row, .ctx-toggle { display:none !important; }'
        '  h2 { color:#1e293b !important; }'
        '}'
    )

    # ----
    # JavaScript (single-quoted strings throughout)
    # ----
    _js = (
        "<script>"
        "function filterGermline(show) {"
        "  document.querySelectorAll('#tbl-germline tr[data-acmg]').forEach(function(tr) {"
        "    var acmg = tr.getAttribute('data-acmg');"
        "    var visible = false;"
        "    if (show === 'all') visible = true;"
        "    else if (show === 'plpvus') visible = ['P','LP','VUS'].indexOf(acmg) >= 0;"
        "    else if (show === 'plp') visible = ['P','LP'].indexOf(acmg) >= 0;"
        "    tr.style.display = visible ? '' : 'none';"
        "  });"
        "  document.querySelectorAll('.germ-filter-btn').forEach(function(b) {"
        "    b.classList.toggle('active', b.dataset.filter === show);"
        "  });"
        "}"
        "function toggleCtx(btn) {"
        "  var tr = btn.closest('tr');"
        "  var ctxRow = tr ? tr.nextElementSibling : null;"
        "  if (!ctxRow || !ctxRow.classList.contains('ctx-row')) return;"
        "  var open = ctxRow.style.display !== 'none';"
        "  ctxRow.style.display = open ? 'none' : '';"
        "  btn.classList.toggle('open', !open);"
        "}"
        "function filterSomatic(show) {"
        "  document.querySelectorAll('#tbl-somatic tr[data-tier]').forEach(function(tr) {"
        "    if (tr.classList.contains('ctx-row')) return;"
        "    var tier = tr.getAttribute('data-tier');"
        "    var visible = false;"
        "    if (show === 'all') visible = true;"
        "    else if (show === 'ti-tii') visible = ['IA','IB','IIC','IID'].indexOf(tier) >= 0;"
        "    else if (show === 'ti') visible = ['IA','IB'].indexOf(tier) >= 0;"
        "    tr.style.display = visible ? '' : 'none';"
        "    // Always collapse ctx-row when filter changes"
        "    var ctxRow = tr.nextElementSibling;"
        "    if (ctxRow && ctxRow.classList.contains('ctx-row')) {"
        "      ctxRow.style.display = 'none';"
        "      var btn = tr.querySelector('.ctx-toggle');"
        "      if (btn) btn.classList.remove('open');"
        "    }"
        "  });"
        "  document.querySelectorAll('.som-filter-btn').forEach(function(b) {"
        "    b.classList.toggle('active', b.dataset.filter === show);"
        "  });"
        "}"
        "document.addEventListener('DOMContentLoaded', function() {"
        "  filterGermline('all');"
        "  filterSomatic('all');"
        "});"
        "</script>"
    )

    # ----
    # Assemble final page
    # ----
    page = (
        '<!doctype html>\n<html lang="en">\n<head>\n'
        '  <meta charset="utf-8"/>\n'
        '  <title>GENUIT Molecular Report \u2014 ' + html.escape(pair) + '</title>\n'
        '  <style>\n    ' + _css + '\n  </style>\n'
        '</head>\n<body>\n<div class="container">\n'
        '  <div style="font-size:11px;letter-spacing:.1em;text-transform:uppercase;opacity:.7;margin-bottom:8px;">GENUIT \u00b7 nf-sarek-clinical</div>\n'
        '  <h1>Molecular Report \u2014 ' + html.escape(patient_id) + '</h1>\n'
        '  <div class="meta">\n'
        '    <div><strong>Run:</strong> ' + html.escape(run_name) + '</div>\n'
        '    <div><strong>Patient:</strong> ' + html.escape(patient_id) + '</div>\n'
        '    <div><strong>Tumor:</strong> ' + html.escape(tumor) + '</div>\n'
        '    <div><strong>Normal:</strong> ' + html.escape(normal) + '</div>\n'
        '    <div><strong>Date:</strong> ' + ids['date'] + '</div>\n'
    + ident_html +
        '    <div><strong>Pipeline:</strong> ' + html.escape(_pipeline_label()) + ' (' + html.escape(_sarek_label()) + ')</div>\n'
        '  </div>\n'
        '  <h2>Dashboard</h2>\n'
        '  <div class="kpi-grid">\n    ' + kpi_cards + '\n  </div>\n'
        + ''.join("  <p class='warn'>&#9888; Copy-number fit quality (ASCAT): " + html.escape(a) + "</p>\n" for a in _ascat_alerts) +
        '  ' + sections + '\n'
        '  <p style="margin-top:40px;font-size:.75rem;opacity:.6;">Generated by nf-sarek-clinical \u00b7 GENUIT</p>\n'
        '</div>\n'
        + _js + '\n'
        '</body>\n</html>\n'
    )
    return page


def text_summary(m):
    """Plain-text summary (PAIR_molecular_report.txt)."""
    ids = m['ids']
    out = ['GENUIT Molecular Report', '======================',
           f"Run:     {ids['run_name']}", f"Pair:    {ids['pair']}", f"Tumor:   {ids['tumor']}",
           f"Normal:  {ids['normal']}", f"Patient: {ids['patient']}", f"Date:    {ids['date']}",
           f"Version: {m['versions']['pipeline']} ({m['versions']['sarek']})", '', 'DASHBOARD METRICS:']
    out += [f'  {k + ":":<20} {v}' for k, v in kpis(m).items()]
    out += ['', 'SECTIONS WITHOUT DATA:']
    out += [f'  {t}' for t in m['unavailable'].values() if t] or ['  none']
    return '\n'.join(out) + '\n'


def main():
    p = argparse.ArgumentParser(description='Laboratory HTML report from the report model.')
    p.add_argument('--model', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--txt', default='', help='plain-text summary')
    a = p.parse_args()
    with open(a.model) as fh:
        m = json.load(fh)
    with open(a.out, 'w') as fh:
        fh.write(render(m))
    if a.txt:
        with open(a.txt, 'w') as fh:
            fh.write(text_summary(m))
    print('HTML report written: ' + a.out, file=sys.stderr)


if __name__ == '__main__':
    main()

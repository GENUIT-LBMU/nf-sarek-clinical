/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINAL HTML REPORT MODULE
    =========================
    Aggregates all postprocessing outputs into a clinical HTML report.
    Auto-detects tumor/normal pair names from ASCAT metrics filename.

    Report sections:
      - KPI dashboard: TMB, MSI, LOH, GI, WGD, ploidy, purity, coverage, HRD
      - Gene coverage plots (interactive bar chart per gene)
      - Germline SNVs table (ACMG classification, filter buttons)
      - Germline CNVs table (GATK)
      - Somatic SNVs table + CIViC annotation (AMP tier, filter buttons)
      - Somatic CNVs (ASCAT/CNVkit) table
      - Somatic CNVs (GATK) table
      - Mutational signatures (SBS96, DBS, ID)
      - Pharmacogenomics table
      - HRD predictions table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FINAL_REPORT {
    tag "$pair"
    label 'process_single'

    publishDir "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(ascat_metrics), path(hrd_predictions), path(snv_germline_tsv), path(snv_somatic_tsv)
    val sarek_outdir
    val run_name

    output:
    tuple val(pair), path("${pair}_molecular_report.html"), emit: html
    tuple val(pair), path("${pair}_molecular_report.txt"),  emit: txt

    stub:
    """
    echo "<html><body><h1>Stub report for ${pair}</h1></body></html>" > ${pair}_molecular_report.html
    echo "Stub report for ${pair}" > ${pair}_molecular_report.txt
    """

    script:
    """
    # Parse pair name: TUMOR_vs_NORMAL
    PAIR="${pair}"
    if [[ "\${PAIR}" =~ ^(.+)_vs_(.+)\$ ]]; then
        TUMOR="\${BASH_REMATCH[1]}"
        NORMAL="\${BASH_REMATCH[2]}"
    else
        TUMOR="\${PAIR}"
        NORMAL="UNKNOWN"
    fi
    PATIENT_ID="\${TUMOR%%_T*}"

    OUT_DIR="${params.outdir}"
    RUN="${run_name}"

    # Collect input paths with fallbacks (|| true prevents pipefail on missing dirs)
    SNV_GERMLINE="\$(find \${OUT_DIR}/germline/snv/extract -name "\${NORMAL}*.germline_snvs.tsv" 2>/dev/null | head -1 || true)"
    SNV_SOMATIC="\$(find \${OUT_DIR}/somatic/snv/extract -name "\${PAIR}*.somatic_snvs.tsv" 2>/dev/null | head -1 || true)"
    CIVIC_TSV="\$(find \${OUT_DIR}/somatic/snv/civic_extract -name "\${PAIR}*.civic.exploded.tsv" 2>/dev/null | head -1 || true)"
    CNV_SOMATIC="\$(find \${OUT_DIR}/somatic/cnv/ascat_cnvkit -name "\${PAIR}*_somatic_cnv_report.tsv" 2>/dev/null | head -1 || true)"
    TMB_FILE="\$(find \${OUT_DIR}/somatic/tmb -name "\${PAIR}*.tmb.tsv" 2>/dev/null | head -1 || true)"
    PHARMACO="\$(find \${OUT_DIR}/germline/pharmaco -name "\${NORMAL}*.pharmaco.tsv" 2>/dev/null | head -1 || true)"
    GENE_COV_N="\$(find \${OUT_DIR}/qc/gene_coverage -name "\${NORMAL}*_gene_coverage.tsv" 2>/dev/null | head -1 || true)"
    GENE_COV_T="\$(find \${OUT_DIR}/qc/gene_coverage -name "\${TUMOR}*_gene_coverage.tsv" 2>/dev/null | head -1 || true)"
    MSI_FILE="\$(find ${sarek_outdir}/variant_calling/msisensorpro -name "\${PAIR}" -type f 2>/dev/null | head -1 || true)"
    COV_N="\$(find ${sarek_outdir}/reports/mosdepth -name "\${NORMAL}*.mosdepth.summary.txt" 2>/dev/null | head -1 || true)"
    COV_T="\$(find ${sarek_outdir}/reports/mosdepth -name "\${TUMOR}*.mosdepth.summary.txt" 2>/dev/null | head -1 || true)"

    # Mutational signature activity files
    SIG_SBS96="\$(find \${OUT_DIR}/somatic/signatures/assignment -name "Assignment_Solution_Activities.txt" -path "*SBS96*" 2>/dev/null | head -1 || true)"
    SIG_DBS="\$(find \${OUT_DIR}/somatic/signatures/assignment -name "Assignment_Solution_Activities.txt" -path "*DBS*" 2>/dev/null | head -1 || true)"
    SIG_ID="\$(find \${OUT_DIR}/somatic/signatures/assignment -name "Assignment_Solution_Activities.txt" -path "*ID*" 2>/dev/null | head -1 || true)"

    # GATK CNV files
    CNV_GERM_GATK="\$(find \${OUT_DIR}/germline/cnv_gatk/annotated -name "\${NORMAL}_segments_annotated.tsv" 2>/dev/null | head -1 || true)"
    CNV_SOM_GATK="\$(find \${OUT_DIR}/somatic/cnv_gatk/annotated -name "\${TUMOR}_annotated.tsv" 2>/dev/null | head -1 || true)"

    # Extract dashboard metrics
    tsv_val() { local f="\$1" col="\$2"; [[ -s "\$f" ]] && awk -F'\\t' -v c="\$col" 'NR==1{for(i=1;i<=NF;i++) if(tolower(\$i)==tolower(c)) col=i} NR==2{print \$col}' "\$f" || echo "NA"; }
    msi_pct()  { local f="\$1"; [[ -s "\$f" ]] && awk 'NR==2 && NF>=3{x=\$3; gsub(/%/,"",x); print x}' "\$f" || echo "NA"; }
    cov_mean() { local f="\$1"; [[ -s "\$f" ]] && awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(tolower(\$i)=="mean") mc=i} \$1=="total_region"{print \$mc; exit}' "\$f" || echo "NA"; }

    TMB_VAL="\$(tsv_val "\$TMB_FILE" TMB_mut_per_mb)"
    MSI_VAL="\$(msi_pct "\$MSI_FILE")"
    LOH_VAL="\$(tsv_val "${ascat_metrics}" LOH)"
    GI_VAL="\$(tsv_val "${ascat_metrics}" GI)"
    WGD_VAL="\$(tsv_val "${ascat_metrics}" WGD)"
    PLOIDY_VAL="\$(tsv_val "${ascat_metrics}" ploidy)"
    PURITY_RAW="\$(tsv_val "${ascat_metrics}" purity)"
    PURITY_PCT="\$(awk -v v="\$PURITY_RAW" 'BEGIN{if(v ~ /^[0-9.]+\$/) printf "%.1f", (v>1.5?v:v*100); else print "NA"}')"
    COV_N_MEAN="\$(cov_mean "\$COV_N")"
    COV_T_MEAN="\$(cov_mean "\$COV_T")"

    # HRD predictions (passed directly as Nextflow input; fallback to find if NO_FILE)
    HRD_FILE="${hrd_predictions}"
    if [[ "\${HRD_FILE}" == "NO_FILE" || ! -s "\${HRD_FILE}" ]]; then
        HRD_FILE="\$(find \${OUT_DIR}/somatic/hrd -name "\${PAIR}*_hrd_predictions*" 2>/dev/null | head -1 || true)"
    fi
    HRD_PROB="\$(tsv_val "\$HRD_FILE" hrd.prob)"
    # HRD status: custom probability threshold >= 0.45 = HRD+ (overrides prediction column)
    HRD_STATUS=\$(awk -v p="\${HRD_PROB}" 'BEGIN{
        if (p == "" || p == "NA" || p == ".") print "NA"
        else if (p+0 >= 0.45) print "HRD+"
        else print "HRD-"
    }')

    # Produce TXT bundle
    {
      echo "GENUIT Molecular Report"
      echo "======================"
      echo "Run:     ${run_name}"
      echo "Pair:    \${PAIR}"
      echo "Tumor:   \${TUMOR}"
      echo "Normal:  \${NORMAL}"
      echo "Patient: \${PATIENT_ID}"
      echo ""
      echo "DASHBOARD METRICS:"
      echo "  TMB:               \${TMB_VAL} mut/Mb"
      echo "  MSI:               \${MSI_VAL} %"
      echo "  LOH:               \${LOH_VAL}"
      echo "  GI:                \${GI_VAL}"
      echo "  WGD:               \${WGD_VAL}"
      echo "  Ploidy:            \${PLOIDY_VAL}"
      echo "  Tumor purity:      \${PURITY_PCT} %"
      echo "  Coverage (Normal): \${COV_N_MEAN} X"
      echo "  Coverage (Tumor):  \${COV_T_MEAN} X"
      echo "  HRD Probability:   \${HRD_PROB}"
      echo "  HRD Status:        \${HRD_STATUS}"
      echo ""
      echo "FILES:"
      for label in "Germline SNVs:\$SNV_GERMLINE" "Somatic SNVs:\$SNV_SOMATIC" \\
                   "CIViC:\$CIVIC_TSV" "Somatic CNVs:\$CNV_SOMATIC" \\
                   "TMB:\$TMB_FILE" "Pharmaco:\$PHARMACO" \\
                   "HRD Predictions:\$HRD_FILE" \\
                   "Signatures SBS96:\$SIG_SBS96" \\
                   "Germline CNVs GATK:\$CNV_GERM_GATK" \\
                   "Somatic CNVs GATK:\$CNV_SOM_GATK"; do
          k="\${label%%:*}"; v="\${label#*:}"
          printf "  %-22s %s\\n" "\$k" "\${v:-NOT FOUND}"
      done
    } > ${pair}_molecular_report.txt

    # Germline-somatic correlation (produces HTML fragment)
    python3 ${projectDir}/bin/germline_somatic_correlate.py \
        --germline      "\${SNV_GERMLINE}" \
        --somatic       "\${SNV_SOMATIC}" \
        --cnv           "\${CNV_SOMATIC}" \
        --cnv-germline  "\${CNV_GERM_GATK}" \
        --germline-panel "${projectDir}/assets/gene_panels/germline_103_cancer_predisposition.txt" \
        --sig           "\${SIG_SBS96}" \
        --tumor         "\${TUMOR}" \
        --tmb           "\${TMB_VAL}" \
        --msi           "\${MSI_VAL}" \
        --hrd           "\${HRD_PROB}" \
        --html-out  correlation_section.html \
        --json-out  correlation_section.json \
        2>&1 || true   # never fail the report if correlation errors

    # Build HTML (embed all data, generate KPI dashboard)
    python3 - <<PYEOF
import os, html, csv

GENE_CTX_TSV = "${projectDir}/assets/gene_annotations/somatic_gene_context.tsv"

def _load_gene_ctx(path):
    ctx = {}
    if not path or not os.path.isfile(path):
        return ctx
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            gene = row.get('GENE','').strip()
            if gene:
                ctx[gene] = {
                    'role':    row.get('ROLE','').strip(),
                    'tier':    row.get('CGC_TIER','').strip(),
                    'pathway': row.get('PATHWAY','').strip(),
                    'summary': row.get('GENE_SUMMARY','').strip(),
                }
    return ctx

gene_ctx = _load_gene_ctx(GENE_CTX_TSV)

def gene_consequence(role, impact, effect):
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
                        'Many oncogenic missense mutations are activating, locking the protein in a '
                        'constitutively active state and driving downstream pathway signaling.')
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

pair       = "${pair}"
tumor      = "\${TUMOR}"
normal     = "\${NORMAL}"
patient_id = "\${PATIENT_ID}"
run_name   = "${run_name}"

_msi_raw = "\${MSI_VAL}"
try:
    _msi_f = float(_msi_raw)
    _msi_status = "MSI detected" if _msi_f >= 5.0 else "MSI not detected (MSS)"
    _msi_disp   = _msi_status + " \u2014 " + _msi_raw + "%"
except ValueError:
    _msi_disp = _msi_raw

metrics = {
    "TMB (mut/Mb)":       "\${TMB_VAL}",
    "MSI":                 _msi_disp,
    "LOH":                 "\${LOH_VAL}",
    "Genomic Instability": "\${GI_VAL}",
    "WGD":                 "\${WGD_VAL}",
    "Ploidy":              "\${PLOIDY_VAL}",
    "Tumor Purity (%)":    "\${PURITY_PCT}",
    "Coverage Normal (X)": "\${COV_N_MEAN}",
    "Coverage Tumor (X)":  "\${COV_T_MEAN}",
    "HRD Probability":     "\${HRD_PROB}",
    "HRD Status":          "\${HRD_STATUS}",
}

sig_sbs96     = "\${SIG_SBS96}"
sig_dbs       = "\${SIG_DBS}"
sig_id        = "\${SIG_ID}"
cnv_germ_gatk = "\${CNV_GERM_GATK}"
cnv_som_gatk  = "\${CNV_SOM_GATK}"

# ----
# CNV harmonisation helpers
# ----
_SOMATIC_38_PANEL = "${projectDir}/assets/gene_panels/somatic_38_oncogenes.txt"
_GERMLINE_PANEL   = "${projectDir}/assets/gene_panels/germline_103_cancer_predisposition.txt"
_AMP_THRESHOLD    = 5

def _load_panel(path):
    genes = set()
    try:
        with open(path) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith('#'):
                    genes.add(g.upper())
    except Exception:
        pass
    return genes

def _harmonize_cnv(ascat_cnvkit_path, gatk_somatic_path):
    # Returns dict: {'homdel': [...], 'amp': [...], 'loh': [...]}
    somatic_panel  = _load_panel(_SOMATIC_38_PANEL)
    germline_panel = _load_panel(_GERMLINE_PANEL)

    a_homdel, a_amp, a_loh = {}, {}, {}
    if ascat_cnvkit_path and os.path.isfile(ascat_cnvkit_path):
        with open(ascat_cnvkit_path) as fh:
            for row in csv.DictReader(fh, delimiter='\\t'):
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
                elif etype == 'AMP' and gene in somatic_panel and cn >= _AMP_THRESHOLD:
                    if gene not in a_amp or cn > a_amp[gene]['cn']:
                        a_amp[gene] = {'cn': cn, 'chrom': chrom, 'log2': log2}
                elif etype == '2nd_hit_LOH' and gene in germline_panel:
                    a_loh[gene] = {'cn': cn, 'loh_type': loh_type, 'chrom': chrom}

    g_homdel, g_amp, g_loh = {}, {}, {}
    if gatk_somatic_path and os.path.isfile(gatk_somatic_path):
        with open(gatk_somatic_path) as fh:
            for row in csv.DictReader(fh, delimiter='\\t'):
                gene  = row.get('GENE', '').strip().upper()
                etype = row.get('EVENT_TYPE', '').strip()
                try:
                    cn = int(float(row.get('CN_ABSOLUTE', -1)))
                except (ValueError, TypeError):
                    cn = -1
                chrom = row.get('CONTIG', '')
                if etype == 'DELETION' and cn == 0 and gene in somatic_panel:
                    g_homdel[gene] = {'cn': 0, 'chrom': chrom}
                elif etype == 'AMPLIFICATION' and cn >= _AMP_THRESHOLD and gene in somatic_panel:
                    if gene not in g_amp or cn > g_amp[gene]['cn']:
                        g_amp[gene] = {'cn': cn, 'chrom': chrom}
                elif etype in ('CN-LOH', 'DEL+LOH') and gene in germline_panel:
                    g_loh[gene] = {
                        'cn': cn,
                        'loh_type': 'cnLOH' if etype == 'CN-LOH' else 'HETDEL+LOH',
                        'chrom': chrom,
                    }

    # Merge
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
    # Sort LOH: HETDEL+LOH first (most severe), then cnLOH, then others
    _loh_order = {'HETDEL+LOH': 0, 'cnLOH': 1}
    loh_list.sort(key=lambda r: (_loh_order.get(r['loh_type'], 2), r['gene']))

    return {'homdel': homdel_list, 'amp': amp_list, 'loh': loh_list}

def _conf_badge(c):
    if c == 'HIGH':
        return "<span class='badge badge-P' style='font-size:10px'>HIGH</span>"
    return "<span class='badge badge-III' style='font-size:10px'>MODERATE</span>"

def _tools_html(t):
    out = ''
    for tool in t.split('+'):
        out += "<span class='badge badge-VUS' style='font-size:10px'>" + html.escape(tool.strip()) + "</span> "
    return out

def _loh_badge(lt):
    colors = {'HETDEL+LOH': 'badge-P', 'cnLOH': 'badge-LP'}
    cls = colors.get(lt, 'badge-III')
    return "<span class='badge " + cls + "' style='font-size:10px'>" + html.escape(lt) + "</span>"

def build_unified_cnv_html(ascat_cnvkit_path, gatk_somatic_path):
    data = _harmonize_cnv(ascat_cnvkit_path, gatk_somatic_path)
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
            "Amplifications (CN &ge; " + str(_AMP_THRESHOLD) + ")</h3>")
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
                + str(_AMP_THRESHOLD) + ") detected in the somatic panel.</p>")

    # LOH
    out += ("<h3 style='color:#eab308;margin-top:16px'>"
            "Loss of Heterozygosity (LOH) &mdash; Germline Panel</h3>")
    if data['loh']:
        out += ("<table><thead><tr>"
                "<th>Gene</th><th>Chr</th><th>LOH Type</th>"
                "<th>CN ASCAT</th><th>CN GATK</th>"
                "<th>Tools</th><th>Confidence</th></tr></thead><tbody>")
        for r in data['loh']:
            out += ("<tr>"
                    "<td><b>" + html.escape(r['gene']) + "</b></td>"
                    "<td>" + html.escape(r['chrom']) + "</td>"
                    "<td>" + _loh_badge(r['loh_type']) + "</td>"
                    "<td>" + (str(r['cn_p']) if r['cn_p'] != '' else '&mdash;') + "</td>"
                    "<td>" + (str(r['cn_s']) if r['cn_s'] != '' else '&mdash;') + "</td>"
                    "<td>" + _tools_html(r['tools']) + "</td>"
                    "<td>" + _conf_badge(r['conf']) + "</td>"
                    "</tr>")
        out += "</tbody></table>"
    else:
        out += "<p class='warn'>No LOH events detected in the germline panel.</p>"

    out += "</div>"
    return out

files = {
    "Germline SNVs":       "\${SNV_GERMLINE}",
    "Somatic SNVs":        "\${SNV_SOMATIC}",
    "CIViC annotation":    "\${CIVIC_TSV}",
    "Somatic CNVs":        "\${CNV_SOMATIC}",
    "TMB":                 "\${TMB_FILE}",
    "Pharmacogenomics":    "\${PHARMACO}",
    "Gene Coverage (N)":   "\${GENE_COV_N}",
    "Gene Coverage (T)":   "\${GENE_COV_T}",
    "HRD Predictions":     "\${HRD_FILE}",
}

# ----
# ACMG classification
# ----
import sys as _sys
_sys.path.insert(0, '${projectDir}/bin')
from final_report_acmg import acmg_criteria_list, acmg_class, _BENIGN_CRIT

# ----
# CIViC / AMP tier
# ----
def load_civic_index(path):
    # Returns dict: (chrom, pos, ref, alt) -> {'level': best_level, 'disease': ..., 'therapy': ..., 'significance': ...}
    index = {}
    if not path or not os.path.isfile(path):
        return index
    level_order = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            key = (r.get('CHROM',''), r.get('POS',''), r.get('REF',''), r.get('ALT',''))
            lv = r.get('CIViC Evidence Level', '').strip().upper()
            if not lv: continue
            existing = index.get(key)
            if existing is None or level_order.get(lv, 9) < level_order.get(existing['level'], 9):
                index[key] = {
                    'level':      lv,
                    'disease':    r.get('CIViC Entity Disease', ''),
                    'therapy':    r.get('CIViC Entity Therapies', ''),
                    'significance': r.get('CIViC Entity Significance', ''),
                    'amp_cat':    r.get('CIViC Assertion AMP Category', ''),
                }
    return index

def amp_tier(row, hdr, civic_index):
    # Returns tier string: 'IA', 'IB', 'IIA', 'IIB', 'IIC', 'III', 'IV'
    def col(name):
        i = hdr.get(name, -1)
        return row[i].strip() if i >= 0 and i < len(row) else '.'

    key   = (col('CHROM'), col('POS'), col('REF'), col('ALT'))
    civic = civic_index.get(key)
    impact = col('ANN[*].IMPACT').upper()

    if civic:
        lv  = civic['level']
        sig = civic['significance'].upper()
        is_sens = 'SENSITIVITY' in sig or 'RESPONSE' in sig
        is_res  = 'RESISTANCE' in sig
        if lv == 'A':
            return 'IA' if (is_sens or is_res) else 'IIA'
        if lv == 'B':
            return 'IB' if (is_sens or is_res) else 'IIB'
        if lv in ('C', 'D'):
            return 'IIC'
        return 'IID'

    if impact in ('HIGH', 'MODERATE'):
        return 'III'
    return 'IV'

# ----
# Load CIViC index once
# ----
civic_index = load_civic_index(files['CIViC annotation'])

# ----
import sys as _sys3; _sys3.path.insert(0, '${projectDir}/bin')
from final_report_sections import (_CPIC_TABLE, _alert_badge, tsv_to_html_pharmaco, fetch_clinical_trials_html, build_methodology_html)

# ----
# Generic TSV -> HTML table
# ----
def tsv_to_html(path, title):
    if not path or not os.path.isfile(path):
        return "<p class='warn'>File not found: " + html.escape(path or "N/A") + "</p>"
    with open(path) as fh:
        lines = [l.rstrip("\\n") for l in fh if l.strip()]
    if not lines:
        return "<p class='warn'>Empty file</p>"
    header = lines[0].split("\\t")
    th = "".join("<th>" + html.escape(c) + "</th>" for c in header)
    tbody = ""
    for line in lines[1:]:
        cells = line.split("\\t")
        td = "".join("<td>" + html.escape(v) + "</td>" for v in cells)
        tbody += "<tr>" + td + "</tr>\\n"
    return ("<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
            "<tbody>" + tbody + "</tbody></table></div>")

# ----
# Germline SNV table with ACMG classification
# ----
def tsv_to_html_germline(path):
    if not path or not os.path.isfile(path):
        return "<p class='warn'>File not found: " + html.escape(path or "N/A") + "</p>"
    with open(path) as fh:
        lines = [l.rstrip("\\n") for l in fh if l.strip()]
    if not lines:
        return "<p class='warn'>Empty file</p>"

    raw_header = lines[0].split("\\t")
    hdr = {c: i for i, c in enumerate(raw_header)}

    # Display column map: display_name -> raw_column_name (or None for computed)
    disp_cols = [
        ('Locus',       None),
        ('Gene',        'ANN[*].GENE'),
        ('HGVS_c',      'ANN[*].HGVS_C'),
        ('HGVS_p',      'ANN[*].HGVS_P'),
        ('Impact',      'ANN[*].IMPACT'),
        ('DP',          None),
        ('VAF',         None),
        ('REVEL',       'dbNSFP_REVEL_score'),
        ('gnomAD_AF',   'dbNSFP_gnomAD_exomes_AF'),
        ('ClinVar',     'dbNSFP_clinvar_clnsig'),
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

    _BENIGN_CRIT = {'BA1', 'BS1', 'BP4', 'BP4_St'}

    def crit_badge(code):
        if code.startswith('ClinVar'):
            cls = 'badge-crit-cv'
        elif code.startswith('BP') or code.startswith('BA') or code.startswith('BS'):
            cls = 'badge-crit-b'
        else:
            cls = 'badge-crit-p'
        return "<span class='badge-crit " + cls + "'>" + html.escape(code) + "</span>"

    def get_col(row, name):
        i = hdr.get(name, -1)
        return row[i].strip() if i >= 0 and i < len(row) else ''

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
    for line in lines[1:]:
        row = line.split("\\t")
        acmg, criteria = acmg_criteria_list(row, hdr)
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
            elif disp_name == 'VAF':
                tds += "<td>" + html.escape(str(vaf_val)) + "</td>"
            else:
                val = get_col(row, raw_col) if raw_col else ''
                tds += "<td>" + html.escape(val) + "</td>"

        tbody += "<tr data-acmg='" + acmg + "'>" + tds + "</tr>\\n"

    return (
        filter_bar +
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody id='tbl-germline'>" + tbody + "</tbody></table></div>"
    )

# ----
# Somatic SNV table with AMP tier classification
# ----
def tsv_to_html_somatic(path, civic_idx):
    if not path or not os.path.isfile(path):
        return "<p class='warn'>File not found: " + html.escape(path or "N/A") + "</p>"
    with open(path) as fh:
        lines = [l.rstrip("\\n") for l in fh if l.strip()]
    if not lines:
        return "<p class='warn'>Empty file</p>"

    raw_header = lines[0].split("\\t")
    hdr = {c: i for i, c in enumerate(raw_header)}

    disp_cols = [
        ('Locus',        None),
        ('Gene',         'ANN[*].GENE'),
        ('HGVS_c',       'ANN[*].HGVS_C'),
        ('HGVS_p',       'ANN[*].HGVS_P'),
        ('Effect',       'ANN[*].EFFECT'),
        ('VAF',          None),
        ('Call',         None),
        ('CIViC Level',  None),
        ('Disease',      None),
        ('Therapy',      None),
        ('AMP Tier',     None),
    ]

    tier_colors = {
        'IA':  'badge-IA',
        'IB':  'badge-IB',
        'IIA': 'badge-IIA',
        'IIB': 'badge-IIB',
        'IIC': 'badge-IIC',
        'IID': 'badge-IIC',
        'III': 'badge-III',
        'IV':  'badge-IV',
    }

    def get_col(row, name):
        i = hdr.get(name, -1)
        return row[i].strip() if i >= 0 and i < len(row) else ''

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
    for line in lines[1:]:
        row = line.split("\\t")
        tier = amp_tier(row, hdr, civic_idx)
        badge_cls = tier_colors.get(tier, 'badge-IV')

        chrom = get_col(row, 'CHROM')
        pos   = get_col(row, 'POS')
        locus = chrom + ':' + pos if chrom and pos else (chrom or pos or '.')

        key   = (chrom, pos, get_col(row, 'REF'), get_col(row, 'ALT'))
        civic = civic_idx.get(key, {})
        civic_level   = civic.get('level', '')
        civic_disease = civic.get('disease', '')
        civic_therapy = civic.get('therapy', '')

        gene_name = get_col(row, 'ANN[*].GENE')
        impact_val = get_col(row, 'ANN[*].IMPACT')
        effect_val = get_col(row, 'ANN[*].EFFECT') or get_col(row, 'ANN[*].ANNOTATION') or ''
        ctx = gene_ctx.get(gene_name, {})
        ctx_role    = ctx.get('role', '')
        ctx_tier    = ctx.get('tier', '')
        ctx_pathway = ctx.get('pathway', '')
        ctx_summary = ctx.get('summary', '')
        conseq_text = gene_consequence(ctx_role, impact_val, effect_val) if ctx_role else ''
        row_id = "ctx-" + str(id(line))

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

        tbody += "<tr data-tier='" + tier + "'>" + tds + "</tr>\\n"

        # Hidden context row
        if ctx_role or ctx_summary:
            role_badge = "<span class='badge badge-role'>" + html.escape(ctx_role or 'Unknown') + "</span>"
            path_badge = "<span class='badge badge-path'>" + html.escape(ctx_pathway) + "</span>" if ctx_pathway else ""
            tier_badge = ("<span class='badge badge-cgc'>CGC Tier " + html.escape(ctx_tier) + "</span>") if ctx_tier else ""
            summary_p  = ("<p class='ctx-p'><strong>Gene function:</strong> " + html.escape(ctx_summary) + "</p>") if ctx_summary else ""
            conseq_p   = ("<p class='ctx-p'><strong>Tumor consequence:</strong> " + html.escape(conseq_text) + "</p>") if conseq_text else ""
            ctx_inner  = ("<div class='ctx-box'>"
                          "<div class='ctx-badges'>" + role_badge + path_badge + tier_badge + "</div>"
                          + summary_p + conseq_p
                          + "</div>")
            tbody += ("<tr class='ctx-row' data-tier='" + tier + "' style='display:none;'>"
                      "<td colspan='11'>" + ctx_inner + "</td></tr>\\n")
        else:
            tbody += ("<tr class='ctx-row' data-tier='" + tier + "' style='display:none;'>"
                      "<td colspan='11'><em style='color:#888;'>No annotation.</em></td></tr>\\n")

    return (
        filter_bar +
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody id='tbl-somatic'>" + tbody + "</tbody></table></div>"
    )

# ----
# Germline CNV GATK table: CN != 2, deduplicated by segment, quality-badged
# ----
def tsv_to_html_germline_cnv(path):
    # One row per exon -> group by (CHROM, POS, END). Only panel genes, QA >= 4.
    MIN_QA = int(os.environ.get('FILT_GERMLINE_CNV_MIN_QA','4'))
    _panel_genes = _load_panel(_GERMLINE_PANEL)
    if not path or not os.path.isfile(path):
        return "<p class='warn'>Germline CNV (GATK) file not found: " + html.escape(path or "N/A") + "</p>"
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        rows = list(reader)
    if not rows:
        return "<p class='warn'>Germline CNV (GATK): empty file</p>"
    # Group rows by segment
    segments = {}
    seg_order = []
    for row in rows:
        cn_raw = row.get('CN', '2').strip()
        try:
            cn_int = int(cn_raw)
        except ValueError:
            cn_int = 2
        if cn_int == 2:
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
            qa_val = int(seg.get('qa', '0') or '0')
        except (ValueError, TypeError):
            qa_val = 0
        if qa_val < MIN_QA:
            continue
        panel_hits = [g for g in seg['genes'] if g in _panel_genes]
        if panel_hits:
            seg['genes'] = panel_hits
            filtered.append(seg)
    if not filtered:
        return (
            "<div class='callout-info'>"
            "<strong>No Reportable Germline CNV Events</strong><br>"
            "The GATK germline copy-number variant analysis did not identify any CN &ne; 2 events "
            "meeting quality threshold in the 103-gene germline panel. "
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
        if cn == 0: return "HOMDEL"
        elif cn == 1: return "DEL (CN=1)"
        elif cn == 3: return "DUP (CN=3)"
        elif cn >= 4: return "AMP (CN=" + str(cn) + ")"
        return alt or "CNV"
    th = ("<th>Locus</th><th>Gene(s)</th><th>Type</th><th>CN</th>"
          "<th>NP (intervals)</th><th>QA</th><th>QS</th><th>Quality</th>")
    tbody = ""
    for seg in filtered:
        locus = seg['chrom'] + ":" + seg['pos'] + "-" + seg['end']
        genes = ", ".join(seg['genes']) if seg['genes'] else "."
        ctype = cn_type(seg['alt'], seg['cn'])
        tbody += (
            "<tr>"
            "<td>" + html.escape(locus) + "</td>"
            "<td><strong>" + html.escape(genes) + "</strong></td>"
            "<td>" + html.escape(ctype) + "</td>"
            "<td>" + html.escape(seg['cn']) + "</td>"
            "<td>" + html.escape(seg['np']) + "</td>"
            "<td>" + html.escape(seg['qa']) + "</td>"
            "<td>" + html.escape(seg['qs']) + "</td>"
            "<td>" + qs_badge(seg['qs']) + "</td>"
            "</tr>\\n"
        )
    note = ("<div style='font-size:.75rem;color:rgba(226,232,240,.6);margin:6px 0 4px;'>"
            + str(len(filtered)) + " segment(s) with CN &ne; 2 (NP &ge; 5 and QA &ge; 10)."
            " Quality: <span class='badge badge-P'>High</span> QS &ge; 30 &nbsp;"
            "<span class='badge badge-LP'>Medium</span> QS 10&ndash;29 &nbsp;"
            "<span class='badge badge-VUS'>Low</span> QS &lt; 10 (manual review recommended)."
            "</div>")
    return (
        note +
        "<div class='tbl-wrap'><table><thead><tr>" + th + "</tr></thead>"
        "<tbody>" + tbody + "</tbody></table></div>"
    )

# ----
# Signature table helper
# ----
SBS_ETIOLOGY_HTML = {
    "SBS1":  ("Clock-like",          "#4e79a7"),
    "SBS2":  ("APOBEC",              "#f28e2b"),
    "SBS3":  ("HRD (BRCA1/2)",       "#e15759"),
    "SBS4":  ("Tobacco smoking",     "#9c755f"),
    "SBS5":  ("Clock-like",          "#4e79a7"),
    "SBS6":  ("MMR deficiency",      "#b07aa1"),
    "SBS7a": ("UV radiation",        "#edc948"),
    "SBS7b": ("UV radiation",        "#edc948"),
    "SBS7c": ("UV radiation",        "#edc948"),
    "SBS7d": ("UV radiation",        "#edc948"),
    "SBS8":  ("HRD-associated",      "#c0392b"),
    "SBS9":  ("Somatic hypermut.",   "#bab0ac"),
    "SBS10a":("POLE/POLD1",          "#ff9da7"),
    "SBS10b":("POLE/POLD1",          "#ff9da7"),
    "SBS10c":("POLE/POLD1",          "#ff9da7"),
    "SBS10d":("POLE/POLD1",          "#ff9da7"),
    "SBS11": ("Alkylating agents",   "#59a14f"),
    "SBS13": ("APOBEC",              "#f28e2b"),
    "SBS14": ("MMR def. + POLE",     "#7b2d8b"),
    "SBS15": ("MMR deficiency",      "#b07aa1"),
    "SBS16": ("Alcohol",             "#ff7043"),
    "SBS18": ("ROS / MUTYH",         "#76b7b2"),
    "SBS20": ("MMR deficiency",      "#b07aa1"),
    "SBS21": ("MMR deficiency",      "#b07aa1"),
    "SBS22a":("Aristolochic acid",   "#8b4513"),
    "SBS22b":("Aristolochic acid",   "#8b4513"),
    "SBS24": ("Aflatoxin",           "#a0522d"),
    "SBS25": ("Chemotherapy",        "#2e8b57"),
    "SBS26": ("MMR deficiency",      "#b07aa1"),
    "SBS28": ("POLE/POLD1",          "#ff9da7"),
    "SBS29": ("Tobacco chewing",     "#6d4c41"),
    "SBS30": ("BER def. (NTHL1)",    "#00acc1"),
    "SBS31": ("Platinum agents",     "#1b5e20"),
    "SBS35": ("Platinum agents",     "#1b5e20"),
    "SBS36": ("BER def. (MUTYH)",    "#0097a7"),
    "SBS38": ("UV (indirect)",       "#f9a825"),
    "SBS44": ("MMR deficiency",      "#b07aa1"),
    "SBS84": ("Somatic hypermut.",   "#bab0ac"),
    "SBS85": ("Somatic hypermut.",   "#bab0ac"),
    "SBS87": ("Thiopurine",          "#388e3c"),
    "SBS88": ("Colibactin",          "#8bc34a"),
    "SBS90": ("Duocarmycin",         "#43a047"),
    "SBS92": ("Tobacco (bladder)",   "#5d4037"),
}
ID_ETIOLOGY_HTML = {
    "ID1":  ("Clock-like / Replication slippage", "#4e79a7"),
    "ID2":  ("Clock-like / Replication slippage", "#4e79a7"),
    "ID3":  ("Tobacco smoking",                   "#9c755f"),
    "ID6":  ("MMR deficiency",                    "#b07aa1"),
    "ID7":  ("MMR deficiency",                    "#b07aa1"),
    "ID11": ("Platinum agents",                   "#1b5e20"),
    "ID13": ("APOBEC",                            "#f28e2b"),
    "ID14": ("POLE/POLD1",                        "#ff9da7"),
}

def sig_to_html(path, sample_name):
    # Read wide-format signature activities; show only non-zero sigs for the sample.
    if not path or not os.path.isfile(path):
        return "<p class='warn'>File not found: " + html.escape(path or "N/A") + "</p>"
    with open(path) as fh:
        lines = [l.rstrip("\\n") for l in fh if l.strip()]
    if len(lines) < 2:
        return "<p class='warn'>No data in signatures file</p>"
    header = lines[0].split("\\t")
    # Detect sig type from header (SBS vs ID)
    sig_type_is_id = any(col.startswith("ID") for col in header[1:])
    etiol_map = ID_ETIOLOGY_HTML if sig_type_is_id else SBS_ETIOLOGY_HTML
    sample_row = None
    for line in lines[1:]:
        parts = line.split("\\t")
        if parts[0] == sample_name or sample_name in parts[0]:
            sample_row = parts
            break
    if sample_row is None:
        return "<p class='warn'>Sample " + html.escape(sample_name) + " not found in signatures file</p>"
    pairs = []
    for i, col in enumerate(header[1:], 1):
        try:
            val = float(sample_row[i]) if i < len(sample_row) else 0.0
        except (ValueError, IndexError):
            val = 0.0
        if val > 0:
            pairs.append((col, int(round(val))))
    if not pairs:
        return "<p class='warn'>No signatures detected (all zero)</p>"
    pairs.sort(key=lambda x: x[1], reverse=True)
    total = sum(c for _, c in pairs)
    th = "<th>Signature</th><th>Mutations Assigned</th><th>Proportion (%)</th><th>COSMIC Aetiology</th>"
    tbody = ""
    for sig, count in pairs:
        pct = "{:.1f}".format(100.0 * count / total) if total > 0 else "0.0"
        etiol, etiol_col = etiol_map.get(sig, ("Unknown", "#888888"))
        badge = ("<span style='background:" + etiol_col + ";color:#fff;padding:2px 8px;"
                 "border-radius:4px;font-size:.72rem;'>" + html.escape(etiol) + "</span>")
        tbody += ("<tr><td><strong>" + html.escape(sig) + "</strong></td>"
                  "<td>" + str(count) + "</td><td>" + pct + "%</td>"
                  "<td>" + badge + "</td></tr>\\n")
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
kpi_cards = ""
for k, v in metrics.items():
    color = hrd_status_color(v) if k == "HRD Status" else "var(--accent)"
    kpi_cards += (
        "<div class='kpi'>"
        "<div class='kpi-label'>" + html.escape(k) + "</div>"
        "<div class='kpi-value' style='color:" + color + "'>" + html.escape(v) + "</div>"
        "</div>"
    )

# ----
# Gene coverage quality table (dual-sample: Normal + Tumor side-by-side)
# ----
def build_coverage_quality_html(path_n, path_t, normal_name, tumor_name):
    # Return HTML for the per-gene coverage quality section.

    def load_cov(path):
        data = {}
        if not path or not os.path.isfile(path):
            return data
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\\t')
            for row in reader:
                g = row.get('gene', '').strip()
                if g:
                    data[g] = row
        return data

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

    cov_n = load_cov(path_n)
    cov_t = load_cov(path_t)
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
            "</tr>\\n"
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

# ----
# Build sections
# ----
sections = ""

# Germline SNVs — specialized table with ACMG
sections += "<section><h2>Germline SNVs</h2>" + tsv_to_html_germline(files['Germline SNVs']) + "</section>\\n"

# Somatic SNVs — specialized table with AMP tier
sections += "<section><h2>Somatic SNVs</h2>" + tsv_to_html_somatic(files['Somatic SNVs'], civic_index) + "</section>\\n"

# Remaining generic sections (skip gene coverage and handled tables)
generic_skip = {'Germline SNVs', 'Somatic SNVs', 'CIViC annotation',
                'Somatic CNVs', 'Pharmacogenomics',
                'Gene Coverage (N)', 'Gene Coverage (T)'}
for title, path in files.items():
    if title in generic_skip:
        continue
    sections += "<section><h2>" + html.escape(title) + "</h2>" + tsv_to_html(path, title) + "</section>\\n"

# Pharmacogenomics -- specialized CPIC-annotated table
sections += ("<section><h2>Pharmacogenomics (PGx) &#8212; CPIC Dose Recommendations</h2>"
             + tsv_to_html_pharmaco(files.get('Pharmacogenomics', ''))
             + "</section>\\n")

# Sequencing Quality — dedicated per-gene coverage page
sections += (
    "<section><h2>Sequencing Quality &#8212; Per-Gene Coverage</h2>"
    + build_coverage_quality_html(
        files.get('Gene Coverage (N)', ''),
        files.get('Gene Coverage (T)', ''),
        normal, tumor)
    + "</section>\\n"
)

# Mutational signature sections — table + profile/assignment plots
sections += ("<section><h2>Mutational Signatures &#8212; SBS96</h2>"
             + sig_to_html(sig_sbs96, pair)
             + sig_plots_to_html(sig_sbs96, 'SBS96')
             + "</section>\\n")
if sig_dbs and os.path.isfile(sig_dbs):
    sections += ("<section><h2>Mutational Signatures &#8212; DBS</h2>"
                 + sig_to_html(sig_dbs, pair)
                 + sig_plots_to_html(sig_dbs, 'DBS')
                 + "</section>\\n")
if sig_id and os.path.isfile(sig_id):
    sections += ("<section><h2>Mutational Signatures &#8212; ID</h2>"
                 + sig_to_html(sig_id, pair)
                 + sig_plots_to_html(sig_id, 'ID')
                 + "</section>\\n")

# Germline-Somatic Correlation section (reads HTML fragment produced by germline_somatic_correlate.py)
_corr_path = 'correlation_section.html'
if os.path.isfile(_corr_path):
    with open(_corr_path) as _cf:
        _corr_html = _cf.read()
else:
    _corr_html = "<p class='warn'>Correlation analysis not available (script output missing)</p>"
sections += '<section><h2>Germline&#8211;Somatic Correlation</h2>' + _corr_html + '</section>\\n'

# CNV sections
# Germline: GATK gCNV — constitutional CN changes in normal sample (unchanged)
sections += "<section><h2>Germline CNVs (GATK) &#8212; CN &ne; 2</h2>" + tsv_to_html_germline_cnv(cnv_germ_gatk) + "</section>\\n"
# Somatic: unified HOMDEL / AMP / LOH from CNVkit+ASCAT+GATK
sections += ("<section><h2>Somatic CNVs &#8212; Reportable Events</h2>"
             + build_unified_cnv_html(files.get('Somatic CNVs', ''), cnv_som_gatk)
             + "</section>\\n")

# Clinical Trials -- ClinicalTrials.gov query for identified findings
_tmb_val_ct = "\${TMB_VAL}"
_msi_val_ct = "\${MSI_VAL}"
_hrd_stat_ct = "\${HRD_STATUS}"
sections += ("<section><h2>Clinical Trials &#8212; Relevant Recruiting Studies</h2>"
             + fetch_clinical_trials_html('correlation_section.json', _tmb_val_ct, _msi_val_ct, _hrd_stat_ct)
             + "</section>\\n")

# Test Methodology -- final section
sections += "<section><h2>Test Methodology</h2>" + build_methodology_html() + "</section>\\n"

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
    '.badge { display:inline-block; padding:2px 8px; border-radius:12px; font-size:.7rem; font-weight:700; }'
    '.badge-P { background:#ef4444; color:#fff; }'
    '.badge-LP { background:#f97316; color:#fff; }'
    '.badge-VUS { background:#eab308; color:#111; }'
    '.badge-LB { background:#3b82f6; color:#fff; }'
    '.badge-B { background:#22c55e; color:#fff; }'
    '.badge-IA, .badge-IB { background:#dc2626; color:#fff; }'
    '.badge-IIA, .badge-IIB { background:#ea580c; color:#fff; }'
    '.badge-IIC { background:#d97706; color:#fff; }'
    '.badge-III { background:#6b7280; color:#fff; }'
    '.badge-IV { background:#374151; color:#888; }'
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
    '.print-btn { position:fixed; bottom:24px; right:24px; padding:10px 20px; background:var(--accent); color:#fff; border:none; border-radius:8px; cursor:pointer; font-size:.9rem; font-weight:600; z-index:999; }'
    '@media print {'
    '  .filter-bar, .print-btn { display:none !important; }'
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
    "    else if (show === 'ti-tii') visible = ['IA','IB','IIA','IIB','IIC'].indexOf(tier) >= 0;"
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
    "function printClinical() {"
    "  filterGermline('plpvus');"
    "  filterSomatic('ti-tii');"
    "  window.print();"
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
    '<!doctype html>\\n<html lang="en">\\n<head>\\n'
    '  <meta charset="utf-8"/>\\n'
    '  <title>GENUIT Molecular Report \\u2014 ' + html.escape(pair) + '</title>\\n'
    '  <style>\\n    ' + _css + '\\n  </style>\\n'
    '</head>\\n<body>\\n<div class="container">\\n'
    '  <div style="font-size:11px;letter-spacing:.1em;text-transform:uppercase;opacity:.7;margin-bottom:8px;">GENUIT \\u00b7 nf-sarek-clinical</div>\\n'
    '  <h1>Molecular Report \\u2014 ' + html.escape(patient_id) + '</h1>\\n'
    '  <div class="meta">\\n'
    '    <div><strong>Run:</strong> ' + html.escape(run_name) + '</div>\\n'
    '    <div><strong>Patient:</strong> ' + html.escape(patient_id) + '</div>\\n'
    '    <div><strong>Tumor:</strong> ' + html.escape(tumor) + '</div>\\n'
    '    <div><strong>Normal:</strong> ' + html.escape(normal) + '</div>\\n'
    '  </div>\\n'
    '  <h2>Dashboard</h2>\\n'
    '  <div class="kpi-grid">\\n    ' + kpi_cards + '\\n  </div>\\n'
    '  ' + sections + '\\n'
    '  <p style="margin-top:40px;font-size:.75rem;opacity:.6;">Generated by nf-sarek-clinical \\u00b7 GENUIT</p>\\n'
    '</div>\\n'
    '<button class="print-btn" onclick="printClinical()">&#128424; Reporte Cl&#237;nico</button>\\n'
    + _js + '\\n'
    '</body>\\n</html>\\n'
)

with open("${pair}_molecular_report.html", "w") as fh:
    fh.write(page)

print("HTML report written: ${pair}_molecular_report.html")
PYEOF
    """
}

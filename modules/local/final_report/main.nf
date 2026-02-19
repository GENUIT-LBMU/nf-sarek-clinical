/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINAL HTML REPORT MODULE
    =========================
    Aggregates all postprocessing outputs into a clinical HTML report.
    Auto-detects tumor/normal pair names from ASCAT metrics filename.

    Report sections:
      - KPI dashboard: TMB, MSI, LOH, GI, WGD, ploidy, purity, coverage
      - Gene coverage plots (interactive bar chart per gene)
      - Germline SNVs table
      - Germline CNVs table
      - Somatic SNVs table + CIViC annotation
      - Somatic CNVs (ASCAT/CNVkit) table
      - Pharmacogenomics table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FINAL_REPORT {
    tag "$pair"
    label 'process_single'

    publishDir "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(ascat_metrics)
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
    // The report generator script reads all output files from the pipeline's outdir
    // It mirrors the logic in report_final_V10_html_gene_coverage.sh but is
    // parameterised via environment variables for portability.
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

    # Collect input paths with fallbacks
    SNV_GERMLINE="\$(find \${OUT_DIR}/germline/snv/extract -name "\${NORMAL}*.germline_snvs.tsv" 2>/dev/null | head -1)"
    SNV_SOMATIC="\$(find \${OUT_DIR}/somatic/snv/extract -name "\${PAIR}*.somatic_snvs.tsv" 2>/dev/null | head -1)"
    CIVIC_TSV="\$(find \${OUT_DIR}/somatic/snv/civic_extract -name "\${PAIR}*.civic.exploded.tsv" 2>/dev/null | head -1)"
    CNV_SOMATIC="\$(find \${OUT_DIR}/somatic/cnv/ascat_cnvkit -name "\${PAIR}*_somatic_cnv_report.tsv" 2>/dev/null | head -1)"
    TMB_FILE="\$(find \${OUT_DIR}/somatic/tmb -name "\${PAIR}*.tmb.tsv" 2>/dev/null | head -1)"
    PHARMACO="\$(find \${OUT_DIR}/germline/pharmaco -name "\${NORMAL}*.pharmaco.tsv" 2>/dev/null | head -1)"
    GENE_COV_N="\$(find \${OUT_DIR}/qc/gene_coverage -name "\${NORMAL}*_gene_coverage.tsv" 2>/dev/null | head -1)"
    GENE_COV_T="\$(find \${OUT_DIR}/qc/gene_coverage -name "\${TUMOR}*_gene_coverage.tsv" 2>/dev/null | head -1)"
    MSI_FILE="\$(find ${sarek_outdir}/variant_calling/msisensorpro -name "\${PAIR}" -o -name "\${PAIR}.tsv" 2>/dev/null | head -1)"
    COV_N="\$(find ${sarek_outdir}/reports/mosdepth -name "\${NORMAL}*.mosdepth.summary.txt" 2>/dev/null | head -1)"
    COV_T="\$(find ${sarek_outdir}/reports/mosdepth -name "\${TUMOR}*.mosdepth.summary.txt" 2>/dev/null | head -1)"

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

    # Produce TXT bundle
    {
      echo "IMAXi Molecular Report"
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
      echo ""
      echo "FILES:"
      for label in "Germline SNVs:\$SNV_GERMLINE" "Somatic SNVs:\$SNV_SOMATIC" \\
                   "CIViC:\$CIVIC_TSV" "Somatic CNVs:\$CNV_SOMATIC" \\
                   "TMB:\$TMB_FILE" "Pharmaco:\$PHARMACO"; do
          k="\${label%%:*}"; v="\${label#*:}"
          printf "  %-20s %s\\n" "\$k" "\${v:-NOT FOUND}"
      done
    } > ${pair}_molecular_report.txt

    # Build HTML (embed all data, generate gauge dashboard)
    python3 - <<PYEOF
import os, html

pair       = "${pair}"
tumor      = "\${TUMOR}"
normal     = "\${NORMAL}"
patient_id = "\${PATIENT_ID}"
run_name   = "${run_name}"

metrics = {
    "TMB (mut/Mb)":       "\${TMB_VAL}",
    "MSI (%)":             "\${MSI_VAL}",
    "LOH":                 "\${LOH_VAL}",
    "Genomic Instability": "\${GI_VAL}",
    "WGD":                 "\${WGD_VAL}",
    "Ploidy":              "\${PLOIDY_VAL}",
    "Tumor Purity (%)":    "\${PURITY_PCT}",
    "Coverage Normal (X)": "\${COV_N_MEAN}",
    "Coverage Tumor (X)":  "\${COV_T_MEAN}",
}

files = {
    "Germline SNVs":       "\${SNV_GERMLINE}",
    "Somatic SNVs":        "\${SNV_SOMATIC}",
    "CIViC annotation":    "\${CIVIC_TSV}",
    "Somatic CNVs":        "\${CNV_SOMATIC}",
    "TMB":                 "\${TMB_FILE}",
    "Pharmacogenomics":    "\${PHARMACO}",
    "Gene Coverage (N)":   "\${GENE_COV_N}",
    "Gene Coverage (T)":   "\${GENE_COV_T}",
}

def tsv_to_html(path, title):
    if not path or not os.path.isfile(path):
        return f"<p class='warn'>File not found: {html.escape(path or 'N/A')}</p>"
    rows = []
    with open(path) as fh:
        lines = [l.rstrip('\\n') for l in fh if l.strip()]
    if not lines:
        return "<p class='warn'>Empty file</p>"
    header = lines[0].split('\\t')
    th = "".join(f"<th>{html.escape(c)}</th>" for c in header)
    tbody = ""
    for line in lines[1:]:
        cells = line.split('\\t')
        td = "".join(f"<td>{html.escape(v)}</td>" for v in cells)
        tbody += f"<tr>{td}</tr>\\n"
    return (f"<div class='tbl-wrap'><table><thead><tr>{th}</tr></thead>"
            f"<tbody>{tbody}</tbody></table></div>")

metric_rows = "".join(
    f"<tr><td>{html.escape(k)}</td><td><strong>{html.escape(v)}</strong></td></tr>"
    for k, v in metrics.items()
)

sections = ""
for title, path in files.items():
    sections += f"<section><h2>{html.escape(title)}</h2>{tsv_to_html(path, title)}</section>\\n"

page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>IMAXi Molecular Report — {html.escape(pair)}</title>
  <style>
    :root {{ --bg:#020617; --fg:#f8fafc; --accent:#f97316; --line:rgba(148,163,184,0.2); }}
    body {{ background:var(--bg); color:var(--fg); font-family:system-ui,sans-serif; margin:0; padding:0; }}
    .container {{ max-width:1180px; margin:0 auto; padding:28px 18px 60px; }}
    h1 {{ font-size:28px; }} h2 {{ font-size:18px; margin-top:28px; border-bottom:1px solid var(--line); padding-bottom:6px; }}
    .meta {{ display:grid; grid-template-columns:1fr 1fr; gap:8px 20px; background:rgba(15,23,42,0.6);
             border:1px solid rgba(249,115,22,0.2); border-radius:12px; padding:14px 18px; margin:14px 0; }}
    .kpi-grid {{ display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin:16px 0; }}
    .kpi {{ background:rgba(15,23,42,0.7); border:1px solid var(--line); border-radius:12px; padding:14px 16px; }}
    .kpi-label {{ font-size:11px; color:rgba(226,232,240,0.7); text-transform:uppercase; letter-spacing:.08em; }}
    .kpi-value {{ font-size:22px; font-weight:700; color:var(--accent); margin-top:4px; }}
    section {{ margin-top:20px; }}
    .tbl-wrap {{ overflow-x:auto; border-radius:10px; border:1px solid var(--line); margin-top:10px; }}
    table {{ border-collapse:collapse; width:100%; font-size:.78rem; }}
    thead {{ background:rgba(249,115,22,0.6); }} th {{ padding:.4rem .6rem; text-align:left; color:#fff; }}
    td {{ padding:.35rem .6rem; border-bottom:1px solid var(--line); color:#e5e7eb; }}
    tr:hover td {{ background:rgba(249,115,22,0.15); }}
    .warn {{ color:#fed7aa; background:rgba(124,45,18,0.3); border-radius:8px; padding:8px 12px; font-size:.85rem; }}
    @media print {{ body {{ background:#fff; color:#111; }} }}
  </style>
</head>
<body>
<div class="container">
  <div style="font-size:11px;letter-spacing:.1em;text-transform:uppercase;opacity:.7;margin-bottom:8px;">IMAXi · nf-sarek-clinical</div>
  <h1>Molecular Report — {html.escape(patient_id)}</h1>
  <div class="meta">
    <div><strong>Run:</strong> {html.escape(run_name)}</div>
    <div><strong>Patient:</strong> {html.escape(patient_id)}</div>
    <div><strong>Tumor:</strong> {html.escape(tumor)}</div>
    <div><strong>Normal:</strong> {html.escape(normal)}</div>
  </div>
  <h2>Dashboard</h2>
  <div class="kpi-grid">
    {"".join(f"<div class='kpi'><div class='kpi-label'>{html.escape(k)}</div><div class='kpi-value'>{html.escape(v)}</div></div>" for k, v in metrics.items())}
  </div>
  {sections}
  <p style="margin-top:40px;font-size:.75rem;opacity:.6;">Generated by nf-sarek-clinical · IMAXi</p>
</div>
</body>
</html>"""

with open("${pair}_molecular_report.html", "w") as fh:
    fh.write(page)

print("HTML report written: ${pair}_molecular_report.html")
PYEOF
    """
}

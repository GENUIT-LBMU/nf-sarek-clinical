/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL PDF REPORT MODULE
    ==========================
    Generates a clinical-grade PDF report (GENUIT Molecular Report) for each
    tumor/normal pair using ReportLab Platypus.

    Report sections:
      - KPI dashboard: TMB, MSI, LOH, GI, WGD, ploidy, purity, coverage, HRD
      - Germline SNVs (ACMG P / LP / VUS classification)
      - Somatic SNVs  (AMP/ASCO/CAP Tier IA – III classification with CIViC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINICAL_PDF_REPORT {
    tag "$pair"
    label 'process_single'

    publishDir "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(ascat_metrics), path(hrd_predictions), path(snv_germline_tsv), path(snv_somatic_tsv)
    path civic_tsv
    val sarek_outdir
    val run_name

    output:
    tuple val(pair), path("${pair}_clinical_report.pdf"), emit: pdf

    stub:
    """
    touch ${pair}_clinical_report.pdf
    """

    script:
    """
    # ── Parse pair name: TUMOR_vs_NORMAL ──────────────────────────────────────
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

    # ── Locate input files ─────────────────────────────────────────────────────
    SNV_GERMLINE="\$(find \${OUT_DIR}/germline/snv/extract   -name "\${NORMAL}*.germline_snvs.tsv" 2>/dev/null | head -1 || true)"
    SNV_SOMATIC="\$(find \${OUT_DIR}/somatic/snv/extract     -name "\${PAIR}*.somatic_snvs.tsv"   2>/dev/null | head -1 || true)"
    CIVIC_TSV="${civic_tsv}"
    TMB_FILE="\$(find \${OUT_DIR}/somatic/tmb                -name "\${PAIR}*.tmb.tsv"            2>/dev/null | head -1 || true)"
    MSI_FILE="\$(find ${sarek_outdir}/variant_calling/msisensorpro -name "\${PAIR}" -type f       2>/dev/null | head -1 || true)"
    COV_N="\$(find ${sarek_outdir}/reports/mosdepth          -name "\${NORMAL}*.mosdepth.summary.txt" 2>/dev/null | head -1 || true)"
    COV_T="\$(find ${sarek_outdir}/reports/mosdepth          -name "\${TUMOR}*.mosdepth.summary.txt"  2>/dev/null | head -1 || true)"
    CNV_SOMATIC="\$(find \${OUT_DIR}/somatic/cnv/ascat_cnvkit  -name "\${PAIR}*_somatic_cnv_report.tsv" 2>/dev/null | head -1 || true)"
    CNV_SOM_GATK="\$(find \${OUT_DIR}/somatic/cnv_gatk/annotated -name "\${TUMOR}_annotated.tsv"          2>/dev/null | head -1 || true)"
    CNV_GERM_GATK="\$(find \${OUT_DIR}/germline/cnv_gatk/annotated -name "\${NORMAL}_segments_annotated.tsv" 2>/dev/null | head -1 || true)"
    PHARMACO="\$(find \${OUT_DIR}/germline/pharmaco -name "\${NORMAL}*.pharmaco.tsv" 2>/dev/null | head -1 || true)"
    GENE_COV_N="\$(find \${OUT_DIR}/qc/gene_coverage -name "\${NORMAL}*_gene_coverage.tsv" 2>/dev/null | head -1 || true)"
    GENE_COV_T="\$(find \${OUT_DIR}/qc/gene_coverage -name "\${TUMOR}*_gene_coverage.tsv"  2>/dev/null | head -1 || true)"

    # HRD predictions: use Nextflow-staged file; fall back to find if placeholder
    HRD_FILE="${hrd_predictions}"
    if [[ "\${HRD_FILE}" == "NO_FILE" || ! -s "\${HRD_FILE}" ]]; then
        HRD_FILE="\$(find \${OUT_DIR}/somatic/hrd -name "\${PAIR}*_hrd_predictions*" 2>/dev/null | head -1 || true)"
    fi

    # ── Metric extraction helpers ───────────────────────────────────────────────
    # tsv_val FILE COLUMN  – extract value from row 2 of a TSV (case-insensitive header match)
    tsv_val() {
        local f="\$1" col="\$2"
        [[ -s "\$f" ]] && awk -F'\\t' -v c="\$col" \
            'NR==1{for(i=1;i<=NF;i++) if(tolower(\$i)==tolower(c)) col=i} NR==2{print \$col}' "\$f" \
            || echo "NA"
    }

    # msi_pct FILE  – extract MSI percentage from MSIsensorPro output (col 3)
    msi_pct() {
        local f="\$1"
        [[ -s "\$f" ]] && awk 'NR==2 && NF>=3{x=\$3; gsub(/%/,"",x); print x}' "\$f" || echo "NA"
    }

    # cov_mean FILE  – extract mean coverage for the total_region from mosdepth summary
    cov_mean() {
        local f="\$1"
        [[ -s "\$f" ]] && awk -F'\\t' \
            'NR==1{for(i=1;i<=NF;i++) if(tolower(\$i)=="mean") mc=i} \$1=="total_region"{print \$mc; exit}' \
            "\$f" || echo "NA"
    }

    # ── Extract dashboard metric values ─────────────────────────────────────────
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

    # ── Correlation pre-computation (produces JSON consumed by PDF script) ─────
    SIG_SBS96="\$(find \${OUT_DIR}/somatic/signatures/assignment -name "Assignment_Solution_Activities.txt" -path "*_SBS96/Assignment_Solution/*" 2>/dev/null | head -1 || true)"
    SIG_ID="\$(find \${OUT_DIR}/somatic/signatures/assignment    -name "Assignment_Solution_Activities.txt" -path "*_ID/Assignment_Solution/*"   2>/dev/null | head -1 || true)"

    python3 ${projectDir}/bin/germline_somatic_correlate.py \\
        --germline      "\${SNV_GERMLINE:-}"   \\
        --somatic       "\${SNV_SOMATIC:-}"    \\
        --cnv           "\${CNV_SOMATIC:-}"    \\
        --cnv-germline  "\${CNV_GERM_GATK:-}" \\
        --germline-panel "${projectDir}/assets/gene_panels/germline_cancer_predisposition.txt" \\
        --sig           "\${SIG_SBS96:-}"      \\
        --sig-id        "\${SIG_ID:-}"         \\
        --tumor         "\${TUMOR}"            \\
        --tmb           "\${TMB_VAL}"          \\
        --msi           "\${MSI_VAL}"          \\
        --hrd           "\$(tsv_val "\$HRD_FILE" hrd.prob)" \\
        --json-out  correlation_section.json \\
        2>&1 || true

    # ── Call Python PDF generator ───────────────────────────────────────────────
    python3 ${projectDir}/bin/clinical_pdf_report.py \\
        --pair        "\${PAIR}"       \\
        --run-name    "\${RUN}"        \\
        --patient     "\${PATIENT_ID}" \\
        --tumor       "\${TUMOR}"      \\
        --normal      "\${NORMAL}"     \\
        --germline    "\${SNV_GERMLINE:-}" \\
        --somatic     "\${SNV_SOMATIC:-}"  \\
        --civic       "\${CIVIC_TSV:-}"    \\
        --tumor-type  "${params.tumor_type ?: ''}"  \\
        --tumor-types "${params.tumor_types ?: ''}" \\
        --tmb         "\${TMB_FILE:-}"     \\
        --ascat       "${ascat_metrics}"   \\
        --hrd         "\${HRD_FILE:-}"     \\
        --cnv-somatic        "\${CNV_SOMATIC:-}"    \\
        --cnv-gatk-somatic   "\${CNV_SOM_GATK:-}"  \\
        --cnv-germline-gatk  "\${CNV_GERM_GATK:-}" \\
        --somatic-cnv-panel  "${projectDir}/assets/gene_panels/somatic_38_oncogenes.txt" \\
        --sig-sbs96          "\${SIG_SBS96:-}"     \\
        --sig-id             "\${SIG_ID:-}"        \\
        --gene-ctx    "${projectDir}/assets/gene_annotations/somatic_gene_context.tsv" \\
        --gene-cov-n       "\${GENE_COV_N:-}" \\
        --gene-cov-t       "\${GENE_COV_T:-}" \\
        --germline-panel   "${projectDir}/assets/gene_panels/germline_cancer_predisposition.txt" \\
        --somatic-panel    "${projectDir}/assets/gene_panels/somatic_snv_panel.txt" \\
        --pharmaco    "\${PHARMACO:-}"  \\
        --corr-json   correlation_section.json \\
        --language    "${params.report_language ?: 'en'}" \\
        --tmb-val     "\${TMB_VAL}"    \\
        --msi-val     "\${MSI_VAL}"    \\
        --loh-val     "\${LOH_VAL}"    \\
        --gi-val      "\${GI_VAL}"     \\
        --wgd-val     "\${WGD_VAL}"    \\
        --ploidy-val  "\${PLOIDY_VAL}" \\
        --purity-val  "\${PURITY_PCT}" \\
        --cov-n-val   "\${COV_N_MEAN}" \\
        --cov-t-val   "\${COV_T_MEAN}" \\
        --out         "${pair}_clinical_report.pdf"
    """
}

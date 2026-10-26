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
    tuple val(pair), path(ascat_metrics), path(hrd_predictions), path(snv_germline_tsv), path(snv_somatic_tsv), path(civic_tsv),
          path(snv_second_hit_tsv), path(snv_tso500_tsv), path(civic_tso500_tsv),
          path(tmb_tsv), path(msi_file), path(cov_n_summary), path(cov_t_summary),
          path(cnv_somatic_tsv), path(cnv_gatk_somatic_tsv), path(cnv_gatk_germline_tsv),
          path(pharmaco_tsv), path(gene_cov_n_tsv), path(gene_cov_t_tsv),
          path(sig_assignment, stageAs: 'sig/assignment'), path(sig_profiles, stageAs: 'sig/matrix_profiles')
    val run_name
    val not_run     // analyses switched off in this run (comma-separated, bin/report_shared.py)
    val report_info // JSON: pipeline version, commit, sarek release, resources (bin/report_shared.py)
    path sarek_versions, stageAs: 'sarek_versions/*'

    output:
    tuple val(pair), path("${pair}_clinical_report.pdf"), emit: pdf

    stub:
    """
    touch ${pair}_clinical_report.pdf
    """

    script:
    """
    # Report inputs are staged per pair (placeholders NO_* when absent)
    staged() { [[ -s "\$1" ]] && echo "\$1" || echo ""; }

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

    RUN="${run_name}"

    # ── Locate input files ─────────────────────────────────────────────────────
    SNV_GERMLINE="\$(staged "${snv_germline_tsv}")"
    SNV_SOMATIC="\$(staged "${snv_somatic_tsv}")"
    CIVIC_TSV="\$(staged "${civic_tsv}")"
    SNV_SECOND_HIT="\$(staged "${snv_second_hit_tsv}")"   # second-hit track: correlation only
    SNV_TSO500="\$(staged "${snv_tso500_tsv}")"           # pan-tumor track (TSO500 outside the core)
    CIVIC_TSO500="\$(staged "${civic_tso500_tsv}")"
    TMB_FILE="\$(staged "${tmb_tsv}")"
    MSI_FILE="\$(staged "${msi_file}")"
    COV_N="\$(staged "${cov_n_summary}")"
    COV_T="\$(staged "${cov_t_summary}")"
    CNV_SOMATIC="\$(staged "${cnv_somatic_tsv}")"
    CNV_SOM_GATK="\$(staged "${cnv_gatk_somatic_tsv}")"
    CNV_GERM_GATK="\$(staged "${cnv_gatk_germline_tsv}")"
    PHARMACO="\$(staged "${pharmaco_tsv}")"
    GENE_COV_N="\$(staged "${gene_cov_n_tsv}")"
    GENE_COV_T="\$(staged "${gene_cov_t_tsv}")"
    # Signature activities inside the staged SigProfilerAssignment directory
    sig_act() { [[ -d sig/assignment ]] && find -L sig/assignment -name "Assignment_Solution_Activities.txt" -path "*_\$1/Assignment_Solution/*" | head -1 || true; }
    SIG_SBS96="\$(sig_act SBS96)"
    SIG_DBS="\$(sig_act DBS)"
    SIG_ID="\$(sig_act ID)"
    export REPORT_NOT_RUN="${not_run}"
    cat > report_info.json <<'INFOEOF'
${report_info}
INFOEOF
    export REPORT_INFO="\$PWD/report_info.json"
    export REPORT_SAREK_VERSIONS="\$(staged "${sarek_versions}")"

    # HRD predictions: the staged file of this pair
    HRD_FILE="\$(staged "${hrd_predictions}")"

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

    # ── Correlation (JSON consumed by the PDF script; a failure fails the task) ─

    python3 ${projectDir}/bin/germline_somatic_correlate.py \\
        --germline      "\${SNV_GERMLINE:-}"   \\
        --somatic       "\${SNV_SOMATIC:-}"    \\
        --somatic-second-hit "\${SNV_SECOND_HIT:-}" \\
        --cnv           "\${CNV_SOMATIC:-}"    \\
        --cnv-germline  "\${CNV_GERM_GATK:-}" \\
        --germline-panel "${projectDir}/assets/gene_panels/germline_cancer_predisposition.txt" \\
        --sig           "\${SIG_SBS96:-}"      \\
        --sig-id        "\${SIG_ID:-}"         \\
        --tumor         "\${TUMOR}"            \\
        --tmb           "\${TMB_VAL}"          \\
        --msi           "\${MSI_VAL}"          \\
        --hrd           "\$(tsv_val "\$HRD_FILE" hrd.prob)" \\
        --json-out  correlation_section.json

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
        --somatic-second-hit "\${SNV_SECOND_HIT:-}" \\
        --somatic-tso500     "\${SNV_TSO500:-}"     \\
        --civic-tso500       "\${CIVIC_TSO500:-}"   \\
        --capture-bed        "${params.panel_bed ?: ''}" \\
        --tumor-type  "${params.tumor_type ?: ''}"  \\
        --tumor-types "${params.tumor_types ?: ''}" \\
        --tmb         "\${TMB_FILE:-}"     \\
        --ascat       "\$(staged "${ascat_metrics}")" \\
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

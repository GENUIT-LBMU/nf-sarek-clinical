/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORT MODEL
    ============
    Germline-somatic correlation and the report model of one pair
    (bin/report_model.py): every value and interpretation both reports show.
    FINAL_REPORT (laboratory HTML) and CLINICAL_PDF_REPORT only render it.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REPORT_MODEL {
    tag "$pair"
    label 'process_single'

    publishDir "${params.outdir}/reporting", mode: params.publish_dir_mode,
        saveAs: { f -> f.endsWith('.report_model.json') ? f : null }

    input:
    tuple val(pair), path(ascat_metrics), path(hrd_predictions), path(snv_germline_tsv), path(snv_somatic_tsv), path(civic_tsv),
          path(snv_second_hit_tsv), path(snv_tso500_tsv), path(civic_tso500_tsv),
          path(tmb_tsv), path(msi_file), path(cov_n_summary), path(cov_t_summary),
          path(cnv_somatic_tsv), path(cnv_gatk_somatic_tsv), path(cnv_gatk_germline_tsv),
          path(pharmaco_tsv), path(gene_cov_n_tsv), path(gene_cov_t_tsv),
          path(sig_assignment, stageAs: 'sig/assignment'), path(sig_profiles, stageAs: 'sig/matrix_profiles')
    val run_name
    val not_run     // analyses switched off in this run (comma-separated, bin/report_shared.py)
    val report_info // JSON: pipeline version, commit, sarek release, resources, code_md5
    path sarek_versions, stageAs: 'sarek_versions/*'
    path samplesheet, stageAs: 'samplesheet/*'   // patient / specimen identification, or NO_SAMPLESHEET

    output:
    tuple val(pair), path("${pair}.report_model.json"), emit: model

    stub:
    """
    echo '{}' > ${pair}.report_model.json
    """

    script:
    """
    # Inputs are staged per pair; placeholders NO_* when absent
    staged() { [[ -s "\$1" ]] && echo "\$1" || echo ""; }

    PAIR="${pair}"
    if [[ "\${PAIR}" =~ ^(.+)_vs_(.+)\$ ]]; then
        TUMOR="\${BASH_REMATCH[1]}"; NORMAL="\${BASH_REMATCH[2]}"
    else
        TUMOR="\${PAIR}"; NORMAL="UNKNOWN"
    fi
    PATIENT_ID="\${TUMOR%%_T*}"

    # Signature activities inside the staged SigProfilerAssignment directory
    sig_act() { [[ -d sig/assignment ]] && find -L sig/assignment -name "Assignment_Solution_Activities.txt" -path "*_\$1/Assignment_Solution/*" | head -1 || true; }
    SIG_SBS96="\$(sig_act SBS96)"; SIG_DBS="\$(sig_act DBS)"; SIG_ID="\$(sig_act ID)"

    export REPORT_NOT_RUN="${not_run}"
    cat > report_info.json <<'INFOEOF'
${report_info}
INFOEOF
    export REPORT_INFO="\$PWD/report_info.json"
    export REPORT_SAREK_VERSIONS="\$(staged "${sarek_versions}")"

    # Metric values the correlation needs (same readers as the model)
    val() { python3 -c "import sys; sys.path.insert(0, '${projectDir}/bin'); import report_model as r; print(\$1)"; }
    TMB_VAL="\$(val "r.tsv_val('\$(staged ${tmb_tsv})', 'TMB_mut_per_mb')")"
    MSI_VAL="\$(val "r.msi_pct('\$(staged ${msi_file})')")"
    HRD_VAL="\$(val "r.tsv_val('\$(staged ${hrd_predictions})', 'hrd.prob')")"

    # Germline-somatic correlation (JSON + HTML fragment; a failure fails the task)
    python3 ${projectDir}/bin/germline_somatic_correlate.py \\
        --germline      "\$(staged ${snv_germline_tsv})" \\
        --somatic       "\$(staged ${snv_somatic_tsv})" \\
        --somatic-second-hit "\$(staged ${snv_second_hit_tsv})" \\
        --cnv           "\$(staged ${cnv_somatic_tsv})" \\
        --cnv-germline  "\$(staged ${cnv_gatk_germline_tsv})" \\
        --germline-panel "${projectDir}/assets/gene_panels/germline_cancer_predisposition.txt" \\
        --sig           "\${SIG_SBS96}" \\
        --sig-id        "\${SIG_ID}" \\
        --tumor         "\${TUMOR}" \\
        --tmb           "\${TMB_VAL}" \\
        --msi           "\${MSI_VAL}" \\
        --hrd           "\${HRD_VAL}" \\
        --html-out      correlation_section.html \\
        --json-out      correlation_section.json

    python3 ${projectDir}/bin/report_model.py \\
        --pair "\${PAIR}" --run-name "${run_name}" --patient "\${PATIENT_ID}" \\
        --tumor "\${TUMOR}" --normal "\${NORMAL}" \\
        --germline           "\$(staged ${snv_germline_tsv})" \\
        --somatic            "\$(staged ${snv_somatic_tsv})" \\
        --civic              "\$(staged ${civic_tsv})" \\
        --somatic-second-hit "\$(staged ${snv_second_hit_tsv})" \\
        --somatic-tso500     "\$(staged ${snv_tso500_tsv})" \\
        --civic-tso500       "\$(staged ${civic_tso500_tsv})" \\
        --capture-bed        "${params.panel_bed ?: ''}" \\
        --tumor-type         "${params.tumor_type ?: ''}" \\
        --tumor-types        "${params.tumor_types ?: ''}" \\
        --tmb                "\$(staged ${tmb_tsv})" \\
        --msi                "\$(staged ${msi_file})" \\
        --ascat              "\$(staged ${ascat_metrics})" \\
        --hrd                "\$(staged ${hrd_predictions})" \\
        --cov-n              "\$(staged ${cov_n_summary})" \\
        --cov-t              "\$(staged ${cov_t_summary})" \\
        --cnv-somatic        "\$(staged ${cnv_somatic_tsv})" \\
        --cnv-gatk-somatic   "\$(staged ${cnv_gatk_somatic_tsv})" \\
        --cnv-germline-gatk  "\$(staged ${cnv_gatk_germline_tsv})" \\
        --somatic-cnv-panel  "${projectDir}/assets/gene_panels/somatic_38_oncogenes.txt" \\
        --sig-sbs96 "\${SIG_SBS96}" --sig-dbs "\${SIG_DBS}" --sig-id "\${SIG_ID}" \\
        --gene-ctx           "${projectDir}/assets/gene_annotations/somatic_gene_context.tsv" \\
        --gene-cov-n         "\$(staged ${gene_cov_n_tsv})" \\
        --gene-cov-t         "\$(staged ${gene_cov_t_tsv})" \\
        --germline-panel     "${projectDir}/assets/gene_panels/germline_cancer_predisposition.txt" \\
        --somatic-panel      "${projectDir}/assets/gene_panels/somatic_snv_panel.txt" \\
        --pharmaco           "\$(staged ${pharmaco_tsv})" \\
        --samplesheet        "\$(staged ${samplesheet})" \\
        --corr-json correlation_section.json --corr-html correlation_section.html \\
        --out "${pair}.report_model.json"
    """
}

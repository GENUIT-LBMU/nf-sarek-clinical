/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINICAL PDF REPORT
    ===================
    Renders the clinical subset of the report model of one pair
    (bin/clinical_pdf_report.py, ReportLab): biomarkers and major findings on
    page 1, germline P/LP/VUS, somatic Tier I-III, CNV / LOH with the events left
    out, signatures, correlation, PGx, trials, coverage, methodology.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINICAL_PDF_REPORT {
    tag "$pair"
    label 'process_single'

    publishDir "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    // model + the signature directories (plots are embedded from them)
    tuple val(pair), path(model), path(sig_assignment, stageAs: 'sig/assignment'), path(sig_profiles, stageAs: 'sig/matrix_profiles')

    output:
    tuple val(pair), path("${pair}_clinical_report.pdf"), emit: pdf

    stub:
    """
    touch ${pair}_clinical_report.pdf
    """

    script:
    """
    python3 ${projectDir}/bin/clinical_pdf_report.py \\
        --model    ${model} \\
        --language "${params.report_language ?: 'en'}" \\
        --out      ${pair}_clinical_report.pdf
    """
}

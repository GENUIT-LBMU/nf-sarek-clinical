/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LABORATORY HTML REPORT
    ======================
    Renders the whole report model of one pair (bin/render_html.py): every
    germline class, every AMP tier, all LOH events, DBS signatures, raw TMB /
    HRD tables, filters and gene context. No classification happens here.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FINAL_REPORT {
    tag "$pair"
    label 'process_single'

    publishDir "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    // model + the signature directories (plots are embedded from them)
    tuple val(pair), path(model), path(sig_assignment, stageAs: 'sig/assignment'), path(sig_profiles, stageAs: 'sig/matrix_profiles')

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
    python3 ${projectDir}/bin/render_html.py \\
        --model ${model} \\
        --out   ${pair}_molecular_report.html \\
        --txt   ${pair}_molecular_report.txt
    """
}

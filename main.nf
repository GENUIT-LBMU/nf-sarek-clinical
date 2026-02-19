#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-sarek-clinical
    =================
    Clinical interpretation postprocessing pipeline for nf-core/sarek WES cancer exomes.

    Usage:
      nextflow run nf-sarek-clinical/main.nf \\
        -profile docker \\
        --sarek_outdir /path/to/sarek/results \\
        --run_name     MY_RUN_2024 \\
        --dbnsfp       /data/dbNSFP4.1a.txt.gz \\
        --snpsift_jar  /opt/snpEff/SnpSift.jar \\
        --panel_bed    /data/hg38/WES_panel.bed \\
        --outdir       /results/MY_RUN_2024_clinical

    GitHub: https://github.com/your-org/nf-sarek-clinical
    Documentation: https://nf-sarek-clinical.readthedocs.io
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK_CLINICAL } from './workflows/sarek_clinical'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    SAREK_CLINICAL()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    def msg    = """
    ====================================================
     nf-sarek-clinical run complete: ${status}
    ====================================================
     Run name  : ${workflow.runName}
     Duration  : ${workflow.duration}
     CPU hours : ${workflow.stats.computeTimeFmt ?: 'N/A'}
     Output    : ${params.outdir}
    ====================================================
    """.stripIndent()
    log.info(msg)
}

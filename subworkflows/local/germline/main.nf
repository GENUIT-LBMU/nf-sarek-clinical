/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GERMLINE POSTPROCESSING SUBWORKFLOW
    =====================================
    Takes Sarek HaplotypeCaller annotated VCFs and runs:
      - Germline SNV: dbNSFP → gene+impact filter → TSV extraction
      - Pharmacogenomics: rsID filter → TSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNPSIFT_DBNSFP_GERMLINE  } from '../../../modules/local/snpsift_germline/main'
include { SNPSIFT_FILTER_GERMLINE  } from '../../../modules/local/snpsift_germline/main'
include { SNPSIFT_EXTRACT_GERMLINE } from '../../../modules/local/snpsift_germline/main'
include { PHARMACO_FILTER          } from '../../../modules/local/pharmaco/main'

workflow GERMLINE_POSTPROCESS {
    take:
    ch_haplotypecaller_vcfs  // channel: [ val(sample), path(vcf.gz) ]
    ch_dbnsfp                // path: dbNSFP database (.txt.gz)
    ch_dbnsfp_tbi            // path: dbNSFP tabix index (.txt.gz.tbi)
    ch_snpsift_jar           // path: SnpSift.jar
    ch_germline_panel        // path: gene panel file
    ch_pharmaco_variants     // path: pharmacogenomics rsID file

    main:
    // ── Germline SNV pipeline ─────────────────────────────────────────────────
    SNPSIFT_DBNSFP_GERMLINE(
        ch_haplotypecaller_vcfs,
        ch_dbnsfp,
        ch_dbnsfp_tbi,
        ch_snpsift_jar
    )

    SNPSIFT_FILTER_GERMLINE(
        SNPSIFT_DBNSFP_GERMLINE.out.vcf,
        ch_germline_panel,
        ch_snpsift_jar,
        params.filters.germline_snv.min_dp
    )

    SNPSIFT_EXTRACT_GERMLINE(
        SNPSIFT_FILTER_GERMLINE.out.vcf,
        ch_snpsift_jar,
        file("${projectDir}/bin/acmg_shared.py")
    )

    // ── Pharmacogenomics ──────────────────────────────────────────────────────
    PHARMACO_FILTER(
        SNPSIFT_DBNSFP_GERMLINE.out.vcf,
        ch_snpsift_jar,
        ch_pharmaco_variants
    )

    emit:
    snv_tsv    = SNPSIFT_EXTRACT_GERMLINE.out.tsv       // [ val(sample), path(*.germline_snvs.tsv) ]
    acmg_tsv   = SNPSIFT_EXTRACT_GERMLINE.out.acmg_tsv  // [ val(sample), path(*.germline_snvs_acmg_annotated.tsv) ]
    pharmaco   = PHARMACO_FILTER.out.tsv                // [ val(sample), path(*.pharmaco.tsv) ]
}

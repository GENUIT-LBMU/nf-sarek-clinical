/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GERMLINE POSTPROCESSING SUBWORKFLOW
    =====================================
    Takes Sarek HaplotypeCaller annotated VCFs and runs:
      - Germline SNV: [snpEff MANE] → [local gnomAD] → [local ClinVar] → dbNSFP → gene+impact filter → TSV extraction
      - Pharmacogenomics: rsID filter → TSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNPEFF_MANE              } from '../../../modules/local/snpeff_mane/main'
include { GNOMAD_ANNOTATE          } from '../../../modules/local/gnomad_annotate/main'
include { CLINVAR_ANNOTATE         } from '../../../modules/local/clinvar_annotate/main'
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
    ch_gnomad                // list: [tsv.gz, tsv.gz.tbi, regions.bed.gz, regions.bed.gz.tbi] or []
    ch_clinvar               // list: [tsv.gz, tsv.gz.tbi] or []
    ch_snpeff_mane           // list: [snpEff data dir, database name] or []

    main:
    // ── Local gnomAD population frequencies (optional) ───────────────────────
    ch_vcfs_for_dbnsfp = ch_haplotypecaller_vcfs

    // ── MANE Select re-annotation (optional) ─────────────────────────────────
    if (ch_snpeff_mane) {
        SNPEFF_MANE(
            ch_vcfs_for_dbnsfp,
            ch_snpeff_mane[0],
            ch_snpeff_mane[1]
        )
        ch_vcfs_for_dbnsfp = SNPEFF_MANE.out.vcf
    }

    if (ch_gnomad) {
        GNOMAD_ANNOTATE(
            ch_vcfs_for_dbnsfp,
            ch_gnomad[0],
            ch_gnomad[1],
            ch_gnomad[2],
            ch_gnomad[3]
        )
        ch_vcfs_for_dbnsfp = GNOMAD_ANNOTATE.out.vcf
    }

    // ── Local ClinVar classifications (optional) ─────────────────────────────
    if (ch_clinvar) {
        CLINVAR_ANNOTATE(
            ch_vcfs_for_dbnsfp,
            ch_clinvar[0],
            ch_clinvar[1]
        )
        ch_vcfs_for_dbnsfp = CLINVAR_ANNOTATE.out.vcf
    }

    // ── Germline SNV pipeline ─────────────────────────────────────────────────
    SNPSIFT_DBNSFP_GERMLINE(
        ch_vcfs_for_dbnsfp,
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

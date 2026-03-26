/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC POSTPROCESSING SUBWORKFLOW
    =====================================
    Takes Sarek Mutect2 annotated VCFs + ASCAT/CNVkit outputs and runs:
      - Somatic SNV: dbNSFP → gene+PASS+impact filter → CIViC → TSV
      - Somatic CNV (ASCAT+CNVkit): Python interpretation → TSV + TXT
      - TMB: PASS+impact+DP filter → count → TMB per Mb
      - Mutational signatures: SigProfilerMatrixGenerator + SigProfilerAssignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNPSIFT_DBNSFP_SOMATIC   } from '../../../modules/local/snpsift_somatic/main'
include { SNPSIFT_FILTER_SOMATIC   } from '../../../modules/local/snpsift_somatic/main'
include { CIVIC_ANNOTATE           } from '../../../modules/local/snpsift_somatic/main'
include { SNPSIFT_EXTRACT_SOMATIC  } from '../../../modules/local/snpsift_somatic/main'
include { CIVIC_EXTRACT            } from '../../../modules/local/snpsift_somatic/main'
include { SOMATIC_CNV_INTERPRET    } from '../../../modules/local/somatic_cnv_interpret/main'
include { TMB_FILTER               } from '../../../modules/local/tmb/main'
include { TMB_CALCULATE            } from '../../../modules/local/tmb/main'
include { MUTSIG_FILTER            } from '../../../modules/local/mutational_signatures/main'
include { MUTSIG_SIGPROFILER       } from '../../../modules/local/mutational_signatures/main'
include { HRD_PROFILER             } from '../../../modules/local/hrd/main'

workflow SOMATIC_POSTPROCESS {
    take:
    ch_mutect2_vcfs     // channel: [ val(pair), path(vcf.gz) ]  snpEff-annotated Mutect2 VCFs
    ch_cnvkit_cns       // channel: [ val(pair), val(tumor), val(normal), path(*.call.cns) ]
    ch_ascat            // channel: [ val(pair), path(segments.txt), path(metrics.txt) ]
    ch_dbnsfp           // path: dbNSFP database (.txt.gz)
    ch_dbnsfp_tbi       // path: dbNSFP tabix index (.txt.gz.tbi)
    ch_snpsift_jar      // path: SnpSift.jar
    ch_somatic_snv_panel // path: somatic SNV gene panel
    ch_somatic_cnv_panel // path: somatic CNV gene panel (38 genes)
    ch_germline_panel   // path: germline gene panel (90 genes, for LOH)
    ch_gene_coords_bed  // path: WES gene coordinates BED
    ch_panel_bed        // path: WES capture BED (for TMB)
    ch_tumor_samples    // channel: [ val(pair), val(tumor_sample) ] for mutsig tumor column

    main:
    // ── Somatic SNV ──────────────────────────────────────────────────────────
    SNPSIFT_DBNSFP_SOMATIC(
        ch_mutect2_vcfs,
        ch_dbnsfp,
        ch_dbnsfp_tbi,
        ch_snpsift_jar
    )

    SNPSIFT_FILTER_SOMATIC(
        SNPSIFT_DBNSFP_SOMATIC.out.vcf,
        ch_somatic_snv_panel,
        ch_snpsift_jar,
        params.filters.somatic_snv.min_dp
    )

    CIVIC_ANNOTATE(
        SNPSIFT_FILTER_SOMATIC.out.vcf,
        params.genome_build
    )

    SNPSIFT_EXTRACT_SOMATIC(
        SNPSIFT_FILTER_SOMATIC.out.vcf,
        ch_snpsift_jar
    )

    CIVIC_EXTRACT(
        CIVIC_ANNOTATE.out.vcf
    )

    // ── Somatic CNV (ASCAT + CNVkit) ─────────────────────────────────────────
    // Join CNVkit .cns with ASCAT segments+metrics per pair
    ch_cnv_joined = ch_cnvkit_cns.join(
        ch_ascat.map { pair, seg, met -> [pair, seg, met] },
        by: 0, remainder: true
    ).map { pair, tumor, normal, cns, seg, met ->
        [pair, tumor, normal, cns,
         seg  ?: file('NO_FILE'),
         met  ?: file('NO_FILE')]
    }

    SOMATIC_CNV_INTERPRET(
        ch_cnv_joined,
        ch_germline_panel,
        ch_somatic_cnv_panel,
        ch_gene_coords_bed ?: file('NO_FILE')
    )

    // ── TMB ───────────────────────────────────────────────────────────────────
    TMB_FILTER(
        ch_mutect2_vcfs,
        ch_snpsift_jar,
        params.filters.somatic_snv.min_dp,
        params.filters.somatic_snv.min_af,
        params.filters.somatic_snv.min_alt_reads
    )

    TMB_CALCULATE(
        TMB_FILTER.out.vcf,
        ch_panel_bed ?: file('NO_FILE')
    )

    // ── Mutational Signatures ────────────────────────────────────────────────
    // Merge tumor sample name with VCF channel
    ch_mutsig_input = ch_mutect2_vcfs.join(ch_tumor_samples, by: 0)
        .map { pair, vcf, tumor_sample -> [pair, tumor_sample, vcf] }

    MUTSIG_FILTER(
        ch_mutsig_input,
        params.filters.somatic_snv.min_alt_reads  // reuse ALT_MIN
    )

    // Collect all filtered VCFs into a staging directory for SigProfiler
    ch_vcf_dir = MUTSIG_FILTER.out.vcf
        .map { pair, vcf -> vcf }
        .collect()
        .map { vcfs ->
            def dir = file("${workDir}/mutsig_vcfs")
            dir.mkdirs()
            vcfs.each { it.copyTo(dir) }
            dir
        }

    MUTSIG_SIGPROFILER(
        params.run_name,
        ch_vcf_dir,
        params.genome_build
    )

    // ── HRD Profiler ─────────────────────────────────────────────────────────
    ch_hrd_tsv = Channel.empty()
    if (!params.skip_hrd) {
        // Join filtered VCF [pair, vcf] with ASCAT segments [pair, seg, met] → [pair, vcf, seg]
        ch_hrd_input = MUTSIG_FILTER.out.vcf
            .join(
                ch_ascat.map { pair, seg, met -> [pair, seg] },
                by: 0
            )
            .filter { pair, vcf, seg -> seg.name != 'NO_FILE' }

        HRD_PROFILER(
            ch_hrd_input,
            params.hrd_organ,
            params.genome_build
        )
        ch_hrd_tsv = HRD_PROFILER.out.predictions
    }

    emit:
    somatic_snv_tsv   = SNPSIFT_EXTRACT_SOMATIC.out.tsv    // [ pair, path(*.somatic_snvs.tsv) ]
    civic_tsv         = CIVIC_EXTRACT.out.tsv               // [ pair, path(*.civic.exploded.tsv) ]
    cnv_tsv           = SOMATIC_CNV_INTERPRET.out.tsv       // [ pair, path(*_somatic_cnv_report.tsv) ]
    cnv_txt           = SOMATIC_CNV_INTERPRET.out.txt       // [ pair, path(*_somatic_cnv_report.txt) ]
    tmb_tsv           = TMB_CALCULATE.out.tsv               // [ pair, path(*.tmb.tsv) ]
    signatures_dir    = MUTSIG_SIGPROFILER.out.assignment   // path(assignment/)
    hrd_tsv           = ch_hrd_tsv                          // [ pair, path(*_hrd_predictions.tsv) ]
}

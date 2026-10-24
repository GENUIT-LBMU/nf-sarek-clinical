/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC POSTPROCESSING SUBWORKFLOW
    =====================================
    Takes Sarek Mutect2 + Strelka2 annotated VCFs + ASCAT/CNVkit outputs and runs:
      - Consensus: Mutect2 ∩ Strelka2 intersection (when Strelka2 available)
      - Somatic SNV: dbNSFP → gene+PASS+impact filter → CIViC → TSV
      - Somatic CNV (ASCAT+CNVkit): Python interpretation → TSV + TXT
      - TMB: PASS+impact+DP filter → count → TMB per Mb (from consensus VCF)
      - Mutational signatures: SigProfilerMatrixGenerator + SigProfilerAssignment (from consensus VCF)
      - HRD: HRProfiler (from consensus VCF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNPSIFT_DBNSFP_SOMATIC   } from '../../../modules/local/snpsift_somatic/main'
include { SNPSIFT_FILTER_SOMATIC   } from '../../../modules/local/snpsift_somatic/main'
include { CIVIC_ANNOTATE           } from '../../../modules/local/snpsift_somatic/main'
include { SNPSIFT_EXTRACT_SOMATIC  } from '../../../modules/local/snpsift_somatic/main'
include { ANNOTATE_CONSENSUS       } from '../../../modules/local/snpsift_somatic/main'
include { CIVIC_EXTRACT            } from '../../../modules/local/snpsift_somatic/main'
include { SOMATIC_CONSENSUS        } from '../../../modules/local/somatic_consensus/main'
include { SOMATIC_CNV_INTERPRET    } from '../../../modules/local/somatic_cnv_interpret/main'
include { TMB_FILTER               } from '../../../modules/local/tmb/main'
include { TMB_CALCULATE            } from '../../../modules/local/tmb/main'
include { MUTSIG_FILTER            } from '../../../modules/local/mutational_signatures/main'
include { MUTSIG_SIGPROFILER       } from '../../../modules/local/mutational_signatures/main'
include { HRD_PROFILER             } from '../../../modules/local/hrd/main'

workflow SOMATIC_POSTPROCESS {
    take:
    ch_mutect2_vcfs     // channel: [ val(pair), path(vcf.gz) ]  snpEff-annotated Mutect2 VCFs
    ch_strelka2_snvs    // channel: [ val(pair), path(vcf.gz) ]  Strelka2 somatic SNVs (may be empty)
    ch_strelka2_indels  // channel: [ val(pair), path(vcf.gz) ]  Strelka2 somatic indels (may be empty)
    ch_cnvkit_cns       // channel: [ val(pair), val(tumor), val(normal), path(*.call.cns) ]
    ch_ascat            // channel: [ val(pair), path(segments.txt), path(metrics.txt) ]
    ch_dbnsfp           // path: dbNSFP database (.txt.gz)
    ch_dbnsfp_tbi       // path: dbNSFP tabix index (.txt.gz.tbi)
    ch_snpsift_jar      // path: SnpSift.jar
    ch_somatic_snv_panel // path: somatic SNV gene panel
    ch_somatic_cnv_panel // path: somatic CNV gene panel (38 genes)
    ch_germline_panel   // path: germline gene panel (105 genes, for LOH)
    ch_gene_coords_bed  // path: WES gene coordinates BED
    ch_panel_bed        // path: WES capture BED (for TMB)
    ch_tumor_samples    // channel: [ val(pair), val(tumor_sample) ] for mutsig tumor column
    ch_fasta            // path: reference FASTA (for bcftools norm)
    ch_fasta_fai        // path: reference FASTA index

    main:

    // ── Consensus VCF (Mutect2 ∩ Strelka2) ─────────────────────────────────
    // Join Mutect2 + Strelka2 SNVs + Strelka2 indels per pair.
    // If Strelka2 is not available (empty channels), fall back to Mutect2 only.
    ch_consensus_input = ch_mutect2_vcfs
        .join(ch_strelka2_snvs, by: 0, remainder: true)
        .join(ch_strelka2_indels, by: 0, remainder: true)

    // Split into pairs that have Strelka2 vs those that don't
    ch_has_strelka = ch_consensus_input
        .filter { pair, mutect2, strelka_snvs, strelka_indels ->
            strelka_snvs != null && strelka_indels != null
        }
    ch_no_strelka = ch_consensus_input
        .filter { pair, mutect2, strelka_snvs, strelka_indels ->
            strelka_snvs == null || strelka_indels == null
        }
        .map { pair, mutect2, strelka_snvs, strelka_indels -> [pair, mutect2] }

    // Run consensus on pairs with both callers
    SOMATIC_CONSENSUS(
        ch_has_strelka,
        ch_fasta,
        ch_fasta_fai,
        params.filters.consensus.min_callers
    )

    // Merge: consensus VCFs + fallback Mutect2-only VCFs
    ch_biomarker_vcfs = SOMATIC_CONSENSUS.out.vcf.mix(ch_no_strelka)

    // ── Somatic SNV (clinical report — always from Mutect2, full annotations) ──
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

    // Annotate somatic TSV with consensus caller column
    // Join extract TSV [pair, tsv] with consensus VCF [pair, vcf]
    // For Mutect2-only pairs (no Strelka2), use an empty VCF placeholder
    ch_extract_with_consensus = SNPSIFT_EXTRACT_SOMATIC.out.tsv
        .join(SOMATIC_CONSENSUS.out.vcf, by: 0, remainder: true)
        .map { pair, tsv, consensus_vcf ->
            [pair, tsv, consensus_vcf ?: file('NO_FILE')]
        }

    ANNOTATE_CONSENSUS(ch_extract_with_consensus)

    CIVIC_EXTRACT(
        CIVIC_ANNOTATE.out.vcf
    )

    // ── Somatic CNV (ASCAT + CNVkit) ─────────────────────────────────────────
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

    // ── TMB (from consensus VCF when available) ─────────────────────────────
    TMB_FILTER(
        ch_biomarker_vcfs,
        ch_snpsift_jar,
        params.filters.somatic_snv.min_dp,
        params.filters.somatic_snv.min_af,
        params.filters.somatic_snv.min_alt_reads
    )

    TMB_CALCULATE(
        TMB_FILTER.out.vcf,
        ch_panel_bed ?: file('NO_FILE')
    )

    // ── Mutational Signatures (from consensus VCF when available) ───────────
    ch_mutsig_input = ch_biomarker_vcfs.join(ch_tumor_samples, by: 0)
        .map { pair, vcf, tumor_sample -> [pair, tumor_sample, vcf] }

    MUTSIG_FILTER(
        ch_mutsig_input,
        params.filters.somatic_snv.min_alt_reads
    )

    ch_vcf_dir = MUTSIG_FILTER.out.vcf
        .map { pair, vcf -> vcf }
        .collect()
        .map { vcfs ->
            def dir = file("${workDir}/mutsig_vcfs")
            dir.mkdirs()
            vcfs.each { it.copyTo(dir) }
            dir
        }

    ch_signatures = Channel.empty()
    if (!params.skip_signatures) {
        MUTSIG_SIGPROFILER(
            params.run_name,
            ch_vcf_dir,
            params.genome_build
        )
        ch_signatures = MUTSIG_SIGPROFILER.out.assignment
    }

    // ── HRD Profiler ─────────────────────────────────────────────────────────
    ch_hrd_tsv = Channel.empty()
    if (!params.skip_hrd) {
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
    somatic_snv_tsv   = ANNOTATE_CONSENSUS.out.tsv
    civic_tsv         = CIVIC_EXTRACT.out.tsv
    cnv_tsv           = SOMATIC_CNV_INTERPRET.out.tsv
    cnv_txt           = SOMATIC_CNV_INTERPRET.out.txt
    tmb_tsv           = TMB_CALCULATE.out.tsv
    signatures_dir    = ch_signatures
    hrd_tsv           = ch_hrd_tsv
    consensus_stats   = SOMATIC_CONSENSUS.out.stats
}

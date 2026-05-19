/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNV SOMATIC GATK SUBWORKFLOW
    ==============================
    Orchestrates the 8-step somatic CNV paired workflow.

    Key design:
      • Steps 1–3 (CollectReadCounts, DenoiseReadCounts, CollectAllelicCounts)
        run independently per sample (tumour AND normal) in parallel.
      • Step 4 (ModelSegments) requires pairing:
        tumour.denoisedCR  +  tumour.allelicCounts  +  normal.allelicCounts
        → joined on patient_id (extracted from sample name by _N[0-9]+/_T[0-9]+).
      • Steps 5–8 run per tumour only.

    Input BAM naming convention (same as Sarek/launch_nt_V6.sh):
      Normal: {PATIENT_ID}_N{XX}   (e.g. 682_N01)
      Tumour: {PATIENT_ID}_T{XX}   (e.g. 682_T01)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK_COLLECT_READ_COUNTS_SOMATIC } from '../../../modules/local/cnv_somatic_gatk/main'
include { GATK_DENOISE_READ_COUNTS         } from '../../../modules/local/cnv_somatic_gatk/main'
include { GATK_COLLECT_ALLELIC_COUNTS      } from '../../../modules/local/cnv_somatic_gatk/main'
include { GATK_MODEL_SEGMENTS              } from '../../../modules/local/cnv_somatic_gatk/main'
include { GATK_CALL_COPY_RATIO_SEGMENTS    } from '../../../modules/local/cnv_somatic_gatk/main'
include { GATK_PLOT_DENOISED_COPY_RATIOS   } from '../../../modules/local/cnv_somatic_gatk/main'
include { GATK_PLOT_MODELED_SEGMENTS       } from '../../../modules/local/cnv_somatic_gatk/main'
include { MANE_ANNOTATE_SOMATIC_CNV        } from '../../../modules/local/cnv_somatic_gatk/main'

workflow CNV_SOMATIC_GATK {
    take:
    ch_bams             // channel: [ val(sample), path(bam), path(bai) ]
    // sample names MUST follow {PATIENT_ID}_N{XX} or {PATIENT_ID}_T{XX} convention
    ch_interval_list    // path: preprocessed.interval_list
    ch_somatic_pon      // path: somatic PoN HDF5
    ch_common_snps      // path: common SNPs VCF (bgzipped)
    ch_common_snps_tbi  // path: common SNPs VCF .tbi
    ch_fasta            // path: reference FASTA (main chroms)
    ch_fasta_fai        // path: reference .fai
    ch_fasta_dict       // path: reference .dict
    ch_mane_txt         // path: exons_mane.txt
    // annotate_cnvs_complete_fixed.R is bundled in bin/ (from FlorPio/gatk-cnv-somatic)

    main:
    // ── Tag each sample with patient_id and status (normal/tumor) ─────────────
    ch_bams_tagged = ch_bams.map { sample, bam, bai ->
        def patient_id = sample.replaceAll(/_N\d*$/, '').replaceAll(/_T\d*$/, '')
        def status     = sample ==~ /.*_N\d*$/ ? 'normal' : 'tumor'
        [sample, patient_id, status, bam, bai]
    }

    // ── Step 1: CollectReadCounts (all samples in parallel) ───────────────────
    GATK_COLLECT_READ_COUNTS_SOMATIC(
        ch_bams_tagged,
        ch_interval_list,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict
    )

    // ── Step 2: DenoiseReadCounts (all samples in parallel) ───────────────────
    GATK_DENOISE_READ_COUNTS(
        GATK_COLLECT_READ_COUNTS_SOMATIC.out.counts,
        ch_somatic_pon
    )

    // ── Step 3: CollectAllelicCounts (all samples in parallel) ────────────────
    GATK_COLLECT_ALLELIC_COUNTS(
        ch_bams_tagged,
        ch_common_snps,
        ch_common_snps_tbi,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict
    )

    // ── Step 4: ModelSegments (per pair) ──────────────────────────────────────
    // We need: [patient_id, tumor_name, normal_name,
    //           tumor.denoisedCR, tumor.allelicCounts, normal.allelicCounts]

    // Separate denoised CR by status
    ch_tumor_denoised = GATK_DENOISE_READ_COUNTS.out.denoised
        .filter { sample, patient_id, status, std_cr, den_cr -> status == 'tumor' }
        .map    { sample, patient_id, status, std_cr, den_cr ->
            [patient_id, sample, den_cr, std_cr]   // [patient, tumor_name, denoisedCR, standardizedCR]
        }

    // Separate allelic counts by status
    ch_tumor_allelic = GATK_COLLECT_ALLELIC_COUNTS.out.allelic_counts
        .filter { sample, patient_id, status, tsv -> status == 'tumor' }
        .map    { sample, patient_id, status, tsv -> [patient_id, sample, tsv] }

    ch_normal_allelic = GATK_COLLECT_ALLELIC_COUNTS.out.allelic_counts
        .filter { sample, patient_id, status, tsv -> status == 'normal' }
        .map    { sample, patient_id, status, tsv -> [patient_id, sample, tsv] }

    // Join tumor denoised + tumor allelic on [patient_id, tumor_name]
    ch_tumor_merged = ch_tumor_denoised
        .join(ch_tumor_allelic, by: [0, 1])
        // → [patient_id, tumor_name, denoisedCR, standardizedCR, tumor_allelicCounts]

    // Join with normal allelic on patient_id only (normal sample name discarded here)
    ch_model_input = ch_tumor_merged
        .join(
            ch_normal_allelic.map { patient_id, normal_name, tsv -> [patient_id, normal_name, tsv] },
            by: 0
        )
        .map { patient_id, tumor, den_cr, std_cr, tumor_allelic, normal_name, normal_allelic ->
            [patient_id, tumor, normal_name, den_cr, tumor_allelic, normal_allelic]
        }
        // → [patient_id, tumor_name, normal_name, denoisedCR, tumor_allelic, normal_allelic]

    GATK_MODEL_SEGMENTS(ch_model_input)

    // ── Step 5: CallCopyRatioSegments (per tumor) ────────────────────────────
    GATK_CALL_COPY_RATIO_SEGMENTS(
        GATK_MODEL_SEGMENTS.out.cr_seg
    )

    // ── Step 6: PlotDenoisedCopyRatios (per tumor) ──────────────────────────
    // Pair standardizedCR + denoisedCR by [patient_id, tumor_name]
    ch_std_den_cr = ch_tumor_denoised
        .map { patient_id, tumor, den_cr, std_cr -> [patient_id, tumor, std_cr, den_cr] }

    GATK_PLOT_DENOISED_COPY_RATIOS(
        ch_std_den_cr,
        ch_fasta_dict,
        params.somatic_min_contig_length
    )

    // ── Step 7: PlotModeledSegments (per tumor) ──────────────────────────────
    // Need: denoised_cr, hets_tsv, model_final_seg → join on [patient_id, tumor]
    ch_plot_modeled_input = ch_tumor_denoised
        .map { patient_id, tumor, den_cr, std_cr -> [patient_id, tumor, den_cr] }
        .join(GATK_MODEL_SEGMENTS.out.hets_tsv,       by: [0, 1])
        .join(GATK_MODEL_SEGMENTS.out.model_final_seg.map { p, t, n, f -> [p, t, f] }, by: [0, 1])
        // → [patient_id, tumor, denoisedCR, hets_tsv, modelFinal.seg]

    GATK_PLOT_MODELED_SEGMENTS(
        ch_plot_modeled_input,
        ch_fasta_dict,
        params.somatic_min_contig_length
    )

    // ── Step 8: MANE Annotation (per tumor) ──────────────────────────────────
    // Gather all 4 ModelSegments outputs + called.seg → join on [patient_id, tumor]
    ch_annotate_input = GATK_MODEL_SEGMENTS.out.model_final_seg
        .map { patient_id, tumor, normal, seg -> [patient_id, tumor, seg] }
        .join(GATK_MODEL_SEGMENTS.out.cr_seg,              by: [0, 1])
        .join(GATK_MODEL_SEGMENTS.out.hets_tsv,            by: [0, 1])
        .join(GATK_CALL_COPY_RATIO_SEGMENTS.out.called_seg, by: [0, 1])
        // → [patient_id, tumor, modelFinal.seg, cr.seg, hets.tsv, called.seg]

    MANE_ANNOTATE_SOMATIC_CNV(
        ch_annotate_input,
        ch_mane_txt
    )

    emit:
    called_seg      = GATK_CALL_COPY_RATIO_SEGMENTS.out.called_seg  // [patient_id, tumor, called.seg]
    model_final_seg = GATK_MODEL_SEGMENTS.out.model_final_seg        // [patient_id, tumor, normal, modelFinal.seg]
    annotated_tsv   = MANE_ANNOTATE_SOMATIC_CNV.out.annotated_tsv   // [patient_id, tumor, annotated.tsv]
    cr_plots        = GATK_PLOT_DENOISED_COPY_RATIOS.out.plots
    modeled_plots   = GATK_PLOT_MODELED_SEGMENTS.out.plots
}

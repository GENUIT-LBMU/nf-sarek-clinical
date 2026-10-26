/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNV GERMLINE GATK SUBWORKFLOW
    ================================
    Orchestrates the 5-step gCNV CASE-mode germline CNV calling pipeline.

    Input:  per-sample BAM (preferred) or CRAM from Sarek preprocessing/recalibrated/
    Output: annotated segments TSV per sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK_COLLECT_READ_COUNTS_GERMLINE      } from '../../../modules/local/cnv_germline_gatk/main'
include { GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY  } from '../../../modules/local/cnv_germline_gatk/main'
include { GATK_GERMLINE_CNV_CALLER               } from '../../../modules/local/cnv_germline_gatk/main'
include { GATK_POSTPROCESS_GERMLINE_CNV_CALLS    } from '../../../modules/local/cnv_germline_gatk/main'
include { MANE_ANNOTATE_GERMLINE_CNV             } from '../../../modules/local/cnv_germline_gatk/main'

workflow CNV_GERMLINE_GATK {
    take:
    ch_bams             // channel: [ val(sample), path(bam|cram), path(bai|crai) ]
    ch_interval_list    // path: preprocessed.interval_list
    ch_ploidy_model     // path: pre-built ploidy model directory
    ch_pon_dir          // path: pre-built gCNV PoN directory
    ch_fasta            // path: reference FASTA
    ch_fasta_fai        // path: reference FASTA .fai
    ch_fasta_dict       // path: reference .dict
    ch_mane_txt         // path: exons_mane.txt
    ch_cram_fasta       // path: FASTA the CRAMs were written with (gcnv_fasta if none)
    ch_cram_fai         // path: its .fai
    ch_cram_dict        // path: its .dict
    ch_gene_panel       // path: germline panel (genes annotated)
    exclude_genes       // val: comma-separated panel genes left out
    ch_lowmap_bed       // path: low-mappability BED (segments flagged, not hidden)
    // exon_annotation.R is bundled in bin/ (from FlorPio/gatk-cnv-somatic)

    main:
    // ── Step 1: CollectReadCounts ─────────────────────────────────────────────
    GATK_COLLECT_READ_COUNTS_GERMLINE(
        ch_bams,
        ch_interval_list,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_cram_fasta,
        ch_cram_fai,
        ch_cram_dict
    )

    // ── Step 2: DetermineGermlineContigPloidy (CASE, per sample) ─────────────
    GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY(
        GATK_COLLECT_READ_COUNTS_GERMLINE.out.counts,
        ch_ploidy_model
    )

    // ── Step 3: GermlineCNVCaller (CASE, per sample) ─────────────────────────
    // Join counts and ploidy on sample name
    def ch_cnvcaller_input = GATK_COLLECT_READ_COUNTS_GERMLINE.out.counts
        .join(GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY.out.ploidy_dir, by: 0)

    GATK_GERMLINE_CNV_CALLER(
        ch_cnvcaller_input,
        ch_pon_dir
    )

    // ── Step 4: PostprocessGermlineCNVCalls ───────────────────────────────────
    // Join cnv_calls_dir and ploidy_dir on sample name
    def ch_postprocess_input = GATK_GERMLINE_CNV_CALLER.out.cnv_calls_dir
        .join(GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY.out.ploidy_dir, by: 0)

    GATK_POSTPROCESS_GERMLINE_CNV_CALLS(
        ch_postprocess_input,
        ch_pon_dir,
        ch_fasta_dict
    )

    // ── Step 5: MANE Annotation ───────────────────────────────────────────────
    MANE_ANNOTATE_GERMLINE_CNV(
        GATK_POSTPROCESS_GERMLINE_CNV_CALLS.out.segments_vcf,
        ch_mane_txt,
        ch_gene_panel,
        exclude_genes,
        ch_lowmap_bed
    )

    emit:
    counts          = GATK_COLLECT_READ_COUNTS_GERMLINE.out.counts
    ploidy          = GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY.out.ploidy_dir
    denoised        = GATK_POSTPROCESS_GERMLINE_CNV_CALLS.out.denoised
    segments_vcf    = GATK_POSTPROCESS_GERMLINE_CNV_CALLS.out.segments_vcf
    annotated_tsv   = MANE_ANNOTATE_GERMLINE_CNV.out.annotated_tsv
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC CNV – GATK CNV SOMATIC PAIRED WORKFLOW
    ================================================
    Runs GATK somatic CNV calling on tumour-normal paired BAMs from Sarek.

    Steps:
      1. CollectReadCounts (T + N)     → sample.counts.hdf5
      2. DenoiseReadCounts (T + N)     → standardizedCR.tsv + denoisedCR.tsv
      3. CollectAllelicCounts (T + N)  → allelicCounts.tsv
      4. ModelSegments (T, with N)     → modelFinal.seg, cr.seg, hets.tsv
      5. CallCopyRatioSegments (T)     → called.seg
      6. PlotDenoisedCopyRatios (T)    → diagnostic PNG/PDF
      7. PlotModeledSegments (T)       → diagnostic PNG/PDF
      8. MANE_ANNOTATE_SOMATIC_CNV (T) → annotated TSV per patient (Rscript)

    Required external assets (not included in the repo):
      --somatic_interval_list   : preprocessed.interval_list
      --somatic_pon_hdf5        : somatic Panel-of-Normals HDF5
      --common_snps_vcf         : common SNPs VCF (bgzipped + tabix indexed)
      --somatic_fasta           : reference FASTA  (main chromosomes only)
      --somatic_dict            : reference .dict
      --cram_fasta              : FASTA (+ .fai, .dict) Sarek wrote the CRAMs with; only for CRAM input
      --gcnv_mane_txt           : MANE transcript exon coordinates
      --somatic_annotate_r      : annotate_cnvs_complete_fixed.R script path
      --somatic_min_contig_length : min contig length for plots (default: 46709983)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ─────────────────────────────────────────────────────────────────────────────
// 1. CollectReadCounts  (shared for tumor and normal – same command)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_COLLECT_READ_COUNTS_SOMATIC {
    tag "$sample"
    label 'process_medium'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/counts", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(patient_id), val(status), path(bam), path(bai)
    // status: 'normal' | 'tumor'
    path interval_list
    path fasta
    path fasta_fai
    path fasta_dict
    path cram_fasta, stageAs: 'cram_ref/*'   // reference of the CRAMs (BAM input: unused)
    path cram_fai,   stageAs: 'cram_ref/*'
    path cram_dict,  stageAs: 'cram_ref/*'

    output:
    tuple val(sample), val(patient_id), val(status),
          path("${sample}.counts.hdf5"), emit: counts

    script:
    // A CRAM is decoded with the reference it was written with; a BAM keeps somatic_fasta
    def ref = bam.name.endsWith('.cram') ? cram_fasta : fasta
    """
    gatk CollectReadCounts \\
        -I ${bam} \\
        -L ${interval_list} \\
        -R ${ref} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --format HDF5 \\
        -O ${sample}.counts.hdf5
    """

    stub:
    """
    touch ${sample}.counts.hdf5
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 2. DenoiseReadCounts  (tumor and normal independently)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_DENOISE_READ_COUNTS {
    tag "$sample"
    label 'process_medium'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/counts", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(patient_id), val(status), path(counts)
    path somatic_pon_hdf5

    output:
    tuple val(sample), val(patient_id), val(status),
          path("${sample}.standardizedCR.tsv"),
          path("${sample}.denoisedCR.tsv"),
          emit: denoised

    script:
    """
    gatk DenoiseReadCounts \\
        -I ${counts} \\
        --count-panel-of-normals ${somatic_pon_hdf5} \\
        --standardized-copy-ratios ${sample}.standardizedCR.tsv \\
        --denoised-copy-ratios     ${sample}.denoisedCR.tsv
    """

    stub:
    """
    touch ${sample}.standardizedCR.tsv ${sample}.denoisedCR.tsv
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 3. CollectAllelicCounts  (tumor and normal independently)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_COLLECT_ALLELIC_COUNTS {
    tag "$sample"
    label 'process_medium'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/counts", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(patient_id), val(status), path(bam), path(bai)
    path common_snps_vcf
    path common_snps_vcf_tbi
    path fasta
    path fasta_fai
    path fasta_dict
    path cram_fasta, stageAs: 'cram_ref/*'   // reference of the CRAMs (BAM input: unused)
    path cram_fai,   stageAs: 'cram_ref/*'
    path cram_dict,  stageAs: 'cram_ref/*'

    output:
    tuple val(sample), val(patient_id), val(status),
          path("${sample}.allelicCounts.tsv"), emit: allelic_counts

    script:
    def ref = bam.name.endsWith('.cram') ? cram_fasta : fasta
    """
    gatk CollectAllelicCounts \\
        -I ${bam} \\
        -L ${common_snps_vcf} \\
        -R ${ref} \\
        -O ${sample}.allelicCounts.tsv
    """

    stub:
    """
    touch ${sample}.allelicCounts.tsv
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 4. ModelSegments  (paired: tumor denoised CR + tumor allelic + normal allelic)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_MODEL_SEGMENTS {
    tag "${patient_id}: ${tumor} vs ${normal}"
    label 'process_medium'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/model_segments/${patient_id}",
               mode: params.publish_dir_mode

    input:
    tuple val(patient_id),
          val(tumor), val(normal),
          path(tumor_denoised_cr),    // ${tumor}.denoisedCR.tsv
          path(tumor_allelic),        // ${tumor}.allelicCounts.tsv
          path(normal_allelic)        // ${normal}.allelicCounts.tsv

    output:
    tuple val(patient_id), val(tumor), val(normal),
          path("${tumor}.modelFinal.seg"), emit: model_final_seg
    tuple val(patient_id), val(tumor),
          path("${tumor}.cr.seg"),         emit: cr_seg
    tuple val(patient_id), val(tumor),
          path("${tumor}.hets.tsv"),       emit: hets_tsv
    tuple val(patient_id), val(tumor),
          path("${tumor}.*"),              emit: all_outputs

    script:
    """
    gatk ModelSegments \\
        --denoised-copy-ratios   ${tumor_denoised_cr} \\
        --allelic-counts         ${tumor_allelic} \\
        --normal-allelic-counts  ${normal_allelic} \\
        --output-prefix          ${tumor} \\
        -O .
    """

    stub:
    """
    touch ${tumor}.modelFinal.seg ${tumor}.cr.seg ${tumor}.hets.tsv
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 5. CallCopyRatioSegments  (tumor only)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_CALL_COPY_RATIO_SEGMENTS {
    tag "${patient_id}: ${tumor}"
    label 'process_single'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/called_segments",
               mode: params.publish_dir_mode

    input:
    tuple val(patient_id), val(tumor), path(cr_seg)

    output:
    tuple val(patient_id), val(tumor), path("${tumor}.called.seg"), emit: called_seg

    script:
    """
    gatk CallCopyRatioSegments \\
        --input  ${cr_seg} \\
        --output ${tumor}.called.seg
    """

    stub:
    """
    touch ${tumor}.called.seg
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 6. PlotDenoisedCopyRatios  (tumor only – diagnostic)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_PLOT_DENOISED_COPY_RATIOS {
    tag "${patient_id}: ${tumor}"
    label 'process_low'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/plots/${patient_id}/CR_plots",
               mode: params.publish_dir_mode

    input:
    tuple val(patient_id), val(tumor),
          path(standardized_cr),   // ${tumor}.standardizedCR.tsv
          path(denoised_cr)        // ${tumor}.denoisedCR.tsv
    path fasta_dict
    val  min_contig_length

    output:
    tuple val(patient_id), val(tumor),
          path("${tumor}*.png"), optional: true, emit: plots
    tuple val(patient_id), val(tumor),
          path("${tumor}*.pdf"), optional: true, emit: plots_pdf

    script:
    """
    gatk PlotDenoisedCopyRatios \\
        --standardized-copy-ratios ${standardized_cr} \\
        --denoised-copy-ratios     ${denoised_cr} \\
        --sequence-dictionary      ${fasta_dict} \\
        --minimum-contig-length    ${min_contig_length} \\
        --output       . \\
        --output-prefix ${tumor}
    """

    stub:
    """
    touch ${tumor}.png
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 7. PlotModeledSegments  (tumor only – diagnostic)
// ─────────────────────────────────────────────────────────────────────────────
process GATK_PLOT_MODELED_SEGMENTS {
    tag "${patient_id}: ${tumor}"
    label 'process_low'
    label 'gatk'

    publishDir "${params.outdir}/somatic/cnv_gatk/plots/${patient_id}/Modeled_segments_plots",
               mode: params.publish_dir_mode

    input:
    tuple val(patient_id), val(tumor),
          path(denoised_cr),    // ${tumor}.denoisedCR.tsv
          path(hets_tsv),       // ${tumor}.hets.tsv
          path(model_final_seg) // ${tumor}.modelFinal.seg
    path fasta_dict
    val  min_contig_length

    output:
    tuple val(patient_id), val(tumor),
          path("${tumor}*.png"), optional: true, emit: plots
    tuple val(patient_id), val(tumor),
          path("${tumor}*.pdf"), optional: true, emit: plots_pdf

    script:
    """
    gatk PlotModeledSegments \\
        --denoised-copy-ratios ${denoised_cr} \\
        --allelic-counts       ${hets_tsv} \\
        --segments             ${model_final_seg} \\
        --sequence-dictionary  ${fasta_dict} \\
        --minimum-contig-length ${min_contig_length} \\
        --output       . \\
        --output-prefix ${tumor}
    """

    stub:
    """
    touch ${tumor}.png
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 8. MANE Annotation  (per patient, uses model_segments directory)
// ─────────────────────────────────────────────────────────────────────────────
process MANE_ANNOTATE_SOMATIC_CNV {
    tag "${patient_id}: ${tumor}"
    label 'process_low'
    label 'r_base'

    publishDir "${params.outdir}/somatic/cnv_gatk/annotated/${patient_id}",
               mode: params.publish_dir_mode

    input:
    tuple val(patient_id), val(tumor),
          path(model_final_seg),   // ${tumor}.modelFinal.seg
          path(cr_seg),            // ${tumor}.cr.seg
          path(hets_tsv),          // ${tumor}.hets.tsv
          path(called_seg)         // ${tumor}.called.seg
    path  mane_txt                 // exons_mane.txt
    path  cnv_genes                // genes to annotate (germline + somatic CNV panels)

    output:
    tuple val(patient_id), val(tumor),
          path("${tumor}_annotated.tsv"),       emit: annotated_tsv
    tuple val(patient_id), val(tumor),
          path("${tumor}_annotated*"),           emit: all_outputs, optional: true

    script:
    // The R script expects --input_dir with ModelSegments outputs, --mane, --output_dir
    """
    mkdir -p seg_input output

    # Stage all ModelSegments files into the expected input directory
    cp ${model_final_seg} ${cr_seg} ${hets_tsv} ${called_seg} seg_input/

    Rscript ${projectDir}/bin/annotate_cnvs_complete_fixed.R \\
        --input_dir seg_input \\
        --mane      ${mane_txt} \\
        --genes     ${cnv_genes} \\
        --output_dir output || echo "WARNING: R annotation script exited non-zero (files may still be present)" >&2

    # Use find to locate annotated output safely (find exits 0 even with no matches,
    # unlike ls glob which exits 1 under set -e when no .tsv files exist)
    # Prefer sample-specific file first, then any annotated file
    annotated="\$(find output/ -maxdepth 1 -type f -name '${tumor}*annotated*' 2>/dev/null | head -1)"
    if [[ -z "\$annotated" ]]; then
        annotated="\$(find output/ -maxdepth 1 -type f \\( -name '*annotated*.tsv' -o -name '*annotated*.txt' \\) 2>/dev/null | sort | head -1)"
    fi
    if [[ -n "\$annotated" ]]; then
        cp "\$annotated" ${tumor}_annotated.tsv
    else
        # Fallback: copy called.seg as TSV if R script produced no annotated output
        echo "WARNING: annotate_cnvs_complete_fixed.R produced no annotated TSV, using fallback" >&2
        cp ${called_seg} ${tumor}_annotated.tsv
    fi

    # Copy any additional outputs (plots, extra TSVs)
    for f in output/*; do
        [[ -f "\$f" ]] && cp "\$f" "${tumor}_annotated_\$(basename \$f)" || true
    done
    """

    stub:
    """
    touch ${tumor}_annotated.tsv
    """
}

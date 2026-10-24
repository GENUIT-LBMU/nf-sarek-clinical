/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GERMLINE CNV – GATK gCNV CASE MODE
    =====================================
    Runs per-sample gCNV calling using pre-built Panel-of-Normals (PoN) and
    ploidy model (both generated offline with COHORT mode on control samples).

    Steps:
      1. CollectReadCounts          → sample.counts.hdf5
      2. DetermineGermlineContigPloidy (CASE) → ploidy calls directory
      3. GermlineCNVCaller (CASE)   → CNV calls directory
      4. PostprocessGermlineCNVCalls → denoised.tsv + intervals.vcf.gz + segments.vcf.gz
      5. MANE_ANNOTATE_GERMLINE_CNV → annotated TSV per sample (Rscript)

    Required external assets (not included in the repo):
      --gcnv_interval_list   : preprocessed.interval_list  (from CreateReadCountPanelOfNormals)
      --gcnv_ploidy_model    : directory with pre-built ploidy model
      --gcnv_pon_dir         : directory with pre-built gCNV PoN model
      --gcnv_fasta           : reference FASTA (hg38.fa)
      --cram_fasta           : FASTA (+ .fai, .dict) Sarek wrote the CRAMs with; only for CRAM input
      --gcnv_dict            : reference .dict (hg38.dict)
      --gcnv_mane_txt        : MANE transcript exon coordinates (exons_mane.txt)
      --gcnv_exon_r          : exon_annotation.R script path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GATK_COLLECT_READ_COUNTS_GERMLINE {
    tag "$sample"
    label 'process_medium'
    label 'gatk'

    publishDir "${params.outdir}/germline/cnv_gatk/counts", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(bam), path(bai)
    path interval_list
    path fasta
    path fasta_fai
    path fasta_dict
    path cram_fasta, stageAs: 'cram_ref/*'   // reference of the CRAMs (BAM input: unused)
    path cram_fai,   stageAs: 'cram_ref/*'
    path cram_dict,  stageAs: 'cram_ref/*'

    output:
    tuple val(sample), path("${sample}.counts.hdf5"), emit: counts

    script:
    // A CRAM is decoded with the reference it was written with; a BAM keeps gcnv_fasta
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


process GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY {
    tag "$sample"
    label 'process_high'
    label 'gatk'

    publishDir "${params.outdir}/germline/cnv_gatk/ploidy", mode: params.publish_dir_mode,
               saveAs: { filename -> filename.startsWith(sample) ? filename : null }

    input:
    tuple val(sample), path(counts)
    path ploidy_model_dir       // pre-built ploidy model directory

    output:
    tuple val(sample), path("${sample}_ploidy"), emit: ploidy_dir
    // ploidy_dir contains: ${sample}-calls/, ${sample}-model/ (if COHORT)

    script:
    """
    export THEANO_FLAGS="base_compiledir=\${PWD}/.theano"
    export HOME=\${PWD}
    mkdir -p ${sample}_ploidy

    gatk DetermineGermlineContigPloidy \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --model         ${ploidy_model_dir} \\
        --input         ${counts} \\
        --output        ${sample}_ploidy \\
        --output-prefix ${sample}
    """

    stub:
    """
    mkdir -p ${sample}_ploidy/${sample}-calls
    touch ${sample}_ploidy/${sample}-calls/contig_ploidy.tsv
    """
}


process GATK_GERMLINE_CNV_CALLER {
    tag "$sample"
    label 'process_high'
    label 'gatk'
    // Note: gCNV CASE mode can be CPU-intensive; PyTensor/gcnvkernel ≥ 0.9 required.
    // Container must use broadinstitute/gatk image with compatible PyTensor version.

    publishDir "${params.outdir}/germline/cnv_gatk/calls", mode: params.publish_dir_mode,
               saveAs: { filename -> filename.startsWith(sample) ? filename : null }

    input:
    tuple val(sample), path(counts), path(ploidy_dir)
    path pon_dir    // pre-built gCNV model (PoN)

    output:
    tuple val(sample), path("${sample}_cnvcalls"), emit: cnv_calls_dir

    script:
    """
    export THEANO_FLAGS="base_compiledir=\${PWD}/.theano"
    export HOME=\${PWD}
    mkdir -p ${sample}_cnvcalls

    gatk GermlineCNVCaller \\
        --run-mode CASE \\
        -I ${counts} \\
        --contig-ploidy-calls ${ploidy_dir}/${sample}-calls \\
        --model               ${pon_dir} \\
        --output              ${sample}_cnvcalls \\
        --output-prefix       ${sample}
    """

    stub:
    """
    mkdir -p ${sample}_cnvcalls/${sample}-calls
    touch ${sample}_cnvcalls/${sample}-calls/SAMPLE_000/log_q_tau_tchunk_0.tsv
    """
}


process GATK_POSTPROCESS_GERMLINE_CNV_CALLS {
    tag "$sample"
    label 'process_medium'
    label 'gatk'

    publishDir "${params.outdir}/germline/cnv_gatk/postprocessed", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(cnv_calls_dir), path(ploidy_dir)
    path pon_dir
    path fasta_dict

    output:
    tuple val(sample), path("${sample}_denoised.tsv"),      emit: denoised
    tuple val(sample), path("${sample}_intervals.vcf.gz"),  emit: intervals_vcf
    tuple val(sample), path("${sample}_intervals.vcf.gz.tbi"), emit: intervals_vcf_tbi
    tuple val(sample), path("${sample}_segments.vcf.gz"),   emit: segments_vcf
    tuple val(sample), path("${sample}_segments.vcf.gz.tbi"), emit: segments_vcf_tbi

    script:
    """
    export THEANO_FLAGS="base_compiledir=\${PWD}/.theano"
    export HOME=\${PWD}
    gatk PostprocessGermlineCNVCalls \\
        --model-shard-path         ${pon_dir} \\
        --calls-shard-path         ${cnv_calls_dir}/${sample}-calls \\
        --contig-ploidy-calls      ${ploidy_dir}/${sample}-calls \\
        --allosomal-contig chrX \\
        --allosomal-contig chrY \\
        --sample-index 0 \\
        --output-denoised-copy-ratios  ${sample}_denoised.tsv \\
        --output-genotyped-intervals   ${sample}_intervals.vcf.gz \\
        --output-genotyped-segments    ${sample}_segments.vcf.gz \\
        --sequence-dictionary          ${fasta_dict}

    tabix -f ${sample}_intervals.vcf.gz
    tabix -f ${sample}_segments.vcf.gz
    """

    stub:
    """
    touch ${sample}_denoised.tsv
    touch ${sample}_intervals.vcf.gz ${sample}_intervals.vcf.gz.tbi
    touch ${sample}_segments.vcf.gz  ${sample}_segments.vcf.gz.tbi
    """
}


process MANE_ANNOTATE_GERMLINE_CNV {
    tag "$sample"
    label 'process_low'
    label 'r_base'

    publishDir "${params.outdir}/germline/cnv_gatk/annotated", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(segments_vcf_gz)
    path mane_txt       // exons_mane.txt
    path gene_panel     // germline panel: first column = gene symbol
    val  exclude_genes  // comma-separated panel genes left out of germline CNV
    path lowmap_bed     // low-mappability regions (scripts/build_lowmap_panel.sh)

    output:
    tuple val(sample), path("${sample}_segments_annotated.tsv"), emit: annotated_tsv
    // exon_annotation.R may produce additional files; capture all with wildcard
    tuple val(sample), path("${sample}_segments*"),              emit: all_outputs

    script:
    """
    # Decompress VCF to a working directory
    mkdir -p vcf_in
    gunzip -c ${segments_vcf_gz} > vcf_in/${sample}_segments.vcf

    # Any R failure fails the task: a missing table must never read as "no CNV"
    Rscript ${projectDir}/bin/exon_annotation.R \\
        --vcf  vcf_in \\
        --MANE ${mane_txt} \\
        --genes ${gene_panel} \\
        --exclude "${exclude_genes}" \\
        --lowmap ${lowmap_bed}

    # exon_annotation.R writes <input>.txt next to each input VCF
    out=vcf_in/${sample}_segments.txt
    if [[ ! -f "\$out" ]]; then
        echo "ERROR: exon_annotation.R did not write \$out" >&2
        exit 1
    fi

    # The report readers need these columns; anything else is a broken table
    header="\$(head -1 "\$out")"
    for col in CHROM POS END CN QA gen LOWMAP LOWMAP_SEG; do
        if ! printf '%s\\n' "\$header" | tr '\\t' '\\n' | grep -qx "\$col"; then
            echo "ERROR: \$out lacks column \$col (header: \$header)" >&2
            exit 1
        fi
    done

    cp "\$out" ${sample}_segments_annotated.tsv
    """

    stub:
    """
    touch ${sample}_segments_annotated.tsv
    """
}

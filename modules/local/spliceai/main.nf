/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLICEAI
    ========
    Splicing predictions (SpliceAI 1.3.1, GENCODE v24 canonical annotation shipped
    with the tool) for the germline candidates that passed the panel + impact
    filter (tens of variants per sample, so it runs on CPU in minutes).

      INFO/SpliceAI  ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
                     distance 50 bp, masked scores (-M 1, recommended by Illumina
                     for variant interpretation)

    SNPSIFT_EXTRACT_GERMLINE turns it into SpliceAI_DS (max of the four deltas),
    used for PP3 / BP4 (ClinGen splicing recommendations, Walker et al. 2023).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SPLICEAI {
    tag "$sample"
    label 'process_medium'
    label 'spliceai'

    input:
    tuple val(sample), path(vcf)
    path fasta
    path fasta_fai

    output:
    tuple val(sample), path("${sample}.spliceai.vcf.gz"), emit: vcf
    path "versions.yml",                                   emit: versions

    script:
    """
    # SpliceAI (pysam) needs ##contig lines; plain VCF in, plain VCF out
    zcat ${vcf} > input.vcf
    if ! grep -q '^##contig=' input.vcf; then
        awk 'BEGIN{OFS=""} {print "##contig=<ID=", \$1, ",length=", \$2, ">"}' ${fasta_fai} > contigs.txt
        awk 'NR==FNR{c=c \$0 "\\n"; next} /^#CHROM/{printf "%s", c} {print}' contigs.txt input.vcf > tmp.vcf
        mv tmp.vcf input.vcf
    fi

    spliceai -I input.vcf -O ${sample}.spliceai.vcf -R ${fasta} -A grch38 -D 50 -M 1
    gzip -c ${sample}.spliceai.vcf > ${sample}.spliceai.vcf.gz
    rm -f input.vcf ${sample}.spliceai.vcf

    n_total=\$(zcat ${sample}.spliceai.vcf.gz | grep -vc '^#' || true)
    n_scored=\$(zcat ${sample}.spliceai.vcf.gz | grep -v '^#' | grep -c 'SpliceAI=' || true)
    echo "[SpliceAI] ${sample}: \${n_total} sites, \${n_scored} scored" >&2

    printf '"%s":\\n    spliceai: %s\\n' "${task.process}" "1.3.1" > versions.yml
    """

    stub:
    """
    touch ${sample}.spliceai.vcf.gz
    touch versions.yml
    """
}

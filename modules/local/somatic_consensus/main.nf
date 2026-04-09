/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC CONSENSUS MODULE (Mutect2 + Strelka2)
    ================================================
    Generates a consensus VCF from Mutect2 and Strelka2 somatic calls.
    Only variants called by BOTH callers (intersection) are retained.

    The consensus VCF preserves Mutect2 FORMAT fields (DP, AD, AF) and snpEff
    annotations via bcftools isec -w1, so all downstream filters (TMB, signatures,
    HRD) work without modification.

    Steps:
      1. Filter PASS from Strelka2 SNVs + indels, merge into one VCF
      2. Filter PASS from Mutect2
      3. Normalize both VCFs (bcftools norm — left-align indels)
      4. Intersect: keep variants called by both callers (bcftools isec -n=2 -w1)
      5. Emit consensus VCF + caller statistics TSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SOMATIC_CONSENSUS {
    tag "$pair"
    label 'process_low'
    label 'bcftools'

    publishDir "${params.outdir}/somatic/consensus", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(mutect2_vcf), path(strelka_snvs), path(strelka_indels)
    path fasta
    path fasta_fai

    output:
    tuple val(pair), path("${pair}.consensus.vcf.gz"),     emit: vcf
    tuple val(pair), path("${pair}.consensus.vcf.gz.tbi"), emit: tbi
    tuple val(pair), path("${pair}.caller_stats.tsv"),     emit: stats

    script:
    """
    # ── 1. Strelka2: filter PASS from SNVs + indels, merge ─────────────────
    bcftools view -f PASS ${strelka_snvs} -Oz -o strelka_snvs_pass.vcf.gz
    tabix -p vcf strelka_snvs_pass.vcf.gz

    bcftools view -f PASS ${strelka_indels} -Oz -o strelka_indels_pass.vcf.gz
    tabix -p vcf strelka_indels_pass.vcf.gz

    bcftools concat -a strelka_snvs_pass.vcf.gz strelka_indels_pass.vcf.gz \\
        | bcftools sort -Oz -o strelka_merged.vcf.gz
    tabix -p vcf strelka_merged.vcf.gz

    # ── 2. Mutect2: filter PASS + haplotype (soft-filter for phased variants) ─
    bcftools view -f PASS,haplotype ${mutect2_vcf} -Oz -o mutect2_pass.vcf.gz
    tabix -p vcf mutect2_pass.vcf.gz

    # ── 3. Normalize both (left-align indels, split multi-allelic) ─────────
    bcftools norm -m-both -f ${fasta} mutect2_pass.vcf.gz -Oz -o mutect2_norm.vcf.gz
    tabix -p vcf mutect2_norm.vcf.gz

    bcftools norm -m-both -f ${fasta} strelka_merged.vcf.gz -Oz -o strelka_norm.vcf.gz
    tabix -p vcf strelka_norm.vcf.gz

    # ── 4. Intersection: variants in BOTH callers ──────────────────────────
    # -n=2  : require exactly 2 files to agree
    # -w1   : write output from file 1 (Mutect2) → preserves FORMAT/INFO fields
    bcftools isec -n=2 -w1 mutect2_norm.vcf.gz strelka_norm.vcf.gz \\
        -Oz -o ${pair}.consensus.vcf.gz
    tabix -p vcf ${pair}.consensus.vcf.gz

    # ── 5. Caller statistics ───────────────────────────────────────────────
    N_MUTECT2=\$(bcftools view -H mutect2_norm.vcf.gz | wc -l)
    N_STRELKA=\$(bcftools view -H strelka_norm.vcf.gz | wc -l)
    N_CONSENSUS=\$(bcftools view -H ${pair}.consensus.vcf.gz | wc -l)

    printf "caller\\tvariants\\n"               >  ${pair}.caller_stats.tsv
    printf "mutect2\\t\${N_MUTECT2}\\n"         >> ${pair}.caller_stats.tsv
    printf "strelka2\\t\${N_STRELKA}\\n"        >> ${pair}.caller_stats.tsv
    printf "consensus\\t\${N_CONSENSUS}\\n"     >> ${pair}.caller_stats.tsv

    echo "[CONSENSUS] ${pair}: mutect2=\${N_MUTECT2} strelka2=\${N_STRELKA} consensus=\${N_CONSENSUS}"
    """

    stub:
    """
    printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n" \\
        | bgzip > ${pair}.consensus.vcf.gz
    tabix -p vcf ${pair}.consensus.vcf.gz
    printf "caller\\tvariants\\nmutect2\\t0\\nstrelka2\\t0\\nconsensus\\t0\\n" > ${pair}.caller_stats.tsv
    """
}

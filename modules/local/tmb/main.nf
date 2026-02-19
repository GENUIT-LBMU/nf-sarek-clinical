/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TUMOR MUTATIONAL BURDEN (TMB) MODULE
    =====================================
    Input:  snpEff-annotated Mutect2 VCFs from Sarek
    Steps:
      1. SnpSift filter: PASS + HIGH/MODERATE + DP >= min_dp + AF >= min_af + AD_alt > min_alt_reads
      2. Count variants per Mb (using panel BED for capture size)
    Output: tmb_per_sample.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TMB_FILTER {
    tag "$pair"
    label 'process_low'
    label 'snpsift'

    publishDir "${params.outdir}/somatic/tmb/filtered", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path snpsift_jar
    val min_dp
    val min_af
    val min_alt_reads

    output:
    tuple val(pair), path("${pair}.tmb_filter.vcf.gz"), emit: vcf

    stub:
    """
    touch ${pair}.tmb_filter.vcf.gz
    """

    script:
    """
    FILTER_EXPR="( FILTER = 'PASS' ) \\
        & (ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE') \\
        & ( GEN[1].DP >= ${min_dp} ) \\
        & ( GEN[1].AF >= ${min_af} ) \\
        & ( GEN[1].AD[1] > ${min_alt_reads} )"

    java -jar ${snpsift_jar} filter "\${FILTER_EXPR}" ${vcf} \\
        | bgzip -c > ${pair}.tmb_filter.vcf.gz
    """
}


process TMB_CALCULATE {
    tag "$pair"
    label 'process_single'
    label 'bcftools'

    publishDir "${params.outdir}/somatic/tmb", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path panel_bed  // optional; pass [] if not available

    output:
    tuple val(pair), path("${pair}.tmb.tsv"), emit: tsv

    stub:
    """
    printf "sample\\tn_variants\\tpanel_bases\\tpanel_mb\\tTMB_mut_per_mb\\n" > ${pair}.tmb.tsv
    printf "${pair}\\t0\\tNA\\tNA\\tNA\\n" >> ${pair}.tmb.tsv
    """

    script:
    def bed_available = panel_bed && panel_bed.name != 'NO_FILE'
    """
    # Index VCF
    tabix -f ${vcf} || true

    # Count variants (split multi-allelic first)
    n_vars=\$(bcftools norm -m -any -Ou ${vcf} | bcftools view -H | wc -l)

    ${bed_available ? """
    # Calculate panel size
    panel_bases=\$(sort -k1,1 -k2,2n ${panel_bed} | bedtools merge -i - | awk 'BEGIN{s=0}{s+=(\$3-\$2)}END{print s}')
    panel_mb=\$(awk -v b="\$panel_bases" 'BEGIN{printf "%.6f", b/1e6}')
    tmb=\$(awk -v n="\$n_vars" -v mb="\$panel_mb" 'BEGIN{if(mb>0) printf "%.4f", n/mb; else print "NA"}')
    """ : """
    panel_bases=NA
    panel_mb=NA
    tmb=NA
    """}

    printf "sample\\tn_variants\\tpanel_bases\\tpanel_mb\\tTMB_mut_per_mb\\n" > ${pair}.tmb.tsv
    printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "${pair}" "\$n_vars" "\$panel_bases" "\$panel_mb" "\$tmb" >> ${pair}.tmb.tsv

    echo "[TMB] ${pair}: n=\$n_vars  panel_mb=\$panel_mb  TMB=\$tmb"
    """
}

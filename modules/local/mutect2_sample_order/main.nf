/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MUTECT2 SAMPLE ORDER MODULE
    ===========================
    Mutect2 does not guarantee the sample order of its VCF (it comes out
    alphabetical). Every somatic step downstream reads GEN[0] as the normal and
    GEN[1] as the tumor, so the VCF is rewritten as normal,tumor using the
    ##normal_sample / ##tumor_sample header lines Mutect2 writes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MUTECT2_SAMPLE_ORDER {
    tag "$pair"
    label 'process_single'
    label 'bcftools'

    input:
    tuple val(pair), path(vcf)

    output:
    tuple val(pair), path("${pair}.mutect2.nt.vcf.gz"), emit: vcf

    stub:
    """
    touch ${pair}.mutect2.nt.vcf.gz
    """

    script:
    """
    NORMAL=\$(bcftools view -h ${vcf} | sed -n 's/^##normal_sample=//p' | head -1)
    TUMOR=\$(bcftools view -h ${vcf} | sed -n 's/^##tumor_sample=//p' | head -1)
    if [ -z "\${NORMAL}" ] || [ -z "\${TUMOR}" ]; then
        echo "ERROR: ${vcf} has no ##normal_sample/##tumor_sample header line (tumor/normal Mutect2 VCF expected)" >&2
        exit 1
    fi
    SAMPLES=\$(bcftools query -l ${vcf} | sort | tr '\\n' ' ')
    EXPECTED=\$(printf '%s\\n%s\\n' "\${NORMAL}" "\${TUMOR}" | sort | tr '\\n' ' ')
    if [ "\${SAMPLES}" != "\${EXPECTED}" ]; then
        echo "ERROR: ${vcf} samples [\${SAMPLES}] do not match normal=\${NORMAL} tumor=\${TUMOR}" >&2
        exit 1
    fi
    bcftools view -s "\${NORMAL},\${TUMOR}" ${vcf} -Oz -o ${pair}.mutect2.nt.vcf.gz
    echo "[MUTECT2_SAMPLE_ORDER] ${pair}: GEN[0]=\${NORMAL} GEN[1]=\${TUMOR}"
    """
}

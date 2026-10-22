/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GNOMAD_ANNOTATE
    ===============
    Adds population frequencies from the local gnomAD panel resource built by
    scripts/build_gnomad_panel.sh (gnomAD joint exomes+genomes, panel regions).

      INFO/GNOMAD_REGION  gene name(s) when the site lies inside the resource
                          regions. Inside + no GNOMAD_AF  => absent from gnomAD.
                          Outside                          => frequency unknown.
      INFO/GNOMAD_*       AF, AC, AN, FAF95 (grpmax), GRPMAX, AF_AMR, AN_AMR,
                          NHOMALT, NHOMALT_XY, FILTER — matched by CHROM/POS/REF/ALT.

    No network access: replaces the Ensembl / gnomAD web-API rescue.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GNOMAD_ANNOTATE {
    tag "$sample"
    label 'process_single'
    label 'bcftools'

    input:
    tuple val(sample), path(vcf)
    path gnomad_tsv
    path gnomad_tsv_tbi
    path gnomad_regions
    path gnomad_regions_tbi

    output:
    tuple val(sample), path("${sample}.gnomad.vcf.gz"), emit: vcf
    path "versions.yml",                                 emit: versions

    script:
    def cols = 'CHROM,POS,REF,ALT,INFO/GNOMAD_FILTER,INFO/GNOMAD_AC,INFO/GNOMAD_AN,INFO/GNOMAD_AF,' +
               'INFO/GNOMAD_NHOMALT,INFO/GNOMAD_FAF95,INFO/GNOMAD_GRPMAX,INFO/GNOMAD_AF_AMR,' +
               'INFO/GNOMAD_AN_AMR,INFO/GNOMAD_NHOMALT_XY'
    """
    cat > gnomad_header.txt <<'HDR'
##INFO=<ID=GNOMAD_REGION,Number=1,Type=String,Description="Site inside the regions of the local gnomAD resource (gene names); outside = frequency unknown">
##INFO=<ID=GNOMAD_FILTER,Number=A,Type=String,Description="gnomAD joint FILTER (PASS, EXOMES_FILTERED, GENOMES_FILTERED, BOTH_FILTERED)">
##INFO=<ID=GNOMAD_AC,Number=A,Type=Integer,Description="gnomAD joint alternate allele count">
##INFO=<ID=GNOMAD_AN,Number=1,Type=Integer,Description="gnomAD joint total number of alleles">
##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD joint (exomes+genomes) allele frequency">
##INFO=<ID=GNOMAD_NHOMALT,Number=A,Type=Integer,Description="gnomAD joint number of homozygous individuals">
##INFO=<ID=GNOMAD_FAF95,Number=A,Type=Float,Description="gnomAD joint maximum filtering allele frequency (Poisson 95% CI) across ancestry groups">
##INFO=<ID=GNOMAD_GRPMAX,Number=A,Type=String,Description="gnomAD joint ancestry group with the maximum allele frequency">
##INFO=<ID=GNOMAD_AF_AMR,Number=A,Type=Float,Description="gnomAD joint allele frequency, admixed American group">
##INFO=<ID=GNOMAD_AN_AMR,Number=1,Type=Integer,Description="gnomAD joint number of alleles, admixed American group">
##INFO=<ID=GNOMAD_NHOMALT_XY,Number=A,Type=Integer,Description="gnomAD joint homozygous/hemizygous count in XY samples">
HDR
    # Record resource version + checksum in the VCF header (read by SNPSIFT_EXTRACT_GERMLINE)
    ver=\$(tabix -H ${gnomad_tsv} | sed -n 's/^## gnomad_version=\\([^ ]*\\).*/\\1/p')
    echo "##gnomAD_resource=<Version=\${ver:-unknown},File=${gnomad_tsv},MD5=\$(md5sum ${gnomad_tsv} | cut -d' ' -f1)>" >> gnomad_header.txt

    bcftools annotate \\
        -h gnomad_header.txt \\
        -a ${gnomad_regions} \\
        -c CHROM,FROM,TO,INFO/GNOMAD_REGION \\
        ${vcf} -Ou \\
    | bcftools annotate \\
        -a ${gnomad_tsv} \\
        -c ${cols} \\
        -Oz -o ${sample}.gnomad.vcf.gz

    n_total=\$(bcftools view -H ${sample}.gnomad.vcf.gz | wc -l)
    n_region=\$(bcftools view -H -i 'INFO/GNOMAD_REGION!="."' ${sample}.gnomad.vcf.gz | wc -l)
    n_found=\$(bcftools view -H -i 'INFO/GNOMAD_AF!="."' ${sample}.gnomad.vcf.gz | wc -l)
    echo "[gnomAD local] ${sample}: \${n_total} sites, \${n_region} in resource regions, \${n_found} with gnomAD AF" >&2

    printf '"%s":\\n    bcftools: %s\\n    gnomad_resource: %s\\n' \\
        "${task.process}" "\$(bcftools --version | head -1 | sed 's/bcftools //')" "\$(basename ${gnomad_tsv})" \\
        > versions.yml
    """

    stub:
    """
    touch ${sample}.gnomad.vcf.gz
    touch versions.yml
    """
}

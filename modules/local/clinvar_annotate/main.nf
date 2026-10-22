/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLINVAR_ANNOTATE
    ================
    Adds ClinVar germline classifications from the local panel resource built by
    scripts/build_clinvar_panel.sh (NCBI ClinVar GRCh38 VCF, panel regions,
    left-normalised). Matched by CHROM/POS/REF/ALT, so indels are covered.

      INFO/CLINVAR_ID       ClinVar VariationID
      INFO/CLINVAR_SIG      aggregate germline classification (CLNSIG)
      INFO/CLINVAR_REVSTAT  review status (CLNREVSTAT) -> stars downstream
      INFO/CLINVAR_CONF     conflicting classifications detail (CLNSIGCONF)
      INFO/CLINVAR_DN       disease name(s) (CLNDN)

    When this resource is used it is the only ClinVar source: a variant absent
    from it is "not in ClinVar" (no fallback to the dbNSFP 4.1a ClinVar fields).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINVAR_ANNOTATE {
    tag "$sample"
    label 'process_single'
    label 'bcftools'

    input:
    tuple val(sample), path(vcf)
    path clinvar_tsv
    path clinvar_tsv_tbi

    output:
    tuple val(sample), path("${sample}.clinvar.vcf.gz"), emit: vcf
    path "versions.yml",                                  emit: versions

    script:
    def cols = 'CHROM,POS,REF,ALT,INFO/CLINVAR_ID,INFO/CLINVAR_SIG,INFO/CLINVAR_REVSTAT,' +
               'INFO/CLINVAR_CONF,INFO/CLINVAR_DN'
    """
    cat > clinvar_header.txt <<'HDR'
##INFO=<ID=CLINVAR_ID,Number=1,Type=String,Description="ClinVar VariationID (local ClinVar resource)">
##INFO=<ID=CLINVAR_SIG,Number=1,Type=String,Description="ClinVar aggregate germline classification (CLNSIG)">
##INFO=<ID=CLINVAR_REVSTAT,Number=1,Type=String,Description="ClinVar review status of the germline classification (CLNREVSTAT)">
##INFO=<ID=CLINVAR_CONF,Number=1,Type=String,Description="ClinVar conflicting germline classifications (CLNSIGCONF)">
##INFO=<ID=CLINVAR_DN,Number=1,Type=String,Description="ClinVar preferred disease name(s) (CLNDN)">
HDR
    # Record resource version + checksum in the VCF header (read by SNPSIFT_EXTRACT_GERMLINE)
    ver=\$(tabix -H ${clinvar_tsv} | sed -n 's/^## clinvar_version=\\([^ ]*\\).*/\\1/p')
    echo "##ClinVar_resource=<Version=\${ver:-unknown},File=${clinvar_tsv},MD5=\$(md5sum ${clinvar_tsv} | cut -d' ' -f1)>" >> clinvar_header.txt

    bcftools annotate \\
        -h clinvar_header.txt \\
        -a ${clinvar_tsv} \\
        -c ${cols} \\
        ${vcf} -Oz -o ${sample}.clinvar.vcf.gz

    n_total=\$(bcftools view -H ${sample}.clinvar.vcf.gz | wc -l)
    n_found=\$(bcftools view -H -i 'INFO/CLINVAR_ID!="."' ${sample}.clinvar.vcf.gz | wc -l)
    echo "[ClinVar local] ${sample}: \${n_total} sites, \${n_found} in ClinVar \${ver:-unknown}" >&2

    printf '"%s":\\n    bcftools: %s\\n    clinvar_resource: %s\\n' \\
        "${task.process}" "\$(bcftools --version | head -1 | sed 's/bcftools //')" "\$(basename ${clinvar_tsv})" \\
        > versions.yml
    """

    stub:
    """
    touch ${sample}.clinvar.vcf.gz
    touch versions.yml
    """
}

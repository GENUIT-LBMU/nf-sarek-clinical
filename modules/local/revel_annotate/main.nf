/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REVEL_ANNOTATE
    ==============
    Adds REVEL v1.3 missense scores from the local panel resource built by
    scripts/build_revel_panel.py (matched by CHROM/POS/REF/ALT).

      INFO/REVEL_SCORE  REVEL score of the MANE Select transcript, or the highest
                        score over the other transcripts when MANE has none
      INFO/REVEL_TX     MANE | other_tx
      INFO/REVEL_AA     score of every amino-acid change ("R>Q:0.123,R>W:0.456");
                        SNPSIFT_EXTRACT_GERMLINE picks the one matching the MANE
                        HGVS_p when REVEL_TX = other_tx

    Replaces the dbNSFP 4.1a REVEL value (first transcript in dbNSFP's list).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REVEL_ANNOTATE {
    tag "$sample"
    label 'process_single'
    label 'bcftools'

    input:
    tuple val(sample), path(vcf)
    path revel_tsv
    path revel_tsv_tbi

    output:
    tuple val(sample), path("${sample}.revel.vcf.gz"), emit: vcf
    path "versions.yml",                                emit: versions

    script:
    """
    cat > revel_header.txt <<'HDR'
##INFO=<ID=REVEL_SCORE,Number=1,Type=Float,Description="REVEL v1.3 score (MANE Select transcript, local REVEL resource)">
##INFO=<ID=REVEL_TX,Number=1,Type=String,Description="Transcript the REVEL score comes from: MANE (MANE Select) or other_tx (highest over other transcripts)">
##INFO=<ID=REVEL_AA,Number=1,Type=String,Description="REVEL v1.3 score per amino-acid change (aaref>aaalt:score, comma-separated)">
HDR
    ver=\$(tabix -H ${revel_tsv} | sed -n 's/^## revel_version=\\([^ ]*\\).*/\\1/p')
    echo "##REVEL_resource=<Version=\${ver:-unknown},File=${revel_tsv},MD5=\$(md5sum ${revel_tsv} | cut -d' ' -f1)>" >> revel_header.txt

    bcftools annotate \\
        -h revel_header.txt \\
        -a ${revel_tsv} \\
        -c CHROM,POS,REF,ALT,INFO/REVEL_SCORE,INFO/REVEL_TX,INFO/REVEL_AA \\
        ${vcf} -Oz -o ${sample}.revel.vcf.gz

    n_total=\$(bcftools view -H ${sample}.revel.vcf.gz | wc -l)
    n_found=\$(bcftools view -H -i 'INFO/REVEL_SCORE!="."' ${sample}.revel.vcf.gz | wc -l)
    echo "[REVEL local] ${sample}: \${n_total} sites, \${n_found} with REVEL v\${ver:-unknown}" >&2

    printf '"%s":\\n    bcftools: %s\\n    revel_resource: %s\\n' \\
        "${task.process}" "\$(bcftools --version | head -1 | sed 's/bcftools //')" "\$(basename ${revel_tsv})" \\
        > versions.yml
    """

    stub:
    """
    touch ${sample}.revel.vcf.gz
    touch versions.yml
    """
}

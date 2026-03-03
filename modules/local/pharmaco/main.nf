/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHARMACOGENOMICS MODULE (DPYD + UGT1A1)
    =========================================
    Filters germline HaplotypeCaller VCFs for known PGx-relevant rsIDs.
    Default panel: DPYD *2A, *13, HapB3, UGT1A1 *28 (fluoropyrimidine/irinotecan risk).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PHARMACO_FILTER {
    tag "$sample"
    label 'process_single'
    label 'snpsift'

    publishDir "${params.outdir}/germline/pharmaco", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf)    // dbNSFP-annotated germline VCF
    path  snpsift_jar
    path  variants_file             // rsID list (format: rsXXXXX per line, with optional columns)

    output:
    tuple val(sample), path("${sample}.pharmaco.tsv"), emit: tsv

    script:
    """
    # Build rsID list from file
    grep -v '^#' ${variants_file} | awk '{print \$1}' | grep '^rs' > _rsids.txt

    if [[ ! -s _rsids.txt ]]; then
        echo "WARNING: No rsIDs found in ${variants_file}" >&2
        printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tDP\\t" > ${sample}.pharmaco.tsv
        printf "ANN[*].GENE\\tANN[*].HGVS_C\\tANN[*].HGVS_P\\tANN[*].IMPACT\\t" >> ${sample}.pharmaco.tsv
        printf "GEN[0].GT\\tGEN[0].AD\\tGEN[0].DP\\tGEN[0].GQ\\t" >> ${sample}.pharmaco.tsv
        printf "dbNSFP_clinvar_id\\tdbNSFP_clinvar_clnsig\\tdbNSFP_clinvar_trait\\t" >> ${sample}.pharmaco.tsv
        printf "dbNSFP_SIFT_pred\\tdbNSFP_Polyphen2_HDIV_pred\\n" >> ${sample}.pharmaco.tsv
        exit 0
    fi

    # Build SnpSift filter expression from rsID file
    python3 - <<'PYEOF'
with open('_rsids.txt') as fh:
    ids = [line.strip() for line in fh if line.strip()]
exprs = ["(ID = '{}')".format(i) for i in ids]
with open('_filter_expr.txt', 'w') as out:
    out.write(' | '.join(exprs))
PYEOF

    java -jar ${snpsift_jar} filter "\$(cat _filter_expr.txt)" ${vcf} \\
        | gzip -c > ${sample}.pharmaco_filtered.vcf.gz

    # Extract fields
    java -jar ${snpsift_jar} extractFields -v -s ";" -e "." ${sample}.pharmaco_filtered.vcf.gz \\
        CHROM POS ID REF ALT QUAL FILTER DP \\
        "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].IMPACT" \\
        "GEN[0].GT" "GEN[0].AD" "GEN[0].DP" "GEN[0].GQ" \\
        "dbNSFP_clinvar_id" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_trait" \\
        "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_pred" \\
        > ${sample}.pharmaco.tsv

    echo "[PHARMACO] ${sample}: \$(wc -l < ${sample}.pharmaco.tsv) variants (incl. header)"
    """

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tdbNSFP_clinvar_clnsig\\n" > ${sample}.pharmaco.tsv
    """
}

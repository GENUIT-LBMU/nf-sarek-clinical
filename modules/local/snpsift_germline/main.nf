/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GERMLINE SNV POSTPROCESSING MODULE
    ===================================
    Steps:
      1. dbNSFP annotation (SnpSift dbnsfp)
      2. Gene panel + impact filter (SnpSift filter)
      3. Field extraction to TSV (SnpSift extractFields)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SNPSIFT_DBNSFP_GERMLINE {
    tag "$sample"
    label 'process_medium'
    label 'snpsift'

    publishDir "${params.outdir}/germline/snv/dbnsfp", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf)
    path dbnsfp
    path snpsift_jar

    output:
    tuple val(sample), path("${sample}.dbnsfp.vcf.gz"), emit: vcf

    stub:
    """
    touch ${sample}.dbnsfp.vcf.gz
    """

    script:
    def dbnsfp_fields = [
        "SIFT_score","SIFT_converted_rankscore","SIFT_pred",
        "SIFT4G_score","SIFT4G_converted_rankscore","SIFT4G_pred",
        "Polyphen2_HDIV_score","Polyphen2_HDIV_rankscore","Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score","Polyphen2_HVAR_rankscore","Polyphen2_HVAR_pred",
        "LRT_score","LRT_converted_rankscore","LRT_pred",
        "MutationTaster_score","MutationTaster_converted_rankscore","MutationTaster_pred",
        "MutationAssessor_score","MutationAssessor_rankscore","MutationAssessor_pred",
        "FATHMM_score","FATHMM_converted_rankscore","FATHMM_pred",
        "PROVEAN_score","PROVEAN_converted_rankscore","PROVEAN_pred",
        "REVEL_score","REVEL_rankscore",
        "MutPred_score","MutPred_rankscore","MutPred_protID","MutPred_AAchange","MutPred_Top5features",
        "PrimateAI_score","PrimateAI_rankscore","PrimateAI_pred",
        "ClinPred_score","ClinPred_rankscore","ClinPred_pred",
        "CADD_raw","CADD_raw_rankscore","CADD_phred",
        "DANN_score","DANN_rankscore",
        "fathmm-MKL_coding_score","fathmm-MKL_coding_rankscore","fathmm-MKL_coding_pred",
        "GERP++_NR","GERP++_RS","GERP++_RS_rankscore",
        "phyloP100way_vertebrate","phyloP100way_vertebrate_rankscore",
        "phastCons100way_vertebrate","phastCons100way_vertebrate_rankscore",
        "gnomAD_exomes_AC","gnomAD_exomes_AN","gnomAD_exomes_AF",
        "gnomAD_genomes_AC","gnomAD_genomes_AN","gnomAD_genomes_AF",
        "clinvar_id","clinvar_clnsig","clinvar_trait","clinvar_review","clinvar_hgvs",
        "clinvar_var_source","clinvar_MedGen_id","clinvar_OMIM_id","clinvar_Orphanet_id",
        "Interpro_domain"
    ].join(",")
    """
    java -jar ${snpsift_jar} dbnsfp \\
        -f ${dbnsfp_fields} \\
        -v \\
        -db ${dbnsfp} \\
        ${vcf} | bgzip -c > ${sample}.dbnsfp.vcf.gz
    """
}


process SNPSIFT_FILTER_GERMLINE {
    tag "$sample"
    label 'process_low'
    label 'snpsift'

    publishDir "${params.outdir}/germline/snv/filtered", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf)
    path gene_panel
    path snpsift_jar
    val min_dp

    output:
    tuple val(sample), path("${sample}.filtered.vcf.gz"), emit: vcf

    stub:
    """
    touch ${sample}.filtered.vcf.gz
    """

    script:
    """
    # Build dynamic filter expression from gene panel file
    genes=\$(grep -v '^#' ${gene_panel} | tr '\\n' '|' | sed 's/|$//')

    python3 - "\$genes" ${min_dp} <<'PYEOF'
import sys, subprocess, os

genes_raw = sys.argv[1]
min_dp    = sys.argv[2]

gene_list = [g.strip() for g in genes_raw.split("|") if g.strip()]
gene_expr  = " | ".join([f"(ANN[*].GENE = '{g}')" for g in gene_list])
filter_expr = f"({gene_expr}) & (ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE') & ( DP >= {min_dp} )"

with open("filter_expr.txt", "w") as fh:
    fh.write(filter_expr)
PYEOF

    java -jar ${snpsift_jar} filter "\$(cat filter_expr.txt)" ${vcf} \\
        | bgzip -c > ${sample}.filtered.vcf.gz
    """
}


process SNPSIFT_EXTRACT_GERMLINE {
    tag "$sample"
    label 'process_single'
    label 'snpsift'

    publishDir "${params.outdir}/germline/snv/extract", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf)
    path snpsift_jar

    output:
    tuple val(sample), path("${sample}.germline_snvs.tsv"), emit: tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\n" > ${sample}.germline_snvs.tsv
    """

    script:
    """
    java -jar ${snpsift_jar} extractFields -v -s ";" -e "." "${vcf}" \\
        CHROM POS ID REF ALT QUAL FILTER DP \\
        "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \\
        "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \\
        "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" \\
        "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" \\
        "GEN[0].GT" "GEN[0].AD" "GEN[0].DP" "GEN[0].GQ" \\
        "dbNSFP_SIFT_pred" "dbNSFP_SIFT4G_pred" \\
        "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HVAR_pred" \\
        "dbNSFP_LRT_pred" "dbNSFP_MutationTaster_pred" \\
        "dbNSFP_MutationAssessor_pred" "dbNSFP_FATHMM_pred" \\
        "dbNSFP_PROVEAN_pred" \\
        "dbNSFP_REVEL_score" "dbNSFP_REVEL_rankscore" \\
        "dbNSFP_MutPred_score" "dbNSFP_MutPred_rankscore" \\
        "dbNSFP_PrimateAI_pred" "dbNSFP_ClinPred_pred" \\
        "dbNSFP_CADD_phred" \\
        "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_genomes_AF" \\
        "dbNSFP_clinvar_id" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_trait" \\
        "dbNSFP_clinvar_review" "dbNSFP_clinvar_hgvs" \\
        "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" \\
        "dbNSFP_Interpro_domain" \\
        > "${sample}.germline_snvs.tsv"
    """
}

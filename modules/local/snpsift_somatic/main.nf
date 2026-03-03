/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC SNV POSTPROCESSING MODULE
    ===================================
    Steps:
      1. dbNSFP annotation (SnpSift dbnsfp)  – same fields as germline
      2. Somatic gene panel + PASS + impact + DP filter
      3. CIViC annotation (civicpy annotate-vcf)
      4. Field extraction to TSV
      5. CIViC field explode to TSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SNPSIFT_DBNSFP_SOMATIC {
    tag "$pair"
    label 'process_medium'
    label 'snpsift'

    publishDir "${params.outdir}/somatic/snv/dbnsfp", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path dbnsfp
    path dbnsfp_tbi
    path snpsift_jar

    output:
    tuple val(pair), path("${pair}.dbnsfp.vcf.gz"), emit: vcf

    stub:
    """
    touch ${pair}.dbnsfp.vcf.gz
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
        "REVEL_score","REVEL_rankscore",
        "MutPred_score","MutPred_rankscore",
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
        "clinvar_MedGen_id","clinvar_OMIM_id","clinvar_Orphanet_id",
        "Interpro_domain"
    ].join(",")
    """
    java -jar ${snpsift_jar} dbnsfp \\
        -f ${dbnsfp_fields} \\
        -v \\
        -db ${dbnsfp} \\
        ${vcf} | gzip -c > ${pair}.dbnsfp.vcf.gz
    """
}


process SNPSIFT_FILTER_SOMATIC {
    tag "$pair"
    label 'process_low'
    label 'snpsift'

    publishDir "${params.outdir}/somatic/snv/filtered", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path gene_panel
    path snpsift_jar
    val min_dp

    output:
    tuple val(pair), path("${pair}.filtered.vcf.gz"), emit: vcf

    stub:
    """
    touch ${pair}.filtered.vcf.gz
    """

    script:
    """
    genes=\$(grep -v '^#' ${gene_panel} | paste -sd '|')

    python3 - "\$genes" ${min_dp} <<'PYEOF'
import sys

genes_raw = sys.argv[1]
min_dp    = sys.argv[2]

gene_list  = [g.strip() for g in genes_raw.split("|") if g.strip()]
gene_expr  = " | ".join([f"(ANN[*].GENE = '{g}')" for g in gene_list])
filter_expr = (
    f"( FILTER = 'PASS' | FILTER = 'haplotype' ) "
    f"& ({gene_expr}) "
    f"& (ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE') "
    f"& ( GEN[1].DP >= {min_dp} )"
)

with open("filter_expr.txt", "w") as fh:
    fh.write(filter_expr)
PYEOF

    java -jar ${snpsift_jar} filter "\$(cat filter_expr.txt)" ${vcf} \\
        | gzip -c > ${pair}.filtered.vcf.gz
    """
}


process CIVIC_ANNOTATE {
    tag "$pair"
    label 'process_low'
    label 'civicpy'

    publishDir "${params.outdir}/somatic/snv/civic", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    val genome_build

    output:
    tuple val(pair), path("${pair}.civic.vcf.gz"), emit: vcf

    stub:
    """
    touch ${pair}.civic.vcf.gz
    """

    script:
    // Reference: GRCh38 → 'GRCh38', GRCh37 → 'GRCh37'
    // Gracefully fall back to passing through the input VCF if CIViC annotation fails
    // (CIVIC_EXTRACT already handles VCFs without CIVIC INFO fields → empty TSV)
    """
    # Ensure CIViC cache is fresh
    CIVIC_CACHE="\${PWD}/civicpy_cache.pkl"
    civicpy update --hard --cache-save-path "\${CIVIC_CACHE}" || true

    if ! CIVICPY_CACHE_FILE="\${CIVIC_CACHE}" \\
         civicpy annotate-vcf \\
             --input-vcf  ${vcf} \\
             --output-vcf ${pair}.civic.vcf.gz \\
             --reference  ${genome_build} \\
             -i submitted 2>&1; then
        echo "WARNING: CIViC annotation failed — passing through unannotated VCF"
        cp ${vcf} ${pair}.civic.vcf.gz
    fi
    """
}


process SNPSIFT_EXTRACT_SOMATIC {
    tag "$pair"
    label 'process_single'
    label 'snpsift'

    publishDir "${params.outdir}/somatic/snv/extract", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path snpsift_jar

    output:
    tuple val(pair), path("${pair}.somatic_snvs.tsv"), emit: tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\n" > ${pair}.somatic_snvs.tsv
    """

    script:
    """
    java -jar ${snpsift_jar} extractFields -v -s ";" -e "." "${vcf}" \\
        CHROM POS ID REF ALT QUAL FILTER DP \\
        "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \\
        "ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" \\
        "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \\
        "ANN[*].AA_POS" "ANN[*].AA_LEN" \\
        "GEN[1].GT" "GEN[1].AD" "GEN[1].AF" "GEN[1].DP" \\
        "dbNSFP_SIFT_pred" "dbNSFP_SIFT4G_pred" \\
        "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HVAR_pred" \\
        "dbNSFP_LRT_pred" "dbNSFP_MutationTaster_pred" \\
        "dbNSFP_MutationAssessor_pred" "dbNSFP_FATHMM_pred" \\
        "dbNSFP_REVEL_score" "dbNSFP_REVEL_rankscore" \\
        "dbNSFP_PrimateAI_pred" "dbNSFP_ClinPred_pred" \\
        "dbNSFP_CADD_phred" \\
        "dbNSFP_gnomAD_exomes_AF" "dbNSFP_gnomAD_genomes_AF" \\
        "dbNSFP_clinvar_id" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_trait" \\
        "dbNSFP_clinvar_review" "dbNSFP_clinvar_hgvs" \\
        "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{print; next} {for(i=1;i<=NF;i++){n=split(\$i,a,";"); \$i=a[1]}; print}' \\
    > "${pair}.somatic_snvs.tsv"
    """
}


process CIVIC_EXTRACT {
    tag "$pair"
    label 'process_single'
    label 'bcftools'

    publishDir "${params.outdir}/somatic/snv/civic_extract", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)

    output:
    tuple val(pair), path("${pair}.civic.exploded.tsv"), emit: tsv

    stub:
    """
    printf "CHROM\\tPOS\\tREF\\tALT\\tCIViC Variant Name\\tCIViC Entity Disease\\n" > ${pair}.civic.exploded.tsv
    """

    script:
    // 39-field CIViC explode (pipe-delimited INFO/CIVIC → TSV columns)
    """
    HEADER="CHROM\\tPOS\\tREF\\tALT\\tAllele\\tConsequence\\tSYMBOL\\tEntrez Gene ID\\t"
    HEADER+=\$(printf "Feature_type\\tFeature\\tHGVSc\\tHGVSp\\t")
    HEADER+=\$(printf "CIViC Variant Name\\tCIViC Variant ID\\tCIViC Variant Aliases\\tCIViC Variant URL\\t")
    HEADER+=\$(printf "CIViC Molecular Profile Name\\tCIViC Molecular Profile ID\\tCIViC Molecular Profile Aliases\\tCIViC Molecular Profile URL\\t")
    HEADER+=\$(printf "CIViC HGVS\\tAllele Registry ID\\tClinVar IDs\\tCIViC Molecular Profile Score\\t")
    HEADER+=\$(printf "CIViC Entity Type\\tCIViC Entity ID\\tCIViC Entity URL\\tCIViC Entity Source\\t")
    HEADER+=\$(printf "CIViC Entity Variant Origin\\tCIViC Entity Status\\tCIViC Entity Significance\\t")
    HEADER+=\$(printf "CIViC Entity Direction\\tCIViC Entity Disease\\tCIViC Entity Therapies\\t")
    HEADER+=\$(printf "CIViC Entity Therapy Interaction Type\\tCIViC Evidence Phenotypes\\t")
    HEADER+=\$(printf "CIViC Evidence Level\\tCIViC Evidence Rating\\t")
    HEADER+=\$(printf "CIViC Assertion ACMG Codes\\tCIViC Assertion AMP Category\\t")
    HEADER+=\$(printf "CIViC Assertion NCCN Guideline\\tCIViC Assertion Regulatory Approval\\tCIViC Assertion FDA Companion Test")

    # Check if CIViC annotation was added
    if ! bcftools view -h ${vcf} | grep -q '^##INFO=<ID=CIVIC'; then
        echo "WARNING: No CIVIC INFO field found in ${pair} — producing empty TSV"
        printf "%s\\n" "\$HEADER" > ${pair}.civic.exploded.tsv
        exit 0
    fi

    printf "%s\\n" "\$HEADER" > ${pair}.civic.exploded.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CIVIC\\n' \\
        -i 'INFO/CIVIC!=""' ${vcf} \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"} {
        chrom=\$1; pos=\$2; ref=\$3; alt=\$4; civ=\$5;
        nE = split(civ, E, /,/);
        for (i=1; i<=nE; i++) {
            nF = split(E[i], F, /\\|/);
            for (k=1; k<=39; k++) if (!(k in F)) F[k]=".";
            print chrom, pos, ref, alt, \\
                  F[1],F[2],F[3],F[4],F[5],F[6],F[7],F[8], \\
                  F[9],F[10],F[11],F[12],F[13],F[14],F[15],F[16], \\
                  F[17],F[18],F[19],F[20],F[21],F[22],F[23],F[24],F[25],F[26], \\
                  F[27],F[28],F[29],F[30],F[31],F[32],F[33],F[34],F[35],F[36],F[37],F[38],F[39];
        }
    }' >> ${pair}.civic.exploded.tsv
    """
}

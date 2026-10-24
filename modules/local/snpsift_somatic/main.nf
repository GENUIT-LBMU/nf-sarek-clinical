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
    path civic_cache    // dated civicpy cache (scripts/build_civic_cache.sh)

    output:
    tuple val(pair), path("${pair}.civic.vcf.gz"), emit: vcf

    stub:
    """
    touch ${pair}.civic.vcf.gz
    """

    script:
    // Reference: GRCh38 → 'GRCh38', GRCh37 → 'GRCh37'
    // The staged cache is used as is: no staleness check (civicpy would
    // re-download a cache older than 7 days) and no remote URL, so the same
    // cache always gives the same annotation. Any failure fails the task:
    // an unannotated VCF must never read as "no CIViC evidence".
    """
    export CIVICPY_CACHE_FILE="\${PWD}/${civic_cache}"
    export CIVICPY_CACHE_TIMEOUT_DAYS=36500
    export CIVICPY_REMOTE_CACHE_URL=file:///nonexistent

    civicpy annotate-vcf \\
        --input-vcf  ${vcf} \\
        --output-vcf ${pair}.civic.vcf.gz \\
        --reference  ${genome_build} \\
        -i accepted -i submitted

    if ! bcftools view -h ${pair}.civic.vcf.gz | grep -q '^##INFO=<ID=CIVIC,'; then
        echo "ERROR: civicpy wrote no CIVIC INFO header in ${pair}.civic.vcf.gz" >&2
        exit 1
    fi
    """
}


process SNPSIFT_EXTRACT_SOMATIC {
    tag "$pair"
    label 'process_single'
    label 'snpsift'

    // No publishDir — ANNOTATE_CONSENSUS publishes the final TSV with Callers column

    input:
    tuple val(pair), path(vcf)
    path snpsift_jar

    output:
    tuple val(pair), path("${pair}.somatic_snvs_raw.tsv"), emit: tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\n" > ${pair}.somatic_snvs_raw.tsv
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
    > "${pair}.somatic_snvs_raw.tsv"
    """
}


process ANNOTATE_CONSENSUS {
    tag "$pair"
    label 'process_single'
    label 'bcftools'

    publishDir "${params.outdir}/somatic/snv/extract", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(tsv), path(consensus_vcf)

    output:
    tuple val(pair), path("${pair}.somatic_snvs.tsv"), emit: tsv

    stub:
    """
    cp ${tsv} ${pair}.somatic_snvs.tsv
    """

    script:
    """
    # Extract consensus positions from VCF (CHROM:POS:REF:ALT)
    if [ -s "${consensus_vcf}" ] && [ "${consensus_vcf}" != "NO_FILE" ]; then
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${consensus_vcf} \\
            | sort -u > consensus_positions.txt
    else
        touch consensus_positions.txt
    fi

    # Annotate TSV with Callers column using awk
    # Header: append "Callers"; data rows: lookup CHROM+POS+REF+ALT in consensus set
    awk -F'\\t' -v OFS='\\t' '
    BEGIN {
        while ((getline line < "consensus_positions.txt") > 0) {
            n = split(line, a, "\\t")
            if (n >= 4) cons[a[1] "\\t" a[2] "\\t" a[3] "\\t" a[4]] = 1
        }
        close("consensus_positions.txt")
        has_cons = length(cons)
    }
    NR == 1 { print \$0, "Callers"; next }
    {
        key = \$1 "\\t" \$2 "\\t" \$4 "\\t" \$5
        if (!has_cons)
            print \$0, "Mutect2"
        else if (key in cons)
            print \$0, "Mutect2+Strelka2"
        else
            print \$0, "Mutect2"
    }
    ' "${tsv}" > "${pair}.somatic_snvs.tsv"
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
    // CIViC explode (pipe-delimited INFO/CIVIC → one TSV row per evidence item).
    // Columns come from the CIVIC header Format string, so they always match
    // the civicpy version that wrote the VCF (5.4.0 added "CIViC Assertion
    // ClinGen Codes"); a record with another field count fails the task.
    """
    # CIVIC_ANNOTATE always writes the header; without it the input is not annotated
    fmt="\$(bcftools view -h ${vcf} | sed -n 's/^##INFO=<ID=CIVIC,.*Format: \\(.*\\)">\$/\\1/p')"
    if [[ -z "\$fmt" ]]; then
        echo "ERROR: no CIVIC INFO field (with Format) in ${vcf}: not a CIViC-annotated VCF" >&2
        exit 1
    fi
    nfields=\$(printf '%s\\n' "\$fmt" | awk -F'|' '{print NF}')
    { printf 'CHROM\\tPOS\\tREF\\tALT\\t'; printf '%s\\n' "\$fmt" | tr '|' '\\t'; } > ${pair}.civic.exploded.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CIVIC\\n' \\
        -i 'INFO/CIVIC!=""' ${vcf} \\
    | awk -F'\\t' -v n="\$nfields" 'BEGIN{OFS="\\t"} {
        nE = split(\$5, E, /,/);
        for (i=1; i<=nE; i++) {
            nF = split(E[i], F, /\\|/);
            if (nF != n) {
                printf "ERROR: %s:%s has %d CIViC fields, header has %d\\n", \$1, \$2, nF, n > "/dev/stderr";
                exit 1;
            }
            line = \$1 OFS \$2 OFS \$3 OFS \$4;
            for (k=1; k<=n; k++) line = line OFS F[k];
            print line;
        }
    }' >> ${pair}.civic.exploded.tsv
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SNPEFF_MANE
    ===========
    Re-annotates the germline VCF with snpEff using a MANE Select database
    (default GRCh38.mane.1.2.refseq: one RefSeq NM transcript per gene).

    Sarek runs snpEff with `-canon`, which keeps the transcript with the longest
    CDS. That is often not the clinical MANE transcript, so HGVS numbering
    differs from ClinVar / literature (e.g. BRCA1 p.Ser1634Gly instead of
    p.Ser1613Gly). The sarek ANN / LOF / NMD fields are removed first and
    replaced by the MANE annotation.

    Database: snpEff `download` of the database, or
      https://snpeff-public.s3.amazonaws.com/databases/v5_2/snpEff_v5_2_<db>.zip
    unzipped so that <snpeff_mane_data>/<snpeff_mane_db>/snpEffectPredictor.bin exists.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SNPEFF_MANE {
    tag "$sample"
    label 'process_low'
    label 'snpeff_mane'

    input:
    tuple val(sample), path(vcf)
    path snpeff_data            // directory containing <snpeff_db>/
    val  snpeff_db              // e.g. GRCh38.mane.1.2.refseq

    output:
    tuple val(sample), path("${sample}.mane.vcf.gz"), emit: vcf
    path "versions.yml",                               emit: versions

    script:
    def mem = task.memory ? Math.max(2, (task.memory.toGiga() as int) - 2) : 4
    """
    # Drop the sarek snpEff annotation (ANN / LOF / NMD) from header and INFO
    zcat ${vcf} | awk 'BEGIN { FS = OFS = "\\t" }
        /^##INFO=<ID=(ANN|LOF|NMD),/ || /^##SnpEff/ { next }
        /^#/ { print; next }
        {
            n = split(\$8, kv, ";"); out = ""
            for (i = 1; i <= n; i++)
                if (kv[i] !~ /^(ANN|LOF|NMD)=/) out = out (out == "" ? "" : ";") kv[i]
            \$8 = (out == "" ? "." : out)
            print
        }' > stripped.vcf

    snpEff -Xmx${mem}g ann \\
        -nodownload -noStats -v \\
        -dataDir \$PWD/${snpeff_data} \\
        ${snpeff_db} stripped.vcf 2> snpeff.log \\
        | gzip -c > ${sample}.mane.vcf.gz
    rm -f stripped.vcf

    n_total=\$(zcat ${sample}.mane.vcf.gz | grep -vc '^#' || true)
    n_ann=\$(zcat ${sample}.mane.vcf.gz | grep -v '^#' | grep -c 'ANN=' || true)
    echo "[snpEff MANE] ${sample}: \${n_total} sites, \${n_ann} annotated with ${snpeff_db}" >&2

    printf '"%s":\\n    snpeff: %s\\n    snpeff_db: %s\\n' \\
        "${task.process}" "\$(snpEff -version 2>&1 | head -1 | cut -f2)" "${snpeff_db}" \\
        > versions.yml
    """

    stub:
    """
    touch ${sample}.mane.vcf.gz
    touch versions.yml
    """
}

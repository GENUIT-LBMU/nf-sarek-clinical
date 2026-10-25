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

    The database also holds the MANE Plus Clinical transcripts (e.g. BRAF
    NM_001374258.1, MUTYH NM_001128425.2) and snpEff may list them before the
    MANE Select one. Their ANN entries are moved after the others, so the first
    matching ANN taken downstream is MANE Select, and their HGVS (with HGVS.p,
    or HIGH impact such as splice sites) is copied
    to INFO/MANE_PLUS_CLINICAL (Gene|Transcript|HGVS.c|HGVS.p, '&'-separated)
    for display next to the MANE Select HGVS.

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
    path plus_clinical          // MANE Plus Clinical transcripts of that release (assets/mane)

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
        ${snpeff_db} stripped.vcf 2> snpeff.log > annotated.vcf
    rm -f stripped.vcf

    # MANE Select first: MANE Plus Clinical ANN entries go last (order kept
    # otherwise); their HGVS (with HGVS.p, or HIGH impact) goes to INFO/MANE_PLUS_CLINICAL
    awk 'BEGIN { FS = OFS = "\\t" }
        NR == FNR { if (\$0 !~ /^#/ && \$1 != "") plus[\$1] = 1; next }
        /^##INFO=<ID=ANN,/ {
            print
            print "##INFO=<ID=MANE_PLUS_CLINICAL,Number=1,Type=String,Description=\\"MANE Plus Clinical HGVS from ${snpeff_db}: Gene|Transcript|HGVS.c|HGVS.p, &-separated. These ANN entries are listed after the MANE Select ones\\">"
            next
        }
        /^#/ { print; next }
        {
            n = split(\$8, kv, ";"); extra = ""; moved = 0
            for (i = 1; i <= n; i++) {
                if (kv[i] !~ /^ANN=/) continue
                m = split(substr(kv[i], 5), a, ","); sel = ""; pl = ""; split("", seen)
                for (j = 1; j <= m; j++) {
                    split(a[j], f, "|")
                    if (f[7] in plus) {
                        pl = pl (pl == "" ? "" : ",") a[j]
                        if (sel == "") moved = 1
                        key = f[4] "|" f[7] "|" f[10] "|" f[11]
                        if (f[10] ~ /^c\\./ && (f[11] != "" || f[3] == "HIGH") && !(key in seen)) {
                            seen[key] = 1
                            extra = extra (extra == "" ? "" : "&") key
                        }
                    } else {
                        sel = sel (sel == "" ? "" : ",") a[j]
                    }
                }
                if (sel == "") moved = 0
                kv[i] = "ANN=" sel (sel != "" && pl != "" ? "," : "") pl
            }
            out = kv[1]
            for (i = 2; i <= n; i++) out = out ";" kv[i]
            if (extra != "") { out = out ";MANE_PLUS_CLINICAL=" extra; nplus++ }
            nmoved += moved
            \$8 = out
            print
        }
        END {
            printf "[snpEff MANE] %d sites with a MANE Plus Clinical HGVS, %d with MANE Select moved first\\n", nplus, nmoved > "/dev/stderr"
        }' ${plus_clinical} annotated.vcf | gzip -c > ${sample}.mane.vcf.gz
    rm -f annotated.vcf

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

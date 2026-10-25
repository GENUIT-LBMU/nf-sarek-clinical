/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC SNV POSTPROCESSING MODULE
    ===================================
    Steps (input: Mutect2 VCF, MANE re-annotated and with local gnomAD
    frequencies when configured; dbNSFP is germline-only):
      1. Somatic gene panel + PASS + impact + DP filter
      2. CIViC annotation (civicpy annotate-vcf)
      3. Field extraction to TSV (+ gnomAD_AF / gnomAD_FAF95 / gnomAD_AF_source)
      4. CIViC field explode to TSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    path gene_panel
    path snpsift_jar
    val max_pop_af              // population AF for AMP Tier IV (filters.somatic_snv.max_pop_af)

    output:
    tuple val(pair), path("${pair}.somatic_snvs_raw.tsv"), emit: tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\n" > ${pair}.somatic_snvs_raw.tsv
    """

    script:
    """
    # Optional INFO fields (SnpSift aborts on fields missing from the header):
    # MANE Plus Clinical HGVS (SNPEFF_MANE) and local gnomAD (GNOMAD_ANNOTATE)
    zcat "${vcf}" | grep '^##' > vcf_header.txt || true
    MANE_FIELDS=""
    if grep -q '^##INFO=<ID=MANE_PLUS_CLINICAL,' vcf_header.txt; then
        MANE_FIELDS="MANE_PLUS_CLINICAL"
    fi
    GNOMAD_FIELDS=""
    GNOMAD_VERSION=""
    if grep -q '^##INFO=<ID=GNOMAD_AF,' vcf_header.txt; then
        GNOMAD_FIELDS="GNOMAD_REGION GNOMAD_FILTER GNOMAD_AF GNOMAD_AC GNOMAD_FAF95 GNOMAD_NHOMALT"
        GNOMAD_VERSION=\$(sed -n 's/^##gnomAD_resource=<Version=\\([^,]*\\),.*/\\1/p' vcf_header.txt)
    fi
    rm -f vcf_header.txt
    export GNOMAD_VERSION

    java -jar ${snpsift_jar} extractFields -v -s ";" -e "." "${vcf}" \\
        CHROM POS ID REF ALT QUAL FILTER DP \\
        "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \\
        "ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" \\
        "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \\
        "ANN[*].AA_POS" "ANN[*].AA_LEN" \\
        "GEN[1].GT" "GEN[1].AD" "GEN[1].AF" "GEN[1].DP" \\
        \${MANE_FIELDS} \${GNOMAD_FIELDS} \\
    > extract_all.tsv

    # One ANN per row: the first one that meets the filter (panel gene AND
    # HIGH/MODERATE), so gene, impact and HGVS all come from the same
    # transcript. Other multi-valued fields keep their first value, except
    # GNOMAD_* (one value per ALT allele, all kept for the max below).
    awk -F'\\t' 'BEGIN{OFS="\\t"}
        FNR==NR { sub(/\\r\$/, ""); gsub(/ /, ""); if (\$0 != "" && \$0 !~ /^#/) panel[\$0] = 1; next }
        FNR==1 {
            for (i = 1; i <= NF; i++) {
                if (index(\$i, "ANN[*].") == 1) ann[i] = 1
                if (index(\$i, "GNOMAD_") == 1) whole[i] = 1
                if (\$i == "ANN[*].GENE")   gc = i
                if (\$i == "ANN[*].IMPACT") ic = i
            }
            print; next
        }
        {
            ng = split(\$gc, g, ";"); split(\$ic, im, ";")
            k = 0
            for (j = 1; j <= ng; j++)
                if ((g[j] in panel) && (im[j] == "HIGH" || im[j] == "MODERATE")) { k = j; break }
            if (k == 0) {
                print "WARNING: " \$1 ":" \$2 " has no ANN with a panel gene and HIGH/MODERATE impact; keeping the first" > "/dev/stderr"
                k = 1
            }
            for (i = 1; i <= NF; i++) {
                if (i in whole) continue
                n = split(\$i, a, ";")
                \$i = ((i in ann) && k <= n) ? a[k] : a[1]
            }
            print
        }' ${gene_panel} extract_all.tsv > picked.tsv
    rm -f extract_all.tsv

    # Population frequency (AMP/ASCO/CAP 2017 Tier IV: MAF >= 1% in the general
    # population). Same sources as the germline gnomAD_AF / gnomAD_AF_source:
    #   gnomAD_v<ver>          local gnomAD resource (max over ALT alleles)
    #   gnomAD_v<ver>_qc_fail  fails QC in gnomAD exomes and genomes (AF still used)
    #   gnomAD_v<ver>_absent   inside the resource regions, not in gnomAD -> 0
    #   unknown                outside the regions, or no --gnomad_tsv
    # gnomAD_POP_AF = FAF95 when available, else AF; POP_COMMON = 1 when
    # gnomAD_POP_AF >= max_pop_af (${max_pop_af}) -> Tier IV in the reports.
    python3 - <<'PYEOF'
import os, sys

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def nums(x):
    return [v for v in (num(t) for t in (x or '').replace(';', ',').split(',')) if v is not None]

tag = 'gnomAD_v' + (os.environ.get('GNOMAD_VERSION') or 'unknown')
max_pop_af = float('${max_pop_af}')
with open('picked.tsv') as fh:
    lines = [l.rstrip('\\n') for l in fh if l.strip()]
hdr = lines[0].split('\\t')
gi = {k: (hdr.index(k) if k in hdr else -1) for k in
      ('GNOMAD_REGION', 'GNOMAD_FILTER', 'GNOMAD_AF', 'GNOMAD_FAF95', 'GNOMAD_NHOMALT')}
def col(row, k):
    i = gi[k]
    return row[i].strip() if 0 <= i < len(row) else '.'
counts = {}
with open('${pair}.somatic_snvs_raw.tsv', 'w') as out:
    out.write('\\t'.join(hdr + ['gnomAD_AF', 'gnomAD_AF_source', 'gnomAD_FAF95', 'gnomAD_NHOMALT',
                               'gnomAD_POP_AF', 'POP_COMMON']) + '\\n')
    for line in lines[1:]:
        row = line.split('\\t') + ['.'] * (len(hdr) - len(line.split('\\t')))
        afs   = col(row, 'GNOMAD_AF').replace(';', ',').split(',')
        filts = col(row, 'GNOMAD_FILTER').replace(';', ',').split(',')
        vals  = [(num(a), filts[k] if k < len(filts) else '.') for k, a in enumerate(afs)]
        vals  = [(a, f) for a, f in vals if a is not None]
        fafs, homs = nums(col(row, 'GNOMAD_FAF95')), nums(col(row, 'GNOMAD_NHOMALT'))
        if vals:
            good = [a for a, f in vals if f != 'BOTH_FILTERED']
            af, src = (max(good), tag) if good else (max(a for a, _f in vals), tag + '_qc_fail')
            faf = max(fafs) if fafs else None
            hom = '%d' % max(homs) if homs else '0'
        elif gi['GNOMAD_REGION'] >= 0 and col(row, 'GNOMAD_REGION') not in ('', '.'):
            af, src, faf, hom = 0.0, tag + '_absent', 0.0, '0'
        else:
            af, src, faf, hom = None, 'unknown', None, '.'
        pop = faf if faf is not None else af
        counts[src] = counts.get(src, 0) + 1
        out.write('\\t'.join(row[:len(hdr)] + [
            repr(af) if af is not None else '.', src,
            repr(faf) if faf is not None else '.', hom,
            repr(pop) if pop is not None else '.',
            '1' if pop is not None and pop >= max_pop_af else '0']) + '\\n')
        if pop is not None and pop >= max_pop_af:
            sys.stderr.write('[gnomAD] ${pair}: %s:%s population AF %.4g >= %s -> Tier IV\\n'
                             % (row[0], row[1], pop, max_pop_af))
sys.stderr.write('[gnomAD] ${pair}: ' + ', '.join(k + '=' + str(v) for k, v in sorted(counts.items())) + '\\n')
PYEOF
    rm -f picked.tsv
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

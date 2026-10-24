/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GERMLINE SNV POSTPROCESSING MODULE
    ===================================
    Steps:
      1. dbNSFP annotation (SnpSift dbnsfp)
      2. Gene panel + impact filter (SnpSift filter): HIGH / MODERATE, plus splice-region,
         synonymous and intronic variants when noncoding_splice is on (Step 2 keeps
         those only with SpliceAI >= 0.2 or ClinVar P/LP)
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
    path dbnsfp_tbi
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
        ${vcf} | gzip -c > ${sample}.dbnsfp.vcf.gz
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
    val noncoding_splice    // true: also LOW / MODIFIER splice-region, synonymous and
                            // intronic variants (kept in Step 2 only if SpliceAI >= 0.2)

    output:
    tuple val(sample), path("${sample}.filtered.vcf.gz"), emit: vcf

    stub:
    """
    touch ${sample}.filtered.vcf.gz
    """

    script:
    """
    # Build dynamic filter expression from gene panel file
    genes=\$(grep -v '^#' ${gene_panel} | paste -sd '|')

    python3 - "\$genes" ${min_dp} ${noncoding_splice} <<'PYEOF'
import sys, subprocess, os

genes_raw = sys.argv[1]
min_dp    = sys.argv[2]
noncoding = sys.argv[3] == 'true'

gene_list = [g.strip() for g in genes_raw.split("|") if g.strip()]
gene_expr  = " | ".join([f"(ANN[*].GENE = '{g}')" for g in gene_list])
impact_expr = "ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE'"
if noncoding:
    # Non-coding splice candidates (GS15): splice region (+3..+8 / -3..-20 and the
    # donor 5th base / polypyrimidine tract), synonymous and intronic; UTR,
    # upstream and downstream stay out. Step 2 keeps them only with SpliceAI >= 0.2.
    effects = ('splice_region_variant', 'splice_donor_5th_base_variant',
               'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant',
               'synonymous_variant', 'intron_variant')
    impact_expr += " | (ANN[*].EFFECT =~ '" + "|".join(effects) + "')"
filter_expr = f"({gene_expr}) & ({impact_expr}) & ( DP >= {min_dp} )"

with open("filter_expr.txt", "w") as fh:
    fh.write(filter_expr)
PYEOF

    java -jar ${snpsift_jar} filter "\$(cat filter_expr.txt)" ${vcf} > filtered.vcf

    # Panel-gene annotations first (stable order): the extraction keeps the first ANN
    # entry, and a non-coding panel variant can be listed after a neighbouring gene
    # (e.g. downstream_gene_variant FPGS before intron_variant ENG)
    python3 - "\$genes" filtered.vcf <<'PYEOF' | gzip -c > ${sample}.filtered.vcf.gz
import sys
panel = set(g.strip() for g in sys.argv[1].split("|") if g.strip())
for line in open(sys.argv[2]):
    if line.startswith("#"):
        sys.stdout.write(line)
        continue
    f = line.rstrip("\\n").split("\\t")
    info = f[7].split(";")
    for i, x in enumerate(info):
        if x.startswith("ANN="):
            ann = x[4:].split(",")
            top = [a for a in ann if len(a.split("|")) > 3 and a.split("|")[3] in panel]
            info[i] = "ANN=" + ",".join(top + [a for a in ann if a not in top])
    f[7] = ";".join(info)
    sys.stdout.write("\\t".join(f) + "\\n")
PYEOF
    rm -f filtered.vcf
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
    path acmg_shared    // bin/acmg_shared.py — staged for in-heredoc Python import
    path protein_index  // ClinVar protein index (PS1 / PM5), or [] (no PS1 / PM5)

    output:
    tuple val(sample), path("${sample}.germline_snvs.tsv"),              emit: tsv
    path("${sample}.germline_all_filtered.tsv"),                          emit: all_filtered_tsv
    tuple val(sample), path("${sample}.germline_snvs_acmg_annotated.tsv"), emit: acmg_tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\tZYGOSITY\\tClinVar_sig\\tClinVar_stars\\tClinVar_source\\tREVEL_score\\tREVEL_source\\tSpliceAI_DS\\tSpliceAI_detail\\tClinVar_PS1\\tClinVar_PM5\\tgnomAD_FAF95\\tgnomAD_AC\\tgnomAD_NHOMALT\\tgnomAD_AF_GRPMAX\\tgnomAD_AC_GRPMAX\\tgnomAD_AN_GRPMAX\\tgnomAD_BA1\\tQC_FLAG\\n" > ${sample}.germline_snvs.tsv
    touch ${sample}.germline_all_filtered.tsv
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\tZYGOSITY\\tClinVar_sig\\tClinVar_stars\\tClinVar_source\\tREVEL_score\\tREVEL_source\\tSpliceAI_DS\\tSpliceAI_detail\\tClinVar_PS1\\tClinVar_PM5\\tQC_FLAG\\tgnomAD_AF\\tgnomAD_AF_source\\tgnomAD_FAF95\\tgnomAD_AC\\tgnomAD_NHOMALT\\tgnomAD_AF_GRPMAX\\tgnomAD_AC_GRPMAX\\tgnomAD_AN_GRPMAX\\tgnomAD_BA1\\tACMG_class\\tACMG_criteria\\tACMG_spec\\tinsil_dam\\tinsil_ben\\tinsil_total\\n" > ${sample}.germline_snvs_acmg_annotated.tsv
    """

    script:
    """
    # Filter thresholds are available as env vars from nextflow.config env {} block

    # Local gnomAD fields exist only when GNOMAD_ANNOTATE ran (--gnomad_tsv);
    # SnpSift aborts on INFO fields missing from the header, so add them conditionally.
    GNOMAD_FIELDS=""
    GNOMAD_VERSION=""
    # (read the whole header: an early-closing pipe would trip pipefail with SIGPIPE)
    zcat "${vcf}" | grep '^##' > vcf_header.txt || true
    if grep -q '^##INFO=<ID=GNOMAD_AF,' vcf_header.txt; then
        GNOMAD_FIELDS="GNOMAD_REGION GNOMAD_FILTER GNOMAD_AF GNOMAD_AC GNOMAD_AN GNOMAD_FAF95 GNOMAD_GRPMAX GNOMAD_AF_AMR GNOMAD_AN_AMR GNOMAD_NHOMALT GNOMAD_NHOMALT_XY"
        GNOMAD_VERSION=\$(sed -n 's/^##gnomAD_resource=<Version=\\([^,]*\\),.*/\\1/p' vcf_header.txt)
        # grpmax AF/AC/AN: resources built with them (ClinGen VCEP frequency rules)
        if grep -q '^##INFO=<ID=GNOMAD_AF_GRPMAX,' vcf_header.txt; then
            GNOMAD_FIELDS="\${GNOMAD_FIELDS} GNOMAD_AF_GRPMAX GNOMAD_AC_GRPMAX GNOMAD_AN_GRPMAX"
        fi
    fi
    # Same for the local ClinVar resource (CLINVAR_ANNOTATE, --clinvar_tsv)
    CLINVAR_FIELDS=""
    CLINVAR_VERSION=""
    if grep -q '^##INFO=<ID=CLINVAR_SIG,' vcf_header.txt; then
        CLINVAR_FIELDS="CLINVAR_ID CLINVAR_SIG CLINVAR_REVSTAT CLINVAR_CONF CLINVAR_DN"
        CLINVAR_VERSION=\$(sed -n 's/^##ClinVar_resource=<Version=\\([^,]*\\),.*/\\1/p' vcf_header.txt)
    fi
    # Local REVEL (REVEL_ANNOTATE, --revel_tsv) and SpliceAI (SPLICEAI, --spliceai)
    REVEL_FIELDS=""
    REVEL_VERSION=""
    if grep -q '^##INFO=<ID=REVEL_SCORE,' vcf_header.txt; then
        REVEL_FIELDS="REVEL_SCORE REVEL_TX REVEL_AA"
        REVEL_VERSION=\$(sed -n 's/^##REVEL_resource=<Version=\\([^,]*\\),.*/\\1/p' vcf_header.txt)
    fi
    SPLICEAI_FIELDS=""
    if grep -q '^##INFO=<ID=SpliceAI,' vcf_header.txt; then
        SPLICEAI_FIELDS="SpliceAI"
    fi
    rm -f vcf_header.txt
    export GNOMAD_VERSION CLINVAR_VERSION REVEL_VERSION

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
        \${GNOMAD_FIELDS} \${CLINVAR_FIELDS} \${REVEL_FIELDS} \${SPLICEAI_FIELDS} \\
        > "${sample}.germline_snvs_raw.tsv"

    # ── Step 1: Normalize multi-annotation fields (keep first value before ;) ──────
    #   e.g. "HIGH;MODIFIER" → "HIGH",  "BRCA1;NBR2" → "BRCA1"
    #   SpliceAI keeps every value (one per gene; Step 1.8 takes the maximum)
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{for(i=1;i<=NF;i++) if(\$i=="SpliceAI") keep=i; print; next} {for(i=1;i<=NF;i++){if(i==keep) continue; n=split(\$i,a,";"); \$i=a[1]}; print}' \\
        "${sample}.germline_snvs_raw.tsv" > "${sample}.germline_all_extracted.tsv"
    rm -f "${sample}.germline_snvs_raw.tsv"

    # ── Step 1.5: population frequency → canonical gnomAD_AF / gnomAD_AF_source ──
    #   LOCAL mode (GNOMAD_* columns present, --gnomad_tsv): no network access.
    #     gnomAD_v<ver>          AF from the local gnomAD joint resource (max over alleles)
    #     gnomAD_v<ver>_absent   inside resource regions, not in gnomAD -> AF 0 (PM2 valid)
    #     gnomAD_v<ver>_qc_fail  fails QC in both gnomAD exomes and genomes: AF kept for
    #                            BA1/BS1, never used for PM2
    #     dbNSFP_v2              outside resource regions -> dbNSFP gnomAD v2.1.1 AF
    #     unknown                outside regions and no dbNSFP AF -> no BA1/BS1/PM2
    #   gnomAD_FAF95: max filtering AF (95% CI) across ancestry groups from the local
    #     resource (0 when absent, "." outside the regions) -> BA1/BS1/PM2 use it
    #     with frequency_metric = 'faf95'; falls back to gnomAD_AF when ".".
    #   gnomAD_AC, gnomAD_NHOMALT, gnomAD_AF_GRPMAX, gnomAD_AC_GRPMAX, gnomAD_AN_GRPMAX
    #     (resource built with the grpmax columns): ClinGen VCEP gene-specific
    #     BA1/BS1/PM2 rules (acmg_shared VCEP_SPECS); 0 when absent, "." outside.
    #   LEGACY mode (no --gnomad_tsv): behaviour below, unchanged.
    #
    # ── Legacy: gnomAD AF rescue (Ensembl VEP + gnomAD GraphQL fallback) ─────────
    #   dbNSFP only covers SNVs — for indels the gnomAD AF column is empty.
    #   This step queries:
    #     1) Ensembl VEP REST API (/vep/human/region) — fast, but often lacks
    #        gnomAD frequencies for indels (no colocated_variants returned).
    #     2) gnomAD GraphQL API (gnomad.broadinstitute.org/api) as fallback
    #        when Ensembl returns no frequency data.  Uses gnomAD v4 dataset.
    #   Adds a gnomAD_AF_source column:  dbNSFP | ensembl | gnomad_api | .
    #   Falls back gracefully if the network is unavailable.
    python3 - <<'PYEOF_GNOMAD'
import sys, json, os, time
try:
    from urllib.request import Request, urlopen
    from urllib.error   import HTTPError, URLError
    _NET = True
except ImportError:
    _NET = False

ENSEMBL   = 'https://rest.ensembl.org'
GNOMAD_GQL = 'https://gnomad.broadinstitute.org/api'
MAX_Q     = int(os.environ.get("FILT_GNOMAD_RESCUE_MAX", "50"))
SLEEP     = 0.35   # seconds between API calls

# ── Source 1: Ensembl VEP ──────────────────────────────────────────────────
def _ensembl_gnomad_af(chrom, pos, ref, alt):
    # Return (af, True) from Ensembl VEP, or (None, False) if no freq data.
    if not _NET:
        return None, False
    c   = chrom.lstrip('chr') or chrom
    end = int(pos) + len(ref) - 1
    url = ('{b}/vep/human/region/{c}:{p}-{e}:1/{a}'
           '?content-type=application/json&gnomadg=1&gnomade=1').format(
           b=ENSEMBL, c=c, p=pos, e=end, a=alt)
    for attempt in range(3):
        try:
            req = Request(url, headers={'Accept': 'application/json'})
            with urlopen(req, timeout=20) as r:
                results = json.loads(r.read().decode())
            best = None
            for entry in results:
                for cv in entry.get('colocated_variants', []):
                    for pop, pdata in cv.get('frequencies', {}).items():
                        af = float(pdata.get('af') or 0)
                        if best is None or af > best:
                            best = af
            return best, (best is not None)
        except HTTPError as e:
            if e.code in (429, 503):
                time.sleep(2 ** (attempt + 1))
                continue
            return None, False
        except (URLError, Exception):
            return None, False
    return None, False

# ── Source 2: gnomAD GraphQL API (v4) ──────────────────────────────────────
def _gnomad_graphql_af(chrom, pos, ref, alt):
    # Return max(exome,genome) AF from gnomAD v4 GraphQL, or None.
    if not _NET:
        return None
    c = chrom.lstrip('chr') or chrom
    vid = c + '-' + str(pos) + '-' + ref + '-' + alt
    query = '{ variant(variantId: "' + vid + '", dataset: gnomad_r4) { exome { af } genome { af } } }'
    body = json.dumps({'query': query}).encode()
    for attempt in range(3):
        try:
            req = Request(GNOMAD_GQL,
                          data=body,
                          headers={'Content-Type': 'application/json',
                                   'Accept': 'application/json'})
            with urlopen(req, timeout=30) as r:
                data = json.loads(r.read().decode())
            v = (data.get('data') or {}).get('variant')
            if v is None:
                return None      # variant not in gnomAD
            best = 0.0
            for subset in ('exome', 'genome'):
                sub = v.get(subset)
                if sub and sub.get('af') is not None:
                    af = float(sub['af'])
                    if af > best:
                        best = af
            return best
        except HTTPError as e:
            if e.code in (429, 502, 503):
                time.sleep(2 ** (attempt + 1))
                continue
            return None
        except (URLError, Exception):
            return None
    return None

# ── Main: iterate variants, rescue missing AFs ────────────────────────────
TSV = '${sample}.germline_all_extracted.tsv'
TMP = '${sample}.germline_all_extracted_rescue_tmp.tsv'

with open(TSV) as fh:
    lines = [l.rstrip('\\n') for l in fh if l.strip()]

if len(lines) <= 1:
    sys.exit(0)

hdr_list = lines[0].split('\\t')
hdr      = {c: i for i, c in enumerate(hdr_list)}
af_i     = hdr.get('dbNSFP_gnomAD_exomes_AF', -1)
chr_i    = hdr.get('CHROM', -1)
pos_i    = hdr.get('POS',   -1)
ref_i    = hdr.get('REF',   -1)
alt_i    = hdr.get('ALT',   -1)

if af_i < 0 or chr_i < 0:
    sys.exit(0)

def _num(v):
    # dbNSFP/gnomAD fields may hold one value per transcript/allele: "0.53,0.53"
    vals = []
    for x in str(v).split(','):
        try:
            vals.append(float(x))
        except (TypeError, ValueError):
            pass
    return max(vals) if vals else None

# ── LOCAL mode: frequencies from the local gnomAD resource, no network ──────
if 'GNOMAD_REGION' in hdr:
    tag = 'gnomAD_v' + (os.environ.get('GNOMAD_VERSION', '') or 'local')
    gi  = {k: hdr.get(k, -1) for k in ('GNOMAD_REGION', 'GNOMAD_AF', 'GNOMAD_FILTER', 'GNOMAD_FAF95',
                                       'GNOMAD_AC', 'GNOMAD_NHOMALT', 'GNOMAD_AF_GRPMAX',
                                       'GNOMAD_AC_GRPMAX', 'GNOMAD_AN_GRPMAX',
                                       'dbNSFP_gnomAD_exomes_AF', 'dbNSFP_gnomAD_genomes_AF')}
    # grpmax columns only exist in resources built with them (ClinGen VCEP rules)
    has_grp = gi['GNOMAD_AF_GRPMAX'] >= 0
    EXTRA = ['GNOMAD_AC', 'GNOMAD_NHOMALT', 'GNOMAD_AF_GRPMAX', 'GNOMAD_AC_GRPMAX', 'GNOMAD_AN_GRPMAX']
    def maxval(row, k):
        # grpmax "." for an allele seen only in bottlenecked groups (ASJ, FIN, ...),
        # which the VCEP rules exclude -> 0
        v = _num(col(row, k))
        if v is None:
            return '0'
        return '%d' % v if k != 'GNOMAD_AF_GRPMAX' else repr(v)
    def col(row, k):
        i = gi[k]
        return row[i].strip() if 0 <= i < len(row) else '.'
    n_base = len(hdr_list)
    hdr_list += ['gnomAD_AF', 'gnomAD_AF_source', 'gnomAD_FAF95']
    if has_grp:
        hdr_list += ['gnomAD_' + k[7:] for k in EXTRA]
    counts = {}
    with open(TMP, 'w') as out:
        out.write('\\t'.join(hdr_list) + '\\n')
        for line in lines[1:]:
            row = line.split('\\t')
            while len(row) < n_base:
                row.append('.')
            afs   = col(row, 'GNOMAD_AF').split(',')
            filts = col(row, 'GNOMAD_FILTER').split(',')
            vals  = [(_num(a), filts[k] if k < len(filts) else '.') for k, a in enumerate(afs)]
            vals  = [(a, f) for a, f in vals if a is not None]
            region = col(row, 'GNOMAD_REGION')
            fafs   = [x for x in (_num(f) for f in col(row, 'GNOMAD_FAF95').split(',')) if x is not None]
            faf    = repr(max(fafs)) if (vals and fafs) else '.'
            if vals:
                good = [a for a, f in vals if f != 'BOTH_FILTERED']
                if good:
                    af, src = repr(max(good)), tag
                else:
                    # Fails QC in exomes AND genomes (paralogs, repeats): the AF is
                    # kept for BA1/BS1 (benign direction) but never supports PM2.
                    af, src = repr(max(a for a, _f in vals)), tag + '_qc_fail'
            elif region not in ('', '.'):
                af, src, faf = '0', tag + '_absent', '0'
            else:
                v2 = [x for x in (_num(col(row, 'dbNSFP_gnomAD_exomes_AF')),
                                  _num(col(row, 'dbNSFP_gnomAD_genomes_AF'))) if x is not None]
                af, src = (repr(max(v2)), 'dbNSFP_v2') if v2 else ('.', 'unknown')
            counts[src] = counts.get(src, 0) + 1
            extra = []
            if has_grp:
                # inside the regions: absent -> 0; outside -> "." (general rules)
                if vals:
                    extra = [maxval(row, k) for k in EXTRA]
                else:
                    extra = ['0' if src.endswith('_absent') else '.'] * len(EXTRA)
            out.write('\\t'.join(row[:n_base] + [af, src, faf] + extra) + '\\n')
    sys.stderr.write('[gnomAD local] ${sample}: ' + ', '.join(
        k + '=' + str(v) for k, v in sorted(counts.items())) + '\\n')
    os.replace(TMP, TSV)
    sys.exit(0)

# ── LEGACY mode ─────────────────────────────────────────────────────────────
# Add gnomAD_AF_source column
hdr_list.append('gnomAD_AF_source')
src_i = len(hdr_list) - 1

n = 0
with open(TMP, 'w') as out:
    out.write('\\t'.join(hdr_list + ['gnomAD_AF']) + '\\n')
    for line in lines[1:]:
        row = line.split('\\t')
        while len(row) < src_i:
            row.append('.')
        row.append('.')          # gnomAD_AF_source placeholder

        af_val = row[af_i].strip()
        needs_rescue = af_val in ('', '.', 'NA', '0', '0.0', '0.00')

        if not needs_rescue:
            row[src_i] = 'dbNSFP'
        elif n < MAX_Q:
            ch  = row[chr_i] if chr_i < len(row) else ''
            pos = row[pos_i] if pos_i < len(row) else ''
            ref = row[ref_i] if ref_i < len(row) else ''
            alt = row[alt_i] if alt_i < len(row) else ''
            if ch and pos and ref not in ('', '.') and alt not in ('', '.'):
                # Try Ensembl VEP first
                af, has_freq = _ensembl_gnomad_af(ch, pos, ref, alt)
                source = None
                if has_freq and af is not None and af > 0:
                    source = 'ensembl'
                # Fallback: gnomAD GraphQL API (also when Ensembl returns 0.0)
                if source is None:
                    time.sleep(SLEEP)
                    af_gql = _gnomad_graphql_af(ch, pos, ref, alt)
                    if af_gql is not None:
                        af = af_gql
                        source = 'gnomad_api'
                    elif has_freq and af is not None:
                        # Ensembl had data (0.0) and gnomAD has no entry
                        source = 'ensembl'
                if source and af is not None:
                    row[af_i] = str(af)
                    row[src_i] = source
                    n += 1
                    sys.stderr.write(
                        '[gnomAD rescue] ${sample} ' + ch + ':' + pos +
                        ' ' + ref + '>' + alt + ' AF=' + str(af) +
                        ' (' + source + ')\\n'
                    )
                time.sleep(SLEEP)
        out.write('\\t'.join(row + [row[af_i]]) + '\\n')

sys.stderr.write('[gnomAD rescue] ${sample}: rescued ' + str(n) + ' variant(s)\\n')
os.replace(TMP, TSV)
PYEOF_GNOMAD

    # ── Step 1.6: zygosity (ZYGOSITY column) ─────────────────────────────────────
    #   From GEN[0].GT (phased "|" treated as "/"):
    #     HET            different alleles (0/1, 1/2)
    #     HOM            same non-ref allele on an autosome (1/1)
    #     HOM/HEMI       1/1 on chrX/chrY outside the PARs (sex not used)
    #     HET(VAF alta)  heterozygous call with VAF > vaf_max (allelic imbalance / artefact?)
    #   No GT → inferred from VAF: >= vaf_max → HOM (HOM/HEMI on non-PAR X/Y), else HET.
    #   No GT and no AD → "."
    #   High VAF is no longer a rejection reason: homozygous / hemizygous variants
    #   in recessive genes (MUTYH, NTHL1, MSH3, MBD4, FANC*, CMMRD) must be reported.
    awk -F'\\t' 'BEGIN { OFS="\\t"; vmax = ENVIRON["FILT_VAF_MAX"]+0 }
    function nonpar(c, p) {
        # GRCh38 pseudoautosomal regions PAR1 / PAR2
        if (c == "chrX") return !((p >= 10001 && p <= 2781479) || (p >= 155701383 && p <= 156030895))
        if (c == "chrY") return !((p >= 10001 && p <= 2781479) || (p >= 56887903 && p <= 57217415))
        return 0
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if (\$i == "CHROM")     chr_col = i
            if (\$i == "POS")       pos_col = i
            if (\$i == "GEN[0].GT") gt_col  = i
            if (\$i == "GEN[0].AD") ad_col  = i
        }
        print \$0, "ZYGOSITY"; next
    }
    {
        vaf = -1
        if (ad_col) {
            n = split(\$(ad_col), arr, ",")
            if (n >= 2) { tot = arr[1] + arr[2]; if (tot > 0) vaf = arr[2] / tot }
        }
        gt = gt_col ? \$(gt_col) : "."
        gsub(/\\|/, "/", gt)
        na = split(gt, al, "/")
        z = "."
        if (na == 2 && al[1] != "." && al[2] != ".") {
            if (al[1] != al[2])      z = (vaf > vmax) ? "HET(VAF alta)" : "HET"
            else if (al[1] != "0")   z = "HOM"
        } else if (vaf >= 0) {
            z = (vaf >= vmax) ? "HOM" : "HET"
        }
        if (z == "HOM" && nonpar(\$(chr_col), \$(pos_col) + 0)) z = "HOM/HEMI"
        print \$0, z
    }' "${sample}.germline_all_extracted.tsv" > "${sample}.germline_all_extracted.zyg.tsv"
    mv "${sample}.germline_all_extracted.zyg.tsv" "${sample}.germline_all_extracted.tsv"

    # ── Step 1.7: canonical ClinVar columns ─────────────────────────────────────
    #   ClinVar_sig, ClinVar_stars, ClinVar_revstat, ClinVar_id, ClinVar_conf,
    #   ClinVar_dn, ClinVar_source. Every downstream consumer reads ClinVar_sig.
    #   LOCAL mode (CLINVAR_* columns, --clinvar_tsv): the local release is the only
    #     source; a variant absent from it is "." (no dbNSFP fallback).
    #     ClinVar_source = ClinVar_<release date>
    #   LEGACY mode: dbNSFP 4.1a ClinVar fields (~2020, SNVs only).
    #     ClinVar_source = dbNSFP_4.1a
    #   Stars from the review status: practice guideline 4, expert panel 3,
    #   multiple submitters no conflicts 2, single submitter / conflicting 1, else 0.
    python3 - <<'PYEOF_CLINVAR'
import os

TSV = '${sample}.germline_all_extracted.tsv'
TMP = '${sample}.germline_all_extracted.clinvar_tmp.tsv'
VER = os.environ.get('CLINVAR_VERSION', '')

def stars(revstat, sig):
    if sig in ('', '.'):
        return '.'
    r = revstat.lower()
    if 'practice_guideline' in r:
        return '4'
    if 'expert_panel' in r:
        return '3'
    if 'multiple_submitters' in r and 'no_conflicts' in r:
        return '2'
    if 'single_submitter' in r or 'conflicting' in r:
        return '1'
    return '0'

with open(TSV) as fin, open(TMP, 'w') as out:
    hdr = fin.readline().rstrip('\\n').split('\\t')
    idx = {c: i for i, c in enumerate(hdr)}
    local = 'CLINVAR_SIG' in idx
    if local:
        src_cols = ('CLINVAR_SIG', 'CLINVAR_REVSTAT', 'CLINVAR_ID', 'CLINVAR_CONF', 'CLINVAR_DN')
        source = 'ClinVar_' + (VER or 'local')
    else:
        src_cols = ('dbNSFP_clinvar_clnsig', 'dbNSFP_clinvar_review', 'dbNSFP_clinvar_id', None, 'dbNSFP_clinvar_trait')
        source = 'dbNSFP_4.1a'
    out.write('\\t'.join(hdr + ['ClinVar_sig', 'ClinVar_stars', 'ClinVar_revstat', 'ClinVar_id',
                                'ClinVar_conf', 'ClinVar_dn', 'ClinVar_source']) + '\\n')
    for line in fin:
        row = line.rstrip('\\n').split('\\t')
        def g(c):
            i = idx.get(c, -1) if c else -1
            v = row[i].strip() if 0 <= i < len(row) else ''
            return v if v else '.'
        sig, rev, cid, conf, dn = (g(c) for c in src_cols)
        if sig == '.':
            rev = conf = dn = '.'
        out.write('\\t'.join(row + [sig, stars(rev, sig), rev, cid, conf, dn, source]) + '\\n')
os.replace(TMP, TSV)
PYEOF_CLINVAR

    # ── Step 1.8: canonical REVEL and SpliceAI columns ──────────────────────────
    #   REVEL_score / REVEL_source:
    #     LOCAL (REVEL_SCORE column, --revel_tsv): REVEL v1.3 score of the MANE Select
    #       transcript (REVEL_v1.3_MANE); when REVEL lacks the MANE transcript ID, the
    #       score of the amino-acid change of the MANE HGVS_p (REVEL_v1.3_MANE_aa), else
    #       the highest over other transcripts (REVEL_v1.3_other_tx); no dbNSFP fallback.
    #     LEGACY: dbNSFP 4.1a REVEL (first transcript), source dbNSFP_4.1a.
    #   SpliceAI_DS: max of DS_AG / DS_AL / DS_DG / DS_DL over the genes scored for the
    #     ALT allele ("." when SpliceAI did not run); SpliceAI_detail: SYMBOL:type=DS@pos.
    python3 - <<'PYEOF_PRED'
import os

TSV = '${sample}.germline_all_extracted.tsv'
TMP = '${sample}.germline_all_extracted.pred_tmp.tsv'
VER = os.environ.get('REVEL_VERSION', '') or '1.3'

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

AA3 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E',
       'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
       'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}

def mane_aa(hgvs_p):
    # p.Arg165Gln -> "R>Q" (missense only)
    h = hgvs_p[2:] if hgvs_p.startswith('p.') else ''
    ref, alt = AA3.get(h[:3]), AA3.get(h[-3:])
    return ref + '>' + alt if ref and alt and h[3:-3].isdigit() else None

with open(TSV) as fin, open(TMP, 'w') as out:
    hdr = fin.readline().rstrip('\\n').split('\\t')
    idx = {c: i for i, c in enumerate(hdr)}
    local_revel = 'REVEL_SCORE' in idx
    out.write('\\t'.join(hdr + ['REVEL_score', 'REVEL_source', 'SpliceAI_DS', 'SpliceAI_detail']) + '\\n')
    for line in fin:
        row = line.rstrip('\\n').split('\\t')
        row += ['.'] * (len(hdr) - len(row))
        g = lambda c: row[idx[c]].strip() if c in idx else '.'
        if local_revel:
            sc, tx = num(g('REVEL_SCORE')), g('REVEL_TX')
            if sc is not None and tx == 'other_tx':
                aa = mane_aa(g('ANN[*].HGVS_P'))
                opts = dict(o.rsplit(':', 1) for o in g('REVEL_AA').split(',') if ':' in o)
                if aa in opts and num(opts[aa]) is not None:
                    sc, tx = num(opts[aa]), 'MANE_aa'
            revel = '%.3f' % sc if sc is not None else '.'
            src = ('REVEL_v' + VER + '_' + tx) if sc is not None else '.'
        else:
            sc = num(g('dbNSFP_REVEL_score').split(',')[0])
            revel = '%.3f' % sc if sc is not None else '.'
            src = 'dbNSFP_4.1a' if sc is not None else '.'
        ds, detail = '.', '.'
        alt = g('ALT')
        for entry in g('SpliceAI').replace(';', ',').split(','):
            f = entry.split('|')
            if len(f) < 10 or f[0] != alt:
                continue
            for k, lab in enumerate(('AG', 'AL', 'DG', 'DL')):
                v = num(f[2 + k])
                if v is not None and (ds == '.' or v > float(ds)):
                    ds = '%.2f' % v
                    detail = '%s:%s=%.2f@%s' % (f[1], lab, v, f[6 + k])
        out.write('\\t'.join(row + [revel, src, ds, detail]) + '\\n')
os.replace(TMP, TSV)
PYEOF_PRED

    # ── Step 1.9: ClinVar same-residue evidence (ACMG PS1 / PM5) ─────────────────
    #   From the ClinVar protein index (scripts/build_clinvar_protein_index.py,
    #   --clinvar_protein_index), for missense variants (MANE HGVS_p):
    #     ClinVar_PS1  P/LP missense with the same amino-acid change, other nucleotide
    #     ClinVar_PM5  P/LP missense with another amino-acid change at the same residue
    #   References need >= clinvar_min_stars stars and are not at exon edges
    #   (SPLICE_REGION, where pathogenicity may come from splicing). "." = none.
    python3 - <<'PYEOF_PS1'
import gzip, os, re

TSV   = '${sample}.germline_all_extracted.tsv'
TMP   = '${sample}.germline_all_extracted.ps1_tmp.tsv'
INDEX = '${protein_index}'
MIN_STARS = int(os.environ.get('ACMG_CLINVAR_MIN_STARS', '2') or 2)
HGVS = re.compile(r'^p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})\$')

ref = {}
if INDEX and os.path.exists(INDEX):
    with gzip.open(INDEX, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\\n').split('\\t')
            if int(f[6]) < MIN_STARS or f[12] == '1':
                continue
            ref.setdefault((f[0], int(f[1])), []).append(f)

with open(TSV) as fin, open(TMP, 'w') as out:
    hdr = fin.readline().rstrip('\\n').split('\\t')
    idx = {c: i for i, c in enumerate(hdr)}
    out.write('\\t'.join(hdr + ['ClinVar_PS1', 'ClinVar_PM5']) + '\\n')
    for line in fin:
        row = line.rstrip('\\n').split('\\t')
        row += ['.'] * (len(hdr) - len(row))
        g = lambda c: row[idx[c]].strip() if c in idx else '.'
        ps1, pm5 = [], []
        m = HGVS.match(g('ANN[*].HGVS_P'))
        if m and 'missense' in g('ANN[*].EFFECT'):
            key = (g('CHROM'), g('POS'), g('REF'), g('ALT'))
            for f in ref.get((g('ANN[*].GENE'), int(m.group(2))), []):
                if (f[8], f[9], f[10], f[11]) == key:
                    continue
                txt = '%s VCV%s %s %s*' % (f[4], f[7], f[5], f[6])
                (ps1 if f[3] == m.group(3) else pm5).append(txt)
        out.write('\\t'.join(row + [', '.join(ps1) or '.', ', '.join(pm5) or '.']) + '\\n')
os.replace(TMP, TSV)
PYEOF_PS1

    # ── Step 1.95: BA1 by the ACMG engine rule set (column gnomAD_BA1) ───────────
    #   acmg_shared.ba1_reason(): ClinGen VCEP gene threshold (VCEP_SPECS) when the
    #   gene has one and the gnomAD grpmax columns are present, else ba1_af on
    #   FAF95 / AF. Text of the Step 2 reason, "." = BA1 not met.
    python3 - <<'PYEOF_BA1'
import os, sys
sys.path.insert(0, '.')
from acmg_shared import ba1_reason, getter_for_row_hdr

TSV = '${sample}.germline_all_extracted.tsv'
TMP = '${sample}.germline_all_extracted.ba1_tmp.tsv'
with open(TSV) as fin, open(TMP, 'w') as out:
    hdr = fin.readline().rstrip('\\n').split('\\t')
    idx = {c: i for i, c in enumerate(hdr)}
    out.write('\\t'.join(hdr + ['gnomAD_BA1']) + '\\n')
    for line in fin:
        row = line.rstrip('\\n').split('\\t')
        row += ['.'] * (len(hdr) - len(row))
        out.write('\\t'.join(row + [ba1_reason(getter_for_row_hdr(row, idx)) or '.']) + '\\n')
os.replace(TMP, TSV)
PYEOF_BA1

    # ── Step 2: gnomAD BA1 + ClinVar + low-VAF filter ────────────────────────────────
    #   Passing variants  → germline_snvs.pre_python.tsv  (stdout)
    #   Rejected variants → germline_all_filtered.tsv  (with FILTER_REASON column)
    #
    #   Rejection reasons:
    #     multiallelic    : more than one ALT allele (e.g. GT 1/2 in the AR CAG/GGC
    #                       repeats); not reportable in a predisposition test
    #     gnomad_BA1      : common polymorphism (ACMG BA1), column gnomAD_BA1 of
    #                       Step 1.95: ClinGen VCEP gene threshold, else gnomAD
    #                       FAF95 (or AF) > ba1_af (5%). Never for ClinVar P/LP
    #                       with >= clinvar_min_stars stars (founder pathogenic
    #                       variants): they reach the ACMG step with their class
    #     clinvar_benign  : ClinVar Benign/Likely_benign (unless also Pathogenic) with
    #                       >= clinvar_min_stars review stars (default 2); B/LB with
    #                       fewer stars stays and is classified by the ACMG score
    #     spliceai_below_threshold : non-coding candidate (LOW / MODIFIER impact:
    #                       splice region, synonymous, intronic) with SpliceAI_DS below
    #                       spliceai_pp3 (0.2) or not scored, unless ClinVar P/LP
    #     vaf_too_low     : VAF below vaf_min_high / vaf_min_moderate (default 0.20),
    #                       only when DP >= vaf_floor_min_dp (default 20) and the
    #                       variant is not ClinVar P/LP (non-conflicting)
    #   High VAF is NOT a rejection reason (homozygous/hemizygous): see ZYGOSITY (Step 1.6)
    #
    #   Adds QC_FLAG ("." or comma-separated):
    #     DP<N      depth below vaf_floor_min_dp: VAF floor not applied
    #     VAF_baja  ClinVar P/LP kept despite VAF below the floor
    #     CNN_SNP / CNN_INDEL / LowQual / <other>   VCF FILTER != PASS (flag only)
    awk -F'\\t' -v FILT="${sample}.germline_all_filtered.tsv" '
    BEGIN { OFS="\\t" }
    NR == 1 {
        print \$0 "\\tQC_FLAG\\tFILTER_REASON" > FILT
        min_dp_floor = ENVIRON["FILT_VAF_FLOOR_MIN_DP"]+0
        cv_min_stars = (ENVIRON["ACMG_CLINVAR_MIN_STARS"] != "") ? ENVIRON["ACMG_CLINVAR_MIN_STARS"]+0 : 2
        use_faf = (ENVIRON["ACMG_FREQ_METRIC"] != "af")
        for (i = 1; i <= NF; i++) {
            if (\$i == "ALT")                        alt_col     = i
            if (\$i == "GEN[0].GT")                  gt_col      = i
            if (\$i == "ANN[*].IMPACT")              impact_col  = i
            if (\$i == "GEN[0].AD")                  ad_col      = i
            if (\$i == "GEN[0].DP")                  dp_col      = i
            if (\$i == "FILTER")                     filter_col  = i
            if (\$i == "ClinVar_sig")                 cs_col      = i
            if (\$i == "ClinVar_stars")               st_col      = i
            if (\$i == "dbNSFP_gnomAD_exomes_AF")     gnomad_v2   = i
            if (\$i == "gnomAD_AF")                   gnomad_can  = i
            if (\$i == "gnomAD_FAF95")                faf_col     = i
            if (\$i == "gnomAD_BA1")                  ba1_col     = i
            if (\$i == "SpliceAI_DS")                 sai_col     = i
        }
        sai_min = (ENVIRON["ACMG_SPLICEAI_PP3"] != "") ? ENVIRON["ACMG_SPLICEAI_PP3"]+0 : 0.2
        gnomad_col = gnomad_can ? gnomad_can : gnomad_v2
        print \$0 "\\tQC_FLAG"; next
    }
    {
        reason = ""
        flag   = ""
        clnsig = (cs_col && cs_col <= NF) ? tolower(\$(cs_col)) : "."
        clinvar_plp = (clnsig ~ /pathogenic/ && clnsig !~ /conflicting/ && clnsig !~ /benign/)
        # Multi-allelic site (ALT "A,B"): not reportable, kept only in the filtered table
        if (alt_col && \$(alt_col) ~ /,/) {
            gt = (gt_col && gt_col <= NF) ? \$(gt_col) : "."
            reason = "multiallelic (GT=" gt "; ALT=" \$(alt_col) ")"
        }
        cv_st0 = (st_col && st_col <= NF) ? \$(st_col) : ""
        clinvar_final_plp = clinvar_plp && ((cv_st0 !~ /^[0-9]+\$/) || (cv_st0 + 0 >= cv_min_stars))
        # BA1: common polymorphism — rule set of the ACMG engine (Step 1.95)
        if (reason == "" && ba1_col && ba1_col <= NF) {
            if (\$(ba1_col) != "." && !clinvar_final_plp) reason = \$(ba1_col)
        } else if (reason == "" && gnomad_col && gnomad_col <= NF) {
            # older tables without gnomAD_BA1: gnomAD FAF95 (or AF) > ba1_af
            af  = \$(gnomad_col) + 0
            faf = (use_faf && faf_col && faf_col <= NF && \$(faf_col) ~ /^[0-9]/) ? \$(faf_col) : ""
            if (faf != "") {
                if (faf + 0 > ENVIRON["ACMG_BA1_AF"]+0)
                    reason = sprintf("gnomad_BA1 (FAF95=%.4f; AF=%.4f)", faf + 0, af)
            } else if (\$(gnomad_col) != "." && af > ENVIRON["ACMG_BA1_AF"]+0) {
                reason = sprintf("gnomad_BA1 (AF=%.4f)", af)
            }
        }
        if (reason == "" && cs_col && cs_col <= NF) {
            # No stars column (older TSVs) or stars not numeric: accept as before
            cv_st = (st_col && st_col <= NF) ? \$(st_col) : ""
            cv_ok = (cv_st !~ /^[0-9]+\$/) || (cv_st + 0 >= cv_min_stars)
            if (clnsig != "." && clnsig ~ /benign/ && clnsig !~ /pathogenic/ && cv_ok) {
                reason = "clinvar_benign (clnsig=" \$(cs_col) (cv_st ~ /^[0-9]+\$/ ? ", " cv_st "*" : "") ")"
            }
        }
        # Non-coding splice candidate (LOW / MODIFIER, GS15): kept only with a SpliceAI
        # prediction >= spliceai_pp3 or ClinVar P/LP (any stars)
        if (reason == "" && impact_col && \$(impact_col) !~ /HIGH|MODERATE/ && !clinvar_plp) {
            ds = (sai_col && sai_col <= NF) ? \$(sai_col) : "."
            if (ds !~ /^[0-9]/ || ds + 0 < sai_min)
                reason = "spliceai_below_threshold (non-coding, DS=" ds ")"
        }
        if (reason == "" && ad_col && ad_col <= NF) {
            n = split(\$(ad_col), arr, ",")
            if (n >= 2) {
                ref = arr[1] + 0; alt = arr[2] + 0; tot = ref + alt
                if (tot > 0) {
                    vaf = alt / tot
                    dp  = (dp_col && \$(dp_col) ~ /^[0-9]+\$/) ? \$(dp_col) + 0 : tot
                    min_vaf = (impact_col && \$(impact_col) ~ /HIGH/) ? ENVIRON["FILT_VAF_MIN_HIGH"]+0 : ENVIRON["FILT_VAF_MIN_MOD"]+0
                    if (dp < min_dp_floor) {
                        flag = "DP<" min_dp_floor
                    } else if (vaf < min_vaf) {
                        if (clinvar_plp) flag = "VAF_baja"
                        else reason = sprintf("vaf_too_low (VAF=%.3f, DP=%d, impact=%s)", vaf, dp, \$(impact_col))
                    }
                }
            }
        }
        # VCF FILTER (GATK CNN tranches): flag only, never a rejection reason
        if (filter_col) {
            fv = \$(filter_col)
            if (fv != "PASS" && fv != "." && fv != "") {
                lab = (fv ~ /^CNN_1D_SNP/) ? "CNN_SNP" : ((fv ~ /^CNN_1D_INDEL/) ? "CNN_INDEL" : fv)
                flag = (flag == "") ? lab : flag "," lab
            }
        }
        if (flag == "") flag = "."
        if (reason != "") {
            print \$0 "\\t" flag "\\t" reason > FILT
            next
        }
        print \$0 "\\t" flag
    }' "${sample}.germline_all_extracted.tsv" > "${sample}.germline_snvs.pre_python.tsv"
    rm -f "${sample}.germline_all_extracted.tsv"

    # ── Step 3: Multi-hit nearby filter (alignment artefacts) ───────────────────
    #   A pair = two variants of the same gene <= multi_hit_distance bp apart (30).
    #   Each pair member is judged on its own and removed when it has
    #     VAF < multi_hit_vaf (0.40)            or
    #     FILTER != PASS  or  GQ < multi_hit_min_gq (20)
    #   A member with normal VAF and good quality is kept.
    #   Exempt: HIGH-impact variants (PVS1) and ClinVar pathogenic/likely_pathogenic.
    #   Rejected variants are appended (with reason and partner) to germline_all_filtered.tsv.
    python3 - <<'PYEOF'
import sys, os
from collections import defaultdict

TSV_IN    = "${sample}.germline_snvs.pre_python.tsv"
TSV_CLEAN = "${sample}.germline_snvs.tsv"
TSV_FILT  = "${sample}.germline_all_filtered.tsv"

MAX_DIST = int(os.environ.get("FILT_MULTI_HIT_DIST", "30"))
MAX_VAF  = float(os.environ.get("FILT_MULTI_HIT_VAF", "0.40"))
MIN_GQ   = float(os.environ.get("FILT_MULTI_HIT_MIN_GQ", "20"))

with open(TSV_IN) as f:
    lines = [l.rstrip("\\n") for l in f if l.strip()]

if len(lines) <= 1:
    with open(TSV_CLEAN, "w") as fc:
        fc.write(lines[0] + "\\n" if lines else "")
    sys.exit(0)

hdr_list = lines[0].split("\\t")
hdr = {c: i for i, c in enumerate(hdr_list)}

def g(row, col):
    i = hdr.get(col, -1)
    return row[i].strip() if 0 <= i < len(row) else "."

def compute_vaf(row):
    try:
        alt = int(g(row, "GEN[0].AD").split(",")[1])
        tot = int(g(row, "GEN[0].DP"))
        return alt / tot if tot > 0 else None
    except Exception:
        return None

def is_exempt(row):
    # Always keep HIGH-impact variants (PVS1 evidence -> ACMG LP minimum)
    if g(row, "ANN[*].IMPACT").upper() == "HIGH":
        return True
    # Always keep ClinVar pathogenic / likely_pathogenic
    clnsig = g(row, "ClinVar_sig").lower()
    if clnsig and clnsig != ".":
        if "pathogenic" in clnsig and "benign" not in clnsig.replace("pathogenic", ""):
            return True
    return False

def bad_reasons(row):
    # Why this pair member looks like an artefact ([] = keep)
    why = []
    v = compute_vaf(row)
    if v is not None and v < MAX_VAF:
        why.append("VAF={:.3f}".format(v))
    filt = g(row, "FILTER")
    if filt not in (".", "", "PASS"):
        why.append("FILTER=" + filt)
    try:
        gq = float(g(row, "GEN[0].GQ"))
        if gq < MIN_GQ:
            why.append("GQ={:g}".format(gq))
    except ValueError:
        pass
    return why

data = [l.split("\\t") for l in lines[1:]]

gene_grp = defaultdict(list)
for i, row in enumerate(data):
    gene_grp[g(row, "ANN[*].GENE")].append(i)

to_filter = {}
for gene, idxs in gene_grp.items():
    if len(idxs) < 2:
        continue
    for a in range(len(idxs)):
        for b in range(a + 1, len(idxs)):
            i1, i2 = idxs[a], idxs[b]
            r1, r2 = data[i1], data[i2]
            try:
                dist = abs(int(g(r1, "POS")) - int(g(r2, "POS")))
            except ValueError:
                continue
            if dist > MAX_DIST:
                continue
            for idx, row, mate in [(i1, r1, r2), (i2, r2, r1)]:
                if idx in to_filter or is_exempt(row):
                    continue
                why = bad_reasons(row)
                if why:
                    to_filter[idx] = (
                        "multi_hit_nearby ({}; gene={}, pair={}:{} {}>{}, dist={} bp)".format(
                            "; ".join(why), gene, g(mate, "CHROM"), g(mate, "POS"),
                            g(mate, "REF"), g(mate, "ALT"), dist)
                    )

with open(TSV_CLEAN, "w") as fc, open(TSV_FILT, "a") as ff:
    fc.write(lines[0] + "\\n")
    for i, row in enumerate(data):
        if i in to_filter:
            ff.write("\\t".join(row) + "\\t" + to_filter[i] + "\\n")
        else:
            fc.write("\\t".join(row) + "\\n")
PYEOF

    rm -f "${sample}.germline_snvs.pre_python.tsv"

    # ── Step 4: ACMG annotation of all candidate variants ────────────────────────
    #   Reads germline_snvs.tsv and appends:
    #     ACMG_class    : P / LP / VUS / LB / B
    #     ACMG_criteria : pipe-separated evidence codes (e.g. PM2_Sup|PP3_Mod)
    #     ACMG_spec     : rule set used (ClinGen VCEP GNxxx vX.Y, or ACMG/SVI general)
    #     insil_dam     : number of in silico tools voting damaging
    #     insil_ben     : number voting benign/tolerated
    #     insil_total   : total tools with a clear verdict
    #   Output: ${sample}.germline_snvs_acmg_annotated.tsv
    python3 - <<'PYEOF'
import sys, os

# ACMG logic from bin/acmg_shared.py (staged into work dir as 'acmg_shared.py')
sys.path.insert(0, '.')
from acmg_shared import (
    count_insil as _count_insil,
    acmg_criteria_list as _acmg_criteria_list_shared,
    acmg_spec as _acmg_spec,
    getter_for_row_hdr,
)

def acmg_criteria_list(row, hdr):
    return _acmg_criteria_list_shared(getter_for_row_hdr(row, hdr))

# --- Annotate all candidate variants ---
TSV_IN  = "${sample}.germline_snvs.tsv"
TSV_OUT = "${sample}.germline_snvs_acmg_annotated.tsv"

with open(TSV_IN) as fh:
    lines = [l.rstrip("\\n") for l in fh if l.strip()]

if not lines:
    open(TSV_OUT, 'w').close()
    sys.exit(0)

raw_header = lines[0].split("\\t")
hdr = {c: i for i, c in enumerate(raw_header)}

def get_col(row, name):
    i = hdr.get(name, -1)
    return row[i].strip() if 0 <= i < len(row) else '.'

new_cols = ['ACMG_class', 'ACMG_criteria', 'ACMG_spec', 'insil_dam', 'insil_ben', 'insil_total']
with open(TSV_OUT, 'w') as out:
    out.write("\\t".join(raw_header + new_cols) + "\\n")
    for line in lines[1:]:
        row = line.split("\\t")
        acmg, criteria = acmg_criteria_list(row, hdr)
        dam, ben = _count_insil(lambda name, r=row: get_col(r, name))
        total = dam + ben
        extra = [acmg, '|'.join(criteria) if criteria else '.',
                 _acmg_spec(getter_for_row_hdr(row, hdr)), str(dam), str(ben), str(total)]
        out.write("\\t".join(row + extra) + "\\n")
PYEOF
    """
}

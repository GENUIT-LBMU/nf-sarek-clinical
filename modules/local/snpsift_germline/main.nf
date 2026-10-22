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
        | gzip -c > ${sample}.filtered.vcf.gz
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

    output:
    tuple val(sample), path("${sample}.germline_snvs.tsv"),              emit: tsv
    path("${sample}.germline_all_filtered.tsv"),                          emit: all_filtered_tsv
    tuple val(sample), path("${sample}.germline_snvs_acmg_annotated.tsv"), emit: acmg_tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\tZYGOSITY\\tQC_FLAG\\n" > ${sample}.germline_snvs.tsv
    touch ${sample}.germline_all_filtered.tsv
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\tZYGOSITY\\tQC_FLAG\\tgnomAD_AF\\tgnomAD_AF_source\\tACMG_class\\tACMG_criteria\\tinsil_dam\\tinsil_ben\\tinsil_total\\n" > ${sample}.germline_snvs_acmg_annotated.tsv
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
    fi
    rm -f vcf_header.txt
    export GNOMAD_VERSION

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
        \${GNOMAD_FIELDS} \\
        > "${sample}.germline_snvs_raw.tsv"

    # ── Step 1: Normalize multi-annotation fields (keep first value before ;) ──────
    #   e.g. "HIGH;MODIFIER" → "HIGH",  "BRCA1;NBR2" → "BRCA1"
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{print; next} {for(i=1;i<=NF;i++){n=split(\$i,a,";"); \$i=a[1]}; print}' \\
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
    gi  = {k: hdr.get(k, -1) for k in ('GNOMAD_REGION', 'GNOMAD_AF', 'GNOMAD_FILTER',
                                       'dbNSFP_gnomAD_exomes_AF', 'dbNSFP_gnomAD_genomes_AF')}
    def col(row, k):
        i = gi[k]
        return row[i].strip() if 0 <= i < len(row) else '.'
    n_base = len(hdr_list)
    hdr_list += ['gnomAD_AF', 'gnomAD_AF_source']
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
            if vals:
                good = [a for a, f in vals if f != 'BOTH_FILTERED']
                if good:
                    af, src = repr(max(good)), tag
                else:
                    # Fails QC in exomes AND genomes (paralogs, repeats): the AF is
                    # kept for BA1/BS1 (benign direction) but never supports PM2.
                    af, src = repr(max(a for a, _f in vals)), tag + '_qc_fail'
            elif region not in ('', '.'):
                af, src = '0', tag + '_absent'
            else:
                v2 = [x for x in (_num(col(row, 'dbNSFP_gnomAD_exomes_AF')),
                                  _num(col(row, 'dbNSFP_gnomAD_genomes_AF'))) if x is not None]
                af, src = (repr(max(v2)), 'dbNSFP_v2') if v2 else ('.', 'unknown')
            counts[src] = counts.get(src, 0) + 1
            out.write('\\t'.join(row[:n_base] + [af, src]) + '\\n')
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

    # ── Step 2: gnomAD BA1 + ClinVar + low-VAF filter ────────────────────────────────
    #   Passing variants  → germline_snvs.pre_python.tsv  (stdout)
    #   Rejected variants → germline_all_filtered.tsv  (with FILTER_REASON column)
    #
    #   Rejection reasons:
    #     gnomad_BA1      : gnomAD AF > 5% (common polymorphism — ACMG BA1), read
    #                       from the canonical gnomAD_AF column ("." = unknown)
    #     clinvar_benign  : ClinVar Benign/Likely_benign (unless also Pathogenic)
    #     vaf_too_low     : VAF below vaf_min_high / vaf_min_moderate (default 0.20),
    #                       only when DP >= vaf_floor_min_dp (default 20) and the
    #                       variant is not ClinVar P/LP (non-conflicting)
    #   High VAF is NOT a rejection reason (homozygous/hemizygous): see ZYGOSITY (Step 1.6)
    #
    #   Adds QC_FLAG ("." or comma-separated):
    #     DP<N      depth below vaf_floor_min_dp: VAF floor not applied
    #     VAF_baja  ClinVar P/LP kept despite VAF below the floor
    awk -F'\\t' -v FILT="${sample}.germline_all_filtered.tsv" '
    BEGIN { OFS="\\t" }
    NR == 1 {
        print \$0 "\\tQC_FLAG\\tFILTER_REASON" > FILT
        min_dp_floor = ENVIRON["FILT_VAF_FLOOR_MIN_DP"]+0
        for (i = 1; i <= NF; i++) {
            if (\$i == "ANN[*].IMPACT")              impact_col  = i
            if (\$i == "GEN[0].AD")                  ad_col      = i
            if (\$i == "GEN[0].DP")                  dp_col      = i
            if (\$i == "dbNSFP_clinvar_clnsig")       cs_col      = i
            if (\$i == "dbNSFP_gnomAD_exomes_AF")     gnomad_v2   = i
            if (\$i == "gnomAD_AF")                   gnomad_can  = i
        }
        gnomad_col = gnomad_can ? gnomad_can : gnomad_v2
        print \$0 "\\tQC_FLAG"; next
    }
    {
        reason = ""
        flag   = ""
        clnsig = (cs_col && cs_col <= NF) ? tolower(\$(cs_col)) : "."
        clinvar_plp = (clnsig ~ /pathogenic/ && clnsig !~ /conflicting/ && clnsig !~ /benign/)
        # BA1: common polymorphism — gnomAD exomes AF > 5%
        if (reason == "" && gnomad_col && gnomad_col <= NF) {
            af = \$(gnomad_col) + 0
            if (\$(gnomad_col) != "." && af > ENVIRON["ACMG_BA1_AF"]+0) {
                reason = sprintf("gnomad_BA1 (AF=%.4f)", af)
            }
        }
        if (reason == "" && cs_col && cs_col <= NF) {
            if (clnsig != "." && clnsig ~ /benign/ && clnsig !~ /pathogenic/) {
                reason = "clinvar_benign (clnsig=" \$(cs_col) ")"
            }
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
        if (flag == "") flag = "."
        if (reason != "") {
            print \$0 "\\t" flag "\\t" reason > FILT
            next
        }
        print \$0 "\\t" flag
    }' "${sample}.germline_all_extracted.tsv" > "${sample}.germline_snvs.pre_python.tsv"
    rm -f "${sample}.germline_all_extracted.tsv"

    # ── Step 3: Multi-hit same-exon low-VAF filter ────────────────────────────────
    #   Removes variants clustered in the same gene + same exon (ANN[*].RANK match
    #   OR positions ≤ 30 bp apart) that also have unbalanced VAF < 0.40.
    #   Exempt: HIGH-impact variants (PVS1) and ClinVar pathogenic/likely_pathogenic.
    #   Rejected variants are appended (with reason) to germline_all_filtered.tsv.
    python3 - <<'PYEOF'
import sys, os
from collections import defaultdict

TSV_IN    = "${sample}.germline_snvs.pre_python.tsv"
TSV_CLEAN = "${sample}.germline_snvs.tsv"
TSV_FILT  = "${sample}.germline_all_filtered.tsv"

INFRAME_DEL_EFFECTS = frozenset([
    "disruptive_inframe_deletion",
    "conservative_inframe_deletion",
])

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
    clnsig = g(row, "dbNSFP_clinvar_clnsig").lower()
    if clnsig and clnsig != ".":
        if "pathogenic" in clnsig and "benign" not in clnsig.replace("pathogenic", ""):
            return True
    return False

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
            rank_same = g(r1, "ANN[*].RANK") == g(r2, "ANN[*].RANK")
            try:
                nearby = abs(int(g(r1, "POS")) - int(g(r2, "POS"))) <= int(os.environ.get("FILT_MULTI_HIT_DIST", "30"))
            except ValueError:
                nearby = False
            if rank_same or nearby:
                for idx, row in [(i1, r1), (i2, r2)]:
                    if idx in to_filter or is_exempt(row):
                        continue
                    effect = g(row, "ANN[*].EFFECT")
                    exon   = g(row, "ANN[*].RANK")
                    v      = compute_vaf(row)
                    if effect in INFRAME_DEL_EFFECTS:
                        to_filter[idx] = (
                            "multi_hit_same_exon_inframe_deletion "
                            "(effect={}, gene={}, exon={})".format(effect, gene, exon)
                        )
                    elif v is not None and v < float(os.environ.get("FILT_MULTI_HIT_VAF", "0.40")):
                        to_filter[idx] = (
                            "multi_hit_same_exon_low_vaf "
                            "(VAF={:.3f}, gene={}, exon={})".format(v, gene, exon)
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
    #     ACMG_criteria : pipe-separated evidence codes (e.g. PM2|BP4(9/9))
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

new_cols = ['ACMG_class', 'ACMG_criteria', 'insil_dam', 'insil_ben', 'insil_total']
with open(TSV_OUT, 'w') as out:
    out.write("\\t".join(raw_header + new_cols) + "\\n")
    for line in lines[1:]:
        row = line.split("\\t")
        acmg, criteria = acmg_criteria_list(row, hdr)
        dam, ben = _count_insil(lambda name, r=row: get_col(r, name))
        total = dam + ben
        extra = [acmg, '|'.join(criteria) if criteria else '.', str(dam), str(ben), str(total)]
        out.write("\\t".join(row + extra) + "\\n")
PYEOF
    """
}

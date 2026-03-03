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

    output:
    tuple val(sample), path("${sample}.germline_snvs.tsv"),              emit: tsv
    path("${sample}.germline_all_filtered.tsv"),                          emit: all_filtered_tsv
    tuple val(sample), path("${sample}.germline_snvs_acmg_annotated.tsv"), emit: acmg_tsv

    stub:
    """
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\n" > ${sample}.germline_snvs.tsv
    touch ${sample}.germline_all_filtered.tsv
    printf "CHROM\\tPOS\\tID\\tREF\\tALT\\tANN[*].GENE\\tANN[*].IMPACT\\tACMG_class\\tACMG_criteria\\tinsil_dam\\tinsil_ben\\tinsil_total\\n" > ${sample}.germline_snvs_acmg_annotated.tsv
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
        > "${sample}.germline_snvs_raw.tsv"

    # ── Step 1: Normalize multi-annotation fields (keep first value before ;) ──────
    #   e.g. "HIGH;MODIFIER" → "HIGH",  "BRCA1;NBR2" → "BRCA1"
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{print; next} {for(i=1;i<=NF;i++){n=split(\$i,a,";"); \$i=a[1]}; print}' \\
        "${sample}.germline_snvs_raw.tsv" > "${sample}.germline_all_extracted.tsv"
    rm -f "${sample}.germline_snvs_raw.tsv"

    # ── Step 1.5: gnomAD AF rescue via Ensembl VEP REST API ──────────────────────
    #   dbNSFP only covers SNVs — for indels the gnomAD AF column is empty.
    #   This step queries the Ensembl VEP REST API (/vep/human/region) for any
    #   variant missing a gnomAD AF, filling the field before the BA1 filter runs.
    #   The Ensembl VEP endpoint returns gnomAD exome+genome frequencies for both
    #   SNVs and indels using standard GRCh38 VCF coordinates.
    #   Falls back gracefully (silent, no change) if the network is unavailable.
    python3 - <<'PYEOF_GNOMAD'
import sys, json, os, time
try:
    from urllib.request import Request, urlopen
    from urllib.error   import HTTPError, URLError
    _NET = True
except ImportError:
    _NET = False

ENSEMBL = 'https://rest.ensembl.org'
MAX_Q   = 50     # safety cap: max variants to rescue per sample
SLEEP   = 0.35   # seconds between API calls (Ensembl limit: ~15 req/s)

def vep_gnomad_af(chrom, pos, ref, alt):
    # Return max gnomAD AF from Ensembl VEP colocated_variants, or None on error.
    if not _NET:
        return None
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
            best = 0.0
            for entry in results:
                for cv in entry.get('colocated_variants', []):
                    for pop, pdata in cv.get('frequencies', {}).items():
                        af = float(pdata.get('af') or 0.0)
                        if af > best:
                            best = af
            return best
        except HTTPError as e:
            if e.code in (429, 503):
                time.sleep(2 ** (attempt + 1))
                continue
            sys.stderr.write('[gnomAD rescue] HTTP ' + str(e.code) + ' for ' +
                             chrom + ':' + pos + ' ' + ref + '>' + alt + '\\n')
            return None
        except (URLError, Exception) as e:
            sys.stderr.write('[gnomAD rescue] ' + chrom + ':' + pos + ' ' +
                             ref + '>' + alt + ': ' + str(e) + '\\n')
            return None
    return None

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
    sys.exit(0)   # column not present — skip silently

n = 0
with open(TMP, 'w') as out:
    out.write(lines[0] + '\\n')
    for line in lines[1:]:
        row = line.split('\\t')
        while len(row) < len(hdr_list):   # guard against truncated rows
            row.append('.')
        af_val = row[af_i].strip()
        if af_val in ('', '.', 'NA') and n < MAX_Q:
            ch  = row[chr_i] if chr_i < len(row) else ''
            pos = row[pos_i] if pos_i < len(row) else ''
            ref = row[ref_i] if ref_i < len(row) else ''
            alt = row[alt_i] if alt_i < len(row) else ''
            if ch and pos and ref not in ('', '.') and alt not in ('', '.'):
                af = vep_gnomad_af(ch, pos, ref, alt)
                if af is not None:
                    row[af_i] = str(af)
                    n += 1
                    sys.stderr.write(
                        '[gnomAD rescue] ${sample} ' + ch + ':' + pos +
                        ' ' + ref + '>' + alt + ' AF=' + str(af) + '\\n'
                    )
                    time.sleep(SLEEP)
        out.write('\\t'.join(row) + '\\n')

sys.stderr.write('[gnomAD rescue] ${sample}: rescued ' + str(n) + ' variant(s)\\n')
os.replace(TMP, TSV)
PYEOF_GNOMAD

    # ── Step 2: gnomAD BA1 + ClinVar + VAF filter ────────────────────────────────
    #   Passing variants  → germline_snvs.pre_python.tsv  (stdout)
    #   Rejected variants → germline_all_filtered.tsv  (with FILTER_REASON column)
    #
    #   Rejection reasons:
    #     gnomad_BA1      : gnomAD exomes AF > 5% (common polymorphism — ACMG BA1)
    #     clinvar_benign  : ClinVar Benign/Likely_benign (unless also Pathogenic)
    #     vaf_too_low     : VAF below threshold (HIGH-impact: 0.25, MODERATE: 0.35)
    #     vaf_too_high    : VAF > 0.75 (likely homozygous or artefact)
    awk -F'\\t' -v FILT="${sample}.germline_all_filtered.tsv" '
    BEGIN { OFS="\\t" }
    NR == 1 {
        print \$0 "\\tFILTER_REASON" > FILT
        for (i = 1; i <= NF; i++) {
            if (\$i == "ANN[*].IMPACT")              impact_col  = i
            if (\$i == "GEN[0].AD")                  ad_col      = i
            if (\$i == "dbNSFP_clinvar_clnsig")       cs_col      = i
            if (\$i == "dbNSFP_gnomAD_exomes_AF")     gnomad_col  = i
        }
        print; next
    }
    {
        reason = ""
        # BA1: common polymorphism — gnomAD exomes AF > 5%
        if (reason == "" && gnomad_col && gnomad_col <= NF) {
            af = \$(gnomad_col) + 0
            if (\$(gnomad_col) != "." && af > 0.05) {
                reason = sprintf("gnomad_BA1 (AF=%.4f)", af)
            }
        }
        if (reason == "" && cs_col && cs_col <= NF) {
            clnsig = tolower(\$(cs_col))
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
                    min_vaf = (impact_col && \$(impact_col) ~ /HIGH/) ? 0.25 : 0.35
                    if (vaf < min_vaf) {
                        reason = sprintf("vaf_too_low (VAF=%.3f, impact=%s)", vaf, \$(impact_col))
                    } else if (vaf > 0.75) {
                        reason = sprintf("vaf_too_high (VAF=%.3f)", vaf)
                    }
                }
            }
        }
        if (reason != "") {
            print \$0 "\\t" reason > FILT
            next
        }
        print
    }' "${sample}.germline_all_extracted.tsv" > "${sample}.germline_snvs.pre_python.tsv"
    rm -f "${sample}.germline_all_extracted.tsv"

    # ── Step 3: Multi-hit same-exon low-VAF filter ────────────────────────────────
    #   Removes variants clustered in the same gene + same exon (ANN[*].RANK match
    #   OR positions ≤ 30 bp apart) that also have unbalanced VAF < 0.40.
    #   Exempt: HIGH-impact variants (PVS1) and ClinVar pathogenic/likely_pathogenic.
    #   Rejected variants are appended (with reason) to germline_all_filtered.tsv.
    python3 - <<'PYEOF'
import sys
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
                nearby = abs(int(g(r1, "POS")) - int(g(r2, "POS"))) <= 30
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
                    elif v is not None and v < 0.40:
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
import sys

# --- Inline ACMG helpers (mirrors bin/final_report_acmg.py) ---
_INSIL_DAMAGING = {
    'dbNSFP_SIFT_pred':             ('D',),
    'dbNSFP_SIFT4G_pred':           ('D',),
    'dbNSFP_Polyphen2_HDIV_pred':   ('D', 'P'),
    'dbNSFP_Polyphen2_HVAR_pred':   ('D', 'P'),
    'dbNSFP_LRT_pred':              ('D',),
    'dbNSFP_MutationTaster_pred':   ('A', 'D'),
    'dbNSFP_MutationAssessor_pred': ('H', 'M'),
    'dbNSFP_FATHMM_pred':           ('D',),
    'dbNSFP_PROVEAN_pred':          ('D',),
}
_INSIL_BENIGN = {
    'dbNSFP_SIFT_pred':             ('T',),
    'dbNSFP_SIFT4G_pred':           ('T',),
    'dbNSFP_Polyphen2_HDIV_pred':   ('B',),
    'dbNSFP_Polyphen2_HVAR_pred':   ('B',),
    'dbNSFP_LRT_pred':              ('N',),
    'dbNSFP_MutationTaster_pred':   ('N', 'P'),
    'dbNSFP_MutationAssessor_pred': ('L', 'N'),
    'dbNSFP_FATHMM_pred':           ('T',),
    'dbNSFP_PROVEAN_pred':          ('N',),
}

def _count_insil(col_fn):
    dam = ben = 0
    for field in _INSIL_DAMAGING:
        # dbNSFP stores multi-transcript predictions comma-separated; take first.
        v = col_fn(field).split(',')[0].strip().upper()
        if v in ('', '.', 'NA'):
            continue
        if v in _INSIL_DAMAGING[field]:
            dam += 1
        elif v in _INSIL_BENIGN.get(field, ()):
            ben += 1
    return dam, ben

def acmg_criteria_list(row, hdr):
    def col(name):
        i = hdr.get(name, -1)
        return row[i].strip() if i >= 0 and i < len(row) else '.'

    clnsig = col('dbNSFP_clinvar_clnsig').lower()
    if clnsig and clnsig != '.':
        if 'conflicting' in clnsig:
            pass
        elif ('pathogenic' in clnsig
              and 'likely_pathogenic' not in clnsig
              and 'likely pathogenic' not in clnsig):
            if 'likely' not in clnsig:
                return 'P', ['ClinVar:P']
        elif 'likely_pathogenic' in clnsig or 'likely pathogenic' in clnsig:
            return 'LP', ['ClinVar:LP']
        elif 'uncertain' in clnsig:
            return 'VUS', ['ClinVar:VUS']
        elif 'likely_benign' in clnsig or 'likely benign' in clnsig:
            return 'LB', ['ClinVar:LB']
        elif 'benign' in clnsig and 'pathogenic' not in clnsig:
            return 'B', ['ClinVar:B']

    try:
        gnomad = float(col('dbNSFP_gnomAD_exomes_AF'))
    except (ValueError, TypeError):
        gnomad = 0.0
    if gnomad > 0.05:
        return 'B', ['BA1']

    score = 0
    criteria = []
    impact = col('ANN[*].IMPACT').upper()
    if impact == 'HIGH':
        score += 8
        criteria.append('PVS1')
    if gnomad < 0.001:
        score += 2
        criteria.append('PM2')

    try:
        revel = float(col('dbNSFP_REVEL_score'))
        if revel >= 0.932:
            score += 4; criteria.append('PP3_St')
        elif revel >= 0.644:
            score += 2; criteria.append('PP3')
        elif revel >= 0.5:
            score += 1; criteria.append('PP3_Sup')
        if revel < 0.15:
            score -= 4; criteria.append('BP4_St')
        elif revel < 0.29:
            score -= 2; criteria.append('BP4')
    except (ValueError, TypeError):
        pass

    if gnomad > 0.01:
        score -= 4
        criteria.append('BS1')

    dam, ben = _count_insil(col)
    total = dam + ben
    if total > 0:
        if ben > dam:
            ratio = ben / total
            if ratio >= 0.75:
                score -= 2; criteria.append('BP4({}/{})'.format(ben, total))
            else:
                score -= 1; criteria.append('BP4_Sup({}/{})'.format(ben, total))
        elif dam > ben:
            ratio = dam / total
            if ratio >= 0.75:
                score += 2; criteria.append('PP3({}/{})'.format(dam, total))
            else:
                score += 1; criteria.append('PP3_Sup({}/{})'.format(dam, total))

    if score >= 10:   return 'P', criteria
    if score >= 6:    return 'LP', criteria
    if score <= -6:   return 'LB', criteria
    return 'VUS', criteria

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

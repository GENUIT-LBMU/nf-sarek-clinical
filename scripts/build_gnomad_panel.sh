#!/usr/bin/env bash
# =============================================================================
# build_gnomad_panel.sh — build the local gnomAD resource used by the
# germline annotation step (GNOMAD_ANNOTATE).
#
# Output (in --outdir):
#   gnomad_v<VER>_panel.tsv.gz(.tbi)   one row per gnomAD allele inside the panel
#                                      regions, a handful of fields (see COLUMNS)
#   gnomad_v<VER>_panel_regions.bed.gz(.tbi)
#                                      the regions that were queried. A variant
#                                      inside these regions and absent from the
#                                      TSV is truly absent from gnomAD; a variant
#                                      outside them has UNKNOWN frequency.
#   gnomad_v<VER>_panel.manifest.txt   version, source URLs, genes, padding, md5
#
# Regions = MANE exons (incl. UTR) of every gene in the --genes files, +/- --padding bp,
# primary chromosomes only (gnomAD has no alt/fix contigs).
# Data are read remotely with tabix from the public gnomAD bucket, so only the
# panel slices are downloaded. Needs bcftools/tabix/bgzip built with libcurl:
#
#   docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -w $PWD \
#       quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
#       bash scripts/build_gnomad_panel.sh \
#           --genes assets/gene_panels/germline_103_cancer_predisposition.txt \
#           --genes assets/gene_panels/somatic_snv_panel.txt \
#           --genes assets/gene_panels/somatic_38_oncogenes.txt \
#           --mane  /path/to/exons_mane.txt --outdir /path/to/refs/gnomad
#
# gnomAD data: https://gnomad.broadinstitute.org (cite Chen et al. 2024, Nature).
# =============================================================================
set -euo pipefail

VERSION="4.1"
PADDING=100
OUTDIR=""
MANE=""
GENE_FILES=()
BASE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release"

usage() { sed -n '2,30p' "$0"; exit 1; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genes)    GENE_FILES+=("$2"); shift 2 ;;
        --mane)     MANE="$2";          shift 2 ;;
        --outdir)   OUTDIR="$2";        shift 2 ;;
        --padding)  PADDING="$2";       shift 2 ;;
        --version)  VERSION="$2";       shift 2 ;;
        --base-url) BASE_URL="$2";      shift 2 ;;
        -h|--help)  usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done
[[ ${#GENE_FILES[@]} -gt 0 && -n "$MANE" && -n "$OUTDIR" ]] || usage
for t in bcftools tabix bgzip; do command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done

mkdir -p "$OUTDIR"
PREFIX="${OUTDIR}/gnomad_v${VERSION}_panel"
TMP="$(cd "$(mktemp -d "${OUTDIR}/.build_XXXXXX")" && pwd)"
trap 'rm -rf "$TMP"' EXIT

# ── 1. Gene list ──────────────────────────────────────────────────────────────
cat "${GENE_FILES[@]}" | grep -v '^#' | awk 'NF{print $1}' | sort -u > "$TMP/genes.txt"
N_GENES=$(wc -l < "$TMP/genes.txt")

# ── 2. Regions: MANE exons +/- padding, merged, BED (0-based) with gene name ──
awk -v pad="$PADDING" 'BEGIN{OFS="\t"}
    NR==FNR {want[$1]=1; next}
    FNR==1 && $1=="chr" {next}
    ($4 in want) && $1 ~ /^chr([0-9]+|X|Y)$/ {s=$2-1-pad; if (s<0) s=0; print $1, s, $3+pad, $4}' \
    "$TMP/genes.txt" "$MANE" \
  | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"}
        $1==c && $2<=e { if ($3>e) e=$3; if (index(","g",", ","$4",")==0) g=g","$4; next }
        { if (c!="") print c,s,e,g; c=$1; s=$2; e=$3; g=$4 }
        END { if (c!="") print c,s,e,g }' > "$TMP/regions.bed"

MISSING=$(comm -23 "$TMP/genes.txt" <(cut -f4 "$TMP/regions.bed" | tr ',' '\n' | sort -u) || true)
if [[ -n "$MISSING" ]]; then
    echo "ERROR: genes not found in MANE file: $(echo $MISSING)" >&2; exit 1
fi
echo "Genes: ${N_GENES}  regions: $(wc -l < "$TMP/regions.bed")  bp: $(awk '{s+=$3-$2}END{print s}' "$TMP/regions.bed")"

# ── 3. Query gnomAD per chromosome ────────────────────────────────────────────
# Columns kept (joint = exomes + genomes). FILTER is kept so that sites failing
# QC in both datasets can be treated as "frequency unknown" downstream.
FMT='%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC_joint\t%INFO/AN_joint\t%INFO/AF_joint\t%INFO/nhomalt_joint\t%INFO/fafmax_faf95_max_joint\t%INFO/grpmax_joint\t%INFO/AF_joint_amr\t%INFO/AN_joint_amr\t%INFO/nhomalt_joint_XY\n'
for chrom in $(cut -f1 "$TMP/regions.bed" | uniq); do
    url="${BASE_URL}/${VERSION}/vcf/joint/gnomad.joint.v${VERSION}.sites.${chrom}.vcf.bgz"
    awk -v c="$chrom" '$1==c' "$TMP/regions.bed" > "$TMP/r.${chrom}.bed"
    echo "  ${chrom}: $(wc -l < "$TMP/r.${chrom}.bed") regions"
    ok=0
    for attempt in 1 2 3; do
        # run inside $TMP: htslib caches the remote .tbi in the working directory
        if (cd "$TMP" && bcftools query -R "r.${chrom}.bed" -f "$FMT" "$url" > "q.${chrom}.tsv" 2> "err.${chrom}"); then
            ok=1; break
        fi
        echo "    retry ${attempt} ($(tail -1 "$TMP/err.${chrom}"))"; sleep $((attempt * 10))
    done
    [[ $ok == 1 ]] || { echo "ERROR: failed to read ${url}" >&2; cat "$TMP/err.${chrom}" >&2; exit 1; }
done

# ── 4. Assemble, deduplicate (overlapping -R regions), compress, index ─────────
{
    printf '## gnomad_version=%s source=%s padding_bp=%s genes=%s\n' "${VERSION}" "gnomad.joint.v${VERSION}.sites" "${PADDING}" "${N_GENES}"
    printf '#CHROM\tPOS\tREF\tALT\tGNOMAD_FILTER\tGNOMAD_AC\tGNOMAD_AN\tGNOMAD_AF\tGNOMAD_NHOMALT\tGNOMAD_FAF95\tGNOMAD_GRPMAX\tGNOMAD_AF_AMR\tGNOMAD_AN_AMR\tGNOMAD_NHOMALT_XY\n'
    for chrom in $(cut -f1 "$TMP/regions.bed" | uniq); do
        sort -k2,2n -k3,3 -k4,4 -u "$TMP/q.${chrom}.tsv"
    done
} | bgzip -c > "${PREFIX}.tsv.gz"
tabix -f -s1 -b2 -e2 -c'#' "${PREFIX}.tsv.gz"

bgzip -c "$TMP/regions.bed" > "${PREFIX}_regions.bed.gz"
tabix -f -p bed "${PREFIX}_regions.bed.gz"

# ── 5. Manifest ───────────────────────────────────────────────────────────────
{
    echo "resource:   gnomAD joint (exomes+genomes) sites, panel subset"
    echo "version:    ${VERSION}"
    echo "source:     ${BASE_URL}/${VERSION}/vcf/joint/"
    echo "built:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "bcftools:   $(bcftools --version | head -1)"
    echo "mane_file:  $(basename "$MANE")"
    echo "padding_bp: ${PADDING}"
    echo "genes:      ${N_GENES} ($(paste -sd, "$TMP/genes.txt"))"
    echo "regions:    $(wc -l < "$TMP/regions.bed")"
    echo "records:    $(zcat "${PREFIX}.tsv.gz" | grep -vc '^#')"
    echo "md5_tsv:    $(md5sum "${PREFIX}.tsv.gz" | cut -d' ' -f1)"
    echo "md5_bed:    $(md5sum "${PREFIX}_regions.bed.gz" | cut -d' ' -f1)"
} > "${PREFIX}.manifest.txt"
cat "${PREFIX}.manifest.txt" | grep -v '^genes'

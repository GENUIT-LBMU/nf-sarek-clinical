#!/usr/bin/env bash
# =============================================================================
# build_lowmap_panel.sh — build the low-mappability BED used to flag germline
# CNV segments (MANE_ANNOTATE_GERMLINE_CNV, --gcnv_lowmap_bed).
#
# Output (in --outdir):
#   GRCh38_lowmappabilityall_panel.bed.gz   GIAB low-mappability regions inside
#                                           the MANE exons of the --genes files,
#                                           +/- --padding bp (BED, 0-based)
#   GRCh38_lowmappabilityall_panel.manifest.txt
#                                           version, source URL, genes, md5
#
# Source: GIAB genome stratifications v3.3, GRCh38_lowmappabilityall (union of
# Umap k100/k250 non-unique regions; covers pseudogenes and paralogs such as
# CHEK2 exons 11-15, PMS2 / PMS2CL, SDHA, FANCD2, NF1). Segmental duplications
# alone are not included (they cover most of NF1 and SDHA).
# Needs curl, awk, sort and bedtools, e.g.:
#
#   docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -w $PWD \
#       quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0 \
#       bash scripts/build_lowmap_panel.sh \
#           --genes assets/gene_panels/germline_cancer_predisposition.txt \
#           --mane  /path/to/exons_mane.txt --outdir /path/to/refs/cnv/lowmap
#
# (the bedtools image has no curl: download the BED first and pass --bed)
# GIAB stratifications: Dwarshuis et al. 2024, Nat Commun 15:9029.
# =============================================================================
set -euo pipefail

VERSION="v3.3"
PADDING=1000
OUTDIR=""
MANE=""
BED=""
GENE_FILES=()
BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications"

usage() { sed -n '2,28p' "$0"; exit 1; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genes)    GENE_FILES+=("$2"); shift 2 ;;
        --mane)     MANE="$2";          shift 2 ;;
        --outdir)   OUTDIR="$2";        shift 2 ;;
        --padding)  PADDING="$2";       shift 2 ;;
        --bed)      BED="$2";           shift 2 ;;
        --version)  VERSION="$2";       shift 2 ;;
        -h|--help)  usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done
[[ ${#GENE_FILES[@]} -gt 0 && -n "$MANE" && -n "$OUTDIR" ]] || usage
for t in bedtools sort awk; do command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done

mkdir -p "$OUTDIR"
PREFIX="${OUTDIR}/GRCh38_lowmappabilityall_panel"
TMP="$(cd "$(mktemp -d "${OUTDIR}/.build_XXXXXX")" && pwd)"
trap 'rm -rf "$TMP"' EXIT
URL="${BASE_URL}/${VERSION}/GRCh38@all/Mappability/GRCh38_lowmappabilityall.bed.gz"

# ── 1. GIAB BED ───────────────────────────────────────────────────────────────
if [[ -z "$BED" ]]; then
    command -v curl >/dev/null || { echo "ERROR: curl not found (or pass --bed)" >&2; exit 1; }
    BED="$TMP/GRCh38_lowmappabilityall.bed.gz"
    curl -sSfL -o "$BED" "$URL"
fi

# ── 2. Regions: MANE exons +/- padding (BED, 0-based) ─────────────────────────
cat "${GENE_FILES[@]}" | grep -v '^#' | awk 'NF{print $1}' | sort -u > "$TMP/genes.txt"
N_GENES=$(wc -l < "$TMP/genes.txt")
awk -v pad="$PADDING" 'BEGIN{OFS="\t"}
    NR==FNR {want[$1]=1; next}
    FNR==1 && $1=="chr" {next}
    ($4 in want) {s=$2-1-pad; if (s<0) s=0; print $1, s, $3+pad, $4}' \
    "$TMP/genes.txt" "$MANE" | sort -k1,1 -k2,2n > "$TMP/regions.bed"
MISSING=$(comm -23 "$TMP/genes.txt" <(cut -f4 "$TMP/regions.bed" | sort -u) || true)
if [[ -n "$MISSING" ]]; then
    echo "ERROR: genes not found in MANE file: $(echo $MISSING)" >&2; exit 1
fi

# ── 3. Intersect, merge, compress ─────────────────────────────────────────────
gzip -dc "$BED" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n > "$TMP/lowmap.bed"
bedtools intersect -a "$TMP/lowmap.bed" -b "$TMP/regions.bed" \
  | sort -k1,1 -k2,2n | bedtools merge > "$TMP/panel.bed"
gzip -nc "$TMP/panel.bed" > "${PREFIX}.bed.gz"

# ── 4. Manifest ───────────────────────────────────────────────────────────────
{
    echo "resource:   GIAB genome stratifications, GRCh38_lowmappabilityall, panel subset"
    echo "version:    ${VERSION}"
    echo "source:     ${URL}"
    echo "md5_source: $(md5sum "$BED" | cut -d' ' -f1)"
    echo "built:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "mane_file:  $(basename "$MANE")"
    echo "padding_bp: ${PADDING}"
    echo "genes:      ${N_GENES} ($(paste -sd, "$TMP/genes.txt"))"
    echo "intervals:  $(wc -l < "$TMP/panel.bed")"
    echo "bp:         $(awk '{s+=$3-$2}END{print s+0}' "$TMP/panel.bed")"
    echo "md5_bed:    $(md5sum "${PREFIX}.bed.gz" | cut -d' ' -f1)"
} > "${PREFIX}.manifest.txt"
grep -v '^genes' "${PREFIX}.manifest.txt"

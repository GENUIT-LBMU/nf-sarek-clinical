#!/usr/bin/env bash
# =============================================================================
# build_clinvar_panel.sh — build the local ClinVar resource used by the
# germline annotation step (CLINVAR_ANNOTATE).
#
# Output (in --outdir):
#   clinvar_<DATE>_panel.tsv.gz(.tbi)  one row per ClinVar allele inside the panel
#                                      regions: CHROM POS REF ALT + ClinVar fields
#   clinvar_<DATE>_panel.manifest.txt  version, source, regions, md5
#
# Input:
#   --clinvar  official NCBI ClinVar GRCh38 VCF (clinvar_YYYYMMDD.vcf.gz, + .tbi)
#              https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
#   --regions  panel regions BED (bgzipped + .tbi); use the one written by
#              scripts/build_gnomad_panel.sh (gnomad_v4.1_panel_regions.bed.gz)
#   --fasta    GRCh38 FASTA with chr-prefixed names (+ .fai), for left-normalisation
#
# Steps: restrict to the regions, rename contigs (1 -> chr1, MT -> chrM),
# left-normalise against --fasta (bcftools norm, records whose REF does not match
# the FASTA are dropped and counted), split multi-allelics, export the fields.
#
#   docker run --rm -u $(id -u):$(id -g) -v $HOME:$HOME -w $PWD \
#       quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
#       bash scripts/build_clinvar_panel.sh \
#           --clinvar /path/to/clinvar_20260822.vcf.gz \
#           --regions /path/to/refs/gnomad/gnomad_v4.1_panel_regions.bed.gz \
#           --fasta   /path/to/hg38.fa --outdir /path/to/refs/clinvar
#
# ClinVar is public domain (NCBI). Cite Landrum et al., Nucleic Acids Res.
# =============================================================================
set -euo pipefail

CLINVAR=""
REGIONS=""
FASTA=""
OUTDIR=""

usage() { sed -n '2,32p' "$0"; exit 1; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --clinvar) CLINVAR="$2"; shift 2 ;;
        --regions) REGIONS="$2"; shift 2 ;;
        --fasta)   FASTA="$2";   shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done
[[ -n "$CLINVAR" && -n "$REGIONS" && -n "$FASTA" && -n "$OUTDIR" ]] || usage
for t in bcftools tabix bgzip; do command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done
for f in "$CLINVAR" "${CLINVAR}.tbi" "$REGIONS" "$FASTA" "${FASTA}.fai"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

# Release date from the VCF header (##fileDate=2026-08-22 -> 20260822)
DATE=$(bcftools view -h "$CLINVAR" | sed -n 's/^##fileDate=\([0-9-]*\).*/\1/p' | tr -d '-')
[[ -n "$DATE" ]] || { echo "ERROR: no ##fileDate in $CLINVAR" >&2; exit 1; }

mkdir -p "$OUTDIR"
PREFIX="${OUTDIR}/clinvar_${DATE}_panel"
TMP="$(cd "$(mktemp -d "${OUTDIR}/.build_XXXXXX")" && pwd)"
trap 'rm -rf "$TMP"' EXIT

# ── 1. Regions in ClinVar contig naming (no chr prefix) ───────────────────────
zcat "$REGIONS" | awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1); if ($1 == "M") $1 = "MT"; print $1, $2, $3}' \
    > "$TMP/regions.nochr.bed"
for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' "$c" "$c"; done > "$TMP/chr_map.txt"
printf 'MT\tchrM\n' >> "$TMP/chr_map.txt"

# ── 2. Subset, rename, normalise, split ───────────────────────────────────────
bcftools view -R "$TMP/regions.nochr.bed" "$CLINVAR" -Ou \
  | bcftools annotate --rename-chrs "$TMP/chr_map.txt" -Ou \
  | bcftools norm -f "$FASTA" -c x -m -any -Ou 2> "$TMP/norm.log" \
  | bcftools sort -T "$TMP" -Oz -o "$TMP/panel.vcf.gz" 2> /dev/null
bcftools index -t "$TMP/panel.vcf.gz"
N_IN=$(bcftools view -H -R "$TMP/regions.nochr.bed" "$CLINVAR" | wc -l)
N_OUT=$(bcftools view -H "$TMP/panel.vcf.gz" | wc -l)

# ── 3. Export fields (already sorted), drop duplicates from overlapping -R regions ──
FMT='%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNSIGCONF\t%INFO/CLNDN\n'
{
    printf '## clinvar_version=%s source=%s\n' "${DATE}" "$(basename "$CLINVAR")"
    printf '#CHROM\tPOS\tREF\tALT\tCLINVAR_ID\tCLINVAR_SIG\tCLINVAR_REVSTAT\tCLINVAR_CONF\tCLINVAR_DN\n'
    bcftools query -f "$FMT" "$TMP/panel.vcf.gz" | awk '!seen[$0]++'
} | bgzip -c > "${PREFIX}.tsv.gz"
tabix -f -s1 -b2 -e2 -c'#' "${PREFIX}.tsv.gz"

# ── 4. Manifest ───────────────────────────────────────────────────────────────
{
    echo "resource:   ClinVar GRCh38 VCF, panel subset"
    echo "version:    ${DATE}"
    echo "source:     $(basename "$CLINVAR") (md5 $(md5sum "$CLINVAR" | cut -d' ' -f1))"
    echo "built:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "bcftools:   $(bcftools --version | head -1)"
    echo "regions:    $(basename "$REGIONS") ($(wc -l < "$TMP/regions.nochr.bed") intervals)"
    echo "fasta:      $(basename "$FASTA")"
    echo "in_regions: ${N_IN}"
    echo "normalised: ${N_OUT} (REF mismatches dropped: $(grep -c 'REF_MISMATCH\|did not match' "$TMP/norm.log" || true))"
    echo "records:    $(zcat "${PREFIX}.tsv.gz" | grep -vc '^#')"
    echo "md5_tsv:    $(md5sum "${PREFIX}.tsv.gz" | cut -d' ' -f1)"
} > "${PREFIX}.manifest.txt"
cat "${PREFIX}.manifest.txt"

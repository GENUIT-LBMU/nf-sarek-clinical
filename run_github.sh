#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# run_github.sh — WES clinical pipeline pulling nf-sarek-clinical from GitHub
# ──────────────────────────────────────────────────────────────────────────────
# Usage:
#   ./run_github.sh --run RUN_NAME --samplesheet /path/to/samplesheet.csv
#   ./run_github.sh --run RUN_NAME --samplesheet /path/to/samplesheet.csv --resume
#   ./run_github.sh --run RUN_NAME --skip-sarek           # clinical report only
#   ./run_github.sh --run RUN_NAME --skip-sarek --resume  # resume clinical report
#
# Arguments:
#   --run          Run identifier (e.g. 20260109_1229). Created at RUNS_BASE/RUN_NAME.
#   --samplesheet  Path to nf-core/sarek CSV samplesheet (Phase 1 only).
#   --resume       Resume both phases using Nextflow cache.
#   --skip-sarek   Skip Phase 1; use sarek results already in <RUN_DIR>/results.
#
# Environment variables:
#   RUNS_BASE      Base directory for run outputs (default: $HOME/runs).
#   CLINICAL_REPO  GitHub org/repo to pull (default: GENUIT-LBMU/nf-sarek-clinical).
#
# Prerequisites:
#   1. Nextflow ≥ 23.04.0 and Docker ≥ 20.10 installed on the host.
#   2. ~/.nextflow/scm configured with GitHub credentials (private repo):
#        providers {
#            github { user = 'USERNAME'; password = 'GITHUB_TOKEN' }
#        }
#   3. lab_params.yml in the same directory as this script (reference file paths).
#   4. sarek.config in the same directory as this script (nf-core/sarek lab params).
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# ── Configuration (override via env vars) ─────────────────────────────────────
RUNS_BASE="${RUNS_BASE:-$HOME/runs}"
CLINICAL_REPO="${CLINICAL_REPO:-GENUIT-LBMU/nf-sarek-clinical}"
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LAB_PARAMS="${SCRIPT_DIR}/lab_params.yml"
SAREK_CONFIG="${SCRIPT_DIR}/conf/sarek.config"

# Validate required config files
[[ -f "${LAB_PARAMS}" ]]    || { echo "ERROR: lab_params.yml not found at ${LAB_PARAMS}"; exit 1; }
[[ -f "${SAREK_CONFIG}" ]]  || { echo "ERROR: sarek.config not found at ${SAREK_CONFIG}";  exit 1; }

# ── Parse arguments ───────────────────────────────────────────────────────────
RUN_NAME=""
SAMPLESHEET=""
RESUME=false
SKIP_SAREK=false

usage() {
    sed -n '3,19p' "${BASH_SOURCE[0]}"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)         RUN_NAME="$2";    shift 2 ;;
        --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
        --resume)      RESUME=true;      shift   ;;
        --skip-sarek)  SKIP_SAREK=true;  shift   ;;
        -h|--help)     usage ;;
        *) echo "ERROR: unknown option '$1'"; usage ;;
    esac
done

[[ -z "$RUN_NAME" ]] && { echo "ERROR: --run RUN_NAME is required"; usage; }
if [[ "$SKIP_SAREK" == false && -z "$SAMPLESHEET" ]]; then
    echo "ERROR: --samplesheet is required unless --skip-sarek is set"
    usage
fi

RUN_DIR="${RUNS_BASE}/${RUN_NAME}"
NF_RESUME=""
[[ "$RESUME" == true ]] && NF_RESUME="-resume"

mkdir -p "${RUN_DIR}"

# ── Phase 1: nf-core/sarek ────────────────────────────────────────────────────
if [[ "$SKIP_SAREK" == false ]]; then

    echo ""
    echo "════════════════════════════════════════════════════════════════════════"
    echo "  Phase 1 — nf-core/sarek  (run: ${RUN_NAME})"
    echo "════════════════════════════════════════════════════════════════════════"

    cat > "${RUN_DIR}/run.config" <<EOF
params {
    run_dir = '${RUN_DIR}'
}
workDir = "${RUN_DIR}/work"
EOF

    nextflow run nf-core/sarek \
        -profile docker \
        -c "${RUN_DIR}/run.config" \
        -c "${SAREK_CONFIG}" \
        -params-file "${LAB_PARAMS}" \
        ${NF_RESUME} \
        --input  "${SAMPLESHEET}" \
        --outdir "${RUN_DIR}/results" \
        --tools  "haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff"

    echo ""
    echo "  Phase 1 complete → ${RUN_DIR}/results"
fi

# ── Verify sarek output exists ─────────────────────────────────────────────────
SAREK_OUTDIR="${RUN_DIR}/results"
if [[ ! -d "${SAREK_OUTDIR}" ]]; then
    echo "ERROR: sarek output not found: ${SAREK_OUTDIR}"
    echo "  Run Phase 1 first (remove --skip-sarek) or check the path."
    exit 1
fi

# ── Phase 2: nf-sarek-clinical (from GitHub) ──────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════════════"
echo "  Phase 2 — ${CLINICAL_REPO}  (run: ${RUN_NAME})"
echo "════════════════════════════════════════════════════════════════════════"

nextflow run "${CLINICAL_REPO}" \
    -profile docker \
    -work-dir "${RUN_DIR}/clinical_work" \
    -params-file "${LAB_PARAMS}" \
    ${NF_RESUME} \
    --sarek_outdir "${SAREK_OUTDIR}" \
    --run_name     "${RUN_NAME}" \
    --outdir       "${RUN_DIR}/clinical_results"

echo ""
echo "════════════════════════════════════════════════════════════════════════"
echo "  Pipeline complete!"
echo "  Clinical reports → ${RUN_DIR}/clinical_results/"
echo "════════════════════════════════════════════════════════════════════════"

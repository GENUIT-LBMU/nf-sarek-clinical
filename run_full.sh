#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# run_full.sh — Full WES clinical pipeline (nf-core/sarek → nf-sarek-clinical)
# ──────────────────────────────────────────────────────────────────────────────
# Usage:
#   ./run_full.sh --run RUN_NAME --samplesheet /path/to/samplesheet.csv
#   ./run_full.sh --run RUN_NAME --samplesheet /path/to/samplesheet.csv --resume
#   ./run_full.sh --run RUN_NAME --skip-sarek          # clinical report only
#   ./run_full.sh --run RUN_NAME --skip-sarek --resume # resume clinical report
#
# Arguments:
#   --run          Run identifier (e.g. 20260109_1229). Folder must exist at
#                  ${RUNS_BASE}/<RUN_NAME> or will be created.
#   --samplesheet  Path to nf-core/sarek CSV samplesheet (Phase 1; with --skip-sarek optional:
#                  its extra columns identify patient and specimens in the reports).
#   --resume       Pass -resume to Nextflow for both phases (use cached tasks).
#   --skip-sarek   Skip Phase 1; go straight to nf-sarek-clinical using the
#                  sarek results already present in <RUN_DIR>/results.
#
# Environment variables:
#   RUNS_BASE      Base directory for run outputs (default: $HOME/runs).
#   SAREK_REVISION nf-core/sarek release to run (default: 3.5.1).
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_BASE="${RUNS_BASE:-$HOME/runs}"
SAREK_REVISION="${SAREK_REVISION:-3.5.1}"   # pinned nf-core/sarek release (env-overridable)

# ── Parse arguments ───────────────────────────────────────────────────────────
RUN_NAME=""
SAMPLESHEET=""
RESUME=false
SKIP_SAREK=false

usage() {
    sed -n '3,18p' "${BASH_SOURCE[0]}"
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

# Use lab_params.yml (real paths) when present; fall back to full_pipeline.params.yml (template)
if [[ -f "${PIPELINE_DIR}/lab_params.yml" ]]; then
    PARAMS_FILE="${PIPELINE_DIR}/lab_params.yml"
else
    PARAMS_FILE="${PIPELINE_DIR}/full_pipeline.params.yml"
fi

# ── Phase 1: nf-core/sarek ────────────────────────────────────────────────────
if [[ "$SKIP_SAREK" == false ]]; then

    echo ""
    echo "════════════════════════════════════════════════════════════════════════"
    echo "  Phase 1 — nf-core/sarek  (run: ${RUN_NAME})"
    echo "════════════════════════════════════════════════════════════════════════"
    echo "  Params file: ${PARAMS_FILE}"

    # Per-run config: tells sarek where to put its work directory
    cat > "${RUN_DIR}/run.config" <<EOF
params {
    run_dir = '${RUN_DIR}'
}
workDir = "${RUN_DIR}/work"
EOF

    # NOTE: --tools is passed explicitly on the CLI (not only in conf/sarek.config).
    # sarek's modules/mutect2.config wraps the MUTECT2_PAIRED process config inside
    # `if (params.tools)`, evaluated at config load time before -c files are applied.
    # If tools were only in conf/sarek.config, MUTECT2_PAIRED would not get
    # --f1r2-tar-gz, causing LearnReadOrientationModel to produce no output and the
    # downstream join (failOnMismatch:true) to throw "Join mismatch" at FilterMutectCalls.
    # Launched from the run directory: a nextflow.config in the launch directory
    # is loaded too, and the one of this repository would reach sarek
    # (resource labels, manifest, parameters).
    SAMPLESHEET="$(realpath "${SAMPLESHEET}")"
    PARAMS_FILE="$(realpath "${PARAMS_FILE}")"
    cd "${RUN_DIR}"
    nextflow run nf-core/sarek \
        -r "${SAREK_REVISION}" \
        -profile docker \
        -c "${RUN_DIR}/run.config" \
        -c "${PIPELINE_DIR}/conf/sarek.config" \
        -params-file "${PARAMS_FILE}" \
        ${NF_RESUME} \
        --input    "${SAMPLESHEET}" \
        --outdir   "${RUN_DIR}/results" \
        --tools    "haplotypecaller,mutect2,strelka,manta,cnvkit,ascat,msisensorpro,snpeff"

    # Release recorded next to sarek's outputs (read by the reports, also with --skip-sarek)
    echo "${SAREK_REVISION}" > "${RUN_DIR}/results/pipeline_info/sarek_release.txt"
    cd "${PIPELINE_DIR}"

    echo ""
    echo "  Phase 1 complete → ${RUN_DIR}/results"
fi

# ── Verify sarek output exists before Phase 2 ─────────────────────────────────
SAREK_OUTDIR="${RUN_DIR}/results"
if [[ ! -d "${SAREK_OUTDIR}" ]]; then
    echo "ERROR: sarek output directory not found: ${SAREK_OUTDIR}"
    echo "  Run Phase 1 first (remove --skip-sarek), or check the path."
    exit 1
fi

# ── Phase 2: nf-sarek-clinical ────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════════════"
echo "  Phase 2 — nf-sarek-clinical  (run: ${RUN_NAME})"
echo "════════════════════════════════════════════════════════════════════════"

cd "${PIPELINE_DIR}"
echo "  Using params file: ${PARAMS_FILE}"

nextflow run main.nf \
    -profile docker \
    -work-dir "${RUN_DIR}/clinical_work" \
    -params-file "${PARAMS_FILE}" \
    ${NF_RESUME} \
    --sarek_outdir  "${SAREK_OUTDIR}" \
    ${SAMPLESHEET:+--samplesheet "$(realpath "${SAMPLESHEET}")"} \
    --run_name      "${RUN_NAME}" \
    --outdir        "${RUN_DIR}/clinical_results"

echo ""
echo "════════════════════════════════════════════════════════════════════════"
echo "  Pipeline complete!"
echo "  Clinical reports → ${RUN_DIR}/clinical_results/"
echo "════════════════════════════════════════════════════════════════════════"

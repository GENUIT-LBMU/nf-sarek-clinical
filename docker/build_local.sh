#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# build_local.sh — Build the 3 custom Docker images required by nf-sarek-clinical
# ──────────────────────────────────────────────────────────────────────────────
# Usage:
#   bash docker/build_local.sh                       # tags: nf-sarek-clinical/<name>:<ver>
#   NAMESPACE=ghcr.io/myorg bash docker/build_local.sh   # custom prefix
#
# After building, run the pipeline normally (default --docker_namespace matches
# the local tags). To use a published registry instead, override at run time:
#   nextflow run ... --docker_namespace ghcr.io/myorg/nf-sarek-clinical
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

NAMESPACE="${NAMESPACE:-nf-sarek-clinical}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Building custom images under namespace: ${NAMESPACE}"
echo ""

docker build -f "${SCRIPT_DIR}/Dockerfile.civicpy"      -t "${NAMESPACE}/civicpy:5.4.0"       "${SCRIPT_DIR}"
docker build -f "${SCRIPT_DIR}/Dockerfile.r_annotation" -t "${NAMESPACE}/r-annotation:1.0.0" "${SCRIPT_DIR}"
docker build -f "${SCRIPT_DIR}/Dockerfile.sigprofiler"  -t "${NAMESPACE}/sigprofiler:0.1.8"  "${SCRIPT_DIR}"

echo ""
echo "Done. Images:"
docker images --filter "reference=${NAMESPACE}/*" --format "  {{.Repository}}:{{.Tag}}  ({{.Size}})"
echo ""
echo "Pipeline default (--docker_namespace nf-sarek-clinical) matches these tags."

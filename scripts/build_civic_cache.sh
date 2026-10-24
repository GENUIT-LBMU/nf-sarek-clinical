#!/usr/bin/env bash
# =============================================================================
# build_civic_cache.sh — freeze a dated CIViC cache for CIVIC_ANNOTATE
# (--civic_cache).
#
# Output (in --outdir):
#   civicpy_cache_<YYYYMMDD>.pkl            CIViC nightly civicpy cache, frozen
#   civicpy_cache_<YYYYMMDD>.manifest.txt   source URL, CIViC snapshot time,
#                                           civicpy version, record counts, md5
#
# The pipeline never downloads CIViC at run time: the same cache gives the same
# annotation, and the report can state the CIViC snapshot date. The cache is a
# Python pickle tied to the civicpy version that reads it, so it is checked by
# loading it with the pipeline's civicpy image:
#
#   docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -v /data/refs:/data/refs \
#       -w $PWD nf-sarek-clinical/civicpy:5.4.0 \
#       bash scripts/build_civic_cache.sh --outdir /data/refs/civic
#
# CIViC data are CC0 (https://civicdb.org; Griffith et al. 2017, Nat Genet).
# =============================================================================
set -euo pipefail

OUTDIR=""
URL="https://civicdb.org/downloads/nightly/nightly-civicpy_cache.pkl"

usage() { sed -n '2,22p' "$0"; exit 1; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir) OUTDIR="$2"; shift 2 ;;
        --url)    URL="$2";    shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done
[[ -n "$OUTDIR" ]] || usage
for t in curl python3 civicpy; do command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done

mkdir -p "$OUTDIR"
TMP="$(mktemp "${OUTDIR}/.civic_XXXXXX.pkl")"
trap 'rm -f "$TMP"' EXIT
curl -sSfL -o "$TMP" "$URL"

# Load with civicpy (no staleness check, no network) and report what is inside
INFO="$(CIVICPY_CACHE_TIMEOUT_DAYS=36500 CIVICPY_REMOTE_CACHE_URL=file:///nonexistent \
    python3 - "$TMP" <<'PYEOF'
import sys
from collections import Counter
import civicpy.civic as civic
if not civic.load_cache(local_cache_path=sys.argv[1], on_stale='ignore'):
    sys.exit('ERROR: civicpy could not load the cache')
snap = civic.CACHE.get('full_cached')
kinds = Counter(type(v).__name__ for k, v in civic.CACHE.items() if not isinstance(k, str))
evid = [v for v in civic.CACHE.values() if type(v).__name__ == 'Evidence']
status = Counter(getattr(e, 'status', '?') for e in evid)
print(snap.strftime('%Y%m%d') if snap else 'unknown')
print(snap.isoformat() if snap else 'unknown')
print(', '.join(f'{k}={n}' for k, n in sorted(kinds.items())))
print(', '.join(f'{k}={n}' for k, n in sorted(status.items())))
PYEOF
)"
SNAP_DAY="$(sed -n 1p <<< "$INFO")"
[[ "$SNAP_DAY" =~ ^[0-9]{8}$ ]] || { echo "ERROR: no CIViC snapshot time in cache" >&2; exit 1; }
PREFIX="${OUTDIR}/civicpy_cache_${SNAP_DAY}"
mv "$TMP" "${PREFIX}.pkl"
trap - EXIT

{
    echo "resource:        CIViC nightly civicpy cache (frozen)"
    echo "source:          ${URL}"
    echo "downloaded:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "civic_snapshot:  $(sed -n 2p <<< "$INFO")"
    echo "civicpy:         $(civicpy --version | awk '{print $NF}')"
    echo "records:         $(sed -n 3p <<< "$INFO")"
    echo "evidence_status: $(sed -n 4p <<< "$INFO")"
    echo "md5:             $(md5sum "${PREFIX}.pkl" | cut -d' ' -f1)"
} > "${PREFIX}.manifest.txt"
cat "${PREFIX}.manifest.txt"

#!/usr/bin/env bash
# run_fit_sfa_tau.sh — fit per-type SFA time constant (tau_a) from raw aisynphys
# long-square sweeps, inside the dedicated venv (with the Zscaler CA for downloads).
#
# One-time deps (beyond the matrix-extraction venv), py3.11:
#   uv pip install --python scripts/aisynphys/.venv/bin/python ipfx h5py matplotlib
# (ipfx 2.x is used via its core SpikeFeatureExtractor; we do NOT import
#  aisynphys.intrinsic_ephys, which pulls the removed ipfx.chirp_features.)
#
# Downloads raw multipatch NWBs (~0.5-2 GB each) into ~/ai_synphys_cache; the
# default --max-experiments 25 keeps this to a few dozen files. Examples:
#   ./run_fit_sfa_tau.sh --list-only              # sampling plan, no downloads
#   ./run_fit_sfa_tau.sh --max-experiments 1      # smoke test on one NWB
#   ./run_fit_sfa_tau.sh                          # pilot (25 experiments)
#   ./run_fit_sfa_tau.sh --n-exp 2                # fast+slow tau (later)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$SCRIPT_DIR/.venv"
PY="$VENV/bin/python"

if [ ! -x "$PY" ]; then
    echo "[ERROR] venv python not found at $PY" >&2
    echo "        Create it and install deps (see the header of this script)." >&2
    exit 1
fi

# Trust the corporate (Zscaler) CA for the NWB / DB downloads, if present.
COMBINED="$VENV/ca-combined.pem"
if [ ! -f "$COMBINED" ] && [ -f "$HOME/zscaler.pem" ]; then
    BASE=""
    for c in /etc/ssl/cert.pem /opt/homebrew/etc/ca-certificates/cert.pem; do
        [ -f "$c" ] && BASE="$c" && break
    done
    { [ -n "$BASE" ] && cat "$BASE"; cat "$HOME/zscaler.pem"; } > "$COMBINED"
fi
if [ -f "$COMBINED" ]; then
    export SSL_CERT_FILE="$COMBINED" REQUESTS_CA_BUNDLE="$COMBINED" CURL_CA_BUNDLE="$COMBINED"
fi

exec "$PY" "$SCRIPT_DIR/fit_sfa_tau.py" "$@"

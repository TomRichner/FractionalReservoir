#!/usr/bin/env bash
# run_extract.sh — run the Campagnola matrix extraction inside its dedicated venv.
#
# Prereqs (one-time, see the header of extract_campagnola_matrices.py and the Phase-2 plan):
#   uv venv --python <pyenv 3.11> scripts/aisynphys/.venv
#   uv pip install --python scripts/aisynphys/.venv/bin/python "numpy<2" scipy pandas \
#       "sqlalchemy<2" statsmodels pyyaml
#   uv pip install --python scripts/aisynphys/.venv/bin/python --no-deps \
#       git+https://github.com/AllenInstitute/neuroanalysis
#   uv pip install --python scripts/aisynphys/.venv/bin/python --no-deps -e <path to aisynphys clone>
#
# NOTE: aisynphys needs Python <= 3.11 (it imports distutils, removed in 3.12).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$SCRIPT_DIR/.venv"
PY="$VENV/bin/python"

if [ ! -x "$PY" ]; then
    echo "[ERROR] venv python not found at $PY" >&2
    echo "        Create it first (see the header of this script / the Phase-2 plan)." >&2
    exit 1
fi

# Trust the corporate (Zscaler) CA for the one-time DB download, if present.
# Combined bundle = system roots + Zscaler; built at venv-creation time.
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

exec "$PY" "$SCRIPT_DIR/extract_campagnola_matrices.py" "$@"

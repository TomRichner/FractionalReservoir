# aisynphys extraction/fit env — how to (re)build the venv

The Python here (`extract_campagnola_matrices.py`, `fit_sfa_tau.py`) pulls Campagnola et al. 2022
data from the Allen `aisynphys` package into `src/connectivity/campagnola/`. It runs in a dedicated,
**gitignored** venv at `scripts/aisynphys/.venv`. This note records exactly how that venv was built,
since it's machine-local and not obvious to reproduce.

## Gotchas (why it's fiddly)

- **Python must be ≤ 3.11.** `aisynphys` imports `distutils`, removed in Python 3.12. We use
  **pyenv 3.11.7** (`~/.pyenv/versions/3.11.7`).
- **Corporate TLS (Zscaler).** Downloads (the ~73 MB SQLite DB, the raw NWBs, pip/git installs) fail
  SSL verification behind Zscaler. Fix: a **combined CA bundle** = system roots + `~/zscaler.pem`,
  exported as `SSL_CERT_FILE`/`REQUESTS_CA_BUNDLE`/`CURL_CA_BUNDLE`. `run_extract.sh` and
  `run_fit_sfa_tau.sh` build `.venv/ca-combined.pem` and export it automatically **if `~/zscaler.pem`
  exists** — so keep that file around. (For `pip`/`uv` installs, also set `PIP_CERT`/`UV_NATIVE_TLS=1`
  or run the install commands with `SSL_CERT_FILE` pointed at the combined bundle.)
- **`aisynphys` + `neuroanalysis` are installed `--no-deps`** (their pinned dep trees don't resolve on
  py3.11); we install their actual runtime deps explicitly.
- **`ipfx` is only needed for the SFA-tau fit** (`fit_sfa_tau.py`), not the matrix extraction. It
  pulls `pynwb`/`h5py`/`matplotlib` (no `allensdk` on modern ipfx). We use **ipfx core** directly
  (`ipfx.feature_extractor.SpikeFeatureExtractor`), *not* `aisynphys.intrinsic_ephys` — that module
  imports `ipfx.chirp_features`, which was removed in ipfx 2.x.

## Rebuild recipe (uv)

Prereqs: `uv` (Homebrew), pyenv 3.11.7 (`pyenv install 3.11.7`), the `aisynphys` repo cloned at
`/Users/richner.thomas/Desktop/local_code/aisynphys`, and `~/zscaler.pem` present.

```bash
cd <repo root>          # the FractionalReservoir checkout
PYENV311="$HOME/.pyenv/versions/3.11.7/bin/python"
VENV=scripts/aisynphys/.venv
AISYNPHYS=/Users/richner.thomas/Desktop/local_code/aisynphys   # local clone of AllenInstitute/aisynphys

# 1. venv on py3.11
uv venv --python "$PYENV311" "$VENV"

# 2. build the combined CA bundle + trust it for the installs (Zscaler)
BASE=/etc/ssl/cert.pem; [ -f /opt/homebrew/etc/ca-certificates/cert.pem ] && BASE=/opt/homebrew/etc/ca-certificates/cert.pem
cat "$BASE" ~/zscaler.pem > "$VENV/ca-combined.pem"
export SSL_CERT_FILE="$VENV/ca-combined.pem" REQUESTS_CA_BUNDLE="$SSL_CERT_FILE" \
       CURL_CA_BUNDLE="$SSL_CERT_FILE" PIP_CERT="$SSL_CERT_FILE" UV_NATIVE_TLS=1

# 3. runtime deps for the DB extraction (matrix pipeline)
uv pip install --python "$VENV/bin/python" "numpy<2" scipy pandas "sqlalchemy<2" statsmodels pyyaml
uv pip install --python "$VENV/bin/python" --no-deps git+https://github.com/AllenInstitute/neuroanalysis
uv pip install --python "$VENV/bin/python" --no-deps -e "$AISYNPHYS"

# 4. EXTRA deps for the raw-sweep SFA-tau fit (fit_sfa_tau.py) — skip if you only need the matrices
uv pip install --python "$VENV/bin/python" ipfx h5py matplotlib
```

### Known-good versions (as built 2026-07)

py 3.11.7 · numpy 1.26.4 · scipy 1.17.1 · pandas 3.0.3 · SQLAlchemy 1.4.54 · statsmodels 0.14.6 ·
neuroanalysis 0.0.7 · aisynphys 0.0.1 · ipfx 2.1.2 · pynwb 3.1.2 · h5py 3.16.0.

## Data cache (NOT in the venv — survives venv/worktree deletion)

`SynphysDatabase.load_current('small')` and raw NWB downloads cache under **`~/ai_synphys_cache/`**
(the SQLite DB + `raw_data_files/<ext_id>/data.nwb`, ~GBs). This lives in `$HOME`, so rebuilding the
venv does **not** re-download it.

## Run

```bash
scripts/aisynphys/run_extract.sh                 # 4x4 Campagnola matrices -> src/connectivity/campagnola/
scripts/aisynphys/run_fit_sfa_tau.sh --list-only # SFA-tau fit: sampling plan (no download)
scripts/aisynphys/run_fit_sfa_tau.sh             # SFA-tau fit pilot (downloads ~dozens of NWBs)
```

Both wrappers auto-load the venv + Zscaler CA. See each script's header for options.

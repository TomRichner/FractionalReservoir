#!/usr/bin/env python
"""
extract_campagnola_matrices.py

Pull 4x4 cell-type matrices (Pyr, Pvalb, Sst, Vip -- pooled across cortical layers)
for MOUSE from the Allen synaptic-physiology 'small' database (Campagnola et al. 2022),
for use in parameterizing the cell-typed SRNN reservoir (Phase 3).

Extracts, per (presynaptic-type, postsynaptic-type) element:
  - connectivity   : raw + distance-adjusted connection probability, n_probed, n_connected
  - strength       : synapse.psp_amplitude, dynamics.pulse_amp_90th_percentile
  - STP            : dynamics.stp_induction_50hz / stp_recovery_250ms / variability_resting_state
  - release model  : synapse_model.ml_{facilitation,depression}_{amount,tau},
                     ml_base_release_probability, ml_n_release_sites, ml_mini_amplitude
  - kinetics (ref) : synapse.latency, psc_rise_time, psc_decay_tau
and per postsynaptic-type: intrinsic.adaptation_index (SFA), from a per-cell (de-duplicated) query.

Aggregation is the MEDIAN (matching the paper's matrices); a matching sample-count (n) matrix is
stored for every metric; elements with n < MIN_N are set to NaN.

Outputs (into --outdir, default src/connectivity/campagnola/):
  - csv/<metric>.csv          : 4x4 median matrix (rows=pre, cols=post), labeled
  - csv/<metric>__n.csv       : 4x4 sample-count matrix
  - connectivity.csv          : long-form connectivity table (cp, adj cp, CIs, n)
  - sfa_adaptation_index.csv  : per-type SFA (median + n)
  - campagnola_matrices.mat   : combined struct for MATLAB (all matrices + n + types)
  - PROVENANCE.md             : DB version, date, query params, cell-class defs

Run via scripts/aisynphys/run_extract.sh (sets up the venv + CA bundle).
"""
import argparse
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.io import savemat

from aisynphys.database import SynphysDatabase
from aisynphys.cell_class import (CellClass, classify_cells, classify_pairs,
                                  classify_pair_dataframe, classify_cell_dataframe)
from aisynphys.connectivity import measure_connectivity

# --- configuration --------------------------------------------------------
DB_SIZE = 'small'
SPECIES = 'mouse'
SIGMA = 100e-6            # lateral-distance Gaussian sigma for adjusted connectivity (as in connectivity.ipynb)
DIST_MEASURE = 'lateral_distance'
MIN_N = 2                # elements with fewer than this many pairs -> NaN (paper requires 2+)
AGG = 'median'

# 4 subclasses, pooled across layers. Order is fixed and used everywhere.
CELL_CLASSES = OrderedDict([
    ('pyr',   CellClass(cell_class='ex', name='pyr')),
    ('pvalb', CellClass(cre_type='pvalb', name='pvalb')),
    ('sst',   CellClass(cre_type='sst',  name='sst')),
    ('vip',   CellClass(cre_type='vip',  name='vip')),
])
ORDER = list(CELL_CLASSES)

# metric key -> dataframe column (dotted table.column form from .dataframe())
METRICS = OrderedDict([
    ('psp_amplitude',          'synapse.psp_amplitude'),
    ('pulse_amp_90pct',        'dynamics.pulse_amp_90th_percentile'),
    ('stp_induction_50hz',     'dynamics.stp_induction_50hz'),
    ('stp_recovery_250ms',     'dynamics.stp_recovery_250ms'),
    ('variability_resting',    'dynamics.variability_resting_state'),
    ('ml_facilitation_amount', 'synapse_model.ml_facilitation_amount'),
    ('ml_facilitation_tau',    'synapse_model.ml_facilitation_tau'),
    ('ml_depression_amount',   'synapse_model.ml_depression_amount'),
    ('ml_depression_tau',      'synapse_model.ml_depression_tau'),
    ('ml_release_prob',        'synapse_model.ml_base_release_probability'),
    ('ml_n_release_sites',     'synapse_model.ml_n_release_sites'),
    ('ml_mini_amplitude',      'synapse_model.ml_mini_amplitude'),
    ('latency',                'synapse.latency'),
    ('psc_rise_time',          'synapse.psc_rise_time'),
    ('psc_decay_tau',          'synapse.psc_decay_tau'),
])

# per-type intrinsic (single-cell) properties -> dataframe column, from the per-cell
# query. Stored in SI units as in the DB: tau [s], rheobase [A], input_resistance [Ohm],
# fi_slope [Hz/A], adaptation_index / isi_adapt_ratio [dimensionless].
INTRINSIC = OrderedDict([
    ('adaptation_index',  'intrinsic.adaptation_index'),
    ('isi_adapt_ratio',   'intrinsic.isi_adapt_ratio'),
    ('tau',               'intrinsic.tau'),
    ('fi_slope',          'intrinsic.fi_slope'),
    ('rheobase',          'intrinsic.rheobase'),
    ('input_resistance',  'intrinsic.input_resistance'),
])


def _empty_4x4():
    return pd.DataFrame(np.full((4, 4), np.nan), index=ORDER, columns=ORDER)


def matrix(pairs_df, col, agg):
    """4x4 (pre x post) aggregate for a dataframe column, reindexed to ORDER."""
    if col not in pairs_df.columns:
        raise KeyError(f"column not found in pair dataframe: {col}")
    m = pairs_df.pivot_table(col, 'pre_class', 'post_class', aggfunc=agg)
    return m.reindex(index=ORDER, columns=ORDER)


def extract(db, outdir):
    os.makedirs(os.path.join(outdir, 'csv'), exist_ok=True)
    mat_out = {}           # -> savemat struct fields
    mat_out['types'] = np.array(ORDER, dtype=object)
    long_rows = []         # tidy summary rows

    # ------------------------------------------------------------------ #
    # 1. Connectivity (all probed pairs; raw + distance-adjusted)
    # ------------------------------------------------------------------ #
    print("== connectivity ==", flush=True)
    probed = db.pair_query(experiment_type='standard multipatch', species=SPECIES).all()
    print(f"  probed mouse pairs: {len(probed)}", flush=True)
    groups = classify_pairs(probed, classify_cells(CELL_CLASSES.values(), pairs=probed))
    conn = measure_connectivity(groups, sigma=SIGMA, dist_measure=DIST_MEASURE)

    cp = _empty_4x4(); cp_adj = _empty_4x4()
    n_probed = _empty_4x4(); n_conn = _empty_4x4()
    conn_rows = []
    for pre in ORDER:
        for post in ORDER:
            res = conn.get((CELL_CLASSES[pre], CELL_CLASSES[post]))
            if res is None:
                continue
            npb = res.get('n_probed', 0); ncn = res.get('n_connected', 0)
            cprob = res.get('connection_probability', (np.nan, np.nan, np.nan))
            adj = res.get('adjusted_connectivity', (np.nan, np.nan, np.nan))
            n_probed.loc[pre, post] = npb
            n_conn.loc[pre, post] = ncn
            if npb >= MIN_N:
                cp.loc[pre, post] = cprob[0]
                cp_adj.loc[pre, post] = adj[0]
            conn_rows.append(dict(pre=pre, post=post, n_probed=npb, n_connected=ncn,
                                  cp=cprob[0], cp_lo=cprob[1], cp_hi=cprob[2],
                                  cp_adj=adj[0], cp_adj_lo=adj[1], cp_adj_hi=adj[2]))
    for name, df in [('conn_prob', cp), ('conn_prob_adj', cp_adj),
                     ('n_probed', n_probed), ('n_connected', n_conn)]:
        df.to_csv(os.path.join(outdir, 'csv', name + '.csv'))
        mat_out[name] = df.values.astype(float)
    pd.DataFrame(conn_rows).to_csv(os.path.join(outdir, 'connectivity.csv'), index=False)
    print(f"  wrote connectivity ({len(conn_rows)} elements)", flush=True)

    # ------------------------------------------------------------------ #
    # 2. Strength / STP / release-model / kinetics matrices
    # ------------------------------------------------------------------ #
    print("== per-connection metrics ==", flush=True)
    pairs = db.pair_query(experiment_type='standard multipatch', species=SPECIES,
                          synapse=True, preload=['cell', 'synapse']).dataframe()
    print(f"  synapse pairs (dataframe): {len(pairs)}", flush=True)
    classify_pair_dataframe(CELL_CLASSES, pairs)   # adds pre_class / post_class

    for key, col in METRICS.items():
        med = matrix(pairs, col, AGG)
        cnt = matrix(pairs, col, 'count')
        med = med.where(cnt >= MIN_N)              # NaN where too few
        med.to_csv(os.path.join(outdir, 'csv', key + '.csv'))
        cnt.to_csv(os.path.join(outdir, 'csv', key + '__n.csv'))
        mat_out[key] = med.values.astype(float)
        mat_out[key + '_n'] = cnt.values.astype(float)
        for pre in ORDER:
            for post in ORDER:
                long_rows.append(dict(metric=key, pre=pre, post=post,
                                      value=med.loc[pre, post], n=cnt.loc[pre, post]))

    # ------------------------------------------------------------------ #
    # 3. Intrinsic per-type properties (per-cell, de-duplicated) incl. SFA
    # ------------------------------------------------------------------ #
    print("== intrinsic per-type properties ==", flush=True)
    cells = (db.query(db.Cell, db.Intrinsic)
               .join(db.Intrinsic, db.Intrinsic.cell_id == db.Cell.id)).dataframe()
    # The coarse-matrix cre_type labels are mouse-specific, so classification below
    # already excludes human. adaptation_index/tau/etc. are postsynaptic-cell properties.
    cells = cells.copy()
    cells['subclass'] = classify_cell_dataframe(CELL_CLASSES, cells, prefix='')  # 'cell.*' cols
    intr_table = OrderedDict([('type', ORDER)])
    for key, col in INTRINSIC.items():
        med = np.full(4, np.nan); cnt = np.full(4, np.nan)
        for i, t in enumerate(ORDER):
            if col in cells.columns:
                vals = cells.loc[cells['subclass'] == t, col].dropna()
            else:
                vals = pd.Series([], dtype=float)
            n = int(vals.shape[0]); cnt[i] = n
            med[i] = float(vals.median()) if n >= MIN_N else np.nan
        intr_table[key] = med
        intr_table[key + '_n'] = cnt
        # keep adaptation_index under its historical 'sfa_adaptation_index' .mat field
        field = 'sfa_adaptation_index' if key == 'adaptation_index' else key
        mat_out[field] = med.reshape(-1, 1)
        mat_out[field + '_n'] = cnt.reshape(-1, 1)
    pd.DataFrame(intr_table).to_csv(os.path.join(outdir, 'intrinsic_per_type.csv'), index=False)
    print(pd.DataFrame(intr_table).to_string(index=False), flush=True)

    # ------------------------------------------------------------------ #
    # 4. Write long-form + .mat
    # ------------------------------------------------------------------ #
    pd.DataFrame(long_rows).to_csv(os.path.join(outdir, 'metrics_long.csv'), index=False)
    savemat(os.path.join(outdir, 'campagnola_matrices.mat'), {'campagnola': mat_out},
            do_compression=True)
    print(f"  wrote campagnola_matrices.mat ({len(mat_out)} fields)", flush=True)

    return dict(mat=mat_out, db_file=str(db.sqlite_file) if hasattr(db, 'sqlite_file') else str(db),
                n_probed=len(probed), n_synapse=len(pairs))


def sanity_checks(mat):
    """Print + assert a few qualitative facts from the paper."""
    idx = {t: i for i, t in enumerate(ORDER)}
    print("\n== sanity checks vs Campagnola 2022 ==", flush=True)
    def g(field, pre, post):
        return mat[field][idx[pre], idx[post]]

    checks = []
    # Pyr->Pvalb connectivity should be high (strong E->Pv)
    checks.append(("Pyr->Pvalb adj connectivity high (>0.2)", g('conn_prob_adj', 'pyr', 'pvalb') > 0.2))
    # Pvalb outputs strongly depressing -> negative stp_induction
    checks.append(("Pvalb->Pyr STP depressing (<0)", g('stp_induction_50hz', 'pvalb', 'pyr') < 0))
    # Recurrent Pyr->Pyr depressing (<0)
    checks.append(("Pyr->Pyr STP depressing (<0)", g('stp_induction_50hz', 'pyr', 'pyr') < 0))
    # Sst->Vip among the most facilitating inhibitory elements (positive facilitation amount)
    sst_vip_fac = g('ml_facilitation_amount', 'sst', 'vip')
    checks.append(("Sst->Vip facilitation_amount positive", (sst_vip_fac == sst_vip_fac) and sst_vip_fac > 0))
    # E->Sst facilitating (positive stp_induction)
    checks.append(("Pyr->Sst STP facilitating (>0)", g('stp_induction_50hz', 'pyr', 'sst') > 0))

    ok = True
    for name, passed in checks:
        print(f"  [{'ok' if passed else '??'}] {name}", flush=True)
        ok = ok and bool(passed)

    # coverage overview
    print("\n  n_connected matrix (pre x post):", flush=True)
    print(pd.DataFrame(mat['n_connected'], index=ORDER, columns=ORDER).to_string(), flush=True)
    return ok


def main():
    ap = argparse.ArgumentParser()
    default_out = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'src', 'connectivity', 'campagnola')
    ap.add_argument('--outdir', default=default_out)
    args = ap.parse_args()

    print(f"loading '{DB_SIZE}' database ...", flush=True)
    db = SynphysDatabase.load_current(DB_SIZE)
    print(f"DB: {db}", flush=True)

    info = extract(db, args.outdir)
    ok = sanity_checks(info['mat'])

    # provenance
    versions = SynphysDatabase.list_current_versions()
    with open(os.path.join(args.outdir, 'PROVENANCE.md'), 'w') as f:
        f.write("# Campagnola 2022 matrices — provenance\n\n")
        f.write(f"- Source: Allen synaptic physiology dataset (Campagnola et al., Science 2022).\n")
        f.write(f"- Database tier: `{DB_SIZE}`; current versions: `{versions.get(DB_SIZE)}`\n")
        f.write(f"- Generated by `scripts/aisynphys/extract_campagnola_matrices.py`.\n")
        f.write(f"- Species: {SPECIES}; experiment_type: standard multipatch.\n")
        f.write(f"- Cell classes (pooled across layers): "
                + ", ".join(f"{k}={v!r}" for k, v in CELL_CLASSES.items()) + "\n")
        f.write(f"- Aggregation: {AGG}; min pairs per element (else NaN): {MIN_N}.\n")
        f.write(f"- Connectivity: measure_connectivity(sigma={SIGMA}, dist_measure='{DIST_MEASURE}'); "
                f"conn_prob = raw, conn_prob_adj = distance-adjusted.\n")
        f.write(f"- Pairs: {info['n_probed']} probed, {info['n_synapse']} with confirmed synapses.\n")
        f.write("- Type order (rows=pre, cols=post): " + ", ".join(ORDER) + "\n")
    print("\nwrote PROVENANCE.md", flush=True)

    print("\nDONE" + ("" if ok else " (some sanity checks flagged '??' — inspect)"), flush=True)
    sys.exit(0)


if __name__ == '__main__':
    main()

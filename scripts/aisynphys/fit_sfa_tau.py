#!/usr/bin/env python
"""
fit_sfa_tau.py

Fit a spike-frequency-adaptation TIME CONSTANT (tau_a) per cortical cell type
(Pyr, Pvalb, Sst, Vip) from the RAW long-square current-clamp sweeps in the Allen
synaptic-physiology dataset (Campagnola et al. 2022). The aisynphys DB only stores
dimensionless adaptation magnitude (adaptation_index / isi_adapt_ratio) and the
passive membrane tau -- NOT an adaptation time constant -- so we go back to the
sweeps IPFX used to compute adaptation_index and fit the firing-rate decay ourselves.

Method (per cell):
  cell.experiment -> expt.data (downloads NWB) ->
  get_intrinsic_recording_dict(data, cell.electrode.device_id)['LP'] ->
  for each suprathreshold long-pulse sweep: MPSweep(rec, -pulse_start) so the step
  begins at t=0, IPFX SpikeFeatureExtractor -> spike times -> instantaneous rate
  f_i = 1/ISI_i at ISI midpoints -> fit  f(t) = f_inf + (f0-f_inf) exp(-t/tau_a).
Aggregation: per-cell median tau over its passing sweeps; per-type median + IQR.

Only ~25 experiments are sampled by default (raw NWBs are ~0.5-2 GB each). Cells
share experiments, so a few dozen NWBs yield many cells across all four types.

Deps beyond the matrix extractor's venv: ipfx, h5py, matplotlib (see run_fit_sfa_tau.sh).
We deliberately do NOT import aisynphys.intrinsic_ephys (it imports ipfx.chirp_features,
gone in ipfx 2.x); we reimplement its small MPSweep adapter and use ipfx core directly.

Outputs (into --outdir, default src/connectivity/campagnola/):
  - sfa_tau_per_type.csv   : committed per-type summary (source of truth for MATLAB)
  - sfa_tau_fits.mat       : per-cell fit table + per-type summary + example curves
  - figures/sfa_tau_fits.png : per-type example exponential fits + tau_a distribution
  - sfa_tau_PROVENANCE.md  : DB version, sample, fit model, QC thresholds

Run via scripts/aisynphys/run_fit_sfa_tau.sh (venv + CA bundle for the downloads).
"""
import argparse
import gc
import os
import sys
import traceback
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.io import savemat

from ipfx.feature_extractor import SpikeFeatureExtractor
from ipfx.sweep import Sweep

from aisynphys.database import SynphysDatabase
from aisynphys.cell_class import CellClass

# --- configuration --------------------------------------------------------
DB_SIZE = 'small'
SPECIES = 'mouse'

# 4 subclasses, pooled across layers -- SAME definitions as extract_campagnola_matrices.py
CELL_CLASSES = OrderedDict([
    ('pyr',   CellClass(cell_class='ex', name='pyr')),
    ('pvalb', CellClass(cre_type='pvalb', name='pvalb')),
    ('sst',   CellClass(cre_type='sst',  name='sst')),
    ('vip',   CellClass(cre_type='vip',  name='vip')),
])
ORDER = list(CELL_CLASSES)

# fit / QC thresholds
MIN_SPIKES = 5          # need >=5 spikes (>=4 ISIs, >=3 fit points for 3 params)
R2_MIN = 0.5            # minimum coefficient of determination for the exp fit
MIN_ADAPT = 0.05       # minimum fractional adaptation (f0-f_inf)/f0 to count as adapting
TAU_BOUNDS = (0.005, 5.0)   # seconds
N_EXAMPLES = 6         # example fits stored/plotted per type


# ------------------------------------------------------------------ #
# Build an IPFX Sweep from a raw recording WITHOUT touching rec.stimulus.items:
# in this ipfx/neuroanalysis combo that property triggers a test-pulse parser that
# raises "Found multiple square pulse in command waveform". So we derive the long
# step timing + holding baseline directly from the piecewise-constant command wave.
# ------------------------------------------------------------------ #
def _pulse_from_command(rec):
    """Return (start, end, holding) of the long square step from the command wave."""
    cmd = rec['command']
    y = np.asarray(cmd.data, float)
    t = np.asarray(cmd.time_values, float)
    I0 = float(np.median(y))                       # holding current (baseline)
    ptp = float(np.ptp(y))
    if ptp <= 0:
        return None
    edges = np.flatnonzero(np.abs(np.diff(y)) > 0.05 * ptp)
    bounds = np.concatenate(([0], edges + 1, [len(y)]))
    best, best_dur = None, 0.0
    for a, b in zip(bounds[:-1], bounds[1:]):       # each segment = one constant level
        if b - a < 2 or abs(np.median(y[a:b]) - I0) < 0.05 * ptp:
            continue                                # too short, or a baseline segment
        dur = t[b - 1] - t[a]
        if dur > best_dur:                          # longest non-baseline level = the step
            best_dur, best = dur, (float(t[a]), float(t[b - 1]))
    if best is None:
        return None
    return best[0], best[1], I0


def build_sweep(rec, k):
    """IPFX Sweep with the step shifted to t=0; returns (sweep, dur) or None."""
    pc = _pulse_from_command(rec)
    if pc is None:
        return None
    start, end, I0 = pc
    dur = end - start
    if not (dur and dur > 0.1):
        return None
    pri = rec['primary']
    cmd = rec['command']
    t = np.asarray(pri.time_values, float) - start          # step begins at t=0
    v = np.asarray(pri.data, float) * 1e3                    # V -> mV
    i = (np.asarray(cmd.data, float) - I0) * 1e12            # A -> pA, holding removed
    ok = (v != 0) & ~np.isnan(v)
    sweep = Sweep(t[ok], v[ok], i[ok], "CurrentClamp", pri.sample_rate, sweep_number=k)
    return sweep, dur


# ------------------------------------------------------------------ #
# adaptation models (single now; double coded for later --n-exp 2)
# ------------------------------------------------------------------ #
def model_single(t, A, f_inf, tau):
    return f_inf + A * np.exp(-t / tau)


def model_double(t, A1, tau1, A2, tau2, f_inf):
    return f_inf + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)


def _r2(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def fit_adaptation(t, f, n_exp=1):
    """Fit instantaneous rate f(t) with a 1- or 2-exponential decay to a floor.

    Returns a dict with the model params, derived tau(s), f0, f_inf, r2 -- or None
    if the fit fails. Extensible: n_exp=2 gives fast/slow taus for the model's n_a=2.
    """
    t = np.asarray(t, float); f = np.asarray(f, float)
    fmax = max(f.max(), 1e-6)
    if n_exp == 1:
        p0 = [max(f[0] - f[-1], 1e-3), max(f[-1], 0.0), 0.1]
        lo = [0.0, 0.0, TAU_BOUNDS[0]]
        hi = [5 * fmax, fmax, TAU_BOUNDS[1]]
        try:
            popt, _ = curve_fit(model_single, t, f, p0=p0, bounds=(lo, hi), maxfev=20000)
        except Exception:
            return None
        A, f_inf, tau = popt
        r2 = _r2(f, model_single(t, *popt))
        return dict(model='single', tau_a=tau, taus=[tau], f0=f_inf + A, f_inf=f_inf,
                    A=A, r2=r2, params=popt.tolist())
    else:
        p0 = [max(f[0] - f[-1], 1e-3) / 2, 0.03, max(f[0] - f[-1], 1e-3) / 2, 0.5, max(f[-1], 0.0)]
        lo = [0.0, TAU_BOUNDS[0], 0.0, TAU_BOUNDS[0], 0.0]
        hi = [5 * fmax, TAU_BOUNDS[1], 5 * fmax, TAU_BOUNDS[1], fmax]
        try:
            popt, _ = curve_fit(model_double, t, f, p0=p0, bounds=(lo, hi), maxfev=40000)
        except Exception:
            return None
        A1, tau1, A2, tau2, f_inf = popt
        taus = sorted([tau1, tau2])
        r2 = _r2(f, model_double(t, *popt))
        return dict(model='double', tau_a=taus[0], tau_fast=taus[0], tau_slow=taus[1],
                    taus=taus, f0=f_inf + A1 + A2, f_inf=f_inf, r2=r2, params=popt.tolist())


def qc_fit(fit):
    """Accept only converged, adapting, well-fit, non-railed fits."""
    if fit is None or not np.isfinite(fit['r2']) or fit['r2'] < R2_MIN:
        return False
    if fit['f0'] <= 0 or (fit['f0'] - fit['f_inf']) / fit['f0'] < MIN_ADAPT:
        return False
    lo, hi = TAU_BOUNDS
    for tau in fit['taus']:
        if not (lo * 1.02 < tau < hi * 0.98):   # railed against a bound -> reject
            return False
    return True


# ------------------------------------------------------------------ #
# per-cell fitting
# ------------------------------------------------------------------ #
def fit_cell(cell, n_exp, rheobase_pA=np.nan):
    """Fit every qualifying suprathreshold long-pulse sweep for one cell.

    Returns a list of per-sweep fit dicts (passing QC). Reuses aisynphys' LP-sweep
    selection + pulse-time parsing; spike times come from IPFX directly.
    """
    from aisynphys.nwb_recordings import get_intrinsic_recording_dict
    out = []
    expt = cell.experiment
    data = expt.data                         # triggers NWB download/cache
    if data is None:
        return out
    dev_id = cell.electrode.device_id
    lp = get_intrinsic_recording_dict(data, dev_id).get('LP', [])
    for k, rec in enumerate(lp):
        try:
            bs = build_sweep(rec, k)
            if bs is None:
                continue
            sweep, dur = bs
            spx = SpikeFeatureExtractor(start=0.0, end=dur)
            sf = spx.process(sweep.t, sweep.v, sweep.i)
            if sf is None or len(sf) < MIN_SPIKES:
                continue
            st = np.asarray(sf['threshold_t'].values, float)
            st = st[(st >= 0) & (st <= dur)]
            if len(st) < MIN_SPIKES:
                continue
            isis = np.diff(st)
            good = isis > 0
            tmid = 0.5 * (st[:-1] + st[1:])[good]
            f = 1.0 / isis[good]
            if len(f) < 3:
                continue
            fit = fit_adaptation(tmid, f, n_exp=n_exp)
            if not qc_fit(fit):
                continue
            if any(tau > 3 * dur for tau in fit['taus']):
                continue                       # tau longer than a short step can constrain
            step = (sweep.t >= 0) & (sweep.t <= dur)
            stim_amp = float(np.median(sweep.i[step]))
            fit.update(dict(n_spikes=int(len(st)), stim_amp_pA=stim_amp,
                            rel_amp_pA=stim_amp - rheobase_pA,
                            t=tmid.tolist(), f=f.tolist(), pulse_dur=float(dur)))
            out.append(fit)
        except Exception:
            continue
    return out


# ------------------------------------------------------------------ #
# sampling
# ------------------------------------------------------------------ #
def load_typed_cells(db):
    """Query mouse cells with intrinsic records; classify each into a subclass."""
    # adaptation_index != None restricts to cells that actually ran the long-square
    # suprathreshold protocol (If_Curve sweeps) and had features extracted -- so the
    # downloadable NWB is guaranteed to contain the sweeps we fit. (A plain Intrinsic
    # join also returns cells with all-null intrinsic rows, whose NWBs hold only
    # synaptic pulse-train stimuli -- no long-square sweeps.)
    q = (db.query(db.Cell, db.Intrinsic, db.Experiment)
           .join(db.Intrinsic, db.Intrinsic.cell_id == db.Cell.id)
           .join(db.Experiment, db.Cell.experiment_id == db.Experiment.id)
           .join(db.Slice, db.Experiment.slice_id == db.Slice.id)
           .filter(db.Slice.species == SPECIES)
           .filter(db.Intrinsic.adaptation_index != None))  # noqa: E711
    rows = q.all()
    cells = []
    for cell, intr, expt in rows:
        t = None
        for name, cls in CELL_CLASSES.items():
            try:
                if cell in cls:
                    t = name
                    break
            except Exception:
                continue
        if t is None:
            continue
        rheo = intr.rheobase
        cells.append(dict(cell=cell, type=t, expt_ext_id=expt.ext_id,
                          rheobase_pA=(rheo * 1e12 if rheo is not None else np.nan),
                          adaptation_index=intr.adaptation_index))
    return cells


def load_dataset(expt, cached_only):
    """Return the experiment's MultiPatchDataset. With cached_only, load ONLY from a
    fully-present local NWB (never trigger a download) and return None if it's missing
    or truncated -- so a preview run touches no network."""
    if not cached_only:
        return expt.data
    import aisynphys.config as config
    p = os.path.join(config.cache_path, 'raw_data_files', expt.ext_id, 'data.nwb')
    if not os.path.exists(p):
        return None
    try:
        from aisynphys.data import MultiPatchDataset
        ds = MultiPatchDataset(p)
        _ = ds.contents                     # force parse; raises on a truncated file
        expt._data = ds                     # so Experiment.data returns it (no download)
        return ds
    except Exception:
        return None


def free_experiment(expt):
    """Release an experiment's loaded NWB (else ~1-2 GB/experiment accumulates -> OOM)."""
    try:
        d = getattr(expt, '_data', None)
        if d is not None:
            close = getattr(d, 'close', None)
            if callable(close):
                try:
                    close()
                except Exception:
                    pass
            expt._data = None
    except Exception:
        pass
    gc.collect()


def select_experiments(cells, max_experiments):
    """Greedily pick experiments to BALANCE per-type coverage within the budget.

    Ranking by raw cell count over-samples pyr-heavy experiments and starves the
    rare types (vip/pvalb). Instead, at each step pick the experiment whose cells
    score highest under weight 1/(1+running_count[type]) -- so once a type is well
    covered its pull drops and experiments rich in the currently-rarest types win.
    Deterministic (alphabetical tie-break).
    """
    by_expt = defaultdict(list)
    for c in cells:
        by_expt[c['expt_ext_id']].append(c)
    if max_experiments <= 0:
        return OrderedDict(sorted(by_expt.items()))
    remaining = dict(by_expt)
    counts = {t: 0 for t in ORDER}
    chosen = OrderedDict()
    for _ in range(min(max_experiments, len(remaining))):
        best, best_score = None, -1.0
        for ext_id in sorted(remaining):
            tc = defaultdict(int)
            for c in remaining[ext_id]:
                tc[c['type']] += 1
            score = sum(n / (1.0 + counts[t]) for t, n in tc.items())
            if score > best_score:
                best_score, best = score, ext_id
        chosen[best] = remaining.pop(best)
        for c in chosen[best]:
            counts[c['type']] += 1
    return chosen


# ------------------------------------------------------------------ #
# aggregation + outputs
# ------------------------------------------------------------------ #
def aggregate(cell_records, n_exp, adapt_index_by_type):
    """cell_records: list of dicts {cell_id, expt, type, cell_tau, sweeps:[...]}.
    adapt_index_by_type: {type: [Campagnola adaptation_index of every sampled cell]}.

    Pure data reduction: report the fitted tau_a (median + IQR) for EVERY type, plus the
    median adaptation_index as context. This makes no modeling decision -- e.g. whether a
    weakly-adapting type (Pvalb, index ~0.002, wide tau_a spread) should use n_a=0 SFA
    timescales is left to the MATLAB Campagnola->parameter mapping layer."""
    per_type = OrderedDict()
    for t in ORDER:
        taus = np.array([r['cell_tau'] for r in cell_records
                         if r['type'] == t and np.isfinite(r['cell_tau'])])
        n_sw = sum(len(r['sweeps']) for r in cell_records if r['type'] == t)
        n_exp_t = len(set(r['expt'] for r in cell_records if r['type'] == t))
        aidx = np.array([a for a in adapt_index_by_type.get(t, []) if a is not None], float)
        aidx_med = float(np.median(aidx)) if len(aidx) else np.nan
        if len(taus):
            tau_med = float(np.median(taus))
            tau_p25 = float(np.percentile(taus, 25))
            tau_p75 = float(np.percentile(taus, 75))
        else:
            tau_med = tau_p25 = tau_p75 = np.nan
        per_type[t] = dict(adaptation_index_median=aidx_med,
                           tau_a_median_s=tau_med, tau_a_p25_s=tau_p25, tau_a_p75_s=tau_p75,
                           n_cells_fit=int(len(taus)), n_sweeps_fit=int(n_sw),
                           n_experiments=int(n_exp_t))
    return per_type


def write_outputs(outdir, cell_records, per_type, n_exp, meta):
    os.makedirs(outdir, exist_ok=True)
    model = 'single' if n_exp == 1 else 'double'

    # --- per-type CSV (committed source of truth) ---
    cols = ['type', 'model', 'tau_a_median_s', 'tau_a_p25_s', 'tau_a_p75_s',
            'adaptation_index_median', 'n_cells_fit', 'n_sweeps_fit', 'n_experiments']
    rows = [dict(type=t, model=model, **per_type[t]) for t in ORDER]
    csv_path = os.path.join(outdir, 'sfa_tau_per_type.csv')
    pd.DataFrame(rows)[cols].to_csv(csv_path, index=False)

    # --- per-cell table + example curves -> .mat ---
    cell_tbl = dict(
        cell_id=np.array([r['cell_id'] for r in cell_records], float),
        type=np.array([r['type'] for r in cell_records], dtype=object),
        experiment=np.array([r['expt'] for r in cell_records], dtype=object),
        tau_a=np.array([r['cell_tau'] for r in cell_records], float),
        n_sweeps=np.array([len(r['sweeps']) for r in cell_records], float),
    )
    examples = {}
    for t in ORDER:
        recs = [r for r in cell_records if r['type'] == t and np.isfinite(r['cell_tau'])]
        recs.sort(key=lambda r: r['cell_tau'])
        picks = recs[:: max(1, len(recs) // N_EXAMPLES)][:N_EXAMPLES] if recs else []
        ex = []
        for r in picks:
            sw = r['sweeps'][len(r['sweeps']) // 2]      # a representative sweep
            tg = np.linspace(0, max(sw['t']), 200)
            if sw['model'] == 'single':
                A, f_inf, tau = sw['params']
                fg = model_single(tg, A, f_inf, tau)
            else:
                fg = model_double(tg, *sw['params'])
            ex.append(dict(t=np.array(sw['t'], float), f=np.array(sw['f'], float),
                           t_fit=tg, f_fit=fg, tau_a=float(sw['tau_a'])))
        examples[t] = ex

    mat = dict(types=np.array(ORDER, dtype=object), model=model,
               per_type={t: per_type[t] for t in ORDER}, cells=cell_tbl,
               examples=examples, meta=meta)
    savemat(os.path.join(outdir, 'sfa_tau_fits.mat'), {'sfa_tau': mat}, do_compression=True)

    make_figure(outdir, cell_records, per_type, examples)
    return csv_path


def make_figure(outdir, cell_records, per_type, examples):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig_dir = os.path.join(outdir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    colors = {'pyr': '#1f77b4', 'pvalb': '#d62728', 'sst': '#2ca02c', 'vip': '#9467bd'}
    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(2, 4, height_ratios=[2, 1.4], hspace=0.35, wspace=0.3)

    # top: per-type example fits
    for k, t in enumerate(ORDER):
        ax = fig.add_subplot(gs[0, k])
        for ex in examples.get(t, []):
            ax.plot(ex['t'], ex['f'], '.', color=colors[t], ms=5, alpha=0.5)
            ax.plot(ex['t_fit'], ex['f_fit'], '-', color=colors[t], lw=1.2, alpha=0.9)
        med = per_type[t]['tau_a_median_s']; nc = per_type[t]['n_cells_fit']
        aidx = per_type[t]['adaptation_index_median']
        sub = (f"median $\\tau_a$={med*1e3:.0f} ms (n={nc})" if np.isfinite(med)
               else f"no reliable fit (n={nc})")
        ax.set_title(f"{t.capitalize()}   idx={aidx:.3f}\n{sub}", fontsize=11)
        ax.set_xlabel('time in step (s)')
        if k == 0:
            ax.set_ylabel('inst. rate (Hz)')
        ax.spines[['top', 'right']].set_visible(False)

    # bottom: tau_a distribution per type (log)
    axd = fig.add_subplot(gs[1, :])
    rng = np.random.default_rng(0)
    for k, t in enumerate(ORDER):
        taus = np.array([r['cell_tau'] for r in cell_records
                         if r['type'] == t and np.isfinite(r['cell_tau'])])
        if len(taus):
            x = k + (rng.random(len(taus)) - 0.5) * 0.3
            axd.scatter(x, taus * 1e3, s=14, color=colors[t], alpha=0.5, edgecolors='none')
            axd.plot([k - 0.25, k + 0.25], [np.median(taus) * 1e3] * 2, color='k', lw=2)
    axd.set_yscale('log')
    axd.set_xticks(range(len(ORDER)))
    axd.set_xticklabels([t.capitalize() for t in ORDER])
    axd.set_ylabel(r'$\tau_a$ (ms)')
    axd.set_title('Per-cell adaptation time constant (black = median)', fontsize=11)
    axd.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Spike-frequency adaptation time constant fit from Campagnola 2022 long-square sweeps',
                 fontsize=13, fontweight='bold')
    out = os.path.join(fig_dir, 'sfa_tau_fits.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"  wrote {out}", flush=True)


def write_provenance(outdir, meta):
    with open(os.path.join(outdir, 'sfa_tau_PROVENANCE.md'), 'w') as f:
        f.write("# SFA time constant (tau_a) fit -- provenance\n\n")
        f.write("- Source: raw long-square current-clamp sweeps, Allen synaptic physiology "
                "(Campagnola et al., Science 2022), fit for a firing-rate adaptation time constant.\n")
        f.write(f"- Database tier: `{meta['db_size']}`; current version: `{meta['db_version']}`.\n")
        f.write(f"- Species: {SPECIES}; cell types: {', '.join(ORDER)} "
                "(same CellClass defs as extract_campagnola_matrices.py).\n")
        f.write(f"- Experiments sampled: {meta['n_experiments']} (ranked by typed-cell count); "
                f"NWB cache: {meta['cache_path']}.\n")
        f.write(f"- Fit model: {meta['model']}  f(t)=f_inf+(f0-f_inf)exp(-t/tau_a), instantaneous "
                "rate 1/ISI at ISI midpoints vs time-in-step.\n")
        f.write(f"- QC: >= {MIN_SPIKES} spikes/sweep, R^2 >= {R2_MIN}, adaptation >= {MIN_ADAPT}, "
                f"tau in {TAU_BOUNDS} s (not railed); per-cell median over passing sweeps.\n")
        f.write("- Pure fit: tau_a is reported for every type (with median adaptation_index as "
                "context). Whether a weakly-adapting type (e.g. Pvalb) uses n_a=0 SFA timescales "
                "is decided later in the MATLAB Campagnola->parameter mapping, not here.\n")
        f.write(f"- Generated by `scripts/aisynphys/fit_sfa_tau.py`.\n")


# ------------------------------------------------------------------ #
def main():
    ap = argparse.ArgumentParser()
    default_out = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'src', 'connectivity', 'campagnola')
    ap.add_argument('--outdir', default=default_out)
    ap.add_argument('--max-experiments', type=int, default=25,
                    help='number of experiments (NWBs) to download+fit; 0 = all (huge)')
    ap.add_argument('--n-exp', type=int, default=1, choices=[1, 2],
                    help='1 = single exponential tau_a; 2 = fast+slow (for model n_a=2)')
    ap.add_argument('--list-only', action='store_true',
                    help='print the sampling plan (per-type cell counts) and exit; no downloads')
    ap.add_argument('--cached-only', action='store_true',
                    help='fit ONLY experiments whose NWB is already fully downloaded; skip the rest (no new downloads)')
    args = ap.parse_args()

    print(f"loading '{DB_SIZE}' database ...", flush=True)
    db = SynphysDatabase.load_current(DB_SIZE)
    print(f"DB: {db}", flush=True)

    print("== querying + classifying mouse intrinsic cells ==", flush=True)
    cells = load_typed_cells(db)
    counts = {t: sum(1 for c in cells if c['type'] == t) for t in ORDER}
    print(f"  typed cells with intrinsic: {len(cells)}  {counts}", flush=True)

    selected = select_experiments(cells, args.max_experiments)
    sel_counts = defaultdict(int)
    for cs in selected.values():
        for c in cs:
            sel_counts[c['type']] += 1
    print(f"  selected {len(selected)} experiments -> per-type cells: {dict(sel_counts)}", flush=True)
    if args.list_only:
        print("list-only: exiting before any download.", flush=True)
        return

    import aisynphys.config as config
    cache_path = getattr(config, 'cache_path', '~/ai_synphys_cache')

    cell_records = []
    for ei, (ext_id, cs) in enumerate(selected.items(), 1):
        print(f"[{ei}/{len(selected)}] experiment {ext_id}  ({len(cs)} typed cells)", flush=True)
        expt_obj = cs[0]['cell'].experiment
        if load_dataset(expt_obj, args.cached_only) is None:
            print(f"    [skip] NWB not fully cached (cached-only mode)", flush=True)
            continue
        try:
            for c in cs:
                try:
                    sweeps = fit_cell(c['cell'], args.n_exp, c['rheobase_pA'])
                except Exception:
                    print(f"    cell {c['cell'].id} FAILED:\n" + traceback.format_exc(), flush=True)
                    continue
                if not sweeps:
                    continue
                cell_tau = float(np.median([s['tau_a'] for s in sweeps]))
                cell_records.append(dict(cell_id=c['cell'].id, type=c['type'], expt=ext_id,
                                         cell_tau=cell_tau, sweeps=sweeps))
        finally:
            free_experiment(expt_obj)        # release the NWB before the next experiment
        n_ok = sum(1 for r in cell_records if r['expt'] == ext_id)
        print(f"    -> {n_ok} cells with a passing tau_a fit so far in this expt", flush=True)

    print(f"\nfit {len(cell_records)} cells total.", flush=True)
    adapt_index_by_type = {t: [] for t in ORDER}
    for cs in selected.values():
        for c in cs:
            adapt_index_by_type[c['type']].append(c['adaptation_index'])
    per_type = aggregate(cell_records, args.n_exp, adapt_index_by_type)
    for t in ORDER:
        pt = per_type[t]
        tau_str = (f"{pt['tau_a_median_s']*1e3:6.1f} ms" if np.isfinite(pt['tau_a_median_s'])
                   else "   n/a")
        print(f"  {t:6s}: tau_a median = {tau_str}  "
              f"(adapt_idx={pt['adaptation_index_median']:.4f}, n_cells={pt['n_cells_fit']}, "
              f"n_sweeps={pt['n_sweeps_fit']})", flush=True)

    versions = SynphysDatabase.list_current_versions()
    meta = dict(db_size=DB_SIZE, db_version=str(versions.get(DB_SIZE)),
                n_experiments=len(selected), cache_path=cache_path,
                model=('single' if args.n_exp == 1 else 'double'))
    csv_path = write_outputs(args.outdir, cell_records, per_type, args.n_exp, meta)
    write_provenance(args.outdir, meta)
    print(f"\nwrote {csv_path} + sfa_tau_fits.mat + figures/sfa_tau_fits.png + PROVENANCE", flush=True)
    print("DONE", flush=True)


if __name__ == '__main__':
    main()

function R = ospec_load_run(run_dir, row)
%OSPEC_LOAD_RUN One run of a run directory's outlier_spectra stage. THIS HEADER IS THE STAGE'S DATA CONTRACT (FD-040).
%
%   R = ospec_load_run(run_dir, row)
%
%   run_dir   a bundle's run directory, or its outlier_spectra folder
%   row       which run: a row NUMBER of the index (ospec_read_index), or
%             one ROW of its table (runs(k, :)), or a struct with fields
%             network, condition and seed (1 x 2)
%   R         the struct ospec_run saved
%
% TYPICAL USE
%   runs = ospec_read_index(run_dir);
%   k = find(strcmp(runs.network, 'EEx1p5') & strcmp(runs.condition, 'sfa3_stdc2') & runs.ok);
%   R = ospec_load_run(run_dir, k(1));
%   [model, params, ws] = ospec_rebuild_model(R);     % for Jacobian blocks; no integration
%   m = ospec_measures(R);                            % the derived measures
%   E = ospec_excursions(R.t, R.rate_top);            % excursions beyond zero, with start, end, cut-off
%
% NOTATION. n neurons, C cell types (2: E, I), N state variables, n_s
% sampled states in the window, n_cal calibration states, n_v states with
% saved eigenvectors, n_f states of the full-Jacobian series, n_m states of
% the followed mode. lambda is an eigenvalue of the frozen-gain effective
% connectivity M, IN UNITS OF W: its rate is (lambda - 1)/tau_d in 1/s.
% Eigenvalues of the full Jacobian are rates in 1/s. Block matrices C x C
% are (postsynaptic, presynaptic). Times are seconds from the start of the
% simulation (the window starts at T_settle).
%
% ---- what the run is --------------------------------------------------
%   R.stage_version    1
%   R.network          char, ospec_network's name ('reference', 'EEx1p5', ...)
%   R.condition        char, the preset's condition
%   R.seed, R.rng_seeds   1 x 2, the model's rng_seeds (the same thing twice)
%   R.preset_name      char
%   R.net              ospec_network's output: name, net_args,
%                      added_structure, target, label
%   R.net_args         cell, the name-value overrides that make the network
%   R.added            [] or what OspecNet added to W
%   R.build            everything ospec_rebuild_model needs (ospec_build_model)
%   R.opts             the options ospec_run was called with
%   R.window_s         measured window after settling (s)
%   R.has_full_series  logical: R.full is not empty
%   R.has_dense        logical: R.dense is not empty
%   R.n, R.N, R.tau_d (s), R.n_per_type (1 x C), R.cell_type_names (1 x C cell)
%   R.S_a              the linear fraction of the piecewise sigmoid
%   R.structural       R, lambda_O, mu_tilde_relative (the class's, theory),
%                      R_realized, K_eig_realized (from the drawn W)
%   R.lya              the class's lya_summary over the first
%                      R.opts.lya_window s of the window (LLE = lambda_1,
%                      1/s); an empty struct when lya_K was 0
%   R.seconds          wall-clock time of each pass
%
% ---- the series, frozen gains (every R.dt) ------------------------------
%   R.t           1 x n_s   sample times (s)
%   R.dt          sampling step (s)
%   R.ev_top      n_top x n_s complex   the rightmost modes of M (Im >= 0), by real part
%   R.rate_top    1 x n_s   (Re ev_top(1, :) - 1)/tau_d, 1/s: THE RIGHTMOST
%                 EIGENVALUE SERIES WITH FROZEN GAINS
%   R.R_eff       1 x n_s   the theoretical bulk radius at each state
%                 (variance profile of M with the state's gains), units of
%                 W. THIS IS THE BULK: unscaled.
%   R.ratio_theory   1 x n_s   Re ev_top(1, :) ./ R_eff: detachment against
%                 the theoretical radius (1 at the radius; NaN where R_eff
%                 is below 1e-6)
%   R.re_keep, R.abs_keep   n_keep x n_s   the largest real parts, the
%                 largest moduli of ALL eigenvalues of M, sorted descending
%                 (outliers are counted against any edge from these)
%   R.n_unstable  1 x n_s   eigenvalues of M with Re > 1 (rate > 0)
%   R.res, R.coh  n_top x n_s   residual of each rightmost eigenpair;
%                 fraction of its eigenvector in the population-uniform span
%   R.K_eig       C x n_s complex, R.K_eff C x C x n_s   the mean part of M
%                 and its eigenvalues (where an outlier sits if outside the bulk)
%   R.dbar, R.drms   C x C x n_s   mean and rms gain by route
%   R.rate, R.frac_on   C x n_s   mean rate by type; fraction of each type
%                 with nonzero slope
%   R.vec_t       1 x n_v   times of the saved eigenvectors (every
%                 R.opts.dt_vectors; a subset of R.t)
%   R.V_top       n x n_top x n_v single complex   unit right eigenvectors
%                 of ev_top at those times (phase arbitrary)
%
% ---- the null-calibrated factors: saved numbers, NOT the edge ----------
%   R.kappa       median over the calibration states of (95th percentile of
%                 the null draws' rightmost real part) / R_eff
%   R.kappa_abs   the same for the largest modulus
%   R.agree_re, R.agree_abs   the same with the null's MEDIAN: how far a
%                 bulk of this size reaches past R_eff with no structure
%
% ---- the calibration states (n_cal spread evenly over the window) -------
%   R.cal.t          1 x n_cal   times (s)
%   R.cal.idx        1 x n_cal   their columns in the series (R.t(R.cal.idx) == R.cal.t)
%   R.cal.S          n_cal x N   THE STATE VECTORS (one row per state; the
%                    layout is params.state_layout of ospec_rebuild_model)
%   R.cal.gain       n x n_cal   per-neuron gain: the slope of the
%                    nonlinearity times the route's output factor (the
%                    available vesicle fraction), averaged over the
%                    neuron's routes with the number of postsynaptic
%                    neurons of each type as weights. Dimensionless: the
%                    factor that scales the neuron's column of W in M.
%   R.cal.gain_route C x n x n_cal   the gain of neuron j on its route to
%                    postsynaptic type a, (a, j, state): M(i, j) =
%                    W(i, j) gain_route(type of i, j). NaN where the neuron
%                    has no synapse onto that type (small networks).
%   R.cal.slope      n x n_cal   slope of the nonlinearity (phi')
%   R.cal.u          n x n_cal   the argument of the nonlinearity relative
%                    to the neuron's own setpoint, u = x - (c/K) sum_k a_k
%                    - S_c_i (same units as x). The rate is the piecewise
%                    sigmoid centered on 0 at u:
%                    SRNNCellTypePairs.piecewiseSigmoid(u, R.S_a, 0).
%   R.cal.rate       n x n_cal   the rate r (0 to 1)
%   R.cal.cell_type  n x 1   cell type index (1 = E, 2 = I; R.cell_type_names)
%   R.cal.S_c        n x 1   the setpoints
%   R.cal.ev_M       n x n_cal complex   THE WHOLE SPECTRUM OF M at each state
%   R.cal.ev_uniform n x n_cal complex   the spectrum with every gain set to
%                    the rms gain (the structure's own outlier: the diamond)
%   R.cal.R_eff      1 x n_cal   the theoretical radius at each state
%   R.null           1 x n_cal struct array, ospec_null at each state:
%                    R_eff, null_re, null_abs (n_null x 1), q95_re, q95_abs,
%                    med_re, med_abs, obs_re, shuf_re (n_shuffle x 1),
%                    obs_pct_in_shuffle, typemean_re, uniform_re,
%                    uniform_mean_re, ev_M (n x 1), ev_null (n x 3: the
%                    first three null draws' spectra), ev_uniform (n x 1)
%
% ---- the full Jacobian ([] when R.has_full_series is false) -------------
%   R.full.t         1 x n_f   times (s), every R.full.dt over the WHOLE window
%   R.full.dt        step (s)
%   R.full.k         eigenvalues asked for per state (k_full)
%   R.full.lam       k x n_f complex   the k rightmost eigenvalues, 1/s, by
%                    real part (NaN where eigs did not converge)
%   R.full.xfrac     k x n_f   fraction of each eigenvector's norm^2 in the
%                    dendritic rows
%   R.full.pr        k x n_f   participation ratio over neurons (1 to n)
%   R.full.coh       k x n_f   population-uniform fraction of the x part
%   R.full.n_converged   1 x n_f
%   R.full.abscissa  1 x n_f   rightmost real part of ANY eigenvalue, 1/s
%   R.full.net       1 x n_f complex   THE RIGHTMOST NETWORK MODE (Im >= 0),
%                    1/s: the rightmost of the k whose eigenvector has
%                    xfrac >= R.full.x_thr. NaN + NaNi where none of the k
%                    is a network mode: it is then to the left of
%                    R.full.bound (a stable network, whose rightmost
%                    eigenvalues are single-neuron slow variables).
%   R.full.net_found 1 x n_f logical
%   R.full.bound     1 x n_f   smallest real part among the converged k:
%                    where net_found is false and n_converged == k, every
%                    network mode has real part below it
%   R.full.kin_edge  1 x n_f   the kinetic edge (rightmost eigenvalue of the
%                    slow block alone), 1/s; -Inf without slow variables
%   R.full.x_thr     the dendritic fraction that makes a network mode (0.1)
%   R.full.seconds   1 x n_f   wall time of each state
%   HOW NETWORK MODES ARE SEPARATED from the single-neuron slow
%   eigenvalues: by the eigenvector, not the eigenvalue
%   (ospec_full_spectrum has the reasons and the numbers). An eigenvalue
%   beyond zero is a network mode under any rule, so EXCURSIONS BEYOND
%   ZERO ARE THE SAME for R.full.abscissa and real(R.full.net):
%     E = ospec_excursions(R.full.t, real(R.full.net))
%   gives their start, end, duration and whether the window cut them off
%   (E.cut_start, E.cut_end). A threshold below zero is valid only above
%   max(R.full.bound(~R.full.net_found)).
%   R.dense          [] or 1 x n_cal struct array, ospec_dense_spectrum at
%                    the calibration states: lam (N x 1 single complex, the
%                    WHOLE spectrum, 1/s, by real part), ov, xfrac
%                    (N x 1 single), n_pos, abscissa, net and net_rank (the
%                    rightmost network mode over the whole spectrum, so
%                    known also in a stable network), kin_edge, cont,
%                    cont_ov, cont_xfrac (the continuation of the rightmost
%                    frozen-gain mode), lam1_M, rate1_M, mode_lam_M,
%                    mode_ov_M, mode_is_rightmost, mode_cont, mode_cont_ov,
%                    mode_cont_xfrac (the same for the intended mode), seconds
%
% ---- the intended mode (for the structure networks, FD-043) -------------
%   R.ref            ospec_mode_reference: the structural eigenvalue lam0,
%                    its right and left eigenvectors r0, l0, the
%                    participation p, the top carriers, ev_W (n x 1: the
%                    whole structural spectrum), ...
%   R.mode           struct of 1 x n_m rows, ospec_snapshot's mode fields
%                    every R.opts.dt_mode: t, lam, ov, is_rightmost,
%                    lam_right, lam_pred, d_rms, d_mean, gain_w, gain_ratio,
%                    shift_obs, shift_pred, on_top, on_rest, gain_top,
%                    gain_rest (only t when dt_mode was 0)
%
% WHAT IS NOT IN THE FILE: the trajectory, the weights (rebuild the model),
% the eigenvectors of M at the calibration states (one eig of
% ospec_effective_matrix gives them).
%
% See also: ospec_read_index, ospec_rebuild_model, ospec_measures,
%           ospec_excursions, ospec_run, run_outlier_spectra_stage

arguments
    run_dir (1,:) char
    row
end
[runs, ~, folder] = ospec_read_index(run_dir);
if istable(row)
    if height(row) ~= 1
        error('ospec_load_run:Row', 'Give one row of the index, not %d.', height(row));
    end
    file = row.file{1};
elseif isstruct(row)
    k = find(strcmp(runs.network, row.network) & strcmp(runs.condition, row.condition) & ...
        runs.seed1 == row.seed(1) & runs.seed2 == row.seed(2));
    if isempty(k)
        error('ospec_load_run:Row', 'The index of %s has no run %s, %s, seed [%d %d].', ...
            folder, row.network, row.condition, row.seed(1), row.seed(2));
    end
    file = runs.file{k};
elseif isnumeric(row) && isscalar(row) && row == round(row) && row >= 1 && row <= height(runs)
    file = runs.file{row};
else
    error('ospec_load_run:Row', 'row must be a row number from 1 to %d, one row of the index table, or a struct with network, condition and seed.', height(runs));
end
path = fullfile(folder, file);
if ~isfile(path)
    error('ospec_load_run:NoFile', 'The run file %s is not there (did that run fail? see the index''s ok and message).', path);
end
L = load(path, 'R');
R = L.R;
end

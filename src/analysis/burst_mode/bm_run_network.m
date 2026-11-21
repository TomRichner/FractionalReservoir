function R = bm_run_network(seed, plan, keep_ev)
%BM_RUN_NETWORK One network of the burst_mode stage: the run, the frozen arm, the spectra along both, and the measures.
%
%   R = bm_run_network(seed, plan, keep_ev)
%
% seed     1 x 2 rng_seeds of the network
% plan     burst_mode_stage_plan's plan
% keep_ev  true: keep every eigenvalue of M at every sample (R.intact.ev,
%          R.frozen.ev; single)
%
% THE RUN (fig_stim_engages_adaptation's): build_from_preset(preset, '')
% with a DC staircase [0 dc_level], hold_dur s per level, the white-noise
% probe noise_intensity, u_ex_scale 1, T_range [0, 2 hold_dur], the
% plan's fs and integrator, the full state stored.
%
% THE FROZEN ARM: from k, the last sample before the step (t < hold_dur),
% the run is integrated on with the class's own integrate() and
% bm_frozen_rhs, so only the dendritic rows move and every adaptive
% variable keeps its value at sample k; the same input. Checked: the
% restarted state equals the run's at k, and the held variables do not
% drift (R.check).
%
% EVERY dt_sample (both arms): M = W diag(d) from the class's Jacobian
% (ospec_effective_matrix, M = I + tau_d J_xx), its eigenvalues, the
% rightmost (Im >= 0) as a rate z1 = (lambda - 1)/tau_d, the coherence of
% its eigenvector (norm^2 in the span of the E and I patterns, ws.U; a
% random vector has 2/n), the predicted bulk radius R_eff (units of W:
% its rate edge is (R_eff - 1)/tau_d), mean gains dbar (post x pre) and
% rates by type; the mean adaptation term (c/K) sum_k a_k and the mean
% depression factor of each depressed route (std_route_factor), by type.
% EVERY dt_full (intact run): the rightmost network mode of the full
% Jacobian (full_jacobian_network_mode, dense).
%
% MEASURES: bursts (bm_detect_bursts on the mean E rate); per window
% (before, during; bm_crossings) for the frozen-gain mode, the full
% Jacobian's network mode and the frozen arm's frozen-gain mode;
% synchrony rho = (chi2 - 1/n)/(1 - 1/n) of the E population's dendritic
% potential at full rate (fig_stim_engages_adaptation's chi2 panel).
% R.summary is one table row (run_burst_mode_stage stacks them).
%
% Ported from scripts/explorations/stim_outlier/so_run_network (pilot,
% 2026-11-21); changed: the bursting figures' stimulus and integrator
% (the pilot's was noise-free rk4 at 400 Hz after 10 s of settling),
% rates throughout instead of eigenvalues of M against Re = 1, the full
% Jacobian's network mode at every dt_full (the pilot: every 100 ms, from
% eigs through ospec_full_spectrum), adaptation and depression traces and
% synchrony for S10.
%
% See also: burst_mode_stage_plan, run_burst_mode_stage, bm_crossings

t_start = tic;
H = plan.hold_dur;
ic = struct('dc_levels', [0 plan.dc_level], 'hold_dur', H, 'ramp_dur', 0, ...
    'noise_intensity', plan.noise_intensity, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
solver_args = {};
if ~isempty(plan.ode_solver); solver_args = {'ode_solver', plan.ode_solver}; end
m = build_from_preset(plan.preset_name, plan.condition, ...
    'input_config', ic, 'u_ex_scale', 1.0, 'fs', plan.fs, 'T_range', [0, 2 * H], ...
    'rng_seeds', seed, 'lya_method', 'none', 'store_full_state', true, 'plot_deci', 1, solver_args{:});
m.run();
if ~ismember(m.ode_solver, [{'rk4'}, stochastic_solver_names()])
    error('bm_run_network:Solver', 'The frozen arm needs a fixed-step integrator; the model runs %s.', m.ode_solver);
end
params = m.get_params();
ws = ospec_weight_stats(params);
t = m.t_out(:);
Sx = m.S_out;
tau = params.tau_d;
step = round(plan.dt_sample * plan.fs);
step_f = round(plan.dt_full / plan.dt_sample);
sel = (1:step:numel(t))';
u_mean = interp1(m.t_ex(:), mean(m.u_ex, 1)', t(sel), 'linear', 'extrap');

R = struct('seed', seed, 'preset_name', plan.preset_name, 'condition', plan.condition, ...
    'tau_d', tau, 'n', params.n, 'n_per_type', params.n_per_type, ...
    'cell_type_names', {params.cell_type_names}, 't_step', H, 'dc_level', plan.dc_level, ...
    'ode_solver', m.ode_solver, 'fs', plan.fs);
R.windows = struct('before', [plan.transient_s, H], 'during', [H + plan.transient_s, 2 * H]);
R.structure = struct('ev_W', eig(ws.Wf), 'R_struct', ws.R_struct, 'K_struct_eig', ws.K_struct_eig);

%% Intact run
R.intact = series(Sx(sel, :), t(sel), params, ws, keep_ev);
R.intact.u = u_mean';
R.bursts = bm_detect_bursts(R.intact.t, R.intact.rate(:, 1), plan.burst);
% The full Jacobian every dt_full, and with keep_ev every dt_sample from
% 1 s before to 0.5 s after the onset of each burst before the step (the
% zoom of Figure 8 B). The measures read the uniform grid only.
coarse = sel(1:step_f:end);
dense = [];
if keep_ev
    on = R.bursts.t_on(R.bursts.t_on < H);
    for b = 1:numel(on)
        dense = [dense; sel(t(sel) >= on(b) - 1 & t(sel) <= on(b) + 0.5)]; %#ok<AGROW>
    end
end
self = unique([coarse; dense]);
R.full = full_series(Sx(self, :), t(self), params, ws, plan.x_thr);
R.full.is_coarse = ismember(self, coarse);

%% Frozen arm
p = m.cached_params;
p.u_interpolant = m.u_interpolant;
x_idx = p.state_layout.x;
held = setdiff(1:p.state_layout.N, x_idx);
k = find(t < H, 1, 'last');
m.build_noise();
[~, S2] = m.integrate(@(tt, s) bm_frozen_rhs(tt, s, p, x_idx), t(k:end), Sx(k, :)');
R.check = struct('t_freeze', t(k), 'restart_mismatch', max(abs(S2(1, :) - Sx(k, :))), ...
    'held_drift', max(abs(S2(:, held) - S2(1, held)), [], 'all'));
fsel = sel(t(sel) >= H - plan.frozen_pre_s);
Sf = Sx(fsel, :);
after = fsel >= k;
Sf(after, :) = S2(fsel(after) - k + 1, :);
R.frozen = series(Sf, t(fsel), params, ws, keep_ev);
R.frozen.t_freeze = t(k);
R.frozen.bursts = bm_detect_bursts(R.frozen.t, R.frozen.rate(:, 1), plan.burst);

%% M check at one baseline state: M = I + tau_d J_xx, read back as W diag(d)
kc = sel(find(t(sel) >= min(20, H / 2), 1));
[Mf, d, J] = ospec_effective_matrix(Sx(kc, :)', params, ws);
Md = zeros(params.n);
for a = 1:params.n_cellTypes
    rows = params.type_indices{a};
    Md(rows, :) = ws.Wf(rows, :) .* d(a, :);
end
R.check.M_vs_I_plus_tau_Jxx = max(abs(Mf - (eye(params.n) + tau * full(J(x_idx, x_idx)))), [], 'all');
R.check.M_vs_W_diag_d = max(abs(Mf - Md), [], 'all');

%% Measures
R.meas = struct();
R.meas.frozen_gain = bm_crossings(R.intact.t, R.intact.z1, R.bursts, R.windows, plan.tol_pre, plan.tol_post);
R.meas.full = bm_crossings(R.full.t(R.full.is_coarse), R.full.net(R.full.is_coarse), R.bursts, R.windows, plan.tol_pre, plan.tol_post);
R.meas.frozen_arm = bm_crossings(R.frozen.t, R.frozen.z1, R.frozen.bursts, ...
    struct('during', R.windows.during), plan.tol_pre, plan.tol_post);
xE = Sx(:, x_idx(params.type_indices{1}));
R.sync = struct();
for w = {'before', 'during'}
    lim = R.windows.(w{1});
    kk = t >= lim(1) & t < lim(2);
    [R.sync.(['chi2_' w{1}]), R.sync.(['rho_' w{1}])] = pop_chi2(xE(kk, :)');
end
t2 = t(k:end);
[~, R.sync.rho_during_frozen] = pop_chi2(S2(t2 >= R.windows.during(1) & t2 < R.windows.during(2), ...
    x_idx(params.type_indices{1}))');
R.seconds = toc(t_start);
R.summary = summary_row(R);
end

%% ------------------------------------------------------------------------
function Q = series(S, ts, params, ws, keep_ev)
% Frozen-gain spectrum, rates, gains, adaptation and depression at each row of S.
ns = size(S, 1); n = params.n; C = params.n_cellTypes;
Q = struct('t', ts(:), 'z1', complex(nan(ns, 1)), 'coh1', nan(ns, 1), 'R_eff', nan(ns, 1), ...
    'dbar', nan(ns, C, C), 'rate', nan(ns, C), 'ev', []);
if keep_ev; Q.ev = complex(nan(ns, n, 'single')); end
for q = 1:ns
    [Mf, ~, ~, pred, pop] = ospec_effective_matrix(S(q, :)', params, ws);
    [V, D] = eig(Mf, 'vector');
    up = find(imag(D) >= 0);
    [~, i] = max(real(D(up)));
    i = up(i);
    v = V(:, i) / norm(V(:, i));
    Q.z1(q) = (D(i) - 1) / params.tau_d;
    Q.coh1(q) = sum(abs(ws.U' * v).^2);
    Q.R_eff(q) = pred.R_eff;
    Q.dbar(q, :, :) = pop.dbar;
    Q.rate(q, :) = pop.rate;
    if keep_ev; Q.ev(q, :) = single(D).'; end
end
[Q.adapt, Q.dep, Q.dep_names] = adaptation_traces(S, params);
end

function F = full_series(S, ts, params, ws, x_thr)
ns = size(S, 1);
F = struct('t', ts(:), 'net', complex(nan(ns, 1)), 'net_coh', nan(ns, 1), 'net_xfrac', nan(ns, 1), ...
    'abscissa', nan(ns, 1), 'n_net_unstable', nan(ns, 1));
for q = 1:ns
    J = SRNNCellTypePairs.compute_Jacobian_fast(S(q, :)', params);
    f = full_jacobian_network_mode(J, ws.is_x, x_thr, ws.U);
    F.net(q) = f.net; F.net_coh(q) = f.net_coh; F.net_xfrac(q) = f.net_xfrac;
    F.abscissa(q) = f.abscissa; F.n_net_unstable(q) = f.n_net_unstable;
end
end

function [adapt, dep, dep_names] = adaptation_traces(S, params)
% Mean adaptation term (c/K) sum_k a_k by type, and mean depression factor
% by depressed route (the model's own factor), at each row of S. As block B
% pilot 2's traces().
[~, a, b] = SRNNCellTypePairs.unpack_and_compute_states(S, params);
names = params.cell_type_names; C = numel(names); nt = size(S, 1);
adapt = nan(nt, C);
for q = 1:C
    if ~isempty(a.(names{q}))
        adapt(:, q) = mean(params.c_eff(q) .* reshape(sum(a.(names{q}), 2), [], nt), 1)';
    end
end
dep = zeros(nt, 0); dep_names = {};
for pre = 1:C
    for post = 1:C
        bb = b.(names{pre}).(names{post});
        if isempty(bb); continue; end
        D = SRNNCellTypePairs.std_route_factor(bb, params, pre, post);
        dep(:, end + 1) = mean(reshape(D, [], nt), 1)'; %#ok<AGROW>
        dep_names{end + 1} = [names{pre} ' to ' names{post}]; %#ok<AGROW>
    end
end
end

function [chi2, rho] = pop_chi2(M)
% Population synchrony of an n x nt matrix (fig_stim_engages_adaptation's
% chi2 panel): chi2 = var_t(mean_i M) / mean_i var_t(M_i), and the
% floor-corrected rho = (chi2 - 1/n)/(1 - 1/n).
den = mean(var(M, 0, 2));
if den < 1e-24
    chi2 = NaN; rho = NaN; return
end
chi2 = var(mean(M, 1)) / den;
n = size(M, 1);
rho = (chi2 - 1 / n) / (1 - 1 / n);
end

function T = summary_row(R)
g = R.meas.frozen_gain; f = R.meas.full; a = R.meas.frozen_arm;
win = @(Q, w) Q.t >= R.windows.(w)(1) & Q.t < R.windows.(w)(2);
kb = win(R.intact, 'before'); kd = win(R.intact, 'during'); kf = win(R.frozen, 'during');
T = table(R.seed(1), R.seed(2), g.before.n_bursts, g.during.n_bursts, ...
    g.before.n_cross, g.before.n_bursts_with_cross, g.before.median_lead, ...
    f.before.n_cross, f.before.n_bursts_with_cross, f.before.median_lead, ...
    g.before.frac_beyond, g.during.frac_beyond, f.before.frac_beyond, f.during.frac_beyond, ...
    g.during.max_re, f.during.max_re, a.during.frac_beyond, a.during.n_bursts, ...
    mean(R.intact.rate(kb, 1)), mean(R.intact.rate(kd, 1)), mean(R.frozen.rate(kf, 1)), ...
    mean(R.intact.dbar(kb, 1, 1)), mean(R.intact.dbar(kd, 1, 1)), ...
    R.sync.rho_before, R.sync.rho_during, R.sync.rho_during_frozen, R.seconds, ...
    'VariableNames', {'seed1', 'seed2', 'n_bursts_before', 'n_bursts_during', ...
    'n_cross_before_frozen_gain', 'n_bursts_with_cross_frozen_gain', 'median_lead_frozen_gain', ...
    'n_cross_before_full', 'n_bursts_with_cross_full', 'median_lead_full', ...
    'frac_beyond_before_frozen_gain', 'frac_beyond_during_frozen_gain', ...
    'frac_beyond_before_full', 'frac_beyond_during_full', ...
    'max_re_during_frozen_gain', 'max_re_during_full', ...
    'frac_beyond_during_frozen_arm', 'n_bursts_during_frozen_arm', ...
    'rate_E_before', 'rate_E_during', 'rate_E_during_frozen_arm', ...
    'gain_EE_before', 'gain_EE_during', 'rho_before', 'rho_during', 'rho_during_frozen_arm', 'seconds'});
end

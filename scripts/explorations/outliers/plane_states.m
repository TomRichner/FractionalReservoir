function [file, Q] = plane_states(network, condition, seed, opts)
%PLANE_STATES Re-simulate one saved outlier run and keep every neuron's operating point on the nonlinearity at the 12 calibration states.
%
%   file = plane_states('reference', 'no_adaptation')
%   [file, Q] = plane_states(network, condition, seed, 'data_dir', ..., 'out_dir', ...)
%
% WHY. outlier_run saves spectra and population summaries of the gains, not
% the states, so the per-neuron arguments of the nonlinearity that
% plane_fig_outliers' inset draws are not in run_*.mat. This simulates the
% same model again (serially, in the client; nothing else of outlier_run is
% repeated), takes the states at the saved run's calibration times
% (R.cal_t: the 12 states whose spectra the complex-plane figure overlays),
% and VERIFIES that it is the same trajectory by comparing with what the
% saved run holds at those states (below). Nothing is written into the
% saved run's folder.
%
% THE MODEL is built as outlier_run built it when the saved runs were made
% (commit 893dbc2, 2026-11-19: build_from_preset on the mu7chain preset,
% 'reference' with nothing changed, 'EEx<a>p<b>' with mu_EE_relative =
% a.b * 7, 'EI3to2' with block means [8.4 -5.6; 8.4 -5.6]; same seeds, run
% length, sampling rate and Lyapunov options, all read from the saved
% run's R.opts). Networks of outlier_network's later grammar (added
% low-rank structure, indegree) are not built here.
%
% WHAT IS KEPT, per neuron (n) and calibration state (n_cal), in Q:
%   u      n x n_cal: the argument of the nonlinearity relative to the
%          neuron's own setpoint, u = x - (c/K) sum_k a_k - S_c_i. The rate
%          is r = phi_0(u) with phi_0 the piecewise sigmoid centered on 0
%          (SRNNCellTypePairs.piecewiseSigmoid(u, S_a, 0)), the same curve
%          for every neuron: the class's phi is a pure translation in S_c.
%   r      n x n_cal: the rate, from the class's own activation function
%          (params.activation_function, which carries the per-neuron S_c)
%   rp     n x n_cal: its slope (params.activation_function_derivative)
%   S_c    n x 1 setpoints (the scalar repeated when there is no spread)
%   S_a    the linear fraction of the piecewise sigmoid
%   type   n x 1: cell type index of each neuron; type_names
%   t      1 x n_cal: the times, equal to the saved run's R.cal_t
%   check  struct, the comparison with the saved run at the same states:
%            ev_M      largest distance from a saved eigenvalue of the
%                      frozen-gain effective connectivity (R.cal.ev_M, all
%                      n at each state) to the nearest re-simulated one
%            dbar, drms  largest absolute difference of the mean and rms
%                      gains by route (R.dbar, R.drms)
%            rate, frac_on   likewise for the mean rate and the fraction
%                      with nonzero slope by type (R.rate, R.frac_on)
%            t         largest difference of the times
%          An exact repeat gives zeros. A chaotic run that differs at
%          round-off is a different trajectory after 60 s; the caller
%          decides what is acceptable (plane_fig_outliers refuses anything
%          above 1e-6).
%
% See also: plane_fig_outliers, outlier_run, outlier_effective_matrix

arguments
    network   (1,:) char
    condition (1,:) char
    seed      (1,2) double = [7 8]
    opts.data_dir (1,:) char = ''
    opts.out_dir  (1,:) char = ''
end

root = fileparts(which('setup_paths'));
data_dir = opts.data_dir; out_dir = opts.out_dir;
if isempty(data_dir); data_dir = fullfile(root, 'data', 'explorations', 'outliers'); end
if isempty(out_dir);  out_dir  = fullfile(root, 'data', 'explorations', 'outliers_plane'); end
if ~isfolder(out_dir); mkdir(out_dir); end

L = load(fullfile(data_dir, outlier_run_name(network, condition, seed, '')), 'R');
R = L.R;
o = R.opts;
if o.n ~= 500
    error('plane_states:Size', 'Only the 500-neuron runs are rebuilt here (saved run has n = %d).', o.n);
end

tok = regexp(network, '^EEx(\d+)p(\d+)$', 'tokens', 'once');
if strcmp(network, 'reference')
    net_args = {};
elseif strcmp(network, 'EI3to2')
    net_args = {'mu_tilde_relative', [8.4 -5.6; 8.4 -5.6]};
elseif ~isempty(tok)
    net_args = {'mu_EE_relative', str2double([tok{1} '.' tok{2}]) * 7};
else
    error('plane_states:Network', 'Network ''%s'' is not one this rebuilds (reference, EEx<a>p<b>, EI3to2).', network);
end
T_end = o.T_settle + o.T_measure;
lya_args = {'lya_method', 'none'};
if o.lya_K > 0
    lya_args = {'lya_method', 'topk', 'lya_K', o.lya_K, 'lya_K_auto', false, ...
        'lya_dt', 0.05, 'lya_T_interval', [o.T_settle T_end], ...
        'lya_warmup', min(5, o.T_settle)};
end
model = build_from_preset(R.preset, condition, net_args{:}, ...
    'rng_seeds', seed, 'T_range', [0 T_end], 'fs', o.fs, ...
    'store_full_state', true, 'verbose', 'near-none', lya_args{:});
t0 = tic;
model.run();
t_sim = toc(t0);

params = model.get_params();
ws = outlier_weight_stats(params);
S_a = model.S_a;
S_c_shared = model.S_c;
step = max(1, round(o.dt_sample * o.fs));
first = find(model.t_out >= o.T_settle, 1, 'first');
idx = first:step:numel(model.t_out);
cal_idx = unique(round(linspace(1, numel(idx), min(o.n_cal, numel(idx)))));
t = reshape(model.t_out(idx(cal_idx)), 1, []);
S_cal = model.S_out(idx(cal_idx), :)';            % N x n_cal
model = [];                                       %#ok<NASGU>  release the trajectory

n = params.n; C = params.n_cellTypes; ti = ws.ti; Lay = params.state_layout;
n_cal = numel(cal_idx);
S_c = params.S_c_vec;
if isempty(S_c); S_c = repmat(S_c_shared, n, 1); end
S_c = S_c(:);
type = zeros(n, 1);
for q = 1:C; type(ti{q}) = q; end

u = zeros(n, n_cal); r = u; rp = u;
chk = struct('ev_M', 0, 'dbar', 0, 'drms', 0, 'rate', 0, 'frac_on', 0, 't', max(abs(t - R.cal_t)));
sample_of = round((R.cal_t - R.t(1)) / R.dt) + 1;   % columns of R.dbar etc. at the calibration states
for k = 1:n_cal
    S = S_cal(:, k);
    x_eff = S(Lay.x);                               % as outlier_effective_matrix and compute_Jacobian_fast form it
    for q = 1:C
        if params.n_a(q) > 0
            a_q = reshape(S(Lay.a{q}), params.n_per_type(q), params.n_a(q));
            x_eff(ti{q}) = x_eff(ti{q}) - params.c_eff(q) .* sum(a_q, 2);
        end
    end
    r(:, k)  = params.activation_function(x_eff);
    rp(:, k) = params.activation_function_derivative(x_eff);
    u(:, k)  = x_eff - S_c;

    [Mf, ~, ~, ~, pop] = outlier_effective_matrix(S, params, ws);
    ev = eig(Mf);
    saved = R.cal(k).ev_M(:);
    chk.ev_M = max(chk.ev_M, max(min(abs(saved - ev.'), [], 2)));
    j = sample_of(k);
    chk.dbar = max(chk.dbar, max(abs(pop.dbar - R.dbar(:, :, j)), [], 'all'));
    chk.drms = max(chk.drms, max(abs(pop.drms - R.drms(:, :, j)), [], 'all'));
    chk.rate = max(chk.rate, max(abs(pop.rate(:) - R.rate(:, j))));
    chk.frac_on = max(chk.frac_on, max(abs(pop.frac_on(:) - R.frac_on(:, j))));
end
% The one curve: the class's rate is phi_0 at u.
curve_err = max(abs(r - SRNNCellTypePairs.piecewiseSigmoid(u, S_a, 0)), [], 'all');
if curve_err > 1e-12
    error('plane_states:Curve', 'The rate is not the centered piecewise sigmoid at u (largest difference %.3g).', curve_err);
end

Q = struct('network', network, 'condition', condition, 'seed', seed, 'preset', R.preset, ...
    'u', u, 'r', r, 'rp', rp, 'S_c', S_c, 'S_a', S_a, 'type', type, ...
    'type_names', {params.cell_type_names}, 't', t, 'check', chk, 'curve_err', curve_err, ...
    'seconds_sim', t_sim, 'source', outlier_run_name(network, condition, seed, ''));
file = fullfile(out_dir, plane_states_name(network, condition, seed));
save(file, 'Q');
fprintf('[plane_states] %s | %s | sim %.0f s | largest difference from the saved run: ev_M %.3g, dbar %.3g, drms %.3g, rate %.3g, frac_on %.3g, t %.3g\n', ...
    network, condition, t_sim, chk.ev_M, chk.dbar, chk.drms, chk.rate, chk.frac_on, chk.t);
end

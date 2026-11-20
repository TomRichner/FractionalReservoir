function [file, R] = outlier_run(network, condition, seed, opts)
%OUTLIER_RUN One network under one condition: simulate, settle, and take the outlier measurements along the trajectory.
%
%   file = outlier_run(network, condition)
%   [file, R] = outlier_run(network, condition, seed, Name, Value, ...)
%
% The network is the paper's (mu7chain preset, 500 neurons, noise 0.025,
% sra1 at 500 Hz) with at most one thing changed:
%   'reference'   as the preset: block means [7 -7; 7 -7], structurally
%                 balanced, NO structural outlier (the mean matrix is
%                 nilpotent)
%   'EEx1p5'      E-to-E mean times 1.5 (the excitation-dominant example of
%                 run_eig_heatmap and the manuscript's imbalance figure):
%                 a real positive structural outlier, about +26 in units of
%                 W against a bulk radius of 5.6
%   'EEx0p5'      E-to-E mean times 0.5: a complex pair with negative real
%                 part (not near the abscissa; here for completeness)
%   'EI3to2'      block means [8.4 -5.6; 8.4 -5.6] on the NOISY parent (the
%                 preset of that name is the noise-free twin; overriding
%                 the means on the parent keeps the noise)
% condition is one of the preset's: 'no_adaptation', 'sfa1_std1' (one
% timescale), 'sfa3_stdc2' (chained multiple timescales).
%
% WHAT IS MEASURED (all after T_settle; see the functions named)
%   fast pass   every dt_sample seconds: the spectrum of the frozen-gain
%               effective connectivity M (outlier_snapshot): the rightmost
%               modes and their eigenvectors, the largest real parts and
%               moduli, the predicted bulk radius and K_eff outliers, the
%               gains and rates by type
%   null        at n_cal states spread over the window: the bulk edge from
%               null draws and the gain controls (outlier_null)
%   full pass   every dt_full seconds, from the start of the window, until
%               full_budget_s of wall time is used: the rightmost
%               eigenvalues of the FULL Jacobian (outlier_full_spectrum),
%               and at n_kin of those states the kinetic edge
%   lambda_1    the top-K Lyapunov exponents over the same window (lya_K;
%               0 skips it)
%
% RUN LENGTH. Runs of different length on one seed are different noise
% realizations (the second Gaussian of sra1 is drawn after the first for
% the whole run), so every run compared must share T_settle + T_measure.
%
% MEMORY. During run() the class holds the whole trajectory: at 80 s,
% 500 Hz and 4000 states (multiple timescales, 500 neurons) 1.3 GB. Only
% the sampled rows are kept afterwards. The fast pass and the null use the
% session's pool when one is up (parfor; n x n matrices only); the full
% pass runs in the client, one state at a time.
%
% Name-value options (defaults)
%   T_settle 60, T_measure 20, fs 500, dt_sample 0.01, n 500,
%   n_top 3, n_keep 30, n_cal 8, n_null 40, n_shuffle 20,
%   dt_full 0.05, k_full 8, full_budget_s 240, n_kin 2, lya_K 3,
%   out_dir '' (data/explorations/outliers under the repository root),
%   tag '' (appended to the file name), verbose true
%
% Output: the .mat file's path (variable R, a struct; see the end of this
% function for its fields) and R itself.
%
% See also: explore_outliers, outlier_measures, fig_outliers

arguments
    network   (1,:) char {mustBeMember(network, {'reference', 'EEx1p5', 'EEx0p5', 'EI3to2'})}
    condition (1,:) char {mustBeMember(condition, {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'})}
    seed      (1,2) double = [7 8]
    opts.T_settle  (1,1) double = 60
    opts.T_measure (1,1) double = 20
    opts.fs        (1,1) double = 500
    opts.dt_sample (1,1) double = 0.01
    opts.n         (1,1) double = 500
    opts.n_top     (1,1) double = 3
    opts.n_keep    (1,1) double = 30
    opts.n_cal     (1,1) double = 8
    opts.n_null    (1,1) double = 40
    opts.n_shuffle (1,1) double = 20
    opts.dt_full   (1,1) double = 0.05
    opts.k_full    (1,1) double = 8
    opts.full_budget_s (1,1) double = 240
    opts.n_kin     (1,1) double = 2
    opts.lya_K     (1,1) double = 3
    opts.out_dir   (1,:) char = ''
    opts.tag       (1,:) char = ''
    opts.verbose   (1,1) logical = true
end

t_all = tic;
preset = ['celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_' ...
    'stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'];
T_end = opts.T_settle + opts.T_measure;

switch network
    case 'reference', net_args = {};
    case 'EEx1p5',    net_args = {'mu_EE_relative', 1.5 * 7};
    case 'EEx0p5',    net_args = {'mu_EE_relative', 0.5 * 7};
    case 'EI3to2',    net_args = {'mu_tilde_relative', [8.4 -5.6; 8.4 -5.6]};
end
size_args = {};
if opts.n ~= 500
    % The size rule of small_network_model: alpha stays 0.2 and the weight
    % unit F is the smaller network's own.
    k_in = max(2, round(0.2 * opts.n));
    size_args = {'n', opts.n, 'indegree', k_in, 'F_ref_n', opts.n, 'F_ref_indegree', k_in};
end
lya_args = {'lya_method', 'none'};
if opts.lya_K > 0
    lya_args = {'lya_method', 'topk', 'lya_K', opts.lya_K, 'lya_K_auto', false, ...
        'lya_dt', 0.05, 'lya_T_interval', [opts.T_settle T_end], ...
        'lya_warmup', min(5, opts.T_settle)};
end

model = build_from_preset(preset, condition, net_args{:}, size_args{:}, ...
    'rng_seeds', seed, 'T_range', [0 T_end], 'fs', opts.fs, ...
    'store_full_state', true, 'verbose', 'near-none', lya_args{:});
say(opts, '[outlier_run] %s | %s | seed [%d %d] | n = %d, %d states, %.0f s at %g Hz (trajectory about %.2f GB)\n', ...
    network, condition, seed(1), seed(2), model.n, model.N_sys_eqs, T_end, opts.fs, ...
    model.N_sys_eqs * T_end * opts.fs * 8 / 1e9);
t0 = tic;
model.run();
t_sim = toc(t0);

params = model.get_params();
ws = outlier_weight_stats(params);
tau_d = params.tau_d;

% The sampled rows, then let go of the trajectory.
step = max(1, round(opts.dt_sample * opts.fs));
first = find(model.t_out >= opts.T_settle, 1, 'first');
idx = first:step:numel(model.t_out);
t = model.t_out(idx);
t = t(:)';
S_sel = model.S_out(idx, :)';                 % N x n_s
lya = struct();
if opts.lya_K > 0
    lya = model.lya_summary();
end
structural = struct('R', model.R, 'lambda_O', model.lambda_O, ...
    'mu_tilde_relative', model.mu_tilde_relative, ...
    'R_realized', ws.R_struct, 'K_eig_realized', ws.K_struct_eig);
N = model.N_sys_eqs;
model = [];                                   %#ok<NASGU>  release the trajectory
n_s = numel(idx);
say(opts, '  simulated in %.0f s; %d sampled states, %.3g s apart\n', t_sim, n_s, step / opts.fs);

n_workers = 0;
pool = gcp('nocreate');
if ~isempty(pool); n_workers = pool.NumWorkers; end

%% Fast pass: the frozen-gain spectrum at every sampled state
t0 = tic;
snap = cell(1, n_s);
n_top = opts.n_top; n_keep = opts.n_keep;
parfor (k = 1:n_s, n_workers)
    snap{k} = outlier_snapshot(S_sel(:, k), params, ws, n_top, n_keep);
end
snap = [snap{:}];
t_fast = toc(t0);
say(opts, '  fast pass: %.0f s (%d workers)\n', t_fast, n_workers);

%% Null draws and gain controls at n_cal states
t0 = tic;
cal_idx = unique(round(linspace(1, n_s, min(opts.n_cal, n_s))));
cal = cell(1, numel(cal_idx));
S_cal = S_sel(:, cal_idx);
n_null = opts.n_null; n_shuffle = opts.n_shuffle;
parfor (k = 1:numel(cal_idx), n_workers)
    cal{k} = outlier_null(S_cal(:, k), params, ws, n_null, n_shuffle);
end
cal = [cal{:}];
t_null = toc(t0);
say(opts, '  null: %.0f s at %d states\n', t_null, numel(cal_idx));

%% Full pass: the rightmost eigenvalues of the full Jacobian, in the client
t0 = tic;
step_full = max(1, round(opts.dt_full / (step / opts.fs)));
full_idx = 1:step_full:n_s;
kin_at = unique(round(linspace(1, numel(full_idx), min(opts.n_kin, numel(full_idx)))));
full = cell(1, numel(full_idx));
n_done = 0;
for k = 1:numel(full_idx)
    if toc(t0) > opts.full_budget_s; break; end
    full{k} = outlier_full_spectrum(S_sel(:, full_idx(k)), params, ws, opts.k_full, false);
    n_done = k;
end
if n_done == 0
    error('outlier_run:NoFullPass', 'full_budget_s = %g left no time for a single state.', opts.full_budget_s);
end
full = [full{1:n_done}];
% The kinetic edge, separately, so the budget above is not spent on it.
kin = struct('edge', nan(1, numel(kin_at)), 'band', nan(numel(kin_at), 3), ...
    't', t(full_idx(kin_at)));
for q = 1:numel(kin_at)
    fq = outlier_full_spectrum(S_sel(:, full_idx(kin_at(q))), params, ws, 2, true);
    kin.edge(q) = fq.kin_edge;
    kin.band(q, :) = fq.kin_band;
end
t_full = toc(t0);
say(opts, '  full pass: %.0f s, %d of %d states done, kinetic edge %s\n', ...
    t_full, n_done, numel(full_idx), mat2str(kin.edge, 3));

%% Assemble
R = struct();
R.network = network; R.condition = condition; R.seed = seed;
R.preset = preset; R.opts = opts;
R.n = params.n; R.N = N; R.tau_d = tau_d; R.n_per_type = params.n_per_type;
R.structural = structural;
R.lya = lya;
R.t = t;                                                   % 1 x n_s
R.dt = step / opts.fs;
R.ev_top   = cat(2, snap.ev_top);                          % n_top x n_s
R.V_top    = cat(3, snap.V_top);                           % n x n_top x n_s, single
R.res      = cat(2, snap.res);
R.coh      = cat(2, snap.coh);
R.re_keep  = cat(2, snap.re_keep);                         % n_keep x n_s
R.abs_keep = cat(2, snap.abs_keep);
R.n_unstable = [snap.n_unstable];
R.R_eff    = [snap.R_eff];
R.K_eig    = cat(2, snap.K_eig);                           % C x n_s
R.K_eff    = cat(3, snap.K_eff);                           % C x C x n_s
R.dbar     = cat(3, snap.dbar);                            % C x C (post, pre) x n_s
R.drms     = cat(3, snap.drms);
R.rate     = cat(1, snap.rate)';                           % C x n_s
R.frac_on  = cat(1, snap.frac_on)';
R.cal = cal;  R.cal_t = t(cal_idx);                        % struct array, see outlier_null
R.full = struct('t', t(full_idx(1:n_done)), 'dt', step_full * step / opts.fs, ...
    'lam', cat(2, full.lam), 'xfrac', cat(2, full.xfrac), 'pr', cat(2, full.pr), ...
    'coh', cat(2, full.coh), 'n_converged', [full.n_converged], ...
    'seconds', [full.seconds], 'n_requested', numel(full_idx));
R.kin = kin;
R.seconds = struct('sim', t_sim, 'fast', t_fast, 'null', t_null, 'full', t_full, 'all', toc(t_all));

out_dir = opts.out_dir;
if isempty(out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'explorations', 'outliers');
end
if ~isfolder(out_dir); mkdir(out_dir); end
file = fullfile(out_dir, outlier_run_name(network, condition, seed, opts.tag));
save(file, 'R', '-v7.3');
say(opts, '  saved %s (%.0f s in all)\n', file, toc(t_all));
end

function say(opts, varargin)
if opts.verbose; fprintf(varargin{:}); end
end

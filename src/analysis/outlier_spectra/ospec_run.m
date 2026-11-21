function [file, R] = ospec_run(network, condition, seed, opts)
%OSPEC_RUN One run of the outlier_spectra stage: one network under one condition on one seed, simulated, settled, and measured along the trajectory.
%
%   [file, R] = ospec_run(network, condition, seed, Name, Value, ...)
%
% Ported from scripts/explorations/outliers/outlier_run.m (FD-040); the
% exploration's header has the question and the hypothesis. R, the struct
% saved (variable R of the .mat file), is documented field by field in the
% header of OSPEC_LOAD_RUN, which is the stage's data contract.
%
% THE NETWORK is the preset's (preset_name: the bundle's noisy
% SRNNCellTypePairs preset, sra1) with what its name changes
% (ospec_network has the whole grammar; tokens are joined by '-'):
%   'reference'     as the preset
%   'EEx1p5'        the E-to-E block mean times 1.5; '<XY>x<a>p<b>' the
%                   block mean mu_XY times a.b, XY (post, pre)
%   'EI3to2', 'ind<k>', 'lrA<a>p<b>W<c>p<d>'   the exploration's other
%                   networks (ospec_network, OspecNet)
% condition is one of the preset's.
%
% WHAT IS MEASURED (all after T_settle)
%   fast pass   every dt_sample seconds: the spectrum of the frozen-gain
%               effective connectivity M (ospec_snapshot): the rightmost
%               modes, the largest real parts and moduli, the predicted
%               bulk radius and K_eff outliers, the gains and rates by
%               type. The eigenvectors of the rightmost modes are kept
%               every dt_vectors seconds, and the INTENDED mode of the
%               network (ospec_mode_reference) is followed every dt_mode
%               seconds, from the same decomposition.
%   null        at n_cal states spread over the window (the CALIBRATION
%               STATES): the bulk edge from null draws and the gain
%               controls (ospec_null), and the whole spectrum of M
%   states      at the same states: the state vector and every neuron's
%               operating point (gain, slope, argument of the nonlinearity
%               relative to its own setpoint, rate), so that nothing
%               downstream re-simulates or recovers gains
%   full pass   (dt_full > 0) every dt_full seconds over the WHOLE window:
%               the k_full rightmost eigenvalues of the FULL Jacobian,
%               which of them are network modes, and the kinetic edge
%               (ospec_full_spectrum). No time budget: the exploration's
%               pass stopped after 200 s of wall time and two of its runs
%               were left incomplete.
%   dense pass  (dense = true) at the calibration states: the WHOLE
%               spectrum of the full Jacobian, the continuation of the
%               frozen-gain outlier and of the intended mode in it, and
%               the rightmost network mode (ospec_dense_spectrum)
%   lambda_1    the top lya_K Lyapunov exponents over the first lya_window
%               seconds of the measured window (0 skips it)
%
% WHAT DIFFERS FROM THE EXPLORATION'S RUN
%   * The trajectory is not held: OspecNet.run_sampled integrates the same
%     grid against the same noise in pieces and keeps the sampled rows, bit
%     for bit those of run(). So a 260 s run fits in a pool worker.
%   * The passes run one state after another in the calling process; the
%     stage runs whole RUNS in parallel (run_outlier_spectra_stage).
%   * No random number is drawn that the seed does not fix: the snapshot
%     takes its eigenvectors from eig, and the global stream is seeded
%     before each state's null (from the run's seed and the state's
%     number). A run is the same in the client and on a pool worker, up to
%     the round-off of the linear algebra (a worker has one thread).
%   * lambda_1 is over at most lya_window seconds, so that the stretch of
%     trajectory held for it stays small in a long run. At the default 20 s
%     and a 20 s window it is the exploration's.
%
% RUN LENGTH. Runs of different length on one seed are different noise
% realizations (the second Gaussian of sra1 is drawn after the first for
% the whole run), so runs compared state by state must share
% T_settle + T_measure.
%
% MEMORY, per process, at 500 neurons and 4000 state variables (multiple
% timescales), 500 Hz: the noise and the external input the class holds
% for the whole run (four n x steps arrays: 0.6 GB at 80 s, 2.1 GB at
% 260 s); one piece of trajectory (chunk_s 10 s: 0.16 GB); the stretch for
% lambda_1 (25 s: 0.4 GB, released before the passes); the sampled states
% (20 s at 10 ms: 64 MB; 200 s at 50 ms: 128 MB; 200 s at 10 ms: 640 MB);
% and in the dense pass one 4000 x 4000 complex eigenvector matrix with
% eig's workspace (about 1 GB while it runs). ESTIMATED, not measured:
% about 2 GB for a 20 s window and 3.5 GB for a 200 s window.
%
% Name-value options (defaults)
%   preset_name (the mu7chain preset), T_settle 60, T_measure 20, fs 500,
%   ode_solver '' (ospec_build_model's rule: sra1 for a noisy preset; the
%   stage hands the bundle's cfg.ode_solver),
%   dt_sample 0.01, n 0 (the preset's; another value builds a smaller
%   network by the rule of small_network_model: indegree 0.2 n and the
%   weight unit the smaller network's own), n_top 3, n_keep 30, n_cal 12,
%   n_null 40, n_shuffle 20, dt_vectors 0.1 (0: every sample), dt_mode 0.1
%   (0: not followed), dt_full 0.1 (0: no full pass), k_full 8,
%   x_thr 0.1 (the dendritic fraction that makes a network mode),
%   dense true, lya_K 3, lya_window 20, chunk_s 10,
%   out_dir '' (nothing is saved and file is ''), verbose true
%
% See also: ospec_load_run, ospec_measures, run_outlier_spectra_stage,
%           outlier_spectra_stage_plan, OspecNet

arguments
    network   (1,:) char
    condition (1,:) char
    seed      (1,2) double = [7 8]
    opts.preset_name (1,:) char = ['celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_' ...
        'stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain']
    opts.T_settle  (1,1) double = 60
    opts.T_measure (1,1) double = 20
    opts.fs        (1,1) double = 500
    opts.ode_solver (1,:) char  = ''     % '' = ospec_build_model's rule (sra1 for a noisy preset); else the bundle's integrator
    opts.dt_sample (1,1) double = 0.01
    opts.n         (1,1) double = 0
    opts.n_top     (1,1) double = 3
    opts.n_keep    (1,1) double = 30
    opts.n_cal     (1,1) double = 12
    opts.n_null    (1,1) double = 40
    opts.n_shuffle (1,1) double = 20
    opts.dt_vectors (1,1) double = 0.1
    opts.dt_mode   (1,1) double = 0.1
    opts.dt_full   (1,1) double = 0.1
    opts.k_full    (1,1) double = 8
    opts.x_thr     (1,1) double = 0.1
    opts.dense     (1,1) logical = true
    opts.lya_K     (1,1) double = 3
    opts.lya_window (1,1) double = 20
    opts.chunk_s   (1,1) double = 10
    opts.out_dir   (1,:) char = ''
    opts.verbose   (1,1) logical = true
end

t_all = tic;
T_end = opts.T_settle + opts.T_measure;

%% The model
p = srnn_param_preset(opts.preset_name);
if ~isfield(p, 'mu_tilde_relative') || ~isequal(size(p.mu_tilde_relative), [2 2])
    error('ospec_run:Preset', ['Preset %s does not state a 2 x 2 mu_tilde_relative: the network names ' ...
        'of this stage scale the block means of a two-type network.'], opts.preset_name);
end
net = ospec_network(network, p.mu_tilde_relative);
size_args = {};
if opts.n > 0 && opts.n ~= p.n
    % The size rule of small_network_model: alpha stays 0.2 and the weight
    % unit F is the smaller network's own.
    k_in = max(2, round(0.2 * opts.n));
    size_args = {'n', opts.n, 'indegree', k_in, 'F_ref_n', opts.n, 'F_ref_indegree', k_in};
end
lya_args = {'lya_method', 'none'};
if opts.lya_K > 0
    lya_args = {'lya_method', 'topk', 'lya_K', opts.lya_K, 'lya_K_auto', false, ...
        'lya_dt', 0.05, 'lya_T_interval', [opts.T_settle, opts.T_settle + min(opts.T_measure, opts.lya_window)], ...
        'lya_warmup', min(5, opts.T_settle)};
end
build = struct('preset_name', opts.preset_name, 'condition', condition, ...
    'net_args', {net.net_args}, 'size_args', {size_args}, 'added_structure', net.added_structure, ...
    'rng_seeds', seed, 'T_range', [0 T_end], 'fs', opts.fs, 'lya_args', {lya_args});
% Only when an integrator is named, so a run that names none saves the
% build record it always did.
if ~isempty(opts.ode_solver); build.ode_solver = opts.ode_solver; end
model = ospec_build_model(build);
say(opts, '[ospec_run] %s | %s | seed [%d %d] | n = %d, %d states, %.0f s at %g Hz\n', ...
    network, condition, seed(1), seed(2), model.n, model.N_sys_eqs, T_end, opts.fs);

params = model.get_params();
ws = ospec_weight_stats(params);
tau_d = params.tau_d;
n = params.n; C = params.n_cellTypes; N = params.N_sys_eqs;
target_value = [];
if ~isempty(model.added); target_value = model.added.target; end
ref = ospec_mode_reference(ws, net.target, target_value);

% The sampled rows of the time grid: every dt_sample from T_settle on.
step = max(1, round(opts.dt_sample * opts.fs));
t_grid = model.t_ex(:)';
first = find(t_grid >= opts.T_settle, 1, 'first');
idx = first:step:numel(t_grid);
t = t_grid(idx);
n_s = numel(idx);

t0 = tic;
[S_sel, lya] = model.run_sampled(idx, opts.chunk_s);          % N x n_s
t_sim = toc(t0);
structural = struct('R', model.R, 'lambda_O', model.lambda_O, ...
    'mu_tilde_relative', model.mu_tilde_relative, ...
    'R_realized', ws.R_struct, 'K_eig_realized', ws.K_struct_eig);
S_c = params.S_c_vec;
if isempty(S_c); S_c = repmat(model.S_c, n, 1); end
S_c = S_c(:);
S_a = model.S_a;
added = model.added;
clear model
say(opts, '  simulated in %.0f s; %d sampled states, %.3g s apart\n', t_sim, n_s, step / opts.fs);

%% Fast pass: the frozen-gain spectrum at every sampled state
t0 = tic;
every = @(dt_x) max(1, round(dt_x / (step / opts.fs)));
vec_idx = 1:every(opts.dt_vectors):n_s;
if opts.dt_vectors <= 0; vec_idx = 1:n_s; end
mode_idx = [];
if opts.dt_mode > 0; mode_idx = 1:every(opts.dt_mode):n_s; end
is_vec = false(1, n_s);  is_vec(vec_idx) = true;
is_mode = false(1, n_s); is_mode(mode_idx) = true;

n_top = min(opts.n_top, floor(n / 2));
n_keep = min(opts.n_keep, n);
F = struct('ev_top', complex(zeros(n_top, n_s)), 'res', zeros(n_top, n_s), 'coh', zeros(n_top, n_s), ...
    're_keep', zeros(n_keep, n_s), 'abs_keep', zeros(n_keep, n_s), 'n_unstable', zeros(1, n_s), ...
    'R_eff', zeros(1, n_s), 'K_eig', complex(zeros(C, n_s)), 'K_eff', zeros(C, C, n_s), ...
    'dbar', zeros(C, C, n_s), 'drms', zeros(C, C, n_s), 'rate', zeros(C, n_s), 'frac_on', zeros(C, n_s));
V_top = complex(zeros(n, n_top, numel(vec_idx), 'single'));
modes = cell(1, numel(mode_idx));
i_vec = 0; i_mode = 0;
for k = 1:n_s
    if is_mode(k)
        s = ospec_snapshot(S_sel(:, k), params, ws, n_top, n_keep, ref);
        i_mode = i_mode + 1;
        modes{i_mode} = s.mode;
    else
        s = ospec_snapshot(S_sel(:, k), params, ws, n_top, n_keep);
    end
    F.ev_top(:, k) = s.ev_top;   F.res(:, k) = s.res;        F.coh(:, k) = s.coh;
    F.re_keep(:, k) = s.re_keep; F.abs_keep(:, k) = s.abs_keep;
    F.n_unstable(k) = s.n_unstable; F.R_eff(k) = s.R_eff;
    F.K_eig(:, k) = s.K_eig;     F.K_eff(:, :, k) = s.K_eff;
    F.dbar(:, :, k) = s.dbar;    F.drms(:, :, k) = s.drms;
    F.rate(:, k) = s.rate(:);    F.frac_on(:, k) = s.frac_on(:);
    if is_vec(k)
        i_vec = i_vec + 1;
        V_top(:, :, i_vec) = s.V_top;
    end
end
mode = struct_of_arrays([modes{:}]);
mode.t = t(mode_idx);
t_fast = toc(t0);
say(opts, '  fast pass: %.0f s (%.3f s a state)\n', t_fast, t_fast / n_s);

%% The calibration states: null draws, gain controls, and every neuron's operating point
t0 = tic;
cal_idx = unique(round(linspace(1, n_s, min(opts.n_cal, n_s))));
n_cal = numel(cal_idx);
S_cal = S_sel(:, cal_idx);
null = cell(1, n_cal);
cal = struct('t', t(cal_idx), 'idx', cal_idx, 'S', S_cal.', ...
    'gain', nan(n, n_cal), 'gain_route', nan(C, n, n_cal), 'slope', zeros(n, n_cal), ...
    'u', zeros(n, n_cal), 'rate', zeros(n, n_cal), 'cell_type', ws.type_of(:), 'S_c', S_c, ...
    'ev_M', complex(zeros(n, n_cal)), 'ev_uniform', complex(zeros(n, n_cal)), 'R_eff', zeros(1, n_cal));
n_per = cellfun(@numel, ws.ti);
has_route = ws.col_sq > 0;                                    % C x n: the neuron has a synapse onto that type
for k = 1:n_cal
    rng(null_seed(seed, k), 'twister');
    null{k} = ospec_null(S_cal(:, k), params, ws, opts.n_null, opts.n_shuffle);
    [~, d, ~, ~, pop] = ospec_effective_matrix(S_cal(:, k), params, ws);
    d(~has_route) = NaN;
    cal.gain_route(:, :, k) = d;
    w = n_per(:) .* has_route;                                % weights of the routes that exist
    d0 = d; d0(~has_route) = 0;
    cal.gain(:, k) = (sum(w .* d0, 1) ./ max(sum(w, 1), realmin)).';
    cal.gain(sum(w, 1) == 0, k) = NaN;
    cal.slope(:, k) = pop.rp;
    cal.rate(:, k)  = pop.r;
    cal.u(:, k)     = pop.x_eff - S_c;
    cal.ev_M(:, k)  = null{k}.ev_M;
    cal.ev_uniform(:, k) = null{k}.ev_uniform;
    cal.R_eff(k)    = null{k}.R_eff;
end
null = [null{:}];
t_null = toc(t0);
say(opts, '  null: %.0f s at %d states\n', t_null, n_cal);

% The null-calibrated edge factors: saved numbers, not the edge. The bulk
% is the unscaled predicted radius R_eff (Tom, 2026-11-20).
ok = isfinite(cal.R_eff) & cal.R_eff > 1e-6;
kap = struct('kappa', NaN, 'kappa_abs', NaN, 'agree_re', NaN, 'agree_abs', NaN);
if any(ok)
    kap.kappa     = median([null(ok).q95_re]  ./ cal.R_eff(ok));
    kap.kappa_abs = median([null(ok).q95_abs] ./ cal.R_eff(ok));
    kap.agree_re  = median([null(ok).med_re]  ./ cal.R_eff(ok));
    kap.agree_abs = median([null(ok).med_abs] ./ cal.R_eff(ok));
end

%% Full pass: the rightmost eigenvalues of the full Jacobian over the whole window
t0 = tic;
full = [];
if opts.dt_full > 0
    full_idx = 1:every(opts.dt_full):n_s;
    n_f = numel(full_idx);
    fs_ = cell(1, n_f);
    for k = 1:n_f
        fs_{k} = ospec_full_spectrum(S_sel(:, full_idx(k)), params, ws, opts.k_full, opts.x_thr);
    end
    fs_ = [fs_{:}];
    full = struct('t', t(full_idx), 'dt', every(opts.dt_full) * step / opts.fs, ...
        'lam', cat(2, fs_.lam), 'xfrac', cat(2, fs_.xfrac), 'pr', cat(2, fs_.pr), ...
        'coh', cat(2, fs_.coh), 'n_converged', [fs_.n_converged], 'seconds', [fs_.seconds], ...
        'abscissa', [fs_.abscissa], 'net', [fs_.net], 'net_found', [fs_.net_found], ...
        'bound', [fs_.bound], 'kin_edge', [fs_.kin_edge], 'x_thr', opts.x_thr, 'k', size(fs_(1).lam, 1));
    say(opts, '  full pass: %.0f s at %d states (%.2f s a state); a network mode among the %d rightmost at %.0f%% of states, beyond zero at %.0f%%\n', ...
        toc(t0), n_f, toc(t0) / n_f, full.k, 100 * mean(full.net_found), 100 * mean(real(full.net) > 0));
end
t_full = toc(t0);

%% Dense pass: the whole full-Jacobian spectrum at the calibration states
t0 = tic;
dense = [];
if opts.dense
    dense = cell(1, n_cal);
    for k = 1:n_cal
        dense{k} = ospec_dense_spectrum(S_cal(:, k), params, ws, ref, opts.x_thr);
    end
    dense = [dense{:}];
    say(opts, '  dense pass: %.0f s at %d states; continuation of the rightmost (median) %s 1/s, rightmost network mode %s 1/s, unstable eigenvalues (median) %g\n', ...
        toc(t0), n_cal, num2str(median([dense.cont]), '%.3g'), num2str(median([dense.net]), '%.3g'), median([dense.n_pos]));
end
t_dense = toc(t0);

%% Assemble (the fields are documented in ospec_load_run)
R = struct();
R.stage_version = 1;
R.network = network; R.condition = condition; R.seed = seed; R.rng_seeds = seed;
R.preset_name = opts.preset_name;
R.net = net; R.net_args = net.net_args; R.added = added;
R.build = build;
R.opts = opts;
R.window_s = opts.T_measure;
R.has_full_series = ~isempty(full);
R.has_dense = ~isempty(dense);
R.n = n; R.N = N; R.tau_d = tau_d; R.n_per_type = params.n_per_type;
R.cell_type_names = params.cell_type_names; R.S_a = S_a;
R.structural = structural;
R.lya = lya;
R.t = t; R.dt = step / opts.fs;
for f = fieldnames(F)'
    R.(f{1}) = F.(f{1});
end
R.rate_top = (real(F.ev_top(1, :)) - 1) ./ tau_d;
R.ratio_theory = nan(1, n_s);
okr = F.R_eff > 1e-6;
R.ratio_theory(okr) = real(F.ev_top(1, okr)) ./ F.R_eff(okr);
R.vec_t = t(vec_idx); R.V_top = V_top;
R.kappa = kap.kappa; R.kappa_abs = kap.kappa_abs; R.agree_re = kap.agree_re; R.agree_abs = kap.agree_abs;
R.cal = cal;
R.null = null;
R.dense = dense;
R.full = full;
R.ref = ref;
R.mode = mode;
R.seconds = struct('sim', t_sim, 'fast', t_fast, 'null', t_null, 'full', t_full, 'dense', t_dense, 'all', toc(t_all));

file = '';
if ~isempty(opts.out_dir)
    if ~isfolder(opts.out_dir); mkdir(opts.out_dir); end
    file = fullfile(opts.out_dir, ospec_run_name(network, condition, seed));
    save(file, 'R', '-v7');
    say(opts, '  saved %s (%.0f s in all)\n', file, toc(t_all));
end
end

%% ------------------------------------------------------------------------
function s = null_seed(seed, k)
% The seed of the global stream before the null of calibration state k: a
% function of the run's seed pair and k only.
s = mod(1000003 * seed(1) + 1009 * seed(2) + k, 2^32 - 1);
end

function A = struct_of_arrays(s)
% A 1 x m struct array with scalar fields as one struct of 1 x m rows (an
% empty struct when s is empty).
A = struct();
if isempty(s); return; end
for f = fieldnames(s)'
    A.(f{1}) = [s.(f{1})];
end
end

function say(opts, varargin)
if opts.verbose; fprintf(varargin{:}); end
end

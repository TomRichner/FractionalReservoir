function out_dir = run_noise_switch_rollouts(cfg)
% RUN_NOISE_SWITCH_ROLLOUTS Copies of one network that share the noise, then stop sharing it.
%
%   out_dir = RUN_NOISE_SWITCH_ROLLOUTS()
%   out_dir = RUN_NOISE_SWITCH_ROLLOUTS('preset_name', p, 'n_copies', 8, 'out_dir', d)
%
% THE INDEPENDENT-NOISE ENSEMBLE (Tom, 2026-11-14). One network (the same W,
% setpoints, time constants and initial condition) is run as n_copies
% copies. Every copy sees the IDENTICAL Wiener path for the first t_shared
% seconds, long enough to leave the initial transient, and from then on
% each copy continues on a noise seed of its own for t_rollout seconds
% (SRNNCellTypePairs.noise_seed_schedule). The copies are therefore one
% trajectory up to the switch and an ensemble after it, and the stage
% measures how fast they drift apart, in which state variables, and how
% that differs between the preset's adaptation conditions.
%
% This is NOT the Lyapunov measurement. There the two trajectories share
% the noise, it cancels in their difference, and what is left is the
% growth of an infinitesimal perturbation. Here the copies stop sharing
% it, so their difference is driven: it measures how the noise spreads the
% state (the covariance the noise maintains), which a stable network has
% too. lambda_1 of one copy over the shared window is stored for the
% comparison (top-K method, K = 1).
%
% WHAT IS MEASURED. Five blocks of per-neuron series, each from the
% model's own plot data: x (dendritic state), r (rate), theta (synaptic
% output, route 1; the preset's routes must carry identical synaptic
% dynamics, SRNNCellTypePairs.routes_identical, so that a neuron has one),
% a (the adaptation term (c/K) sum_k a_k that enters the rate) and b (the
% route's depression factor from SRNNCellTypePairs.std_route_factor: the
% releasable pool b_1 on a chain route, the product of the pools on a
% product route). A block holds the neurons that have the variable; a
% condition without SFA or without depression has no a or no b block.
% Per block X (neurons x time x copies), at every stored time:
%   rms_diff    root of the mean over copy pairs and neurons of
%               (X_i - X_j)^2. Computed as sqrt(2 variance): the mean
%               over pairs of a squared difference is exactly twice the
%               unbiased across-copy variance.
%   variance    across-copy variance (unbiased), mean over neurons
%   mean_norm   RMS over neurons of the ensemble mean
%   mean_fluct_norm   the same after subtracting each neuron's temporal
%               mean over the reference window: what is left of the
%               common fluctuation once the copies average it away
%   corr        Pearson correlation across neurons between two copies'
%               deviations from that temporal mean, mean over copy pairs
% and, so that blocks can be compared, the first four divided by the
% block's own temporal standard deviation `scale` (variance by its
% square): the root of the mean over neurons of the temporal variance over
% the reference window, the second half of the shared window (the first
% half is left to the initial transient). They carry the suffix _rel. Two
% independent draws of a stationary process differ by sqrt(2) of its
% standard deviation, so rms_diff_rel near 1.41 means the copies have
% become independent samples of the block's fluctuation.
%
% THE FIT (per block and condition; `fit`), on rms_diff_rel against time
% since the switch:
%   saturation   mean over the last fifth of the rollout
%   saturated    true when that differs from the mean over the fifth
%                before it by under 10 %; false says the rollout was too
%                short for the block (slow adaptation variables)
%   t_half       time to first reach half the saturation (interpolated)
%   growth_rate  least-squares slope of log(rms_diff_rel) against time
%                (1/s), from the first sample after the switch through
%                the first at or above half the saturation
%   rise_exponent  slope of the same samples in log-log: about 0.5 when
%                the early rise is the noise's own diffusion, larger when
%                the rise is exponential; says how to read growth_rate
% growth_rate and rise_exponent are NaN with fewer than three samples.
%
% THE RUN. Each copy is built separately from the same seeds
% (build_from_preset; the copies run in the session's pool when it holds
% one, pool_ensure, serially otherwise) at the integrator and rate the
% sweeps use for the preset at cfg.run_mode
% (analysis_run_config('param_space', ...): sra1 at 400 Hz for a noisy
% preset at 'medium', the default). The stage ASSERTS that the copies are
% bit-identical up to and including the switch (every stored sample of
% every block, and the full state vector at the switch time) and that
% they differ afterwards.
%
% THE FILE, <out_dir>/noise_switch_rollouts_data.mat: settings, cond_names
% and results, one element per condition (name, title, lambda_1, t,
% blocks, block_labels, stats.<block>, example). Statistics are stored at
% cfg.store_fs (default 200 Hz: every other sample of a 400 Hz run, plain
% subsampling, no filtering). example holds the series of two copies
% (1 and 2) for the plotted neurons only, the first n_show_E and
% n_show_I of each type as fig_example_timeseries chooses them, in single
% precision at cfg.example_fs (default 100 Hz, subsampled again).
% git_provenance.txt is written beside it.
%
% Options not described above: conditions (names, in order; empty = all
% the preset states; select_conditions), noise_seed (the shared path's
% seed), rollout_seeds (one per copy; empty = noise_seed + (1:n_copies)),
% rng_seeds (network and stimulus seeds), model_args (name-value overrides
% of the preset, for tests on a small network).
%
% See also: fig_noise_switch_rollouts, wiener_increments, SRNNCellTypePairs,
%           noise_switch_rollouts_example

arguments
    cfg.preset_name   (1,:) char   = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'
    cfg.conditions                 = {}
    cfg.t_shared      (1,1) double {mustBePositive} = 10
    cfg.t_rollout     (1,1) double {mustBePositive} = 5
    cfg.n_copies      (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(cfg.n_copies, 2)} = 8
    cfg.noise_seed    (1,1) double {mustBeInteger, mustBeNonnegative} = 101
    cfg.rollout_seeds (1,:) double {mustBeInteger, mustBeNonnegative} = []
    cfg.rng_seeds     (1,2) double = [1 2]
    cfg.run_mode      (1,:) char   = 'medium'
    cfg.store_fs      (1,1) double {mustBePositive} = 200
    cfg.example_fs    (1,1) double {mustBePositive} = 100
    cfg.n_show_E      (1,1) double {mustBeInteger, mustBeNonnegative} = 25
    cfg.n_show_I      (1,1) double {mustBeInteger, mustBeNonnegative} = 25
    cfg.out_dir       (1,:) char   = ''
    cfg.verbose                    = 'minimal'
    cfg.model_args    (1,:) cell   = {}
end

setup_paths();
project_root = fileparts(which('setup_paths'));

%% Settings, resolved before any compute
if ~ismember(cfg.run_mode, run_mode_names())
    error('run_noise_switch_rollouts:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
        cfg.run_mode, strjoin(run_mode_names(), ', '));
end
[preset, model_class, conditions] = srnn_param_preset(cfg.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_noise_switch_rollouts:WrongClass', ...
        'Preset ''%s'' is written for %s; this stage runs SRNNCellTypePairs only.', ...
        cfg.preset_name, model_class);
end
cond_names = select_conditions(cellfun(@(c) c.name, conditions, 'UniformOutput', false), cfg.conditions);
n_cond = numel(cond_names);

M = cfg.n_copies;
rollout_seeds = cfg.rollout_seeds;
if isempty(rollout_seeds), rollout_seeds = cfg.noise_seed + (1:M); end
if numel(rollout_seeds) ~= M || numel(unique(rollout_seeds)) ~= M
    error('run_noise_switch_rollouts:BadSeeds', ...
        'rollout_seeds must hold one distinct seed per copy (%d).', M);
end

run_cfg = analysis_run_config('param_space', cfg.run_mode, preset);
T = cfg.t_shared + cfg.t_rollout;
P = struct();
P.preset_name   = cfg.preset_name;
P.verbose       = cfg.verbose;
P.model_args    = cfg.model_args;
P.rng_seeds     = cfg.rng_seeds;
P.fs            = run_cfg.model.fs;
P.ode_solver    = run_cfg.model.ode_solver;
P.T_range       = [0, T];
P.t_shared      = cfg.t_shared;
P.noise_seed    = cfg.noise_seed;
P.rollout_seeds = rollout_seeds;
P.store_deci    = whole_ratio(P.fs, cfg.store_fs, 'store_fs');
P.n_show        = [cfg.n_show_E, cfg.n_show_I];
% lambda_1 of copy 1 over the second half of the shared window, the
% window the temporal scale is taken over; warm-up as fig_example_timeseries.
P.lya_T_interval = [cfg.t_shared / 2, cfg.t_shared];
P.lya_warmup     = cfg.t_shared / 4;
P.lya_dt         = 0.05;
store_fs     = P.fs / P.store_deci;
example_deci = whole_ratio(store_fs, cfg.example_fs, 'example_fs');

if isempty(cfg.out_dir)
    out_dir = fullfile(project_root, 'data', 'noise_switch_rollouts');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir), mkdir(out_dir); end

vprintf(cfg.verbose, 'minimal', ['[noise_switch_rollouts] preset=%s: %d condition(s) x %d copies, ' ...
    'shared %g s then %g s on own seeds, %s at %g Hz\n'], ...
    cfg.preset_name, n_cond, M, cfg.t_shared, cfg.t_rollout, P.ode_solver, P.fs);

pool = pool_ensure('verbose', cfg.verbose);
n_pool = 0;                              % no pool: the copies run serially in the client
if ~isempty(pool), n_pool = pool.NumWorkers; end

%% Run: per condition, the copies in parallel, reduced in the client
t_stage = tic;
results = struct('name', cond_names, 'title', titles_by_name(cond_names, cond_names));
for i = 1:n_cond
    t0 = tic;
    copies = cell(1, M);
    cname = cond_names{i};
    parfor m = 1:M
        copies{m} = rollout_copy(P, cname, m);
    end
    c1 = copies{1};
    assert_shared_then_split(copies, c1.t, cfg.t_shared, cname);

    results(i).lambda_1       = c1.lambda_1;
    results(i).lya_T_interval = P.lya_T_interval;
    results(i).n              = c1.n;
    results(i).N_sys_eqs      = c1.N_sys_eqs;
    results(i).type_names     = c1.type_names;
    results(i).n_per_type     = c1.n_per_type;
    results(i).t              = c1.t(:)';
    results(i).blocks         = fieldnames(c1.series)';
    results(i).block_labels   = c1.labels;
    results(i).identical_before_switch = true;
    results(i).stats   = struct();
    results(i).example = example_series(copies(1:2), c1, example_deci);
    for k = 1:numel(results(i).blocks)
        blk = results(i).blocks{k};
        X = cellfun(@(c) c.series.(blk), copies, 'UniformOutput', false);
        X = cat(3, X{:});                               % neurons x time x copies
        results(i).stats.(blk) = ensemble_statistics(X, results(i).t, cfg.t_shared);
        results(i).stats.(blk).neurons = c1.neurons.(blk);
    end
    sx = results(i).stats.x.fit;
    vprintf(cfg.verbose, 'verbose', ['  %-14s %d copies in %4.0f s | lambda_1 %+.3f /s | x: growth %.3g /s, ' ...
        't_half %.3g s, saturation %.3g\n'], cname, M, toc(t0), c1.lambda_1, ...
        sx.growth_rate, sx.t_half, sx.saturation);
end

%% Save
sigma = c1.sigma_u_noise;
settings = struct('preset_name', cfg.preset_name, 'model_class', model_class, ...
    'run_mode', cfg.run_mode, 'ode_solver', P.ode_solver, 'fs', P.fs, ...
    'store_fs', store_fs, 'example_fs', store_fs / example_deci, ...
    't_shared', cfg.t_shared, 't_rollout', cfg.t_rollout, 'T_range', P.T_range, ...
    'reference_window', [cfg.t_shared / 2, cfg.t_shared], ...
    'n_copies', M, 'noise_seed', cfg.noise_seed, 'rollout_seeds', rollout_seeds, ...
    'rng_seeds', cfg.rng_seeds, 'sigma_u_noise', sigma, 'n', c1.n, ...
    'example_copies', [1 2], 'model_args', {cfg.model_args}, ...
    'lya_method', 'topk', 'lya_K', 1, 'lya_dt', P.lya_dt, 'lya_warmup', P.lya_warmup, ...
    'n_workers', n_pool, 'minutes', toc(t_stage) / 60, 'created', char(datetime('now')));
if verbose_level(cfg.verbose) >= verbose_level('verbose')
    settings.git = capture_git_provenance(out_dir, project_root);
else
    [~, settings.git] = evalc('capture_git_provenance(out_dir, project_root)');
end
mat_file = fullfile(out_dir, 'noise_switch_rollouts_data.mat');
save(mat_file, 'results', 'cond_names', 'settings');
vprintf(cfg.verbose, 'minimal', '[noise_switch_rollouts] %.1f min -> %s\n', settings.minutes, mat_file);
end

%% ------------------------------------------------------------------------
function c = rollout_copy(P, cname, m)
% One copy: build from the shared seeds, run with the shared noise seed and
% this copy's seed from the switch on, and reduce to the block series.
rng(0, 'twister');                   % workers default to threefry; a seed must mean one network
args = {'rng_seeds', P.rng_seeds, 'T_range', P.T_range, 'fs', P.fs, ...
    'ode_solver', P.ode_solver, 'noise_seed', P.noise_seed, ...
    'noise_seed_schedule', [P.t_shared, P.rollout_seeds(m)], ...
    'plot_deci', P.store_deci, 'store_full_state', true, 'verbose', P.verbose};
if m == 1
    args = [args, {'lya_method', 'topk', 'lya_K', 1, 'lya_K_auto', false, ...
        'lya_dt', P.lya_dt, 'lya_T_interval', P.lya_T_interval, 'lya_warmup', P.lya_warmup}];
else
    args = [args, {'lya_method', 'none'}];
end
model = build_from_preset(P.preset_name, cname, args{:}, P.model_args{:});
if model.sigma_u_noise == 0
    error('run_noise_switch_rollouts:NoNoise', ...
        'Preset ''%s'' has sigma_u_noise = 0: there is no noise to share or to stop sharing.', ...
        P.preset_name);
end
params = model.get_params();
if ~SRNNCellTypePairs.routes_identical(params)
    error('run_noise_switch_rollouts:RoutesDiffer', ...
        ['Preset ''%s'' (%s) carries different synaptic dynamics on different routes; ' ...
         'the theta and b blocks are one series per neuron and need them identical.'], ...
        P.preset_name, cname);
end
model.run();

c = struct();
c.t = model.plot_data.t;
[c.series, c.neurons, c.labels] = block_series(model.plot_data, params, model);
i_switch = round((P.t_shared - P.T_range(1)) * P.fs) + 1;
c.S_switch = model.S_out(i_switch, :);
c.lambda_1 = NaN;
if m == 1, c.lambda_1 = model.lya_results.LLE; end
c.n = model.n;
c.N_sys_eqs = model.N_sys_eqs;
c.type_names = model.cell_type_names;
c.n_per_type = model.n_per_type;
c.type_indices = model.type_indices;
c.sigma_u_noise = model.sigma_u_noise;
c.shown = cell(1, numel(c.type_names));      % network indices of the plotted neurons
for q = 1:numel(c.type_names)
    want = P.n_show(min(q, 2));              % E count for type 1, I count for the rest
    c.shown{q} = c.type_indices{q}(1:min(want, c.n_per_type(q)));
end
end

%% ------------------------------------------------------------------------
function [series, neurons, labels] = block_series(pd, params, model)
% The five blocks as neurons x time matrices, rows in network order, with
% the network index of each row. Every quantity is the model's own: the
% synaptic output and the route arrays from plot_data, the depression
% factor from std_route_factor. A block no neuron has is left out.
names = model.cell_type_names;
idx   = model.type_indices;
C     = numel(names);
series = struct(); neurons = struct(); labels = struct();

X = cell(C, 1); R = cell(C, 1); TH = cell(C, 1);
A = cell(C, 1); B = cell(C, 1); A_idx = cell(1, C); B_idx = cell(1, C);
is_chain = false; n_b_max = 0;
for q = 1:C
    X{q} = pd.x.(names{q});
    R{q} = pd.r.(names{q});
    posts = fieldnames(pd.synaptic_output.(names{q}));
    TH{q} = pd.synaptic_output.(names{q}).(posts{1});       % route 1 == every route
    a = pd.a.(names{q});                                    % n_q x n_a x n_t, or []
    if ~isempty(a)
        A{q} = params.c_eff(q) * reshape(sum(a, 2), size(a, 1), size(a, 3));
        A_idx{q} = idx{q};
    end
    b_posts = fieldnames(pd.b.(names{q}));
    b = pd.b.(names{q}).(b_posts{1});                       % n_q x n_b x n_t, or []
    if ~isempty(b)
        post_q = find(strcmp(names, b_posts{1}), 1);
        D = SRNNCellTypePairs.std_route_factor(b, params, q, post_q);
        B{q} = reshape(D, size(b, 1), size(b, 3));
        B_idx{q} = idx{q};
        is_chain = is_chain || params.std_chain(q, post_q);
        n_b_max = max(n_b_max, size(b, 2));
    end
end
all_idx = [idx{:}];
series.x = vertcat(X{:});      neurons.x = all_idx;     labels.x = 'x';
series.r = vertcat(R{:});      neurons.r = all_idx;     labels.r = 'r';
series.theta = vertcat(TH{:}); neurons.theta = all_idx; labels.theta = 'synaptic output \theta';
if ~isempty([A_idx{:}])
    series.a = vertcat(A{:}); neurons.a = [A_idx{:}];
    labels.a = 'SFA  (c/K)\Sigma_k a_k';
end
if ~isempty([B_idx{:}])
    series.b = vertcat(B{:}); neurons.b = [B_idx{:}];
    if is_chain
        labels.b = 'b_1 (releasable pool)';
    elseif n_b_max > 1
        labels.b = 'depression \Pi_m b_m';
    else
        labels.b = 'depression b';
    end
end
end

%% ------------------------------------------------------------------------
function assert_shared_then_split(copies, t, t_shared, cname)
% The premise of the stage: one trajectory up to the switch, an ensemble
% after it. Bit for bit on every stored sample of every block and on the
% full state vector at the switch time; and the copies must differ at the
% end, or the schedule did nothing.
c1 = copies{1};
before = t <= t_shared + 1e-9;
blocks = fieldnames(c1.series);
for m = 2:numel(copies)
    cm = copies{m};
    same = isequal(cm.S_switch, c1.S_switch);
    for k = 1:numel(blocks)
        same = same && isequal(cm.series.(blocks{k})(:, before), c1.series.(blocks{k})(:, before));
    end
    if ~same
        error('run_noise_switch_rollouts:CopiesDiffer', ...
            '%s: copy %d is not bit-identical to copy 1 before the switch at %g s.', cname, m, t_shared);
    end
    if isequal(cm.series.x(:, end), c1.series.x(:, end))
        error('run_noise_switch_rollouts:CopiesIdentical', ...
            '%s: copy %d equals copy 1 at the end of the rollout; the seed switch had no effect.', cname, m);
    end
end
end

%% ------------------------------------------------------------------------
function S = ensemble_statistics(X, t, t_shared)
% The per-time statistics of one block, X = neurons x time x copies, and
% the fit to its normalized RMS difference. See the header for each field.
M = size(X, 3);
S = struct();
S.n_neurons = size(X, 1);

% Variance of the deviations from copy 1 (the variance does not depend on
% the shift): identical copies then give exactly 0, where the variance of
% the values themselves leaves the round-off of their mean.
S.variance  = mean(var(X - X(:, :, 1), 0, 3), 1);
S.rms_diff  = sqrt(2 * S.variance);
mu          = mean(X, 3);
S.mean_norm = sqrt(mean(mu .^ 2, 1));

ref      = t >= t_shared / 2 & t <= t_shared + 1e-9;
X_ref    = X(:, ref, 1);                     % every copy is this before the switch
baseline = mean(X_ref, 2);
S.scale  = sqrt(mean(var(X_ref, 0, 2)));
S.mean_fluct_norm = sqrt(mean((mu - baseline) .^ 2, 1));

% Mean over copy pairs of the correlation across neurons. With unit-norm
% centered columns z_i, sum_{i<j} z_i'z_j = (|sum_i z_i|^2 - M) / 2.
Z = X - baseline;
Z = Z - mean(Z, 1);
Z = Z ./ sqrt(sum(Z .^ 2, 1));
S.corr = (sum(sum(Z, 3) .^ 2, 1) - M) / (M * (M - 1));

scale = S.scale;
if scale == 0, scale = NaN; end              % a constant block has no scale
S.rms_diff_rel        = S.rms_diff / scale;
S.variance_rel        = S.variance / scale ^ 2;
S.mean_norm_rel       = S.mean_norm / scale;
S.mean_fluct_norm_rel = S.mean_fluct_norm / scale;

after = t > t_shared + 1e-9;
S.fit = rise_fit(t(after) - t_shared, S.rms_diff_rel(after));
end

%% ------------------------------------------------------------------------
function F = rise_fit(t, y)
% Saturation, half time and early growth of y(t), t the time since the
% switch (first sample after it onward). See the header.
F = struct('saturation', NaN, 'saturated', false, 't_half', NaN, ...
    'growth_rate', NaN, 'rise_exponent', NaN, 'fit_window', [NaN NaN], 'n_fit', 0);
n = numel(t);
if n < 10 || any(~isfinite(y)), return; end
fifth = max(1, floor(n / 5));
last  = mean(y(n - fifth + 1:n));
prev  = mean(y(n - 2 * fifth + 1:n - fifth));
F.saturation = last;
F.saturated  = abs(last - prev) < 0.1 * abs(last);

k = find(y >= 0.5 * last, 1);
if isempty(k), return; end
if k == 1
    F.t_half = t(1) * (0.5 * last) / y(1);               % y is 0 at the switch
else
    F.t_half = interp1(y(k-1:k), t(k-1:k), 0.5 * last);
end
if k >= 3 && all(y(1:k) > 0)
    p = polyfit(t(1:k), log(y(1:k)), 1);
    q = polyfit(log(t(1:k)), log(y(1:k)), 1);
    F.growth_rate   = p(1);
    F.rise_exponent = q(1);
    F.fit_window    = [t(1), t(k)];
    F.n_fit         = k;
end
end

%% ------------------------------------------------------------------------
function E = example_series(two, c1, deci)
% The series of two copies for the plotted neurons, single precision,
% subsampled by deci. Per block: n_shown x time x 2, with the cell type
% of each row.
keep = 1:deci:numel(c1.t);
shown = [c1.shown{:}];
type_of = zeros(1, c1.n);
for q = 1:numel(c1.type_indices), type_of(c1.type_indices{q}) = q; end
E = struct('t', c1.t(keep)', 'copies', [1 2], 'series', struct(), ...
    'neurons', struct(), 'type', struct());
blocks = fieldnames(c1.series);
for k = 1:numel(blocks)
    blk = blocks{k};
    [has, rows] = ismember(shown, c1.neurons.(blk));
    rows = rows(has);
    E.series.(blk)  = single(cat(3, two{1}.series.(blk)(rows, keep), two{2}.series.(blk)(rows, keep)));
    E.neurons.(blk) = shown(has);
    E.type.(blk)    = type_of(shown(has));
end
E.type_names = c1.type_names;
end

%% ------------------------------------------------------------------------
function k = whole_ratio(fs_from, fs_to, what)
% The subsampling factor from fs_from down to fs_to; it must be a whole number.
k = fs_from / fs_to;
if abs(k - round(k)) > 1e-9 || k < 1
    error('run_noise_switch_rollouts:BadRate', ...
        '%s = %g Hz is not a whole fraction of %g Hz.', what, fs_to, fs_from);
end
k = round(k);
end

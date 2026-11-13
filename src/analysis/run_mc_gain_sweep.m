function [out_dir, results_file] = run_mc_gain_sweep(opts)
% RUN_MC_GAIN_SWEEP Memory capacity and lambda_1 at every point of a gain-by-bias sweep.
%
%   out_dir = RUN_MC_GAIN_SWEEP()                          % 'essential': a first look
%   out_dir = RUN_MC_GAIN_SWEEP('run_mode', 'fast', 'output_dir', d)
%   [out_dir, results_file] = RUN_MC_GAIN_SWEEP('gain_levels', g, 'bias_levels', b, ...
%       'seeds', s, 'conditions', {'no_adaptation', 'sfa3_stdc2'})
%
% Block C (FD-010): is the memory-capacity result a property of adaptation
% or of one well-tuned network? For every adaptation condition of the
% preset, on several network seeds, this measures total memory capacity and
% the largest Lyapunov exponent over a grid of
%   gain   level_of_chaos, the multiplier on every connection weight
%   bias   S_c, the centre of the per-neuron setpoint draw (the per-neuron
%          spread sigma_S_c stays as the preset has it)
% run_memory_capacity measures the same task at the preset network only.
%
% ONE JOB is one (gain, bias, seed). It builds every condition on the same
% seeds -- the same W, input weights and input sequence, only the adaptation
% differs -- so the conditions are PAIRED exactly as run_memory_capacity
% pairs its trials (verify_shared_build asserts it once, on the first job).
% The stimulus seed follows the network seed as it does there
% (seed + seed_stim_base - seed_net_base). The input weights use the
% reservoir's fixed rng_seed_input, so they are the same in every job.
% Per condition:
%   memory capacity   by SRNN_ESN_reservoir.run_memory_capacity_lean: the
%                     task, protocol (mc_protocol: 10 s washout, 600 s
%                     training, 150 s test, delays to 15 s, fs 200, 0.3 s
%                     holds, synaptic readout) and scoring of
%                     run_memory_capacity, integrated segment by segment so
%                     no state trajectory and no n x T input table is held.
%                     Bit-identical to the existing path
%                     (test_mc_lean_path), at tens of MB instead of GB.
%   lambda_1          by a separate short run() of the same reservoir with
%                     the top-K estimator (lya_K exponents, no retry;
%                     lya_T_sec seconds, accumulated over the last
%                     lya_acc_sec): once with the memory-task input on
%                     (lle_input, the value plotted: it is the first
%                     lya_T_sec of the same input sequence) and once with no
%                     input (lle_noinput).
%   rates             each neuron's time-mean rate over the training and
%                     test windows of the memory task: the mean over E and
%                     over I neurons, and the fractions of neurons pinned
%                     (above 0.9) and silent (below 0.01), so a low memory
%                     capacity can be read as saturated, silent or chaotic.
%
% SIZES, by run_mode (explicit gain_levels / bias_levels / seeds override):
%   essential    9 gains x  4 biases x  1 seed    gains as fast, biases
%                                                [-0.6 -0.2 0.2 0.6]
%   fast         9 gains x  7 biases x  4 seeds   gains [0.1 0.15 0.22
%                                                0.33 0.5 0.7 1.0 1.4 2.0],
%                                                biases -0.6:0.2:0.6
%   medium(2)   13 gains x  9 biases x  5 seeds   biases -0.6:0.15:0.6
%   production  17 gains x 13 biases x 10 seeds   biases -0.6:0.1:0.6
% Gains run from 0.1 to 2 and always include 1.0. In medium and production
% they are log-spaced in two pieces that meet at 1 (nine and twelve equal
% log steps from 0.1 to 1, then three and four from 1 to 2: ratios 1.29 and
% 1.26 per step in medium, 1.21 and 1.19 in production), rounded to three
% significant figures:
%   medium      0.1 0.129 0.167 0.215 0.278 0.359 0.464 0.599 0.774
%               1 1.26 1.59 2
%   production  0.1 0.121 0.147 0.178 0.215 0.261 0.316 0.383 0.464 0.562
%               0.681 0.825 1 1.19 1.41 1.68 2
% Seeds are seed_net_base + 1, 2, ... The protocol is the same in every mode.
% WHY THIS GRID (FD-010, the lead's decisions 14 and 16, from the pilot of
% 2026-11-12, which ran gains 0.1 to 3 and biases -0.2 to 0.8): bias 0.8
% silenced the network and gain 3 gave no memory in any condition, so both
% ends are dropped; the no-adaptation peak lay on the pilot's low-bias edge
% (-0.2), so the bias axis is extended down to -0.6; and the preset's own
% point (gain 1, bias 0.2) is a grid point in every mode. The bias levels
% are MATLAB colon ranges, so "0.2" is the colon's 0.19999999999999996 (as
% in the main run of 2026-11-12, whose job files carry those values); the
% figures and statistics take the nearest level to the preset's bias.
%
% RESUMABLE. Each finished job writes <out_dir>/jobs/job_*.mat; a rerun into
% the same output directory skips the jobs that are there (and reruns the
% ones recorded as failed, unless retry_failed is false). A job file made
% with a different preset, protocol, integrator or model override is
% refused, not mixed in. A condition that throws is recorded as failed with
% its message, and the other conditions and jobs carry on.
%
% OUTPUT, in <output_dir>/mc_gain_sweep/ (default output_dir:
% data/mc_gain_sweep_<timestamp>/, so the default lands in data/ and the
% figures' run_dir is that dated folder):
%   mc_gain_sweep_results.mat   T (table) and settings (struct), below
%   mc_gain_sweep_results.csv   T without the R^2 curve
%   jobs/job_*.mat              one per job
%   git_provenance.txt          (+ working_changes.patch when dirty)
% T has one row per (gain, bias, seed, condition):
%   gain, bias, seed, seed_stim, condition, mc_total, horizon (s),
%   R2_delay (1 x n_delay, R^2 against delay), lle_input, lle_noinput (1/s),
%   rate_E, rate_I, frac_pinned, frac_silent, seconds (wall time of the
%   condition), failed (logical), message (the error text of a failed row).
% settings: conditions, delay_s, preset_name, preset_gain, preset_bias,
%   run_mode, gain_levels, bias_levels, seeds, protocol, ode_solver, the
%   Lyapunov settings, model_args, n, n_jobs, n_failed, minutes, git.
%
% Options beyond the grid:
%   preset_name   a SRNNCellTypePairs preset with noise off (default the
%                 mu7chain bundle's noise-free memory-capacity preset)
%   conditions    names of the preset's conditions to run (default all)
%   ode_solver    default 'sra1', as the bundle's memory capacity uses
%   model_args    name-value overrides of the preset (cell), applied after
%                 it and before the job's gain and bias
%   protocol      a struct overriding fields of mc_protocol (T_wash_sec,
%                 T_train_sec, ...). FOR TESTS: a shorter training window
%                 makes the readout under-determined at the paper's n.
%   lya_K, lya_T_sec, lya_acc_sec   the Lyapunov run (5, 30 s, 20 s)
%   segment_len   samples per integration segment of the lean path
%   retry_failed  rerun the jobs an earlier run recorded as failed (true)
%   use_parallel  false runs the jobs serially in the client (true)
%
% The pool is the session's allocation (pool_set); with none the jobs run
% serially.
%
% See also: SRNN_ESN_reservoir.run_memory_capacity_lean, mc_protocol,
%           run_memory_capacity, read_mc_gain_sweep, verify_shared_build

arguments
    opts.preset_name  (1,:) char   = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain'
    opts.run_mode     (1,:) char   = 'essential'
    opts.gain_levels  (1,:) double = []
    opts.bias_levels  (1,:) double = []
    opts.seeds        (1,:) double = []
    opts.conditions   (1,:) cell   = {}
    opts.output_dir   (1,:) char   = ''
    opts.ode_solver   (1,:) char   = 'sra1'
    opts.verbose                   = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (or a logical); see verbose_level
    opts.model_args   (1,:) cell   = {}
    opts.protocol     (1,1) struct = struct()
    opts.lya_K        (1,1) double = 5
    opts.lya_T_sec    (1,1) double = 30
    opts.lya_acc_sec  (1,1) double = 20
    opts.segment_len  (1,1) double = 600
    opts.retry_failed (1,1) logical = true
    opts.use_parallel (1,1) logical = true
end

setup_paths();
t_start = tic;

%% Sizes -- BEFORE any compute, so a bad mode fails in milliseconds
P = mc_protocol();
% The grid by mode; see SIZES in the header for the levels and why.
fast_gains = [0.1 0.15 0.22 0.33 0.5 0.7 1.0 1.4 2.0];
switch opts.run_mode
    case 'essential'
        gain_levels = fast_gains;
        bias_levels = [-0.6 -0.2 0.2 0.6];
        n_seeds = 1;
    case 'fast'
        gain_levels = fast_gains;
        bias_levels = -0.6:0.2:0.6;
        n_seeds = 4;
    case {'medium', 'medium2'}
        gain_levels = [0.1 0.129 0.167 0.215 0.278 0.359 0.464 0.599 0.774 1 1.26 1.59 2];
        bias_levels = -0.6:0.15:0.6;
        n_seeds = 5;
    case 'production'
        gain_levels = [0.1 0.121 0.147 0.178 0.215 0.261 0.316 0.383 0.464 0.562 ...
                       0.681 0.825 1 1.19 1.41 1.68 2];
        bias_levels = -0.6:0.1:0.6;
        n_seeds = 10;
    otherwise
        error('run_mc_gain_sweep:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
            opts.run_mode, strjoin(run_mode_names(), ', '));
end
seeds       = P.seed_net_base + (1:n_seeds);
if ~isempty(opts.gain_levels); gain_levels = opts.gain_levels; end
if ~isempty(opts.bias_levels); bias_levels = opts.bias_levels; end
if ~isempty(opts.seeds);       seeds       = opts.seeds;       end

%% Protocol: mc_protocol, with the (test-only) overrides
over = fieldnames(opts.protocol);
unknown = setdiff(over, fieldnames(P));
if ~isempty(unknown)
    error('run_mc_gain_sweep:BadProtocolField', ...
        'protocol has field(s) mc_protocol does not: %s.', strjoin(unknown', ', '));
end
for k = 1:numel(over)
    P.(over{k}) = opts.protocol.(over{k});
end
if opts.lya_acc_sec <= 0 || opts.lya_acc_sec >= opts.lya_T_sec
    error('run_mc_gain_sweep:BadLyapunovWindow', ...
        'lya_acc_sec must be positive and shorter than lya_T_sec.');
end

%% The preset: the network, and its conditions
[preset, model_class, conditions] = srnn_param_preset(opts.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_mc_gain_sweep:BadModelClass', ...
        'Preset ''%s'' is written for %s; SRNN_ESN_reservoir subclasses SRNNCellTypePairs.', ...
        opts.preset_name, model_class);
end
if isfield(preset, 'mu_S_c') && ~isempty(preset.mu_S_c)
    error('run_mc_gain_sweep:PinnedSetpointMean', ...
        ['Preset ''%s'' sets mu_S_c, which overrides S_c as the centre of the setpoint ' ...
         'draw, so the bias axis (S_c) would do nothing.'], opts.preset_name);
end
all_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
if isempty(opts.conditions)
    cond_names = all_names;
else
    cond_names = opts.conditions;
    bad = setdiff(cond_names, all_names);
    if ~isempty(bad)
        error('run_mc_gain_sweep:NoSuchCondition', ...
            'Preset ''%s'' has no condition %s (it has: %s).', ...
            opts.preset_name, strjoin(bad, ', '), strjoin(all_names, ', '));
    end
end
[~, where] = ismember(cond_names, all_names);
cond_args = cellfun(@(c) struct2namevalue(rmfield(c, 'name')), conditions(where), ...
    'UniformOutput', false);
n_cond = numel(cond_names);

% Noise is off for block C, and a stochastic preset would need the solver
% checked against it: the model class's own check, before any compute.
sigma = 0;
if isfield(preset, 'sigma_u_noise'); sigma = preset.sigma_u_noise; end
check_noise_settings(sigma, opts.ode_solver, 'run_mc_gain_sweep');

%% What every reservoir of the run shares
fs = P.fs;
hold_len  = max(1, round(P.T_hold * fs));
d_max     = round(P.d_max_sec * fs);
d_max_eff = max(1, floor(d_max / hold_len));
delay_s   = (1:d_max_eff) * (hold_len / fs);

J = struct();                       % the job recipe, sent to the workers
J.base_args = [struct2namevalue(preset), opts.model_args, ...
    {'fs', fs, 'ode_solver', opts.ode_solver, ...
     'input_type', P.input_type, 'u_f_cutoff', P.u_f_cutoff, 'u_alpha', P.u_alpha, ...
     'T_hold', P.T_hold, 'verbose', 'near-none'}];
J.mc_args = {'T_wash', round(P.T_wash_sec * fs), 'T_train', round(P.T_train_sec * fs), ...
    'T_test', round(P.T_test_sec * fs), 'd_max', d_max, ...
    'lya_method', 'none', 'lean_input', true};
% The Lyapunov run: lya_T_sec of the same task (the sequence's first
% seconds), top-K with a fixed small K, accumulated over the last
% lya_acc_sec. lya_T_interval is set explicitly: its default comes from
% T_range, which the reservoir's input does not follow.
lya_n = round(opts.lya_T_sec * fs);
lya_t_end = (lya_n - 1) / fs;
lya_t0 = opts.lya_T_sec - opts.lya_acc_sec;
J.lya_args = {'T_wash', round(lya_t0 * fs), 'T_train', lya_n - round(lya_t0 * fs), ...
    'T_test', 0, 'd_max', d_max, ...
    'lya_method', 'topk', 'lya_K', opts.lya_K, 'lya_K_auto', false, 'lya_dt', 0.05, ...
    'lya_T_interval', [lya_t0, lya_t_end], 'lya_warmup', min(5, lya_t0), ...
    'store_decimated_state', false, 'store_full_state', false};
J.cond_names = cond_names;
J.cond_args = cond_args;
J.readout_signal = P.readout_signal;
J.R2_threshold = P.R2_threshold_for_horizon;
J.delay_s = delay_s;
J.segment_len = opts.segment_len;
J.seed_offset = P.seed_stim_base - P.seed_net_base;
% What makes two runs the same experiment; a job file with another
% signature is refused. Run-size choices (grid, seeds) are not part of it.
J.signature = struct('preset_name', opts.preset_name, 'protocol', P, ...
    'ode_solver', opts.ode_solver, 'model_args', {opts.model_args}, ...
    'lya', [opts.lya_K, opts.lya_T_sec, opts.lya_acc_sec], 'conditions', {cond_names});

%% Output directory
project_root = fileparts(which('setup_paths'));
if isempty(opts.output_dir)
    run_root = fullfile(project_root, 'data', ...
        ['mc_gain_sweep_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
else
    run_root = opts.output_dir;
end
out_dir = fullfile(run_root, 'mc_gain_sweep');
job_dir = fullfile(out_dir, 'jobs');
if ~exist(job_dir, 'dir'); mkdir(job_dir); end

%% The job table
[G, B, S] = ndgrid(gain_levels, bias_levels, seeds);
jobs = struct('gain', num2cell(G(:)), 'bias', num2cell(B(:)), 'seed', num2cell(S(:)));
n_jobs = numel(jobs);
job_files = cell(n_jobs, 1);
todo = true(n_jobs, 1);
for k = 1:n_jobs
    job_files{k} = fullfile(job_dir, job_file_name(jobs(k)));
    if exist(job_files{k}, 'file')
        old = load(job_files{k}, 'job');
        if ~isequaln(old.job.signature, J.signature)
            error('run_mc_gain_sweep:DifferentRun', ...
                ['%s was made with a different preset, protocol, integrator, conditions or ' ...
                 'model override. Use another output_dir; results are not mixed.'], job_files{k});
        end
        todo(k) = opts.retry_failed && any(old.job.failed);
    end
end
vprintf(opts.verbose, 'minimal', ...
    '[mc_gain_sweep] preset=%s run_mode=%s: %d gains x %d biases x %d seeds = %d jobs x %d conditions (%d already done)\n', ...
    opts.preset_name, opts.run_mode, numel(gain_levels), numel(bias_levels), numel(seeds), ...
    n_jobs, n_cond, nnz(~todo));

%% Fairness check, once: the first job's conditions share their build
if any(todo) && n_cond > 1
    chk = cell(1, n_cond);
    for c = 1:n_cond
        a = job_args(J, jobs(1), c);
        chk{c} = SRNN_ESN_reservoir(a{:}, J.mc_args{:});
        chk{c}.build();
    end
    % tau_a and synapse_config are what the conditions vary; the network,
    % the input weights and the input sequence must be shared.
    verify_shared_build(chk, {'tau_a', 'synapse_config'}, {'W', 'W_in', 'u_scalar', 't_ex'});
    clear chk
end

%% Run
pool = [];
if opts.use_parallel
    pool = pool_ensure('verbose', opts.verbose);
end
n_par = 0;
if ~isempty(pool); n_par = pool.NumWorkers; end
todo_idx = find(todo);
todo_jobs = jobs(todo_idx);
todo_files = job_files(todo_idx);
n_todo = numel(todo_idx);
verbose = opts.verbose;
parfor (k = 1:n_todo, n_par)
    job = run_one_job(J, todo_jobs(k), todo_files{k});
    vprintf(verbose, 'minimal', '  job %d/%d gain %.4g bias %.4g seed %d: %.0f s, MC %s%s\n', ...
        k, n_todo, job.gain, job.bias, job.seed, sum(job.seconds), ...
        strtrim(sprintf('%.2f ', job.mc_total)), failed_text(job));
    if any(job.failed)
        first = find(job.failed, 1);
        vfail('[mc_gain_sweep] job gain %.4g bias %.4g seed %d, %s FAILED: %s\n', ...
            job.gain, job.bias, job.seed, job.conditions{first}, job.message{first});
    end
end

%% Consolidate from the job files (this run's and earlier ones')
rows = cell(n_jobs, 1);
for k = 1:n_jobs
    L = load(job_files{k}, 'job');
    rows{k} = job_rows(L.job);
end
T = vertcat(rows{:});

settings = struct();
settings.conditions  = cond_names;
settings.delay_s     = delay_s;
settings.preset_name = opts.preset_name;
% The preset's own gain and bias (with model_args applied), read off a
% reservoir constructed without a job's values.
probe = SRNN_ESN_reservoir(J.base_args{:}, J.cond_args{1}{:});
settings.preset_gain = probe.level_of_chaos;
settings.preset_bias = probe.S_c;
settings.run_mode    = opts.run_mode;
settings.gain_levels = gain_levels;
settings.bias_levels = bias_levels;
settings.seeds       = seeds;
settings.protocol    = P;
settings.ode_solver  = opts.ode_solver;
settings.model_args  = opts.model_args;
settings.lya_K       = opts.lya_K;
settings.lya_T_sec   = opts.lya_T_sec;
settings.lya_acc_sec = opts.lya_acc_sec;
settings.lya_T_interval = [lya_t0, lya_t_end];
settings.segment_len = opts.segment_len;
settings.n          = probe.n;
settings.eta         = probe.eta;               % ridge parameter (the class's)
settings.rng_seed_input = probe.rng_seed_input; % input weights: the same in every job
settings.n_jobs      = n_jobs;
settings.n_run       = n_todo;
settings.n_failed    = nnz(T.failed);
settings.n_workers   = n_par;
settings.minutes     = toc(t_start) / 60;
if verbose_level(opts.verbose) >= 2
    settings.git = capture_git_provenance(out_dir, project_root);
else
    [~, settings.git] = evalc('capture_git_provenance(out_dir, project_root)');
end

results_file = fullfile(out_dir, 'mc_gain_sweep_results.mat');
save(results_file, 'T', 'settings');
writetable(removevars(T, 'R2_delay'), fullfile(out_dir, 'mc_gain_sweep_results.csv'));

vprintf(opts.verbose, 'minimal', '[mc_gain_sweep] %d jobs (%d run now, %d condition rows failed), %.1f min -> %s\n', ...
    n_jobs, n_todo, settings.n_failed, settings.minutes, results_file);
end

%% ======================================================================
function job = run_one_job(J, spec, file)
% One (gain, bias, seed): every condition on the same seeds. Writes its own
% file, so a crash costs this job and not the run.
%
% The generator is named before anything is built: a pool worker's default
% is not the client's, and a seed must mean the same network on both.
rng(0, 'twister');
n_cond = numel(J.cond_names);
job = struct();
job.gain = spec.gain;
job.bias = spec.bias;
job.seed = spec.seed;
job.seed_stim = spec.seed + J.seed_offset;
job.conditions = J.cond_names;
job.signature = J.signature;
job.mc_total    = nan(1, n_cond);
job.horizon     = nan(1, n_cond);
job.R2_delay    = nan(n_cond, numel(J.delay_s));
job.lle_input   = nan(1, n_cond);
job.lle_noinput = nan(1, n_cond);
job.rate_E      = nan(1, n_cond);
job.rate_I      = nan(1, n_cond);
job.frac_pinned = nan(1, n_cond);
job.frac_silent = nan(1, n_cond);
job.seconds     = nan(1, n_cond);
job.failed      = false(1, n_cond);
job.message     = repmat({''}, 1, n_cond);

for c = 1:n_cond
    t_c = tic;
    try
        args = job_args(J, spec, c);

        % Memory capacity, lean path.
        esn = SRNN_ESN_reservoir(args{:}, J.mc_args{:});
        esn.build();
        [mc, r2, res] = esn.run_memory_capacity_lean('readout_signal', J.readout_signal, ...
            'segment_len', J.segment_len, 'verbose', false);
        r2 = r2(:)';
        if numel(r2) ~= numel(J.delay_s)
            error('run_mc_gain_sweep:DelayCount', 'Expected %d delays, got %d.', ...
                numel(J.delay_s), numel(r2));
        end
        % The horizon, as run_memory_capacity defines it: the longest delay
        % whose R^2 is above the threshold.
        idx = find(r2 > J.R2_threshold, 1, 'last');
        if isempty(idx); horizon = 0; else; horizon = J.delay_s(idx); end

        % Rates: per-neuron time means over the training and test windows.
        rate = res.rate_mean;
        if ~all(isfinite(rate)) || ~all(isfinite(r2))
            % compute_R2 clamps with max(0, min(1, .)), which turns a NaN
            % into 1: a trajectory that left the finite numbers would
            % otherwise score a perfect memory. It is a failed condition.
            error('run_mc_gain_sweep:NonFinite', ...
                'The trajectory or the R^2 curve is not finite (gain %g, bias %g).', ...
                spec.gain, spec.bias);
        end
        params = esn.cached_params;
        by_type = mat2cell(rate, params.n_per_type(:), 1);
        names = params.cell_type_names;
        rate_E = NaN; rate_I = NaN;
        if any(strcmp(names, 'E')); rate_E = mean(by_type{strcmp(names, 'E')}); end
        if any(strcmp(names, 'I')); rate_I = mean(by_type{strcmp(names, 'I')}); end
        esn = []; res = []; %#ok<NASGU>  release before the Lyapunov runs

        % lambda_1, input on and off: a short run() of the same reservoir.
        lle = nan(1, 2);
        extra = {{}, {'u_scale', 0, 'u_offset', 0}};
        for v = 1:2
            m = SRNN_ESN_reservoir(args{:}, J.lya_args{:}, extra{v}{:});
            m.build();
            m.run();
            lle(v) = m.lya_results.LLE;
            m = []; %#ok<NASGU>
        end
        if ~all(isfinite(lle))
            % The estimator returns NaN when the trajectory it follows is
            % not finite. (The rate can hide that: a clipping nonlinearity
            % maps a NaN potential to a finite rate.)
            error('run_mc_gain_sweep:NonFiniteLyapunov', ...
                'lambda_1 is not finite (gain %g, bias %g): the trajectory diverged.', ...
                spec.gain, spec.bias);
        end

        job.mc_total(c)    = mc;
        job.horizon(c)     = horizon;
        job.R2_delay(c, :) = r2;
        job.lle_input(c)   = lle(1);
        job.lle_noinput(c) = lle(2);
        job.rate_E(c)      = rate_E;
        job.rate_I(c)      = rate_I;
        job.frac_pinned(c) = mean(rate > 0.9);
        job.frac_silent(c) = mean(rate < 0.01);
    catch ME
        job.failed(c)  = true;
        job.message{c} = sprintf('%s: %s', ME.identifier, ME.message);
    end
    job.seconds(c) = toc(t_c);
end
save_job(file, job);
end

function args = job_args(J, spec, c)
% Constructor arguments of one job's condition c: seeds, the preset and the
% run's shared settings, the condition, then the job's gain and bias LAST
% (the constructor is last-write-wins, so the swept values override the
% preset's).
args = [{'rng_seeds', [spec.seed, spec.seed + J.seed_offset]}, J.base_args, ...
    J.cond_args{c}, {'level_of_chaos', spec.gain, 'S_c', spec.bias}];
end

function save_job(file, job)
% Written to a temporary name and renamed, so an interrupted save never
% leaves a file a rerun would take for a finished job.
tmp = [file '.tmp'];
save(tmp, 'job', '-mat');
movefile(tmp, file, 'f');
end

function name = job_file_name(spec)
% Named by the job's values, not its position, so a rerun with a different
% grid still finds the jobs it shares with an earlier one.
name = sprintf('job_g%s_b%s_s%d.mat', number_tag(spec.gain), number_tag(spec.bias), spec.seed);
end

function s = number_tag(v)
% A number as a file-name token that reads back exactly: %.17g, with the
% characters a file name should not carry replaced.
s = sprintf('%.17g', v);
s = strrep(strrep(strrep(s, '-', 'm'), '.', 'p'), '+', '');
end

function T = job_rows(job)
% One table row per condition of a job.
n = numel(job.conditions);
T = table(repmat(job.gain, n, 1), repmat(job.bias, n, 1), repmat(job.seed, n, 1), ...
    repmat(job.seed_stim, n, 1), job.conditions(:), job.mc_total(:), job.horizon(:), ...
    job.R2_delay, job.lle_input(:), job.lle_noinput(:), job.rate_E(:), job.rate_I(:), ...
    job.frac_pinned(:), job.frac_silent(:), job.seconds(:), job.failed(:), job.message(:), ...
    'VariableNames', {'gain', 'bias', 'seed', 'seed_stim', 'condition', 'mc_total', 'horizon', ...
    'R2_delay', 'lle_input', 'lle_noinput', 'rate_E', 'rate_I', 'frac_pinned', 'frac_silent', ...
    'seconds', 'failed', 'message'});
end

function s = failed_text(job)
s = '';
if any(job.failed)
    s = sprintf('  (%d condition(s) FAILED)', nnz(job.failed));
end
end

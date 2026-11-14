function [out_dirs, results_files] = run_mc_gain_sweep_set(arg_lists)
% RUN_MC_GAIN_SWEEP_SET Several gain sweeps of memory capacity, their jobs run as one pool.
%
%   [out_dirs, results_files] = RUN_MC_GAIN_SWEEP_SET({args_1, args_2, ...})
%
%   arg_lists      a cell with one entry per run: the name-value options of
%                  run_mc_gain_sweep for that run, as a cell (see its header
%                  and mc_gain_sweep_options)
%   out_dirs       each run's mc_gain_sweep folder, in the order given
%   results_files  each run's mc_gain_sweep_results.mat
%
% A RUN is one network preset swept over gain against one second axis: what
% run_mc_gain_sweep does, which is a set of one. A bundle's block C stage
% has several (a gain x bias grid, a gain x self-coupling grid, and variant
% networks along gain only; mc_gain_sweep_stage_plan), and they are run
% here together for two reasons:
%   * EVERY RUN IS PREPARED BEFORE ANY JOB COMPUTES: its options, preset,
%     conditions, second axis and existing job files are checked first, so a
%     mistake in the last run costs seconds and not the hours of the runs
%     before it.
%   * THE JOBS OF ALL RUNS GO THROUGH ONE PARFOR. A run with fewer jobs than
%     the pool has workers would otherwise leave workers idle for as long as
%     one job takes (four small runs one after another take four job-times;
%     pooled, 24 jobs on 12 workers take two).
% Each run still writes its own folder, job files and results table, exactly
% as it would alone, so it is read and resumed alone too.
%
% The runs must write to different output directories, and agree on
% use_parallel (there is one loop). Everything else is per run.
%
% ERRORS carry run_mc_gain_sweep's identifiers (run_mc_gain_sweep:...): the
% checks are that function's, wherever its runs are made.
%
% See also: run_mc_gain_sweep, mc_gain_sweep_options, mc_gain_sweep_grid,
%           run_mc_gain_sweep_stage, mc_gain_sweep_job_rows

arguments
    arg_lists (1,:) cell
end

setup_paths();
t_start = tic;

%% Prepare every run -- BEFORE any compute, so a bad run fails in seconds
n_runs = numel(arg_lists);
runs = cell(1, n_runs);
for r = 1:n_runs
    runs{r} = prepare_run(mc_gain_sweep_options(arg_lists{r}{:}));
end
out_dirs = cellfun(@(R) R.out_dir, runs, 'UniformOutput', false);
results_files = cell(1, n_runs);
if n_runs == 0; return; end
if numel(unique(out_dirs)) < n_runs
    error('run_mc_gain_sweep:SameOutputDir', ...
        'Two runs of the set write to one folder; each needs its own output_dir.');
end
use_parallel = cellfun(@(R) R.opts.use_parallel, runs);
if any(use_parallel ~= use_parallel(1))
    error('run_mc_gain_sweep:MixedParallel', ...
        'The runs of a set share one loop: use_parallel must be the same for all.');
end
verbose = runs{1}.opts.verbose;      % of the lines about the set as a whole

%% One job list: every run's jobs still to do
which_run = zeros(0, 1);
specs = struct('gain', {}, 'second', {}, 'seed', {});
files = cell(0, 1);
for r = 1:n_runs
    idx = find(runs{r}.todo);
    which_run = [which_run; repmat(r, numel(idx), 1)];   %#ok<AGROW>
    specs = [specs; runs{r}.jobs(idx)];                  %#ok<AGROW>
    files = [files; runs{r}.job_files(idx)];             %#ok<AGROW>
end
% Each job with the recipe of its own run (small: names and settings), so
% the loop slices the list and no worker is sent the recipes of other runs.
recipes = cellfun(@(R) R.J, runs, 'UniformOutput', false);
recipes = recipes(which_run);
n_todo = numel(which_run);

%% Run
pool = [];
if use_parallel(1)
    pool = pool_ensure('verbose', verbose);
end
n_par = 0;
if ~isempty(pool); n_par = pool.NumWorkers; end
if n_runs > 1
    vprintf(verbose, 'minimal', '[mc_gain_sweep] %d runs, %d jobs to run as one pool on %d workers\n', ...
        n_runs, n_todo, n_par);
end
parfor (k = 1:n_todo, n_par)
    J = recipes{k};
    job = run_one_job(J, specs(k), files{k});
    column = J.axis2.column;      % the job field and console name of the second axis
    vprintf(J.verbose, 'minimal', '  job %d/%d %sgain %.4g %s %.4g seed %d: %.0f s, MC %s%s\n', ...
        k, n_todo, J.console_tag, job.gain, column, job.(column), job.seed, sum(job.seconds), ...
        strtrim(sprintf('%.2f ', job.mc_total)), failed_text(job));
    if any(job.failed)
        first = find(job.failed, 1);
        vfail('[mc_gain_sweep] job %sgain %.4g %s %.4g seed %d, %s FAILED: %s\n', ...
            J.console_tag, job.gain, column, job.(column), job.seed, job.conditions{first}, job.message{first});
    end
end

%% Consolidate each run from its job files (this call's and earlier ones')
for r = 1:n_runs
    results_files{r} = consolidate_run(runs{r}, n_par, toc(t_start) / 60);
end
end

%% ======================================================================
function R = prepare_run(opts)
% Everything one run needs before its jobs compute: its sizes, protocol,
% preset and conditions checked, the job recipe J (sent to the workers),
% the output folder, and the job table with the jobs already on disk
% marked done. Returns the run as a struct; nothing is simulated beyond
% the one build of the fairness check.

%% Sizes -- first, so a bad mode fails in milliseconds
P = mc_protocol();
axis2 = mc_sweep_second_axis(opts.second_axis);
G = mc_gain_sweep_grid(opts.run_mode, 'second_axis', opts.second_axis, ...
    'gain_levels', opts.gain_levels, 'second_levels', given_second_levels(axis2, opts), ...
    'seeds', opts.seeds);
gain_levels = G.gain_levels;
second_levels = G.second_levels;
seeds = G.seeds;

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
if opts.lya_acc_sec >= opts.lya_T_sec
    error('run_mc_gain_sweep:BadLyapunovWindow', ...
        'lya_acc_sec must be shorter than lya_T_sec.');
end

%% The preset: the network, and its conditions
[preset, model_class, conditions] = srnn_param_preset(opts.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_mc_gain_sweep:BadModelClass', ...
        'Preset ''%s'' is written for %s; SRNN_ESN_reservoir subclasses SRNNCellTypePairs.', ...
        opts.preset_name, model_class);
end
check_second_axis(axis2, preset, opts.preset_name);
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
[delay_s, ~, d_max] = mc_delay_axis(P);

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
LYA_DT = 0.05;          % s between re-orthonormalisations, as the sweeps use
LYA_WARMUP_MAX = 5;     % s the basis aligns before accumulating (the class default)
lya_n = round(opts.lya_T_sec * fs);
lya_t_end = (lya_n - 1) / fs;
lya_t0 = opts.lya_T_sec - opts.lya_acc_sec;
J.lya_args = {'T_wash', round(lya_t0 * fs), 'T_train', lya_n - round(lya_t0 * fs), ...
    'T_test', 0, 'd_max', d_max, ...
    'lya_method', 'topk', 'lya_K', opts.lya_K, 'lya_K_auto', false, 'lya_dt', LYA_DT, ...
    'lya_T_interval', [lya_t0, lya_t_end], 'lya_warmup', min(LYA_WARMUP_MAX, lya_t0), ...
    'store_decimated_state', false, 'store_full_state', false};
J.cond_names = cond_names;
J.cond_args = cond_args;
J.readout_signal = P.readout_signal;
J.R2_threshold = P.R2_threshold_for_horizon;
J.delay_s = delay_s;
J.segment_len = opts.segment_len;
J.seed_offset = P.seed_stim_base - P.seed_net_base;
J.pinned_rate = 0.9;    % a neuron whose mean rate is above this counts as pinned
J.silent_rate = 0.01;   % ... and below this as silent
% What makes two runs the same experiment; a job file with another
% signature is refused. Run-size choices (grid, seeds) are not part of it.
% The second axis is part of it, as a field ONLY when it is not the bias:
% the signature of a bias run is then the one its job files have carried
% since before the axis could be chosen, so they still resume.
J.signature = struct('preset_name', opts.preset_name, 'protocol', P, ...
    'ode_solver', opts.ode_solver, 'model_args', {opts.model_args}, ...
    'lya', [opts.lya_K, opts.lya_T_sec, opts.lya_acc_sec], 'conditions', {cond_names});
if ~axis2.is_default
    J.signature.second_axis = axis2.property;
end
J.axis2 = axis2;
J.verbose = opts.verbose;

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
% How a job's console line names its run: the folder the run writes under.
[~, run_name] = fileparts(run_root);
J.console_tag = sprintf('[%s] ', run_name);

%% The job table
[Gg, Bb, Ss] = ndgrid(gain_levels, second_levels, seeds);
jobs = struct('gain', num2cell(Gg(:)), 'second', num2cell(Bb(:)), 'seed', num2cell(Ss(:)));
n_jobs = numel(jobs);
job_files = cell(n_jobs, 1);
todo = true(n_jobs, 1);
for k = 1:n_jobs
    job_files{k} = fullfile(job_dir, job_file_name(jobs(k), axis2));
    if exist(job_files{k}, 'file')
        old = load(job_files{k}, 'job');
        if ~isequaln(old.job.signature, J.signature)
            error('run_mc_gain_sweep:DifferentRun', ...
                ['%s was made with a different preset, protocol, integrator, conditions, ' ...
                 'second axis or model override. Use another output_dir; results are not mixed.'], job_files{k});
        end
        todo(k) = opts.retry_failed && any(old.job.failed);
    end
end
vprintf(opts.verbose, 'minimal', ...
    '[mc_gain_sweep] preset=%s run_mode=%s: %d gains x %d %s levels x %d seeds = %d jobs x %d conditions (%d already done)\n', ...
    opts.preset_name, opts.run_mode, numel(gain_levels), numel(second_levels), axis2.column, numel(seeds), ...
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
end

R = struct('opts', opts, 'P', P, 'axis2', axis2, 'J', J, 'out_dir', out_dir, ...
    'project_root', project_root, 'gain_levels', gain_levels, 'second_levels', second_levels, ...
    'seeds', seeds, 'cond_names', {cond_names}, 'delay_s', delay_s, ...
    'lya_T_interval', [lya_t0, lya_t_end], 'jobs', jobs, 'job_files', {job_files}, 'todo', todo);
end

function results_file = consolidate_run(R, n_workers, minutes)
% One run's table and settings, from every job file of its table (this
% call's and earlier ones'), written beside them.
opts = R.opts;
axis2 = R.axis2;
n_jobs = numel(R.jobs);
rows = cell(n_jobs, 1);
for k = 1:n_jobs
    L = load(R.job_files{k}, 'job');
    rows{k} = mc_gain_sweep_job_rows(L.job);
end
T = vertcat(rows{:});

settings = struct();
settings.conditions  = R.cond_names;
settings.delay_s     = R.delay_s;
settings.preset_name = opts.preset_name;
% The second axis's settings are named by its column: preset_bias and
% bias_levels for S_c, preset_<property> and <property>_levels otherwise.
settings.second_axis = axis2.property;
[settings.preset_gain, settings.(['preset_' axis2.column])] = ...
    mc_gain_sweep_preset_point(opts.preset_name, opts.model_args, axis2.property);
settings.run_mode    = opts.run_mode;
settings.gain_levels = R.gain_levels;
settings.([axis2.column '_levels']) = R.second_levels;
settings.seeds       = R.seeds;
settings.protocol    = R.P;
settings.ode_solver  = opts.ode_solver;
settings.model_args  = opts.model_args;
settings.lya_K       = opts.lya_K;
settings.lya_T_sec   = opts.lya_T_sec;
settings.lya_acc_sec = opts.lya_acc_sec;
settings.lya_T_interval = R.lya_T_interval;
settings.segment_len = opts.segment_len;
% What every job's reservoir shares, read off one constructed without a
% job's values.
probe = SRNN_ESN_reservoir(R.J.base_args{:}, R.J.cond_args{1}{:});
settings.n           = probe.n;
settings.eta         = probe.eta;               % ridge parameter (the class's)
settings.rng_seed_input = probe.rng_seed_input; % input weights: the same in every job
settings.n_jobs      = n_jobs;
settings.n_run       = nnz(R.todo);
settings.n_failed    = nnz(T.failed);
settings.n_workers   = n_workers;
settings.minutes     = minutes;
if verbose_level(opts.verbose) >= 2
    settings.git = capture_git_provenance(R.out_dir, R.project_root);
else
    [~, settings.git] = evalc('capture_git_provenance(R.out_dir, R.project_root)');
end

results_file = fullfile(R.out_dir, 'mc_gain_sweep_results.mat');
save(results_file, 'T', 'settings');
writetable(removevars(T, 'R2_delay'), fullfile(R.out_dir, 'mc_gain_sweep_results.csv'));

vprintf(opts.verbose, 'minimal', '[mc_gain_sweep] %d jobs (%d run now, %d condition rows failed), %.1f min -> %s\n', ...
    n_jobs, settings.n_run, settings.n_failed, minutes, results_file);
end

function second_levels = given_second_levels(axis2, opts)
% The levels of the second axis the caller gave: second_levels, or for the
% bias its own option bias_levels. Empty when neither (the size table's).
if ~isempty(opts.second_levels) && ~isempty(opts.bias_levels)
    error('run_mc_gain_sweep:TwoLevelOptions', ...
        'Give second_levels or bias_levels (its name for the default axis), not both.');
end
if ~axis2.is_default && ~isempty(opts.bias_levels)
    error('run_mc_gain_sweep:BiasLevelsOnOtherAxis', ...
        'bias_levels are levels of S_c; the second axis is %s. Use second_levels.', axis2.property);
end
second_levels = [opts.second_levels, opts.bias_levels];   % at most one is non-empty
end

function check_second_axis(axis2, preset, preset_name)
% The second axis must be a property the reservoir lets a caller set, other
% than the ones a job already sets; and the preset must not disable it.
RESERVED = {'level_of_chaos', 'rng_seeds'};
meta = ?SRNN_ESN_reservoir;
prop = meta.PropertyList(strcmp({meta.PropertyList.Name}, axis2.property));
if isempty(prop) || ~isequal(prop.SetAccess, 'public') || prop.Constant || ...
        (prop.Dependent && isempty(prop.SetMethod)) || ismember(axis2.property, RESERVED)
    error('run_mc_gain_sweep:BadSecondAxis', ...
        ['second_axis ''%s'' is not a publicly settable property of SRNN_ESN_reservoir ' ...
         '(or is one a job sets itself: %s).'], axis2.property, strjoin(RESERVED, ', '));
end
if strcmp(axis2.property, 'S_c') && isfield(preset, 'mu_S_c') && ~isempty(preset.mu_S_c)
    error('run_mc_gain_sweep:PinnedSetpointMean', ...
        ['Preset ''%s'' sets mu_S_c, which overrides S_c as the centre of the setpoint ' ...
         'draw, so the bias axis (S_c) would do nothing.'], preset_name);
end
end

function job = run_one_job(J, spec, file)
% One (gain, second-axis value, seed): every condition on the same seeds.
% The job struct carries the second-axis value under the axis's column name
% (job.bias for S_c; mc_sweep_second_axis). Writes its own
% file, so a crash costs this job and not the run.
%
% The generator is named before anything is built: a pool worker's default
% is not the client's, and a seed must mean the same network on both.
rng(0, 'twister');
n_cond = numel(J.cond_names);
job = struct();
column = J.axis2.column;
job.gain = spec.gain;
job.(column) = spec.second;
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
                'The trajectory or the R^2 curve is not finite (gain %g, %s %g).', ...
                spec.gain, column, spec.second);
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
                'lambda_1 is not finite (gain %g, %s %g): the trajectory diverged.', ...
                spec.gain, column, spec.second);
        end

        job.mc_total(c)    = mc;
        job.horizon(c)     = horizon;
        job.R2_delay(c, :) = r2;
        job.lle_input(c)   = lle(1);
        job.lle_noinput(c) = lle(2);
        job.rate_E(c)      = rate_E;
        job.rate_I(c)      = rate_I;
        job.frac_pinned(c) = mean(rate > J.pinned_rate);
        job.frac_silent(c) = mean(rate < J.silent_rate);
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
% run's shared settings, the condition, then the job's gain and its value on
% the second axis LAST (the constructor is last-write-wins, so the swept
% values override the preset's).
args = [{'rng_seeds', [spec.seed, spec.seed + J.seed_offset]}, J.base_args, ...
    J.cond_args{c}, {'level_of_chaos', spec.gain, J.axis2.property, spec.second}];
end

function save_job(file, job)
% Written to a temporary name and renamed, so an interrupted save never
% leaves a file a rerun would take for a finished job.
tmp = [file '.tmp'];
save(tmp, 'job', '-mat');
movefile(tmp, file, 'f');
end

function name = job_file_name(spec, axis2)
% Named by the job's values, not its position, so a rerun with a different
% grid still finds the jobs it shares with an earlier one. The second axis's
% token is 'b' for the bias and the property's name otherwise.
name = sprintf('job_g%s_%s%s_s%d.mat', number_tag(spec.gain), axis2.file_token, ...
    number_tag(spec.second), spec.seed);
end

function s = number_tag(v)
% A number as a file-name token that reads back exactly: %.17g, with the
% characters a file name should not carry replaced.
s = sprintf('%.17g', v);
s = strrep(strrep(strrep(s, '-', 'm'), '.', 'p'), '+', '');
end

function s = failed_text(job)
s = '';
if any(job.failed)
    s = sprintf('  (%d condition(s) FAILED)', nnz(job.failed));
end
end

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
% A RUN is one network preset swept over a grid of two axes (the gain
% against a second axis by default; mc_sweep_axis): what run_mc_gain_sweep
% does, which is a set of one. A bundle's block C stage
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
specs = struct('first', {}, 'second', {}, 'seed', {});
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
    point = point_text(J, specs(k));
    vprintf(J.verbose, 'minimal', '  job %d/%d %s%s seed %d: %.0f s, MC %s%s\n', ...
        k, n_todo, J.console_tag, point, job.seed, sum(job.seconds), ...
        strtrim(sprintf('%.2f ', job.mc_total)), failed_text(job));
    if any(job.failed)
        first = find(job.failed, 1);
        vfail('[mc_gain_sweep] job %s%s seed %d, %s FAILED: %s\n', ...
            J.console_tag, point, job.seed, job.conditions{first}, job.message{first});
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
G = mc_gain_sweep_grid(opts.run_mode, 'first_axis', opts.first_axis, 'second_axis', opts.second_axis, ...
    'first_levels', opts.first_levels, 'gain_levels', opts.gain_levels, ...
    'second_levels', opts.second_levels, 'bias_levels', opts.bias_levels, 'seeds', opts.seeds);
axis1 = G.first_axis;
axis2 = G.second_axis;
first_levels = G.first_levels;
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
% The readout: mc_protocol's unless the run names another. It is a field of
% the protocol, and so of the settings and of every job's signature: a run
% on the default readout carries the protocol, and the signature, it always
% has (its job files still resume), and a run on another readout is a
% different experiment that is never mixed with it.
if ~isempty(opts.readout_signal)
    P.readout_signal = opts.readout_signal;
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
check_axis(axis2, axis1, preset, opts.preset_name);   % the second first: naming the gain there is
check_axis(axis1, axis2, preset, opts.preset_name);   % the likelier mistake, and its error says so
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
% Each axis is part of it, as a field (first_axis, second_axis) ONLY when it
% is not its default: the signature of a gain x bias run is then the one
% its job files have carried since before the axes could be chosen, so they
% still resume.
J.signature = struct('preset_name', opts.preset_name, 'protocol', P, ...
    'ode_solver', opts.ode_solver, 'model_args', {opts.model_args}, ...
    'lya', [opts.lya_K, opts.lya_T_sec, opts.lya_acc_sec], 'conditions', {cond_names});
for A = [axis1, axis2]
    if ~A.is_default
        J.signature.(A.signature_field) = A.property;
    end
end
J.axis1 = axis1;
J.axis2 = axis2;
J.verbose = opts.verbose;

%% Every model of the run is valid -- before any compute or any folder
% Each grid point's model, in every condition, constructed with the
% point's own values on both axes: the constructor validates the model
% (among other things its noise against the integrator), so an axis that
% changes what is legal -- a noise level, say -- is checked at every level
% it takes, not once at the preset's. Seed-independent, so the first seed
% stands for all; a few ms a model.
[P1, P2] = ndgrid(first_levels, second_levels);
check_job_models(J, struct('first', num2cell(P1(:)), 'second', num2cell(P2(:)), 'seed', seeds(1)));

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
[F1, F2, Ss] = ndgrid(first_levels, second_levels, seeds);
jobs = struct('first', num2cell(F1(:)), 'second', num2cell(F2(:)), 'seed', num2cell(Ss(:)));
n_jobs = numel(jobs);
job_files = cell(n_jobs, 1);
todo = true(n_jobs, 1);
for k = 1:n_jobs
    job_files{k} = fullfile(job_dir, job_file_name(jobs(k), axis1, axis2));
    if exist(job_files{k}, 'file')
        old = load(job_files{k}, 'job');
        if ~isequaln(old.job.signature, J.signature)
            error('run_mc_gain_sweep:DifferentRun', ...
                ['%s was made with a different preset, protocol, integrator, conditions, ' ...
                 'axes or model override. Use another output_dir; results are not mixed.'], job_files{k});
        end
        todo(k) = opts.retry_failed && any(old.job.failed);
    end
end
vprintf(opts.verbose, 'minimal', ...
    '[mc_gain_sweep] preset=%s run_mode=%s: %d %s x %d %s levels x %d seeds = %d jobs x %d conditions (%d already done)\n', ...
    opts.preset_name, opts.run_mode, numel(first_levels), axis1.column, numel(second_levels), axis2.column, ...
    numel(seeds), n_jobs, n_cond, nnz(~todo));

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

R = struct('opts', opts, 'P', P, 'axis1', axis1, 'axis2', axis2, 'J', J, 'out_dir', out_dir, ...
    'project_root', project_root, 'first_levels', first_levels, 'second_levels', second_levels, ...
    'seeds', seeds, 'cond_names', {cond_names}, 'delay_s', delay_s, ...
    'lya_T_interval', [lya_t0, lya_t_end], 'jobs', jobs, 'job_files', {job_files}, 'todo', todo);
end

function results_file = consolidate_run(R, n_workers, minutes)
% One run's table and settings, from every job file of its table (this
% call's and earlier ones'), written beside them.
opts = R.opts;
axis1 = R.axis1;
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
% Each axis's settings are named by its column (mc_sweep_axis): preset_gain
% and gain_levels for the gain, preset_bias and bias_levels for S_c,
% preset_<property> and <property>_levels otherwise.
settings.first_axis  = axis1.property;
settings.second_axis = axis2.property;
[settings.(['preset_' axis1.column]), settings.(['preset_' axis2.column])] = ...
    mc_gain_sweep_preset_point(opts.preset_name, opts.model_args, axis2.property, axis1.property);
settings.run_mode    = opts.run_mode;
settings.([axis1.column '_levels']) = R.first_levels;
settings.([axis2.column '_levels']) = R.second_levels;
settings.seeds       = R.seeds;
settings.protocol    = R.P;
settings.readout_signal = R.P.readout_signal;   % what the memory readout saw (also in protocol)
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

function check_axis(A, other, preset, preset_name)
% An axis must be a property the reservoir lets a caller set, other than
% the one a job sets itself (rng_seeds) and the other axis; and the preset
% must not disable it.
meta = ?SRNN_ESN_reservoir;
prop = meta.PropertyList(strcmp({meta.PropertyList.Name}, A.property));
reserved = {'rng_seeds', other.property};
if isempty(prop) || ~isequal(prop.SetAccess, 'public') || prop.Constant || ...
        (prop.Dependent && isempty(prop.SetMethod)) || ismember(A.property, reserved)
    error(['run_mc_gain_sweep:Bad' capitalized(A.position) 'Axis'], ...
        ['%s ''%s'' is not a publicly settable property of SRNN_ESN_reservoir ' ...
         '(or is one a job sets itself: %s).'], A.signature_field, A.property, strjoin(reserved, ', '));
end
if strcmp(A.property, 'S_c') && isfield(preset, 'mu_S_c') && ~isempty(preset.mu_S_c)
    error('run_mc_gain_sweep:PinnedSetpointMean', ...
        ['Preset ''%s'' sets mu_S_c, which overrides S_c as the centre of the setpoint ' ...
         'draw, so the bias axis (S_c) would do nothing.'], preset_name);
end
end

function check_job_models(J, points)
% Construct the model of every grid point in every condition, as its job
% will, and turn the first constructor that refuses into an error of the
% stage that says which point it is.
for k = 1:numel(points)
    for c = 1:numel(J.cond_names)
        args = job_args(J, points(k), c);
        try
            SRNN_ESN_reservoir(args{:}, J.mc_args{:});
        catch ME
            error('run_mc_gain_sweep:BadJobModel', ...
                'The model at %s, condition %s, is refused: %s (%s)', ...
                point_text(J, points(k)), J.cond_names{c}, ME.message, ME.identifier);
        end
    end
end
end

function s = point_text(J, spec)
% A grid point in words, by the axes' column names: 'gain 0.5 bias 0.2'.
s = sprintf('%s %.4g %s %.4g', J.axis1.column, spec.first, J.axis2.column, spec.second);
end

function s = capitalized(s)
s(1) = upper(s(1));
end

function job = run_one_job(J, spec, file)
% One grid point (first-axis value, second-axis value) and seed: every
% condition on the same seeds. The job struct carries each axis's value
% under the axis's column name (job.gain and job.bias for the default axes;
% mc_sweep_axis). Writes its own file, so a crash costs this job and not
% the run.
%
% The generator is named before anything is built: a pool worker's default
% is not the client's, and a seed must mean the same network on both.
rng(0, 'twister');
n_cond = numel(J.cond_names);
job = struct();
job.(J.axis1.column) = spec.first;
job.(J.axis2.column) = spec.second;
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
                'The trajectory or the R^2 curve is not finite (%s).', point_text(J, spec));
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
                'lambda_1 is not finite (%s): the trajectory diverged.', point_text(J, spec));
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
% run's shared settings, the condition, then the job's values on the two
% axes LAST (the constructor is last-write-wins, so the swept values
% override the preset's). A property on neither axis, the gain included
% when it is not one, keeps the preset's value.
args = [{'rng_seeds', [spec.seed, spec.seed + J.seed_offset]}, J.base_args, ...
    J.cond_args{c}, {J.axis1.property, spec.first, J.axis2.property, spec.second}];
end

function save_job(file, job)
% Written to a temporary name and renamed, so an interrupted save never
% leaves a file a rerun would take for a finished job.
tmp = [file '.tmp'];
save(tmp, 'job', '-mat');
movefile(tmp, file, 'f');
end

function name = job_file_name(spec, axis1, axis2)
% Named by the job's values, not its position, so a rerun with a different
% grid still finds the jobs it shares with an earlier one. Each axis's
% token is its file_token (mc_sweep_axis): 'g' for the gain and 'b' for the
% bias, so a gain x bias job is job_g<gain>_b<bias>_s<seed>.mat as always,
% and the property's name for any other axis.
name = sprintf('job_%s%s_%s%s_s%d.mat', axis1.file_token, number_tag(spec.first), ...
    axis2.file_token, number_tag(spec.second), spec.seed);
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

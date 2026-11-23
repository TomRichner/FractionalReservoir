function detail = run_outlier_spectra_stage(cfg, run_dir, plan)
% RUN_OUTLIER_SPECTRA_STAGE A bundle's outlier_spectra stage: one file per (network, condition, seed) and an index, into <run_dir>/outlier_spectra.
%
%   detail = RUN_OUTLIER_SPECTRA_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (unused beyond what the plan already holds)
%   run_dir  the bundle's run directory
%   plan     from outlier_spectra_stage_plan(cfg)
%   detail   the folder written, and what was run, on one line
%
% The run half of the opt-in stage 'outlier_spectra' (paper_opt_in_stages;
% FD-040). Every run of the plan goes through ospec_run into
% <run_dir>/<plan.folder>/ospec_<network>_<condition>_s<a>-<b>.mat, and
% ospec_index.mat lists them (table runs; ospec_read_index documents its
% columns, ospec_load_run the contents of a run's file).
%
% PARALLEL OVER RUNS. Whole runs are spread over the session's pool
% (pool_ensure), one run to a worker at a time, the long ones handed out
% first; inside a run everything is serial. With no pool (allocation 0)
% the runs go one after another in the client. A run is the same either
% way up to the round-off of its linear algebra (a worker has one thread):
% the model is built on the twister whatever the worker's generator is, and
% the null draws are seeded per state (ospec_run).
%
% A RUN THAT FAILS does not stop the others: its error is kept in the
% index (ok false, message), the index is written, and THEN the stage
% errors, naming the runs that failed, so the entry point records the stage
% as failed while the finished runs stay readable.
%
% A RERUN INTO THE SAME FOLDER, like the other stages, starts over: the old
% index is removed first, every run of the plan is computed again and its
% file overwritten. Files of an earlier plan that the new one does not
% name are left on disk and are not in the index. (The bundle's entry
% point refuses a non-empty run directory before any stage is reached.)
%
% MEMORY per worker and TIME: see ospec_run and FD-040's report. Estimated,
% not measured at full size: about 3.5 GB for a 200 s run with multiple
% timescales.
%
% See also: outlier_spectra_stage_plan, ospec_run, ospec_read_index,
%           ospec_load_run, run_all_paper_analyses

arguments
    cfg     (1,1) struct %#ok<INUSA>
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = plan.verbose;
out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
index_file = fullfile(out_dir, 'ospec_index.mat');
if isfile(index_file); delete(index_file); end

runs = plan.runs;
n_runs = numel(runs);
args = cell(1, n_runs);
for k = 1:n_runs
    r = runs(k);
    a = {'preset_name', plan.preset_name, 'T_settle', plan.T_settle, 'T_measure', r.window_s, ...
        'fs', plan.fs, 'dt_sample', r.dt_sample, 'n', plan.n_override, ...
        'n_top', plan.n_top, 'n_keep', plan.n_keep, 'n_cal', plan.n_cal, ...
        'n_null', plan.n_null, 'n_shuffle', plan.n_shuffle, ...
        'dt_vectors', plan.dt_vectors, 'dt_mode', plan.dt_mode, 'k_full', plan.k_full, ...
        'x_thr', plan.x_thr, 'lya_K', plan.lya_K, 'lya_window', plan.lya_window_s, ...
        'chunk_s', plan.chunk_s, 'out_dir', out_dir, ...
        'verbose', verbose_level(verbose) >= verbose_level('verbose')};
    % The bundle's integrator (plan.ode_solver; '' or absent = ospec_build_model's rule).
    if isfield(plan, 'ode_solver') && ~isempty(plan.ode_solver)
        a = [a, {'ode_solver', plan.ode_solver}]; %#ok<AGROW>
    end
    if r.full
        a = [a, {'dt_full', plan.dt_full, 'dense', true}]; %#ok<AGROW>
    else
        a = [a, {'dt_full', 0, 'dense', false}]; %#ok<AGROW>
    end
    args{k} = a;
end

pool = pool_ensure();
n_workers = 0;
if ~isempty(pool); n_workers = pool.NumWorkers; end
% The windows as the runs have them (only the long networks differ from
% plan.window_s, FD-057), and the networks with the full-Jacobian passes.
other = runs([runs.window_s] ~= plan.window_s);
window_txt = '';
if ~isempty(other)
    [w_other, i_first] = unique([other.window_s], 'stable');
    parts = arrayfun(@(k) sprintf('%g s for %s', w_other(k), strjoin(unique({other([other.window_s] == w_other(k)).network}, ...
        'stable'), ', ')), 1:numel(i_first), 'UniformOutput', false);
    window_txt = [' (' strjoin(parts, '; ') ')'];
end
vprintf(verbose, 'minimal', ...
    '  outlier_spectra: %d run(s) = %d network(s) x %d condition(s) x %d seed(s); window %g s%s; full Jacobian for %s; %d worker(s)\n', ...
    n_runs, numel(plan.networks), numel(plan.conditions), plan.n_seeds, plan.window_s, window_txt, ...
    strjoin_or(unique({runs([runs.full]).network}, 'stable'), 'none'), n_workers);

% The long runs first: with one run handed out at a time, a worker that
% finishes a short run takes the next one, and no long run is left for last.
cost = [runs.window_s] ./ [runs.dt_sample] + 50 * [runs.full] .* [runs.window_s];
[~, order] = sort(cost, 'descend');

ok = false(1, n_runs); seconds = nan(1, n_runs); message = repmat({''}, 1, n_runs);
t_stage = tic;
if n_workers == 0
    for q = 1:n_runs
        k = order(q);
        [ok(k), seconds(k), message{k}] = one_run(runs(k), args{k}, verbose);
    end
else
    runs_o = runs(order); args_o = args(order);
    ok_o = false(1, n_runs); seconds_o = nan(1, n_runs); message_o = cell(1, n_runs);
    one_at_a_time = parforOptions(pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
    parfor (q = 1:n_runs, one_at_a_time)
        [ok_o(q), seconds_o(q), message_o{q}] = one_run(runs_o(q), args_o{q}, verbose);
    end
    ok(order) = ok_o; seconds(order) = seconds_o; message(order) = message_o;
end

%% The index: what the plan asked for, and what is on disk
seed = cat(1, runs.seed);
runs_table = table({runs.network}', {runs.condition}', seed(:, 1), seed(:, 2), {runs.file}', ...
    [runs.full]', [runs.full]', [runs.window_s]', [runs.dt_sample]', ok', seconds', message', ...
    'VariableNames', {'network', 'condition', 'seed1', 'seed2', 'file', 'has_full_series', 'has_dense', ...
    'window_s', 'dt_sample', 'ok', 'seconds', 'message'});
settings = rmfield(plan, 'runs');
settings.n_workers = n_workers;
settings.minutes = toc(t_stage) / 60;
save_index(index_file, runs_table, settings);

if ~all(ok)
    bad = find(~ok);
    error('run_outlier_spectra_stage:RunsFailed', ...
        '%d of %d run(s) of the outlier_spectra stage failed (the index %s lists them all). The first, %s: %s', ...
        numel(bad), n_runs, index_file, runs(bad(1)).file, message{bad(1)});
end
detail = sprintf('%s (%d runs: %d networks x %d conditions x %d seeds; %.1f min)', ...
    out_dir, n_runs, numel(plan.networks), numel(plan.conditions), plan.n_seeds, settings.minutes);
end

%% ------------------------------------------------------------------------
function [ok, seconds, message] = one_run(r, args, verbose)
% One run, with its error caught so that the others go on.
t0 = tic;
ok = false; message = '';
try
    ospec_run(r.network, r.condition, r.seed, args{:});
    ok = true;
catch err
    message = sprintf('%s: %s', err.identifier, err.message);
    fprintf(2, '  outlier_spectra: %s FAILED: %s\n', r.file, message);
end
seconds = toc(t0);
if ok
    vprintf(verbose, 'minimal', '  outlier_spectra: %s (%.0f s)\n', r.file, seconds);
end
end

function save_index(index_file, runs, settings)
save(index_file, 'runs', 'settings');
end

function s = strjoin_or(c, if_empty)
if isempty(c); s = if_empty; else; s = strjoin(c, ', '); end
end

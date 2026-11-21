function detail = run_burst_mode_stage(cfg, run_dir, plan)
%RUN_BURST_MODE_STAGE The opt-in stage burst_mode: one file per network and an index, into <run_dir>/burst_mode.
%
%   detail = RUN_BURST_MODE_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (unused beyond what the plan holds)
%   run_dir  the bundle's run directory
%   plan     from burst_mode_stage_plan(cfg)
%   detail   the folder written and what was run, on one line
%
% Every network of the plan goes through bm_run_network into
% <run_dir>/burst_mode/burst_mode_s<a>-<b>.mat (variable R; bm_run_network
% documents its fields), and burst_mode_index.mat holds `runs` (one row
% per network: bm_run_network's summary row plus file, ok, message) and
% `settings` (the plan without its runs, the number of workers, minutes).
% A rerun into the same folder starts over: the index is removed first and
% every file is written again.
%
% PARALLEL OVER NETWORKS on the session's pool (pool_ensure), one network
% to a worker at a time; serial in the client with no pool. A network that
% fails does not stop the others: the index is written, then the stage
% errors naming the first failure.
%
% COST (workstation, 2026-11-21, 2 pool workers): about 3 to 4 minutes a
% network at the default size (100 s simulated at 500 Hz with sra1;
% 10,000 frozen-gain spectra of the 50 x 50 M along the run and 2,600 along
% the frozen arm; 2,000 dense eigendecompositions of the 225 x 225 full
% Jacobian, and about 1,500 more around the bursts of the first network).
% The measured times are in the index (runs.seconds, settings.minutes).
%
% See also: burst_mode_stage_plan, bm_run_network, bm_load_runs

arguments
    cfg     (1,1) struct %#ok<INUSA>
    run_dir (1,:) char
    plan    (1,1) struct
end

out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
index_file = fullfile(out_dir, 'burst_mode_index.mat');
if isfile(index_file); delete(index_file); end

runs = plan.runs;
n_runs = numel(runs);
pool = pool_ensure();
n_workers = 0;
if ~isempty(pool); n_workers = pool.NumWorkers; end
vprintf(plan.verbose, 'minimal', '  burst_mode: %d network(s) of %s (%s), step %g at %g s; %d worker(s)\n', ...
    n_runs, plan.preset_name, plan.condition, plan.dc_level, plan.hold_dur, n_workers);

rows = cell(1, n_runs); ok = false(1, n_runs); message = repmat({''}, 1, n_runs);
t_stage = tic;
if n_workers == 0
    for k = 1:n_runs
        [rows{k}, ok(k), message{k}] = one_run(runs(k), plan, out_dir);
    end
else
    one_at_a_time = parforOptions(pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
    parfor (k = 1:n_runs, one_at_a_time)
        [rows{k}, ok(k), message{k}] = one_run(runs(k), plan, out_dir);
    end
end

template = rows(ok);
if isempty(template)
    runs_table = table();
else
    blank = template{1};
    blank{1, :} = NaN;
    for k = find(~ok); rows{k} = blank; end
    runs_table = vertcat(rows{:});
end
runs_table.file = {runs.file}';
runs_table.keep_ev = [runs.keep_ev]';
runs_table.ok = ok';
runs_table.message = message';
settings = rmfield(plan, 'runs');
settings.n_workers = n_workers;
settings.minutes = toc(t_stage) / 60;
runs = runs_table;
save(index_file, 'runs', 'settings');

if ~all(ok)
    bad = find(~ok);
    error('run_burst_mode_stage:RunsFailed', '%d of %d network(s) of the burst_mode stage failed (%s). The first, %s: %s', ...
        numel(bad), n_runs, index_file, plan.runs(bad(1)).file, message{bad(1)});
end
detail = sprintf('%s (%d networks; %.1f min)', out_dir, n_runs, settings.minutes);
end

%% ------------------------------------------------------------------------
function [row, ok, message] = one_run(r, plan, out_dir)
row = []; ok = false; message = '';
try
    R = bm_run_network(r.seed, plan, r.keep_ev);
    save_run(fullfile(out_dir, r.file), R);
    row = R.summary;
    ok = true;
    vprintf(plan.verbose, 'minimal', '  burst_mode: %s (%.0f s)\n', r.file, R.seconds);
catch err
    message = sprintf('%s: %s', err.identifier, err.message);
    fprintf(2, '  burst_mode: %s FAILED: %s\n', r.file, message);
end
end

function save_run(file, R)
save(file, 'R', '-v7.3');
end

function detail = run_ied_network_mode_stage(cfg, plan)
%RUN_IED_NETWORK_MODE_STAGE The ied_stim bundle's network_mode stage: one file per network and an index, into <run_dir>/network_mode.
%
%   detail = RUN_IED_NETWORK_MODE_STAGE(cfg)
%   detail = RUN_IED_NETWORK_MODE_STAGE(cfg, plan)    plan = ied_network_mode_stage_plan(cfg)
%
% cfg is ied_stim_config(mode). Every network of the plan goes through
% inm_run_network into <cfg.run_dir>/network_mode/ied_netmode_seed<s>.mat
% (variable N), and ied_netmode_index.mat holds `runs` (inm_run_network's
% summary rows plus file, ok, message) and `settings` (the plan without its
% runs, workers, minutes). A rerun starts over. Parallel over networks on
% the session's pool (pool_ensure); serial without one. A network that
% fails does not stop the others; the stage errors at the end naming it.
%
% It does not need the bundle's simulation step to have run: it simulates
% its two blocks itself. When the bundle's block files are in cfg.run_dir
% the re-simulation is compared with them (runs.repro_max_diff).
%
% COST, measured on the workstation (2026-11-21, 2 pool workers, 500 Hz):
% about 1.2 to 1.5 s per dense eigendecomposition of the 1190 x 1190
% Jacobian on a one-thread worker (0.43 s on the 14-thread client) and
% 0.31 s per simulated second. 'essential' (1 network, 20 s windows,
% 3 discharges): 2 x 20 + 2 x 3 x 17 = 142 eigendecompositions and two
% 80-s blocks, a few minutes; the measured time is in the index.
%
% See also: ied_network_mode_stage_plan, inm_run_network, fig_ied_network_mode

arguments
    cfg (1,1) struct
    plan struct = struct([])
end
if isempty(plan); plan = ied_network_mode_stage_plan(cfg); end
out_dir = fullfile(plan.run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
index_file = fullfile(out_dir, 'ied_netmode_index.mat');
if isfile(index_file); delete(index_file); end
runs = plan.runs; n_runs = numel(runs);
pool = pool_ensure();
n_workers = 0;
if ~isempty(pool); n_workers = pool.NumWorkers; end
fprintf('[network_mode] %d network(s), %g s windows, %d discharges each; %d worker(s)\n', ...
    n_runs, plan.window_s, plan.n_events, n_workers);
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
if any(ok)
    blank = rows{find(ok, 1)};
    blank{1, :} = NaN;
    for k = find(~ok); rows{k} = blank; end
    T = vertcat(rows{:});
else
    T = table();
end
T.file = {runs.file}'; T.ok = ok'; T.message = message';
settings = rmfield(plan, 'runs');
settings.n_workers = n_workers;
settings.minutes = toc(t_stage) / 60;
runs = T; %#ok<NASGU>
save(index_file, 'runs', 'settings');
if ~all(ok)
    bad = find(~ok);
    error('run_ied_network_mode_stage:RunsFailed', '%d of %d network(s) failed (%s). The first, %s: %s', ...
        numel(bad), n_runs, index_file, plan.runs(bad(1)).file, message{bad(1)});
end
detail = sprintf('%s (%d networks; %.1f min)', out_dir, n_runs, settings.minutes);
end

function [row, ok, message] = one_run(r, plan, out_dir)
row = []; ok = false; message = '';
try
    N = inm_run_network(r.seed, plan);
    save_one(fullfile(out_dir, r.file), N);
    row = N.summary; ok = true;
    fprintf('[network_mode] %s (%.0f s)\n', r.file, N.seconds);
catch err
    message = sprintf('%s: %s', err.identifier, err.message);
    fprintf(2, '[network_mode] %s FAILED: %s\n', r.file, message);
end
end

function save_one(file, N) %#ok<INUSD>
save(file, 'N', '-v7.3');
end

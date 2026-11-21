function detail = run_outlier_structures_stage(cfg, run_dir, plan)
% RUN_OUTLIER_STRUCTURES_STAGE A bundle's outlier_structures stage: the structural check, then one file per run and the indexes, into <run_dir>/outlier_structures.
%
%   detail = RUN_OUTLIER_STRUCTURES_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (passed on; the plan holds what is used)
%   run_dir  the bundle's run directory
%   plan     from outlier_structures_stage_plan(cfg)
%   detail   the folder written, and what was run and skipped, on one line
%
% The run half of the opt-in stage 'outlier_structures' (FD-043). In order:
%
%   1. THE STRUCTURAL CHECK, in the client, before any simulation: every
%      weight draw of plan.checks goes through ostruct_structural_check
%      (one build and one eig of W, no integration). A draw whose intended
%      eigenvalue is not where its name says is recorded with the reason,
%      and its runs are skipped, not run blindly.
%   2. THE STAGE'S OWN INDEX, ostruct_index.mat, written BEFORE any run:
%      table runs (every run of the plan: key, group, base, network,
%      condition, seed1, seed2, realization, noise_seed, file, checked: the
%      run's draw passed the check), table checks (one row per draw,
%      ostruct_structural_check's fields), settings (the plan without its
%      lists), networks (plan.networks). So a stage that fails part way
%      still says what was meant and what was skipped.
%   3. THE RUNS, through run_outlier_spectra_stage with this plan (the runs
%      of the draws that passed): the same parallel machinery (one run to
%      a pool worker at a time through pool_ensure, serial in the client
%      when the session holds no pool), ospec_run with the intended mode
%      followed and the dense pass, one file per run in the same folder,
%      and ospec_index.mat beside ostruct_index.mat with each run's ok,
%      seconds and message (ospec_read_index reads it). A run that fails
%      does not stop the others; the stage then errors
%      (run_outlier_spectra_stage:RunsFailed) after both indexes are on
%      disk. Its progress lines say "outlier_spectra": it is that code.
%
% If every draw fails the check, nothing is run and the stage errors
% (run_outlier_structures_stage:NothingToRun) after writing its index.
% A rerun into the same folder starts over, as the other stages do.
%
% ostruct_load reads the result (both indexes, the checks, every run and
% its measures).
%
% See also: outlier_structures_stage_plan, ostruct_structural_check,
%           run_outlier_spectra_stage, ostruct_load

arguments
    cfg     (1,1) struct
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = plan.verbose;
out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
index_file = fullfile(out_dir, 'ostruct_index.mat');
if isfile(index_file); delete(index_file); end

%% 1. The structural check
n_chk = numel(plan.checks);
rows = cell(1, n_chk);
t0 = tic;
for k = 1:n_chk
    rows{k} = ostruct_structural_check(plan, plan.checks(k).key, plan.checks(k).seed);
    if ~rows{k}.pass
        fprintf(2, '  outlier_structures: %s on seed [%d %d] FAILS the structural check and is skipped: %s\n', ...
            rows{k}.key, rows{k}.seed1, rows{k}.seed2, rows{k}.reason);
    end
end
checks = struct2table([rows{:}], 'AsArray', true);
vprintf(verbose, 'minimal', '  outlier_structures: structural check of %d weight draw(s) in %.0f s, %d passed\n', ...
    n_chk, toc(t0), nnz(checks.pass));

%% 2. The stage's index, before any run
runs = plan.runs;
seed = cat(1, runs.seed);
checked = false(numel(runs), 1);
for k = 1:numel(runs)
    q = strcmp(checks.key, runs(k).key) & checks.seed1 == seed(k, 1) & checks.seed2 == seed(k, 2);
    checked(k) = any(q) && checks.pass(find(q, 1));
end
runs_table = table({runs.key}', {runs.group}', {runs.base}', {runs.network}', {runs.condition}', ...
    seed(:, 1), seed(:, 2), [runs.realization]', [runs.noise_seed]', {runs.file}', checked, ...
    'VariableNames', {'key', 'group', 'base', 'network', 'condition', 'seed1', 'seed2', 'realization', ...
    'noise_seed', 'file', 'checked'});
settings = rmfield(plan, {'runs', 'checks', 'networks'});
networks = plan.networks;
save(index_file, 'runs_table', 'checks', 'settings', 'networks');

%% 3. The runs, through the outlier_spectra stage's machinery
if ~any(checked)
    error('run_outlier_structures_stage:NothingToRun', ...
        'No weight draw of the outlier_structures stage passed the structural check (%s lists why); nothing was run.', index_file);
end
inner = plan;
inner.runs = rmfield(runs(checked), {'key', 'group', 'base', 'realization', 'noise_seed'});
inner.networks = unique({inner.runs.network}, 'stable');
inner.full_networks = inner.networks;
run_outlier_spectra_stage(cfg, run_dir, inner);

n_skip = nnz(~checked);
detail = sprintf('%s (%d runs: %d networks x %d conditions x %d seeds; %d run(s) of %d draw(s) skipped by the structural check)', ...
    out_dir, nnz(checked), numel(plan.networks), numel(plan.conditions), plan.n_seeds, n_skip, nnz(~checks.pass));
end

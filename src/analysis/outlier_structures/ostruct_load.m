function S = ostruct_load(run_dir, opts)
%OSTRUCT_LOAD The outlier_structures stage of one run directory: its indexes, its structural checks, and the measures of every run. THE ONE READER of the stage.
%
%   S = ostruct_load(run_dir)
%   S = ostruct_load(run_dir, 'keep_runs', true, 'keys', {'EEx1p5'}, 'realization', 1)
%
%   run_dir   a bundle's run directory, or its outlier_structures folder
%
% Options
%   keep_runs    false: also return each run's struct R (ospec_load_run;
%                a few MB each at 500 neurons)
%   measures     true: compute each run's measures (ospec_measures and
%                ostruct_mode_measures); false reads the indexes only
%   keys         {} (all): read only these networks (keys of the plan)
%   realization  [] (all): read only these seed rows (1 = the first seed)
%
% S
%   folder     the stage's folder, absolute
%   runs       table, one row per run of the PLAN (ostruct_index.mat's
%              runs_table: key, group, base, network, condition, seed1,
%              seed2, realization, noise_seed, file, checked) joined with
%              the outlier_spectra machinery's index (ok, seconds,
%              message). A run skipped by the structural check has ok
%              false and the message 'skipped: ...'.
%   checks     table, one row per weight draw (ostruct_structural_check)
%   settings   the plan without its lists
%   networks   the plan's networks (key, name, group, target, weight_seed,
%              expect, label)
%   T          table, one row per FINISHED run read (ok and selected):
%              key, group, base, realization, noise_seed, then the fields
%              of ostruct_mode_measures' row. Empty when measures is false.
%   MM         cell, the series of each row of T (ostruct_mode_measures' mm)
%   M          cell, ospec_measures of each row of T
%   Rs         cell, the runs (only when keep_runs)
%
% Errors (ostruct_load:NoStage) when the folder holds no ostruct_index.mat.
% The run files themselves are read through ospec_read_index and
% ospec_load_run, whose header documents a run's contents.
%
% See also: run_outlier_structures_stage, ostruct_mode_measures,
%           ospec_load_run, ospec_read_index

arguments
    run_dir (1,:) char
    opts.keep_runs (1,1) logical = false
    opts.measures (1,1) logical = true
    opts.keys = {}
    opts.realization double = []
end
folder = run_dir;
if ~isfile(fullfile(folder, 'ostruct_index.mat'))
    folder = fullfile(run_dir, 'outlier_structures');
end
index_file = fullfile(folder, 'ostruct_index.mat');
if ~isfile(index_file)
    error('ostruct_load:NoStage', ...
        'No ostruct_index.mat in %s or its outlier_structures folder: the outlier_structures stage has not run there.', run_dir);
end
L = load(index_file, 'runs_table', 'checks', 'settings', 'networks');
runs = L.runs_table;
n = height(runs);
runs.ok = false(n, 1);
runs.seconds = nan(n, 1);
runs.message = repmat({''}, n, 1);
runs.message(~runs.checked) = {'skipped: the weight draw failed the structural check'};
if isfile(fullfile(folder, 'ospec_index.mat'))
    ospec = ospec_read_index(folder);
    [found, at] = ismember(runs.file, ospec.file);
    runs.ok(found) = ospec.ok(at(found));
    runs.seconds(found) = ospec.seconds(at(found));
    runs.message(found) = ospec.message(at(found));
end
S = struct('folder', folder, 'runs', runs, 'checks', L.checks, 'settings', L.settings, ...
    'networks', L.networks, 'T', table(), 'MM', {{}}, 'M', {{}}, 'Rs', {{}});

pick = runs.ok;
if ~isempty(opts.keys)
    keys = cellstr(opts.keys);
    unknown = setdiff(keys, runs.key);
    if ~isempty(unknown)
        error('ostruct_load:Key', 'The stage has no network %s; it has %s.', strjoin(unknown, ', '), ...
            strjoin(unique(runs.key, 'stable')', ', '));
    end
    pick = pick & ismember(runs.key, keys);
end
if ~isempty(opts.realization)
    pick = pick & ismember(runs.realization, opts.realization);
end
if ~opts.measures && ~opts.keep_runs; return; end

idx = find(pick);
rows = cell(1, numel(idx)); MM = cell(1, numel(idx)); M = cell(1, numel(idx)); Rs = cell(1, numel(idx));
for q = 1:numel(idx)
    k = idx(q);
    R = ospec_load_run(folder, struct('network', runs.network{k}, 'condition', runs.condition{k}, ...
        'seed', [runs.seed1(k), runs.seed2(k)]));
    if opts.measures
        M{q} = ospec_measures(R);
        [row, MM{q}] = ostruct_mode_measures(R, M{q});
        meta = struct('key', runs.key{k}, 'group', runs.group{k}, 'base', runs.base{k}, ...
            'realization', runs.realization(k), 'noise_seed', runs.noise_seed(k));
        rows{q} = merge_rows(meta, row);
    end
    if opts.keep_runs; Rs{q} = R; end
end
if opts.measures && ~isempty(rows)
    T = struct2table([rows{:}], 'AsArray', true);
    for c = {'key', 'group', 'base', 'network', 'condition', 'target'}
        if ischar(T.(c{1})); T.(c{1}) = {T.(c{1})}; end      % a single run gives char columns
    end
    S.T = T;
end
S.MM = MM; S.M = M;
if opts.keep_runs; S.Rs = Rs; end
end

%% ------------------------------------------------------------------------
function s = merge_rows(a, b)
s = a;
for f = fieldnames(b)'
    s.(f{1}) = b.(f{1});
end
end

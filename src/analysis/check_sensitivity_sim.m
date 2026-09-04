function model = check_sensitivity_sim(opts)
% CHECK_SENSITIVITY_SIM Re-simulate one grid point of a finished sensitivity run.
%
%   model = CHECK_SENSITIVITY_SIM()
%   model = CHECK_SENSITIVITY_SIM('target_level_of_chaos', 1.5, 'target_rep', 3)
%   model = CHECK_SENSITIVITY_SIM('psa_dir', d)
%
% Reloads the parameters of a single grid point from a completed 1-D sensitivity
% analysis, re-simulates that exact model, and plots its time series with
% model.plot().
%
% Motivation: at high level_of_chaos the network stays highly chaotic even under
% two-timescale STD. This lets you inspect the actual dynamics of one such run.
% It rebuilds the model exactly as ParamSpaceAnalysis2 did in run_single_job --
% same network seed, same adaptation condition, same model_defaults -- so the
% reconstructed trajectory matches the sweep.
%
% ARGUMENTS (all optional; these were base-workspace variables read with
% exist() when this was a script)
%   psa_dir                directory holding psa_object.mat. '' auto-finds the
%                          most recent level_of_chaos sensitivity run under
%                          data/param_space (searched recursively).
%   target_level_of_chaos  grid level to reload.
%   target_condition       condition name, see psa.conditions.
%   target_rep             which 'reps' grid value to use, if reps is an axis.
%   sim_T_range            override the swept window.
%   overrides              struct of extra property overrides, applied last.
%
% EITHER MODEL CLASS, since 2026-09-02. It used to construct SRNNModel2
% directly and read the condition through the n_a_E / n_b_E vocabulary, which
% meant it could no longer inspect ANY sweep the paper produces -- every one is
% SRNNCellTypePairs now. The reconstruction is delegated to
% ParamSpaceAnalysis2.rebuild_model, which ends in feval(obj.model_class, ...)
% and applies run_single_job's own precedence, so this file no longer needs to
% know either class's vocabulary. Most of what was here was deleted, not
% rewritten.
%
% ONE BEHAVIOUR CHANGE, unavoidable and arguably a fix: the old sim_c_E argument
% forced c_E = 1/3 on EVERY call, which contradicted this function's stated
% purpose of reproducing the sweep exactly. c_E does not exist on
% SRNNCellTypePairs (the equivalent is the 1 x C row `c`), so it could not be
% ported as-is. It is replaced by the general `overrides` struct, which defaults
% to empty -- so the default is now a faithful reconstruction, and exploring a
% different adaptation strength is an explicit act:
%
%   check_sensitivity_sim('overrides', struct('c', [0.5 0]))       % Pairs
%   check_sensitivity_sim('overrides', struct('c_E', 1/3))         % SRNNModel2
%
% See also: run_sensitivity_analysis, ParamSpaceAnalysis2/rebuild_model

arguments
    opts.psa_dir               (1,:) char   = ''
    opts.target_level_of_chaos (1,1) double = 2.0
    % '' -> the most-adapted condition of whatever run is loaded, resolved by
    % full_adaptation_condition. Naming a literal here only worked while every
    % preset used the same name for that regime; they no longer do. A saved run
    % from before the rename still works -- the fallback below tries the old
    % 'sfa_and_std' when no name matches the sfaX_stdY pattern.
    opts.target_condition      (1,:) char   = ''
    opts.target_rep            (1,1) double = 1
    opts.sim_T_range           (1,2) double = [0, 60]
    opts.overrides             (1,1) struct = struct()
end

setup_paths();
project_root = fileparts(which('setup_paths'));

%% Locate and load the PSA object
psa_dir = opts.psa_dir;
if isempty(psa_dir)
    search_root = fullfile(project_root, 'data', 'param_space');
    hits = dir(fullfile(search_root, '**', 'psa_object.mat'));
    % Keep only level_of_chaos sensitivity runs (folder name carries the param).
    if ~isempty(hits)
        hits = hits(contains({hits.folder}, 'level_of_chaos'));
    end
    assert(~isempty(hits), ['No level_of_chaos psa_object.mat found under %s. ', ...
        'Run run_sensitivity_analysis first, or pass psa_dir.'], search_root);
    [~, newest] = max([hits.datenum]);
    psa_dir = hits(newest).folder;
end
fprintf('Loading PSA from:\n  %s\n', psa_dir);
psa = ParamSpaceAnalysis2.from_dir(psa_dir);

%% Sanity checks on the loaded sweep
assert(ismember('level_of_chaos', psa.grid_params), ...
    'This PSA does not sweep level_of_chaos (grid_params: %s).', ...
    strjoin(psa.grid_params, ', '));
assert(~isempty(psa.all_configs), ...
    ['Loaded PSA has no all_configs (run may be incomplete). ', ...
     'Point psa_dir at a completed sensitivity run.']);
fprintf('Model class: %s\n', psa.model_class);

cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
if isempty(opts.target_condition)
    % The run's most-adapted regime. Falls back to the pre-rename name so a
    % saved run of any age can be reloaded.
    try
        opts.target_condition = full_adaptation_condition(psa.conditions);
    catch
        opts.target_condition = 'sfa_and_std';
    end
    fprintf('Condition: %s (most adapted; pass target_condition to override)\n', ...
        opts.target_condition);
end
ci = find(strcmp(cond_names, opts.target_condition), 1);
assert(~isempty(ci), 'Condition ''%s'' not found. Available: %s.', ...
    opts.target_condition, strjoin(cond_names, ', '));
condition = psa.conditions{ci};

%% Pick the config: level_of_chaos closest to target, matching rep
has_reps = ismember('reps', psa.grid_params);
n_cfg = numel(psa.all_configs);
loc_vals = nan(n_cfg, 1);
rep_vals = nan(n_cfg, 1);
for i = 1:n_cfg
    loc_vals(i) = psa.all_configs{i}.level_of_chaos;
    if has_reps, rep_vals(i) = psa.all_configs{i}.reps; end
end

cand = true(n_cfg, 1);
if has_reps
    cand = (rep_vals == opts.target_rep);
    assert(any(cand), 'No configs with reps == %g (available reps: %s).', ...
        opts.target_rep, num2str(unique(rep_vals(~isnan(rep_vals)))'));
end
idxs = find(cand);
[~, k] = min(abs(loc_vals(idxs) - opts.target_level_of_chaos));
config_idx = idxs(k);
config = psa.all_configs{config_idx};

fprintf('Selected config_idx=%d: level_of_chaos=%g', config_idx, config.level_of_chaos);
if has_reps, fprintf(', reps=%g', config.reps); end
fprintf('  (target level_of_chaos=%g)\n', opts.target_level_of_chaos);

%% The saved result for that grid point, which is what rebuild_model needs
% rebuild_model reconstructs from the RESULT, not the config: the network is
% pinned by res.network_seed, which run_single_job derived from the grid
% POSITION (config_idx*100 + offset). Reading it back rather than recomputing
% the formula is what keeps this correct across a subsetted run, where a point's
% position is unchanged by how many neighbours were skipped.
assert(isfield(psa.results, opts.target_condition), ...
    'This run has no results for condition ''%s''. Available: %s.', ...
    opts.target_condition, strjoin(fieldnames(psa.results)', ', '));
res_cell = psa.results.(opts.target_condition);
assert(numel(res_cell) >= config_idx && ~isempty(res_cell{config_idx}), ...
    ['Condition ''%s'' has no stored result at config_idx %d -- that grid ' ...
     'point was skipped or the run is incomplete. Pick another ' ...
     'target_level_of_chaos, or point psa_dir at a finished run.'], ...
    opts.target_condition, config_idx);
res = res_cell{config_idx};
fprintf('network_seed = %d (offset %g)\n', res.network_seed, psa.network_seed_offset);

%% Rebuild, exactly as the sweep did
% Delegated rather than reimplemented. This file used to assemble model_args by
% hand -- condition fields enumerated as {'n_a_E','n_a_I','n_b_E','n_b_I'}, grid
% params, then model_defaults minus both -- which is run_single_job's precedence
% written a second time, in SRNNModel2's vocabulary, free to drift from it.
% rebuild_model is the same logic kept beside the driver it mirrors, and it is
% class-generic.
model = psa.rebuild_model(res);

%% Overrides, applied after construction
% T_range and an emptied lya_T_interval so the LLE spans the whole window --
% this is a tool for LOOKING at the dynamics, and the swept window is usually
% too short to see them. Anything else the caller names goes on top.
%
% Assigned rather than passed as name-value pairs because rebuild_model returns
% a constructed model. Both classes are handle classes with public properties
% and no setters, so direct assignment is the supported route.
model.T_range = opts.sim_T_range;
model.lya_T_interval = [];

ov = fieldnames(opts.overrides);
for k = 1:numel(ov)
    if ~isprop(model, ov{k})
        error('check_sensitivity_sim:UnknownOverride', ...
            '''%s'' is not a property of %s. Its properties include: %s.', ...
            ov{k}, class(model), strjoin(sort(properties(model))', ', '));
    end
    model.(ov{k}) = opts.overrides.(ov{k});
end

%% Report, build, run, plot
% The condition is printed generically -- whatever fields it happens to carry --
% rather than assuming n_a_E / n_b_E. On SRNNCellTypePairs a condition carries
% tau_a and a synapse_config struct instead.
fprintf('Condition ''%s'': %s\n', condition.name, describe_condition(condition));
fprintf('Overrides: T_range = [%g, %g] s, lya_T_interval = []', ...
    opts.sim_T_range(1), opts.sim_T_range(2));
if ~isempty(ov); fprintf(', %s', strjoin(ov', ', ')); end
fprintf('\n\n');

model.build();
model.run();
model.plot();

fprintf('\nReconstructed run complete.');
if isprop(model, 'lya_results') && isfield(model.lya_results, 'LLE')
    fprintf(' LLE = %.4f\n', model.lya_results.LLE);
else
    fprintf('\n');
end
end

%% ------------------------------------------------------------------------
function s = describe_condition(condition)
% One line naming whatever a condition actually sets, whichever class it is for.
%
% SRNNModel2 conditions carry n_a_E / n_b_E counts; SRNNCellTypePairs ones carry
% a tau_a cell and a whole synapse_config struct. Enumerating either by name
% here would put a third copy of that vocabulary in the repo, so this prints
% what is there.
f = setdiff(fieldnames(condition)', {'name'});
if isempty(f); s = '(no parameters)'; return; end
parts = cell(1, numel(f));
for k = 1:numel(f)
    parts{k} = sprintf('%s=%s', f{k}, summarize(condition.(f{k})));
end
s = strjoin(parts, ', ');
end

function s = summarize(v)
% Compact rendering: numbers as-is, a synapse_config as its route list.
if isnumeric(v) || islogical(v)
    s = mat2str(v, 3);
elseif iscell(v)
    s = ['{' strjoin(cellfun(@(x) mat2str(x, 3), v, 'UniformOutput', false), ' ') '}'];
elseif isstruct(v)
    pres = fieldnames(v)';
    routes = {};
    for a = 1:numel(pres)
        posts = fieldnames(v.(pres{a}))';
        for b = 1:numel(posts)
            routes{end+1} = [pres{a} '->' posts{b}]; %#ok<AGROW>
        end
    end
    if isempty(routes); s = '(none)'; else; s = strjoin(routes, ','); end
else
    s = class(v);
end
end

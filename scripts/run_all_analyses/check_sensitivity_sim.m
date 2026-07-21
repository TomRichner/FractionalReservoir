% check_sensitivity_sim.m
% Reload the parameters of a single level_of_chaos = 2.0 run from a completed
% 1D sensitivity analysis, re-simulate that exact model, and plot its time
% series with model.plot().
%
% Motivation: with level_of_chaos = 2.0 the network stays highly chaotic even
% under two-timescale STD. This script lets you inspect the actual dynamics of
% one such run. It rebuilds the model exactly as ParamSpaceAnalysis2 did in
% run_single_job -- same network seed (config_idx*100 + network_seed_offset),
% same adaptation condition (n_a_E / n_b_E), and same model_defaults (including
% tau_b_E_rec) -- so the reconstructed trajectory matches the sweep.
%
% Configuration (override any of these in the base workspace before running):
%   psa_dir               - directory holding psa_object.mat. '' auto-finds the
%                           most recent level_of_chaos sensitivity run under
%                           data/param_space (searched recursively).
%   target_level_of_chaos - grid level to reload (default 2.0).
%   target_condition      - condition name, see psa.conditions (default
%                           'sfa_and_std').
%   target_rep            - which 'reps' grid value to use, if reps is an axis
%                           (default 1).
%
% See also: run_sensitivity_analysis, ParamSpaceAnalysis2, SRNNModel2

close all; clc;
setup_paths();

%% ---- Configuration (base-workspace overrides honored) ----
if ~exist('psa_dir', 'var');               psa_dir = '';                      end
if ~exist('target_level_of_chaos', 'var'); target_level_of_chaos = 2.0;       end
if ~exist('target_condition', 'var');      target_condition = 'sfa_and_std';  end
if ~exist('target_rep', 'var');            target_rep = 1;                    end
if ~exist('sim_T_range', 'var');           sim_T_range = [0, 60];             end  % override the swept window
if ~exist('sim_c_E', 'var');               sim_c_E = 1/3;                     end  % override SFA strength

project_root = fileparts(fileparts(which('setup_paths')));

%% ---- Locate and load the PSA object ----
if isempty(psa_dir)
    search_root = fullfile(project_root, 'data', 'param_space');
    hits = dir(fullfile(search_root, '**', 'psa_object.mat'));
    % Keep only level_of_chaos sensitivity runs (folder name carries the param).
    if ~isempty(hits)
        hits = hits(contains({hits.folder}, 'level_of_chaos'));
    end
    assert(~isempty(hits), ['No level_of_chaos psa_object.mat found under %s. ', ...
        'Run run_sensitivity_analysis first, or set psa_dir manually.'], search_root);
    [~, newest] = max([hits.datenum]);
    psa_dir = hits(newest).folder;
end
psa_file = fullfile(psa_dir, 'psa_object.mat');
assert(isfile(psa_file), 'psa_object.mat not found in %s', psa_dir);
fprintf('Loading PSA from:\n  %s\n', psa_file);
loaded = load(psa_file, 'psa');
psa = loaded.psa;

%% ---- Sanity checks on the loaded sweep ----
assert(ismember('level_of_chaos', psa.grid_params), ...
    'This PSA does not sweep level_of_chaos (grid_params: %s).', strjoin(psa.grid_params, ', '));
assert(~isempty(psa.all_configs), ...
    ['Loaded PSA has no all_configs (run may be incomplete). ', ...
     'Point psa_dir at a completed sensitivity run.']);

cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
ci = find(strcmp(cond_names, target_condition), 1);
assert(~isempty(ci), 'Condition ''%s'' not found. Available: %s.', ...
    target_condition, strjoin(cond_names, ', '));
condition = psa.conditions{ci};

%% ---- Pick the config: level_of_chaos closest to target, matching rep ----
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
    cand = (rep_vals == target_rep);
    assert(any(cand), 'No configs with reps == %g (available reps: %s).', ...
        target_rep, num2str(unique(rep_vals(~isnan(rep_vals)))'));
end
idxs = find(cand);
[~, k] = min(abs(loc_vals(idxs) - target_level_of_chaos));
config_idx = idxs(k);
config = psa.all_configs{config_idx};

fprintf('Selected config_idx=%d: level_of_chaos=%g', config_idx, config.level_of_chaos);
if has_reps, fprintf(', reps=%g', config.reps); end
fprintf('  (target level_of_chaos=%g)\n', target_level_of_chaos);

%% ---- Network seed (matches run_batched_simulation) ----
% Prefer the seed recorded in the saved result; fall back to the PSA formula.
network_seed = config_idx * 100 + psa.network_seed_offset;
if isfield(psa.results, target_condition)
    res_cell = psa.results.(target_condition);
    if numel(res_cell) >= config_idx && ~isempty(res_cell{config_idx}) ...
            && isfield(res_cell{config_idx}, 'network_seed')
        network_seed = res_cell{config_idx}.network_seed;
    end
end
fprintf('network_seed = %d (offset %g)\n', network_seed, psa.network_seed_offset);

%% ---- Assemble model_args exactly as ParamSpaceAnalysis2.run_single_job ----
model_args = {'rng_seeds', [network_seed, network_seed + 1]};

% Condition parameters (adaptation counts).
cond_fields = {'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'};
for f = cond_fields
    if isfield(condition, f{1})
        model_args = [model_args, {f{1}, condition.(f{1})}]; %#ok<AGROW>
    end
end

% Grid parameters (handle vector params via the saved lookup, like the driver).
for p = 1:numel(psa.grid_params)
    pname = psa.grid_params{p};
    if isfield(psa.vector_param_lookup, pname)
        model_args = [model_args, {pname, psa.vector_param_lookup.(pname){config.(pname)}}]; %#ok<AGROW>
    else
        model_args = [model_args, {pname, config.(pname)}]; %#ok<AGROW>
    end
end

% Model defaults, skipping grid params and condition-controlled counts.
df = fieldnames(psa.model_defaults);
for d = 1:numel(df)
    fn = df{d};
    if ~ismember(fn, psa.grid_params) && ~ismember(fn, cond_fields)
        model_args = [model_args, {fn, psa.model_defaults.(fn)}]; %#ok<AGROW>
    end
end

% Overrides (appended last so they win over anything from model_defaults):
%   - T_range / lya_T_interval: run the full window and let the LLE span it.
%   - c_E: SFA strength for E neurons (differs from the swept run's default).
model_args = [model_args, {'T_range', sim_T_range, 'lya_T_interval', [], 'c_E', sim_c_E}];

%% ---- Report, build, run, plot ----
fprintf('Condition ''%s'': n_a_E=%d, n_b_E=%d', condition.name, ...
    local_getfield(condition, 'n_a_E', 0), local_getfield(condition, 'n_b_E', 0));
tau_rec = local_getfield(psa.model_defaults, 'tau_b_E_rec', []);
if ~isempty(tau_rec), fprintf(', tau_b_E_rec=[%s]', num2str(tau_rec)); end
fprintf('\nOverrides: T_range = [%g, %g] s, c_E = %g\n\n', sim_T_range(1), sim_T_range(2), sim_c_E);

model = SRNNModel2(model_args{:});
model.build();
model.run();
model.plot();

fprintf('\nReconstructed run complete.');
if isprop(model, 'lya_results') && isfield(model.lya_results, 'LLE')
    fprintf(' LLE = %.4f\n', model.lya_results.LLE);
else
    fprintf('\n');
end

%% ---- Local helper ----
function v = local_getfield(s, f, default_val)
% LOCAL_GETFIELD Return s.(f) if present, else default_val.
if isstruct(s) && isfield(s, f)
    v = s.(f);
else
    v = default_val;
end
end

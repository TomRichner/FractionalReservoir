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
%   sim_c_E                override SFA strength.
%
% SRNNModel2 ONLY. It constructs SRNNModel2 directly and reads the condition
% through the n_a_E / n_b_E vocabulary, neither of which SRNNCellTypePairs
% shares. Pointed at a Pairs run it will fail on the constructor rather than
% silently reconstruct the wrong model.
%
% See also: run_sensitivity_analysis, ParamSpaceAnalysis2, SRNNModel2

arguments
    opts.psa_dir               (1,:) char   = ''
    opts.target_level_of_chaos (1,1) double = 2.0
    opts.target_condition      (1,:) char   = 'sfa_and_std'
    opts.target_rep            (1,1) double = 1
    opts.sim_T_range           (1,2) double = [0, 60]
    opts.sim_c_E               (1,1) double = 1/3
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
assert(strcmp(psa.model_class, 'SRNNModel2'), ...
    ['check_sensitivity_sim reconstructs an SRNNModel2, but this run used %s. ', ...
     'Its condition vocabulary (n_a_E / n_b_E) does not apply.'], psa.model_class);

cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
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

%% Network seed (matches run_batched_simulation)
% Prefer the seed recorded in the saved result; fall back to the PSA formula.
network_seed = config_idx * 100 + psa.network_seed_offset;
if isfield(psa.results, opts.target_condition)
    res_cell = psa.results.(opts.target_condition);
    if numel(res_cell) >= config_idx && ~isempty(res_cell{config_idx}) ...
            && isfield(res_cell{config_idx}, 'network_seed')
        network_seed = res_cell{config_idx}.network_seed;
    end
end
fprintf('network_seed = %d (offset %g)\n', network_seed, psa.network_seed_offset);

%% Assemble model_args exactly as ParamSpaceAnalysis2.run_single_job does
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
model_args = [model_args, {'T_range', opts.sim_T_range, ...
    'lya_T_interval', [], 'c_E', opts.sim_c_E}];

%% Report, build, run, plot
fprintf('Condition ''%s'': n_a_E=%d, n_b_E=%d', condition.name, ...
    local_getfield(condition, 'n_a_E', 0), local_getfield(condition, 'n_b_E', 0));
tau_rec = local_getfield(psa.model_defaults, 'tau_b_E_rec', []);
if ~isempty(tau_rec), fprintf(', tau_b_E_rec=[%s]', num2str(tau_rec)); end
fprintf('\nOverrides: T_range = [%g, %g] s, c_E = %g\n\n', ...
    opts.sim_T_range(1), opts.sim_T_range(2), opts.sim_c_E);

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
end

%% ------------------------------------------------------------------------
function v = local_getfield(s, f, default_val)
% Return s.(f) if present, else default_val.
if isstruct(s) && isfield(s, f)
    v = s.(f);
else
    v = default_val;
end
end

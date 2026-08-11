% test_psa_saveload.m
% Verify saveobj/loadobj round-trip for ParamSpaceAnalysis2, including the
% resolved_defaults parameter snapshot, and that effective_param answers
% identically before and after a save/load cycle.
%
% Runs a deliberately tiny sweep (serial, small n, short T_range) so it
% finishes in seconds.
%
% See also: ParamSpaceAnalysis2, SRNNModel2

close all
clc;
fprintf('=== Testing ParamSpaceAnalysis2 saveobj/loadobj ===\n\n');

%% Create and configure a minimal PSA
fprintf('Creating ParamSpaceAnalysis2 object...\n');
psa = ParamSpaceAnalysis2(...
    'n_levels', 2, ...
    'batch_size', 10, ...
    'note', 'saveload_test', ...
    'verbose', true ...
    );
psa.use_parallel = false;
psa.output_dir = fullfile(tempdir, 'psa_saveload_test');

% 2 params x 2 levels = 4 combinations, x 2 conditions = 8 simulations
psa.add_grid_parameter('level_of_chaos', [1, 1.5]);
psa.add_grid_parameter('f', [0.4, 0.6]);

psa.set_conditions({ ...
    struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
    struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1) ...
    });

% Fast model defaults
psa.model_defaults.n = 20;
psa.model_defaults.indegree = 10;
psa.model_defaults.T_range = [0, 10];
psa.model_defaults.fs = 200;
psa.model_defaults.ode_solver = @ode_rk4;
psa.model_defaults.lya_T_interval = [5, 10];
psa.model_defaults.lya_method = 'benettin';
psa.store_local_lya = true;

%% Run
fprintf('\nRunning analysis (8 simulations)...\n');
psa.run();

%% Save, clear, reload
save_file = fullfile(psa.output_dir, 'psa_object_test.mat');
fprintf('\nSaving PSA object to: %s\n', save_file);
save(save_file, 'psa');

fprintf('Clearing and reloading...\n');
psa_original = psa;
clear psa;
loaded = load(save_file);
psa_loaded = loaded.psa;

%% Compare properties
fprintf('\n=== Property Comparison ===\n');
all_passed = true;

props = {'grid_params', 'param_ranges', 'n_levels', 'conditions', ...
    'model_defaults', 'verbose', ...
    'batch_size', 'output_dir', 'note', 'vector_param_lookup'};
for i = 1:numel(props)
    all_passed = compare_property(props{i}, ...
        psa_original.(props{i}), psa_loaded.(props{i})) && all_passed;
end

% resolved_defaults holds anonymous function handles (the activation and its
% derivative), and MATLAB does not consider an anonymous handle isequal to
% itself after a save/load cycle. Compare handles by their string form, the
% same way same_config does.
all_passed = check('resolved_defaults', ...
    structs_equivalent(psa_original.resolved_defaults, psa_loaded.resolved_defaults)) && all_passed;

%% resolved_defaults content
fprintf('\n=== resolved_defaults ===\n');
rd = psa_original.resolved_defaults;
all_passed = check('is non-empty after run()', ~isempty(fieldnames(rd))) && all_passed;
all_passed = check('captures an explicit override (fs = 200)', ...
    isfield(rd, 'fs') && isequal(rd.fs, 200)) && all_passed;
all_passed = check('captures a set_defaults side effect (activation_function)', ...
    isfield(rd, 'activation_function') && ...
    isa(rd.activation_function, 'function_handle')) && all_passed;
all_passed = check('captures an untouched class default (tau_d)', ...
    isfield(rd, 'tau_d') && ...
    isequal(rd.tau_d, ParamSpaceAnalysis2.class_default('tau_d'))) && all_passed;
all_passed = check('excludes grid parameters (f, level_of_chaos)', ...
    ~isfield(rd, 'f') && ~isfield(rd, 'level_of_chaos')) && all_passed;
all_passed = check('excludes condition fields (n_a_E, n_b_E)', ...
    ~isfield(rd, 'n_a_E') && ~isfield(rd, 'n_b_E')) && all_passed;

%% effective_param survives the round trip
fprintf('\n=== effective_param round trip ===\n');
cond_names = fieldnames(psa_original.results);
first_cond = cond_names{1};
res = psa_original.results.(first_cond){1};
ep_names = {'f', 'level_of_chaos', 'n_a_E', 'tau_d', 'fs', 'T_range'};
ep_ok = true;
for i = 1:numel(ep_names)
    a = psa_original.effective_param(res, ep_names{i});
    b = psa_loaded.effective_param(psa_loaded.results.(first_cond){1}, ep_names{i});
    if ~isequaln(a, b); ep_ok = false; end
end
all_passed = check('identical answers before/after load', ep_ok) && all_passed;
all_passed = check('condition field resolves from condition_name', ...
    isequal(psa_original.effective_param(res, 'n_a_E'), ...
    condition_value(psa_original, res.condition_name, 'n_a_E'))) && all_passed;

%% same_config: a run is poolable with itself, and with its reload
fprintf('\n=== same_config ===\n');
[tf, reason] = psa_original.same_config(psa_loaded);
all_passed = check(sprintf('self-poolable (%s)', reason), tf) && all_passed;

% A run differing only in a parameter NOT named in model_defaults of either
% run would previously have compared equal; resolved_defaults catches it.
fprintf('\nRunning a second sweep with tau_d = 0.2 ...\n');
psa_diff = build_twin(psa_original);
psa_diff.model_defaults.tau_d = 0.2;
psa_diff.run();
[tf_diff, reason_diff] = psa_original.same_config(psa_diff);
all_passed = check(sprintf('refuses a differing tau_d (%s)', reason_diff), ...
    ~tf_diff && contains(reason_diff, 'tau_d')) && all_passed;

% Legacy runs (saved before resolved_defaults existed) are refused by default
% and poolable only under allow_legacy. Uses a real pre-refactor run if one is
% present in data/; those directories are gitignored, so skip when absent.
legacy_file = find_legacy_psa();
if isempty(legacy_file)
    fprintf('  legacy-run checks: SKIP (no pre-refactor psa_object.mat in data/)\n');
else
    lg = load(legacy_file);
    fn = fieldnames(lg);
    psa_legacy = lg.(fn{1});
    all_passed = check('legacy run loads with empty resolved_defaults', ...
        isempty(fieldnames(psa_legacy.resolved_defaults))) && all_passed;
    [tf_leg, reason_leg] = psa_legacy.same_config(psa_legacy);
    all_passed = check(sprintf('refused by default (%s)', reason_leg), ...
        ~tf_leg && contains(reason_leg, 'predate')) && all_passed;
    tf_leg_ok = psa_legacy.same_config(psa_legacy, 'allow_legacy', true);
    all_passed = check('poolable with allow_legacy', tf_leg_ok) && all_passed;
    all_passed = check('effective_param still works on a legacy run', ...
        ~isempty(psa_legacy.effective_param([], 'tau_d'))) && all_passed;
end

%% Results survived
fprintf('\n=== Results ===\n');
orig_LLEs = cellfun(@(r) r.LLE, psa_original.results.(first_cond));
loaded_LLEs = cellfun(@(r) r.LLE, psa_loaded.results.(first_cond));
all_passed = check(sprintf('results.%s.LLE', first_cond), ...
    isequaln(orig_LLEs, loaded_LLEs)) && all_passed;

%% Final result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');
fprintf('\nTest output directory: %s\n', psa_original.output_dir);

%% Helper functions
function tf = structs_equivalent(a, b)
% Field-by-field equality, comparing function handles by their string form.
tf = false;
fa = sort(fieldnames(a));
fb = sort(fieldnames(b));
if ~isequal(fa, fb); return; end
for i = 1:numel(fa)
    va = a.(fa{i});
    vb = b.(fa{i});
    if isa(va, 'function_handle') || isa(vb, 'function_handle')
        if ~(isa(va, 'function_handle') && isa(vb, 'function_handle') ...
                && strcmp(func2str(va), func2str(vb)))
            fprintf('    (handle mismatch on %s)\n', fa{i});
            return;
        end
    elseif ~isequaln(va, vb)
        fprintf('    (value mismatch on %s)\n', fa{i});
        return;
    end
end
tf = true;
end

function twin = build_twin(src)
% A second PSA with the same grid and conditions as SRC, so same_config can
% only differ on model parameters.
twin = ParamSpaceAnalysis2('n_levels', src.n_levels, 'batch_size', 10, ...
    'note', 'saveload_test_twin', 'verbose', false);
twin.use_parallel = false;
twin.output_dir = fullfile(tempdir, 'psa_saveload_test_twin');
for i = 1:numel(src.grid_params)
    p = src.grid_params{i};
    twin.add_grid_parameter(p, src.param_ranges.(p));
end
twin.set_conditions(src.conditions);
twin.model_defaults = src.model_defaults;
twin.store_local_lya = src.store_local_lya;
end

function f = find_legacy_psa()
% First psa_object.mat under data/ that predates resolved_defaults, or ''.
f = '';
root = fileparts(which('setup_paths'));
hits = dir(fullfile(root, 'data', 'param_space', '**', 'psa_object.mat'));
for i = 1:numel(hits)
    candidate = fullfile(hits(i).folder, hits(i).name);
    try
        L = load(candidate);
        fn = fieldnames(L);
        p = L.(fn{1});
        if isa(p, 'ParamSpaceAnalysis2') && isempty(fieldnames(p.resolved_defaults))
            f = candidate;
            return;
        end
    catch
        % Unreadable or a different class -- keep looking.
    end
end
end

function val = condition_value(psa, cond_name, field)
val = [];
for i = 1:numel(psa.conditions)
    if strcmp(psa.conditions{i}.name, cond_name)
        val = psa.conditions{i}.(field);
        return;
    end
end
end

function passed = compare_property(name, orig, loaded)
passed = check(name, isequaln(orig, loaded));
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end

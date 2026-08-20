% test_psa_loaders.m
% Verify that a run can be read back off disk as the object that produced it.
%
% This file exists because of two RESOLVED loader bugs, both fixed 2026-08-14.
% `tau_a_E` is a VECTOR parameter, so it cannot be a grid coordinate directly --
% add_vector_parameter pre-builds the concrete vectors into vector_param_lookup
% and the grid axis carries a LEVEL INDEX, which run_single_job decodes before
% constructing the model. The runs were always correct; the loader was not.
% param_space_summary.mat never stored the lookup, so the old load_results left
% effective_param handing back 1..n_levels dressed up as seconds, silently. The
% same loader also dropped model_class, so a loaded SRNNCellTypePairs run
% reported SRNNModel2.
%
% Both are gone because psa_object.mat is now the single authoritative artifact:
% run() writes it before batching (so an interrupted run keeps its config) and
% again on completion, and ParamSpaceAnalysis2.from_dir is the only way in.
%
% The vector parameter is the whole point of this file -- it is the only kind
% that exposed the bug, and test_psa_saveload never exercised a loader path.
%
% See also: ParamSpaceAnalysis2.from_dir, save_object, test_psa_saveload

fprintf('=== Testing ParamSpaceAnalysis2 loaders ===\n\n');
all_passed = true;

root = fullfile(tempdir, 'psa_loader_test');
if exist(root, 'dir'); rmdir(root, 's'); end
mkdir(root);
cleanup_root = onCleanup(@() rmdir(root, 's'));

%% A tiny sweep over a VECTOR parameter
fprintf('Running a small tau_a_E sweep...\n');
psa = ParamSpaceAnalysis2('n_levels', 3, 'verbose', false, ...
    'use_parallel', false, 'batch_size', 2, 'output_dir', root, ...
    'note', 'loaders', 'randomize_order', false);
psa.model_defaults = struct('n', 12, 'indegree', 6, 'fs', 200, ...
    'T_range', [0 2], 'lya_T_interval', [1 2], 'lya_warmup', 1, ...
    'ode_solver', 'rk4', 'lya_method', 'benettin');
psa.set_conditions({struct('name', 'sfa_only', 'n_a_E', 3, 'n_b_E', 0)});
psa.add_vector_parameter('tau_a_E', ...
    'vary_element', 'last', 'fixed_value', 0.25, 'vary_range', [1 30], ...
    'n_elements', 3, 'spacing', 'log', 'level_spacing', 'linear');
psa.add_grid_parameter('level_of_chaos', [0.8, 1.2]);   % a scalar axis alongside
evalc('psa.run();');
run_dir = psa.output_dir;
res_ref = psa.results.sfa_only{1};
tau_ref = psa.effective_param(res_ref, 'tau_a_E');

all_passed = check('the sweep produced a 3-element tau_a_E', ...
    numel(tau_ref) == 3) && all_passed;

%% run() writes psa_object.mat itself -- no script has to
all_passed = check('run() wrote psa_object.mat', ...
    isfile(fullfile(run_dir, 'psa_object.mat'))) && all_passed;
S = load(fullfile(run_dir, 'psa_object.mat'));
all_passed = check('...under the canonical variable name ''psa''', ...
    isfield(S, 'psa')) && all_passed;

%% THE BUG: a vector parameter must decode to a value, not a level index
q = ParamSpaceAnalysis2.from_dir(run_dir);
tau_loaded = q.effective_param(q.results.sfa_only{1}, 'tau_a_E');
all_passed = check('effective_param returns the VECTOR, not a level index', ...
    numel(tau_loaded) == 3 && isequal(tau_loaded, tau_ref)) && all_passed;
all_passed = check('...and every level decodes to a distinct vector', ...
    numel(unique(cellfun(@(r) max(q.effective_param(r, 'tau_a_E')), ...
        q.results.sfa_only))) == 3) && all_passed;

%% model_class survives the round trip
all_passed = check('model_class survives from_dir', ...
    strcmp(q.model_class, psa.model_class)) && all_passed;
% The mechanism, without paying for a SRNNCellTypePairs sweep.
mc_dir = fullfile(root, 'model_class_only');
mkdir(mc_dir);
tmp = ParamSpaceAnalysis2('verbose', false, 'output_dir', mc_dir);
tmp.model_class = 'SRNNCellTypePairs';
tmp.save_object();
all_passed = check('a non-default model_class round trips', ...
    strcmp(ParamSpaceAnalysis2.from_dir(mc_dir).model_class, ...
        'SRNNCellTypePairs')) && all_passed;

%% Config and results both survive
all_passed = check('grid_params survive', ...
    isequal(q.grid_params, psa.grid_params)) && all_passed;
all_passed = check('resolved_defaults survive', ...
    isequal(q.resolved_defaults, psa.resolved_defaults)) && all_passed;
all_passed = check('every result survives', ...
    numel(q.results.sfa_only) == numel(psa.results.sfa_only)) && all_passed;

%% The early write must carry everything needed to interpret a crashed run
% run() saves the object before batching (ParamSpaceAnalysis2.run, just after
% generate_grid), so a run that dies part-way leaves its full configuration next
% to the temp_batches/ that hold its results. The decode table is the part that
% matters: without it the level indices in those results are uninterpretable.
%
% The crash state itself is not constructible here -- consolidate deletes
% temp_batches/ on success, and results/has_run are SetAccess=private -- so the
% recovery path is verified against the real interrupted runs on disk
% (data/param_space/run_all_aug_13_26_23_37/FAILED_OOM_tau_*), not synthesised.
all_passed = check('the saved object carries the vector decode table', ...
    isfield(q.vector_param_lookup, 'tau_a_E') && ...
    numel(q.vector_param_lookup.tau_a_E) == psa.n_levels) && all_passed;
all_passed = check('the saved object carries the full grid', ...
    numel(q.all_configs) == numel(psa.all_configs)) && all_passed;
all_passed = check('per-condition result files exist alongside it', ...
    isfile(fullfile(run_dir, 'sfa_only', ...
        'param_space_results_sfa_only.mat'))) && all_passed;

%% Legacy variable names still load -- selection is by CLASS, not by name
legacy_dir = fullfile(root, 'legacy_name');
mkdir(legacy_dir);
psa_tau_a = q;
save(fullfile(legacy_dir, 'psa_object.mat'), 'psa_tau_a');
all_passed = check('an object saved as ''psa_tau_a'' still loads', ...
    isa(ParamSpaceAnalysis2.from_dir(legacy_dir), 'ParamSpaceAnalysis2')) && all_passed;

%% Fail closed: a vector parameter whose lookup is not populated
% Reachable through the public API: add_vector_parameter registers the parameter
% immediately, but the lookup is only built when the grid is generated.
unbuilt = ParamSpaceAnalysis2('verbose', false);
unbuilt.add_vector_parameter('tau_a_E', ...
    'vary_element', 'last', 'fixed_value', 0.25, 'vary_range', [1 30], ...
    'n_elements', 3, 'spacing', 'log', 'level_spacing', 'linear');
[threw, err] = capture(@() unbuilt.effective_param( ...
    struct('config', struct('tau_a_E', 3)), 'tau_a_E'));
all_passed = check('an undecodable vector parameter ERRORS, not returns an index', ...
    threw && strcmp(err.identifier, 'ParamSpaceAnalysis2:MissingVectorLookup')) ...
    && all_passed;
all_passed = check('...and the message says it is a level index', ...
    threw && contains(err.message, 'LEVEL INDEX')) && all_passed;

%% A scalar grid axis is unaffected by that guard -- config.(name) IS its value
loc_ref = psa.effective_param(res_ref, 'level_of_chaos');
all_passed = check('a scalar grid axis still reads back as its value', ...
    isscalar(loc_ref) && ...
    isequal(q.effective_param(q.results.sfa_only{1}, 'level_of_chaos'), loc_ref)) ...
    && all_passed;

%% A directory with no run in it says so
empty_dir = fullfile(root, 'empty');
mkdir(empty_dir);
[threw, err] = capture(@() ParamSpaceAnalysis2.from_dir(empty_dir));
all_passed = check('a directory with no psa_object.mat errors clearly', ...
    threw && strcmp(err.identifier, 'ParamSpaceAnalysis2:NoPsaObject')) && all_passed;
[threw, ~] = capture(@() ParamSpaceAnalysis2.from_dir( ...
    fullfile(root, 'does_not_exist')));
all_passed = check('a missing directory errors', threw) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function [threw, err] = capture(fcn)
threw = false; err = [];
try
    fcn();
catch err
    threw = true;
end
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

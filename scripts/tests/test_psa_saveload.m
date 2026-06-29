% test_psa_saveload.m
% Test script to verify saveobj/loadobj round-trip for ParamSpaceAnalysis
%
% This script creates a minimal PSA run with fast parameters, saves the
% object, loads it back, and verifies all properties match.

clear all
close all
clc;
fprintf('=== Testing ParamSpaceAnalysis saveobj/loadobj ===\n\n');

%% Setup paths
setup_paths();

%% Create and configure PSA with fast parameters
fprintf('Creating ParamSpaceAnalysis object...\n');
psa = ParamSpaceAnalysis(...
    'n_levels', 2, ...
    'batch_size', 10, ...
    'note', 'saveload_test', ...
    'verbose', true ...
    );

% Add minimal grid - 2 params x 2 levels = 4 combinations x 4 conditions = 16 sims
psa.add_grid_parameter('level_of_chaos', [1, 1.5]);
psa.add_grid_parameter('EI_imbalance', [0.8, 1.2]);

% Configure fast model defaults
psa.model_defaults.n = 20;
psa.model_defaults.T_range = [-5, 10];
psa.model_defaults.fs = 100;
psa.model_defaults.lya_method = 'benettin';
psa.store_local_lya = true;

%% Run the analysis
fprintf('\nRunning analysis (16 simulations)...\n');
psa.run();

%% Save the object
save_file = fullfile(psa.output_dir, 'psa_object_test.mat');
fprintf('\nSaving PSA object to: %s\n', save_file);
save(save_file, 'psa');

%% Clear and reload
fprintf('Clearing and reloading...\n');
psa_original = psa;
clear psa;
loaded = load(save_file);
psa_loaded = loaded.psa;

%% Compare all properties
fprintf('\n=== Property Comparison ===\n');
all_passed = true;

% Helper function for comparison
compare_prop = @(name, orig, loaded) compare_property(name, orig, loaded);

% Configuration Properties
all_passed = all_passed && compare_prop('grid_params', psa_original.grid_params, psa_loaded.grid_params);
all_passed = all_passed && compare_prop('param_ranges', psa_original.param_ranges, psa_loaded.param_ranges);
all_passed = all_passed && compare_prop('n_levels', psa_original.n_levels, psa_loaded.n_levels);
all_passed = all_passed && compare_prop('conditions', psa_original.conditions, psa_loaded.conditions);

% Model Default Properties
all_passed = all_passed && compare_prop('model_defaults', psa_original.model_defaults, psa_loaded.model_defaults);
all_passed = all_passed && compare_prop('verbose', psa_original.verbose, psa_loaded.verbose);

% Execution Properties
all_passed = all_passed && compare_prop('batch_size', psa_original.batch_size, psa_loaded.batch_size);
all_passed = all_passed && compare_prop('output_dir', psa_original.output_dir, psa_loaded.output_dir);
all_passed = all_passed && compare_prop('note', psa_original.note, psa_loaded.note);

% Check results exist and have correct structure
has_results_orig = ~isempty(fieldnames(psa_original.results));
has_results_loaded = ~isempty(fieldnames(psa_loaded.results));
if has_results_orig && has_results_loaded
    % Compare LLE from first condition
    cond_names = fieldnames(psa_original.results);
    first_cond = cond_names{1};
    orig_LLEs = cellfun(@(r) r.LLE, psa_original.results.(first_cond));
    loaded_LLEs = cellfun(@(r) r.LLE, psa_loaded.results.(first_cond));
    if isequal(orig_LLEs, loaded_LLEs)
        fprintf('  results.%s.LLE: PASS\n', first_cond);
    else
        fprintf('  results.%s.LLE: FAIL\n', first_cond);
        all_passed = false;
    end
else
    fprintf('  results: FAIL (empty)\n');
    all_passed = false;
end

%% Final result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Plot distributions for visual comparison
fprintf('\nGenerating comparison plots...\n');

% Original object plots
fprintf('Plotting original PSA distributions...\n');
psa_original.plot('metric', 'LLE');
set(gcf, 'Name', 'Original PSA - LLE Distribution');

psa_original.plot('metric', 'mean_rate');
set(gcf, 'Name', 'Original PSA - Mean Rate Distribution');

% Loaded object plots
fprintf('Plotting loaded PSA distributions...\n');
psa_loaded.plot('metric', 'LLE');
set(gcf, 'Name', 'Loaded PSA - LLE Distribution');

psa_loaded.plot('metric', 'mean_rate');
set(gcf, 'Name', 'Loaded PSA - Mean Rate Distribution');

fprintf('\nDone! Compare the 4 figure windows to verify distributions match.\n');

%% Helper function
function passed = compare_property(name, orig, loaded)
if isequal(orig, loaded)
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end

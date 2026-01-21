function psa = load_and_plot_param_space_analysis(results_dir)
% LOAD_AND_PLOT_PARAM_SPACE_ANALYSIS Load a previous PSA run and plot results
%
% Usage:
%   psa = load_and_plot_param_space_analysis('/path/to/param_space_...')
%
% This function:
%   1. Loads the PSA object from psa_object.mat if available
%   2. Otherwise, creates a new object and loads results
%   3. Checks if consolidation is needed (temp_batches exists)
%   4. Consolidates if necessary
%   5. Plots LLE and mean_rate distributions
%
% Input:
%   results_dir - Path to a param_space_* output directory
%
% Output:
%   psa - The loaded ParamSpaceAnalysis object
%
% See also: ParamSpaceAnalysis, run_param_space_analysis

if nargin < 1 || isempty(results_dir)
    error('load_and_plot_param_space_analysis:NoInput', ...
        'Please provide the path to a param_space_* results directory.');
end

if ~exist(results_dir, 'dir')
    error('load_and_plot_param_space_analysis:NotFound', ...
        'Directory not found: %s', results_dir);
end

fprintf('=== Loading Parameter Space Analysis ===\n');
fprintf('Directory: %s\n\n', results_dir);

%% Try to load full PSA object first
psa_file = fullfile(results_dir, 'psa_object.mat');
if exist(psa_file, 'file')
    fprintf('Loading PSA object from: %s\n', psa_file);
    loaded = load(psa_file);
    psa = loaded.psa;
    fprintf('Loaded PSA object successfully.\n');
else
    % Fall back to creating new object and loading results
    fprintf('No psa_object.mat found. Creating new object and loading results...\n');
    psa = ParamSpaceAnalysis();
    psa.output_dir = results_dir;
end

%% Check if consolidation is needed
temp_dir = fullfile(results_dir, 'temp_batches');
if exist(temp_dir, 'dir')
    fprintf('\nFound unconsolidated batch files in temp_batches/\n');
    fprintf('Running consolidation...\n');
    psa.consolidate();
    fprintf('Consolidation complete.\n');
elseif ~psa.has_run
    % Need to load results if object doesn't have them
    fprintf('Loading results from per-condition MAT files...\n');
    psa.load_results(results_dir);
end

%% Display summary
fprintf('\n=== Analysis Summary ===\n');
fprintf('Output directory: %s\n', psa.output_dir);
if ~isempty(psa.grid_params)
    fprintf('Grid parameters: %s\n', strjoin(psa.grid_params, ', '));
    fprintf('Levels per parameter: %d\n', psa.n_levels);
end
if ~isempty(psa.conditions)
    cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
    fprintf('Conditions: %s\n', strjoin(cond_names, ', '));
end

%% Count successful results
if ~isempty(fieldnames(psa.results))
    for i = 1:length(psa.conditions)
        cond_name = psa.conditions{i}.name;
        if isfield(psa.results, cond_name)
            results = psa.results.(cond_name);
            n_success = sum(cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, results));
            fprintf('  %s: %d/%d successful\n', cond_name, n_success, length(results));
        end
    end
end

%% Plot results
fprintf('\nGenerating plots...\n');
psa.plot('metric', 'LLE');
psa.plot('metric', 'mean_rate');

fprintf('\nDone!\n');

end

% Fig_2_fraction_excitatory_load_and_plot.m
% Load saved PSA results and regenerate plots + stats without re-running.
%
% Usage:
%   1. Edit 'results_dir' below to point at your param_space_* directory
%   2. Run this script (setup_paths must be called first)
%
% Ported from ConnectivityAdaptation to use ParamSpaceAnalysis2 / SRNNModel2.
%
% See also: Fig_2_fraction_excitatory_analysis, ParamSpaceAnalysis2

%% ========== USER CONFIGURATION ==========

% Point this at your param_space_* output directory
results_dir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
    'data', 'param_space', 'param_space_fig2_frac_exc_nLevs_5_EDIT_ME');

% Plotting options
transient_skip = 3;
periods_to_plot = [0 1 1];
metrics_to_hist = {'br', 'lle'};

% Output control
save_figs = false;
save_stats = true;

% ==========================================

% NOTE: save_figs is set in the USER CONFIGURATION block above. A
% `master_save_figs` override used to be read out of the base workspace here;
% that protocol is gone (no script in this repo passes settings to another
% through workspace variables -- see CLAUDE.md), and the block silently did
% nothing whenever the variable was absent.

%% Setup paths
setup_paths();

% Derive project_root from setup_paths.m (at the repo root) so this tolerates living
% in a subdirectory such as scripts/EI_balance/.
project_root = fileparts(which('setup_paths'));
figs_root = fullfile(project_root, 'figs');

%% Validate and load PSA object
if ~exist(results_dir, 'dir')
    error('Results directory not found:\n  %s', results_dir);
end

fprintf('=== Loading PSA Object ===\n');
fprintf('Directory: %s\n', results_dir);
% from_dir reads psa_object.mat, fills in the per-condition results when the
% object does not already carry them, and errors clearly if the run is absent.
psa = ParamSpaceAnalysis2.from_dir(results_dir);
fprintf('Loaded PSA object: %d combinations, %d conditions\n\n', ...
    psa.num_combinations, length(psa.conditions));

%% Plot results
[~, figs_hist] = load_and_make_unit_histograms(psa.output_dir, 'Metrics', metrics_to_hist);
fig_paired_swarm = load_and_plot_lle_by_stim_period(psa.output_dir, ...
    'transient_skip', transient_skip, 'periods_to_plot', periods_to_plot);

% Combine into 3x4 layout
fig_combined = concatenate_figs([figs_hist, fig_paired_swarm], 'vertical', ...
    'HideTitlesAfterFirstRow', true);

%% Statistical Analysis: Stim vs No-Stim (Wilcoxon signed-rank)
condition_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
num_conditions = length(condition_names);

condition_titles_map = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

if isfield(psa.model_defaults, 'input_config')
    n_steps = psa.model_defaults.input_config.n_steps;
    no_stim_pattern = psa.model_defaults.input_config.no_stim_pattern;
else
    n_steps = 3;
    no_stim_pattern = false(1, 3);
    no_stim_pattern(1:2:end) = true;
end

if isfield(psa.model_defaults, 'T_range')
    T_stim = psa.model_defaults.T_range(2);
else
    T_stim = 50;
end
step_period = T_stim / n_steps;

stats_lines = {};
stats_lines{end+1} = '=== Statistical Analysis: Stim vs No-Stim ===';
stats_lines{end+1} = sprintf('Data directory: %s', results_dir);
stats_lines{end+1} = sprintf('Test: Wilcoxon signed-rank (paired, non-parametric)');
stats_lines{end+1} = sprintf('Transient skip: %.2f s', transient_skip);
stats_lines{end+1} = sprintf('Step period: %.2f s, No-stim pattern: [%s]', ...
    step_period, strjoin(string(no_stim_pattern), ', '));
stats_lines{end+1} = '';

for c_idx = 1:num_conditions
    cond_name = condition_names{c_idx};
    results_cell = psa.results.(cond_name);

    n_valid = 0;
    stim_means_all = [];
    no_stim_means_all = [];

    for k = 1:length(results_cell)
        res = results_cell{k};
        if ~isstruct(res) || ~isfield(res, 'success') || ~res.success
            continue;
        end
        if ~isfield(res, 'local_lya') || isempty(res.local_lya)
            continue;
        end

        t_lya = res.t_lya;
        local_lya = res.local_lya;

        valid_mask = t_lya >= 0;
        t_lya = t_lya(valid_mask);
        local_lya = local_lya(valid_mask);

        step_means = NaN(1, n_steps);
        for step_idx = 1:n_steps
            step_start = (step_idx - 1) * step_period + transient_skip;
            step_end = step_idx * step_period;
            step_mask = t_lya >= step_start & t_lya < step_end;
            if any(step_mask)
                step_means(step_idx) = mean(local_lya(step_mask), 'omitnan');
            end
        end

        stim_mean = mean(step_means(~no_stim_pattern), 'omitnan');
        no_stim_mean = mean(step_means(no_stim_pattern), 'omitnan');

        if ~isnan(stim_mean) && ~isnan(no_stim_mean)
            n_valid = n_valid + 1;
            stim_means_all(end+1) = stim_mean; %#ok<SAGROW>
            no_stim_means_all(end+1) = no_stim_mean; %#ok<SAGROW>
        end
    end

    if condition_titles_map.isKey(cond_name)
        display_name = condition_titles_map(cond_name);
    else
        display_name = strrep(cond_name, '_', ' ');
    end

    if n_valid >= 2
        stim_means_all = stim_means_all(:);
        no_stim_means_all = no_stim_means_all(:);

        [p_value, ~, ~] = signrank(stim_means_all, no_stim_means_all);

        differences = stim_means_all - no_stim_means_all;
        cohens_d = mean(differences) / std(differences);

        stats_lines{end+1} = sprintf('%s (n=%d pairs):', display_name, n_valid); %#ok<SAGROW>
        stats_lines{end+1} = sprintf('  Wilcoxon signed-rank p-value: %.4g', p_value); %#ok<SAGROW>
        stats_lines{end+1} = sprintf('  Cohen''s d: %.4f', cohens_d); %#ok<SAGROW>
        stats_lines{end+1} = sprintf('  Median LLE difference (stim - no-stim): %.4f', ...
            median(stim_means_all) - median(no_stim_means_all)); %#ok<SAGROW>
        stats_lines{end+1} = ''; %#ok<SAGROW>
    else
        stats_lines{end+1} = sprintf('%s: Insufficient valid pairs (n=%d)', ...
            display_name, n_valid); %#ok<SAGROW>
        stats_lines{end+1} = ''; %#ok<SAGROW>
    end
end

% Print stats
fprintf('\n');
for i = 1:length(stats_lines)
    fprintf('%s\n', stats_lines{i});
end

% Write stats to text file
if save_stats
    if save_figs
        stats_dir = fullfile(figs_root, 'fraction_excitatory_analysis');
    else
        stats_dir = results_dir;
    end
    if ~exist(stats_dir, 'dir')
        mkdir(stats_dir);
    end
    stats_file = fullfile(stats_dir, 'stats_results.txt');
    fid = fopen(stats_file, 'w');
    if fid == -1
        warning('Could not open stats file for writing: %s', stats_file);
    else
        fprintf(fid, 'Generated: %s\n\n', string(datetime('now')));
        for i = 1:length(stats_lines)
            fprintf(fid, '%s\n', stats_lines{i});
        end
        fclose(fid);
        fprintf('Stats written to: %s\n', stats_file);
    end
end

%% Panel letters
try
    drawnow
    AddLetters2Plots(fig_combined, ...
        {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)'}, ...
        'FontSize', 14, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.04);
catch ME
    warning('AddLetters2Plots failed (non-critical): %s', ME.message);
end

%% Save figures
if save_figs
    save_dir = fullfile(figs_root, 'fraction_excitatory_analysis');
    save_some_figs_to_folder_2(save_dir, 'fraction_excitatory', [], {'fig', 'svg', 'png', 'jp2'});
    fprintf('Figures saved to %s\n', save_dir);
end

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('Output directory: %s\n', psa.output_dir);
fprintf('Grid parameters: %s\n', strjoin(psa.grid_params, ', '));
fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false), ', '));
fprintf('\nDone! Plots and stats generated from: %s\n', psa.output_dir);

function load_and_plot_lle_by_stim_period(results_dir, options)
% LOAD_AND_PLOT_LLE_BY_STIM_PERIOD Plot mean local LLE per stim step
%
% Usage:
%   load_and_plot_lle_by_stim_period('/path/to/param_space_...')
%   load_and_plot_lle_by_stim_period(results_dir, 'transient_skip', 1.0)
%
% This function:
%   1. Loads the PSA object from psa_object.mat
%   2. Validates that local_lya time series exist
%   3. Computes mean local LLE for each step (using t >= 0)
%   4. Plots lines for each simulation with 'o' markers
%
% Input:
%   results_dir     - Path to a param_space_* output directory
%   options         - Name-value pairs:
%       'transient_skip' - Seconds to skip at start of each step (default: 0)
%
% See also: ParamSpaceAnalysis, load_and_plot_param_space_analysis

arguments
    results_dir (1,:) char
    options.transient_skip (1,1) double = 0
end

if ~exist(results_dir, 'dir')
    error('load_and_plot_lle_by_stim_period:NotFound', ...
        'Directory not found: %s', results_dir);
end

fprintf('=== Loading PSA for LLE by Stim Period ===\n');
fprintf('Directory: %s\n', results_dir);
fprintf('Transient skip: %.2f s\n\n', options.transient_skip);

%% Load PSA object
psa_file = fullfile(results_dir, 'psa_object.mat');
if ~exist(psa_file, 'file')
    % Also check for test file name
    psa_file = fullfile(results_dir, 'psa_object_test.mat');
end
if ~exist(psa_file, 'file')
    error('load_and_plot_lle_by_stim_period:NoObject', ...
        'psa_object.mat not found in: %s', results_dir);
end

loaded = load(psa_file);
psa = loaded.psa;

%% Validate local_lya exists
% Check first condition's first successful result
cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
has_local_lya = false;

for c_idx = 1:length(cond_names)
    cond_name = cond_names{c_idx};
    if isfield(psa.results, cond_name)
        for k = 1:length(psa.results.(cond_name))
            res = psa.results.(cond_name){k};
            if isstruct(res) && isfield(res, 'local_lya') && ~isempty(res.local_lya)
                has_local_lya = true;
                break;
            end
        end
    end
    if has_local_lya, break; end
end

if ~has_local_lya
    error('load_and_plot_lle_by_stim_period:NoLocalLya', ...
        ['Results do not contain local_lya time series.\n' ...
        'Re-run PSA with store_local_lya = true.']);
end

%% Extract step configuration from model_defaults
if isfield(psa.model_defaults, 'input_config')
    n_steps = psa.model_defaults.input_config.n_steps;
    no_stim_pattern = psa.model_defaults.input_config.no_stim_pattern;
else
    % Use SRNNModel defaults
    n_steps = 3;
    no_stim_pattern = false(1, 3);
    no_stim_pattern(1:2:end) = true;
end

% Get T_stim (positive time portion)
if isfield(psa.model_defaults, 'T_range')
    T_stim = psa.model_defaults.T_range(2);
else
    T_stim = 50;  % Default
end

step_period = T_stim / n_steps;
fprintf('Step configuration: %d steps, period = %.2f s\n', n_steps, step_period);
fprintf('No-stim pattern: [%s]\n\n', strjoin(string(no_stim_pattern), ', '));

%% Compute mean LLE per step for each simulation
num_conditions = length(cond_names);

% Pre-allocate storage
all_means = cell(num_conditions, 1);

for c_idx = 1:num_conditions
    cond_name = cond_names{c_idx};
    results_cell = psa.results.(cond_name);

    % Count successful results with local_lya
    n_valid = 0;
    for k = 1:length(results_cell)
        res = results_cell{k};
        if isstruct(res) && isfield(res, 'success') && res.success && ...
                isfield(res, 'local_lya') && ~isempty(res.local_lya)
            n_valid = n_valid + 1;
        end
    end

    means_matrix = NaN(n_valid, n_steps);
    sim_idx = 0;

    for k = 1:length(results_cell)
        res = results_cell{k};
        if ~isstruct(res) || ~isfield(res, 'success') || ~res.success
            continue;
        end
        if ~isfield(res, 'local_lya') || isempty(res.local_lya)
            continue;
        end

        sim_idx = sim_idx + 1;
        t_lya = res.t_lya;
        local_lya = res.local_lya;

        % Truncate to t >= 0
        valid_mask = t_lya >= 0;
        t_lya = t_lya(valid_mask);
        local_lya = local_lya(valid_mask);

        % Compute mean for each step
        for step_idx = 1:n_steps
            step_start = (step_idx - 1) * step_period + options.transient_skip;
            step_end = step_idx * step_period;

            step_mask = t_lya >= step_start & t_lya < step_end;

            if any(step_mask)
                means_matrix(sim_idx, step_idx) = mean(local_lya(step_mask), 'omitnan');
            end
        end
    end

    all_means{c_idx} = means_matrix;
    fprintf('Condition %s: %d valid simulations\n', cond_name, sim_idx);
end

%% Create figure
condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

fig = figure('Name', 'Mean Local LLE by Step', ...
    'Position', [100, 100, 300 * num_conditions, 350]);

% Get lines colormap
cmap_lines = lines;

for c_idx = 1:num_conditions
    cond_name = cond_names{c_idx};
    means_matrix = all_means{c_idx};

    ax = subplot(1, num_conditions, c_idx);

    n_sims = size(means_matrix, 1);

    % Build color matrix for each simulation
    sim_colors = zeros(n_sims, 3);
    for sim_idx = 1:n_sims
        color_idx = mod(sim_idx - 1, size(cmap_lines, 1)) + 1;
        sim_colors(sim_idx, :) = cmap_lines(color_idx, :);
    end

    % Build labels for steps
    step_labels = cell(1, n_steps);
    for step_idx = 1:n_steps
        if no_stim_pattern(step_idx)
            step_labels{step_idx} = 'no-stim';
        else
            step_labels{step_idx} = 'stim';
        end
    end

    % Use paired_beeswarm for beeswarm chart with connecting lines
    paired_beeswarm(means_matrix, 'Axes', ax, ...
        'Colors', sim_colors, ...
        'MarkerSize', 0.8, ...
        'LineWidth', 0.5, ...
        'Labels', step_labels, ...
        'Alpha', 0.7, ...
        'SortStyle', 'hex', ...
        'ShowYAxis', c_idx == 1);

    % Labels
    if condition_titles.isKey(cond_name)
        title(condition_titles(cond_name), 'FontWeight', 'normal');
    else
        title(strrep(cond_name, '_', ' '), 'FontWeight', 'normal');
    end

    if c_idx == 1
        ylabel('Mean Local LLE');
    end

    xtickangle(ax, 45);
    box(ax, 'off');

end

% Link y-axes and set tight ylims
drawnow;
ax_handles = findobj(fig, 'Type', 'Axes');
linkaxes(ax_handles, 'y');
axis(ax_handles, 'tight');
% Restore x-limits since axis tight affects both x and y
for i = 1:length(ax_handles)
    xlim(ax_handles(i), [0.5, n_steps + 0.5]);
end

%% Save figure
fig_dir = fullfile(results_dir, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
saveas(fig, fullfile(fig_dir, 'lle_by_stim_period.png'));
saveas(fig, fullfile(fig_dir, 'lle_by_stim_period.fig'));
fprintf('\nFigure saved to: %s\n', fig_dir);

fprintf('\nDone!\n');
end

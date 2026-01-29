function load_and_plot_lle_by_stim_period(results_dir, options)
% LOAD_AND_PLOT_LLE_BY_STIM_PERIOD Plot mean local LLE per stim step
%
% Usage:
%   load_and_plot_lle_by_stim_period('/path/to/param_space_...')
%   load_and_plot_lle_by_stim_period(results_dir, 'transient_skip', 1.0)
%   load_and_plot_lle_by_stim_period(results_dir, 'periods_to_plot', [1 1 0])
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
%       'transient_skip'  - Seconds to skip at start of each step (default: 0)
%       'periods_to_plot' - Logical vector selecting which periods to plot
%                           e.g., [1 1 0] plots first 2 periods but not 3rd
%                           (default: all periods)
%
% See also: ParamSpaceAnalysis, load_and_plot_param_space_analysis

arguments
    results_dir (1,:) char
    options.transient_skip (1,1) double = 0
    options.periods_to_plot (1,:) logical = logical([])
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
all_f_values = cell(num_conditions, 1);  % Store f values for coloring

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
    f_values = NaN(n_valid, 1);  % Store f values for each simulation
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

        % Extract f value if available
        if isfield(res, 'config') && isfield(res.config, 'f')
            f_values(sim_idx) = res.config.f;
        elseif isfield(psa.model_defaults, 'f')
            f_values(sim_idx) = psa.model_defaults.f;
        end

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
    all_f_values{c_idx} = f_values;
    fprintf('Condition %s: %d valid simulations\n', cond_name, sim_idx);
end

%% Determine global f value range for consistent coloring
% Collect all f values across conditions
all_f_combined = [];
for c_idx = 1:num_conditions
    all_f_combined = [all_f_combined; all_f_values{c_idx}];
end
all_f_combined = all_f_combined(~isnan(all_f_combined));

if ~isempty(all_f_combined)
    f_min = min(all_f_combined);
    f_max = max(all_f_combined);
    has_f_variation = (f_max > f_min);
else
    f_min = 0;
    f_max = 1;
    has_f_variation = false;
end

if has_f_variation
    fprintf('Coloring by f value: [%.3f, %.3f]\n', f_min, f_max);
end

%% Create figure
condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

fig = figure('Name', 'Mean Local LLE by Step', ...
    'Position', [100, 100, 300 * num_conditions, 350]);

% Get colormap for f-value coloring (blue=low f, red=high f, gray middle)
cmap_f = blue_gray_red_colormap(256);
n_colors = size(cmap_f, 1);

for c_idx = 1:num_conditions
    cond_name = cond_names{c_idx};
    means_matrix = all_means{c_idx};
    f_values = all_f_values{c_idx};

    ax = subplot(1, num_conditions, c_idx);

    n_sims = size(means_matrix, 1);

    % Build color matrix for each simulation based on f value
    sim_colors = zeros(n_sims, 3);
    for sim_idx = 1:n_sims
        if has_f_variation && ~isnan(f_values(sim_idx))
            % Map f value to colormap index
            f_normalized = (f_values(sim_idx) - f_min) / (f_max - f_min);
            f_normalized = max(0, min(1, f_normalized));  % Clamp to [0, 1]
            color_idx = round(f_normalized * (n_colors - 1)) + 1;
            sim_colors(sim_idx, :) = cmap_f(color_idx, :);
        else
            % Default gray if no f variation or missing f value
            sim_colors(sim_idx, :) = [0.5 0.5 0.5];
        end
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

    % Apply periods_to_plot filter
    if isempty(options.periods_to_plot)
        % Default: plot all periods
        periods_mask = true(1, n_steps);
    else
        if length(options.periods_to_plot) ~= n_steps
            error('load_and_plot_lle_by_stim_period:InvalidPeriodsToPlot', ...
                'periods_to_plot length (%d) must match n_steps (%d)', ...
                length(options.periods_to_plot), n_steps);
        end
        periods_mask = options.periods_to_plot;
    end

    % Filter data to selected periods
    means_matrix_filtered = means_matrix(:, periods_mask);
    step_labels_filtered = step_labels(periods_mask);
    n_steps_plot = sum(periods_mask);

    % Add dashed black line at y=0 (stability reference) - plot first so it's behind data
    line(ax, [0.3, n_steps_plot + 0.7], [0, 0], ...
        'Color', [0 0 0 0.5], ...  % Black with 50% opacity
        'LineStyle', '--', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off');  % Don't show in legend
    hold(ax, 'on');

    % Use paired_beeswarm for beeswarm chart with connecting lines
    % Only show y-axis for first subplot since all are linked
    paired_beeswarm(means_matrix_filtered, 'Axes', ax, ...
        'Colors', sim_colors, ...
        'MarkerSize', 0.9, ...
        'LineWidth', 0.75, ...
        'Labels', step_labels_filtered, ...
        'Alpha', 0.7, ...
        'SortStyle', 'hex', ...
        'ShowYAxis', (c_idx == 1));

    % Labels
    if condition_titles.isKey(cond_name)
        title(condition_titles(cond_name), 'FontWeight', 'normal');
    else
        title(strrep(cond_name, '_', ' '), 'FontWeight', 'normal');
    end

    if c_idx == 1
        ylabel('Mean Local LLE');
    else
        % Hide y-axis for subplots 2-4 (same scale as first)
        ax.YAxis.Visible = 'off';
    end

    hold(ax, 'off');

    xtickangle(ax, 45);
    box(ax, 'off');

end

% Link y-axes for ALL conditions (same scale across all subplots)
drawnow;
ax_handles = findobj(fig, 'Type', 'Axes');

% Link ALL axes together
if length(ax_handles) >= 2
    linkaxes(ax_handles, 'y');
end

% Set x-limits for all axes (use n_steps_plot for filtered count)
for i = 1:length(ax_handles)
    xlim(ax_handles(i), [0.5, n_steps_plot + 0.5]);
end

% Compute global y-limits across ALL conditions
global_ymin = inf;
global_ymax = -inf;
for i = 1:length(ax_handles)
    children = ax_handles(i).Children;
    for j = 1:length(children)
        if isprop(children(j), 'YData')
            ydata = children(j).YData;
            ydata = ydata(isfinite(ydata));
            if ~isempty(ydata)
                global_ymin = min(global_ymin, min(ydata));
                global_ymax = max(global_ymax, max(ydata));
            end
        end
    end
end

% Compute padding and limits
global_range = global_ymax - global_ymin;
if global_range > 0
    global_padding = 0.05 * global_range;
else
    global_padding = 0.1;
end

% Ensure upper limit is at least 0.05 so 0 is visible
global_upper = max(global_ymax + global_padding, 0.05);
global_lower = global_ymin - global_padding;

% Set ylim on one of the linked axes (will apply to all linked)
if ~isempty(ax_handles)
    ylim(ax_handles(1), [global_lower, global_upper]);
end

%% Save figure
fig_dir = fullfile(results_dir, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
saveas(fig, fullfile(fig_dir, 'lle_by_stim_period.png'));
saveas(fig, fullfile(fig_dir, 'lle_by_stim_period.fig'));
fprintf('\nFigure saved to: %s\n', fig_dir);

%% Create separate colorbar figure for f values
if has_f_variation
    fig_cb = figure('Name', 'f Value Colorbar', ...
        'Position', [500, 100, 150, 350]);

    % Create a dummy image to get a colorbar
    ax_cb = axes(fig_cb);

    % Create gradient image (low values at bottom, high at top)
    n_colors = 256;
    gradient_img = repmat(linspace(0, 1, n_colors)', 1, 2);
    imagesc(ax_cb, [0 1], [f_min f_max], gradient_img);
    colormap(ax_cb, cmap_f);

    % Configure axes
    ax_cb.XTick = [];
    ax_cb.YDir = 'normal';
    ylabel(ax_cb, 'f (Fraction Excitatory)', 'FontSize', 12);
    box(ax_cb, 'off');

    % Set aspect ratio
    pbaspect(ax_cb, [0.3 1 1]);

    % Save colorbar figure
    saveas(fig_cb, fullfile(fig_dir, 'lle_by_stim_period_colorbar.png'));
    saveas(fig_cb, fullfile(fig_dir, 'lle_by_stim_period_colorbar.fig'));
    fprintf('Colorbar figure saved to: %s\n', fig_dir);
end

fprintf('\nDone!\n');
end


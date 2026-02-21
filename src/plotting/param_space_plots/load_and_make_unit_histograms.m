function [psa, fig_handles] = load_and_make_unit_histograms(results_dir, options)
% LOAD_AND_MAKE_UNIT_HISTOGRAMS Load PSA results and create unit histograms colored by f
%
% Usage:
%   psa = load_and_make_unit_histograms('/path/to/param_space_...')
%   psa = load_and_make_unit_histograms('/path/to/param_space_...', 'NormalizeMode', 'probability')
%   [psa, figs] = load_and_make_unit_histograms(..., 'Metrics', {'lle', 'r', 'br'})
%
% This function:
%   1. Loads the PSA object from psa_object.mat if available
%   2. Otherwise, creates a new object and loads results
%   3. Checks if consolidation is needed (temp_batches exists)
%   4. Consolidates if necessary
%   5. Plots distributions using unit_histogram_patch with f-value coloring
%
% Input:
%   results_dir   - Path to a param_space_* output directory
%
% Options:
%   NormalizeMode - 'count' (default) or 'probability'
%   Metrics       - Cell array of metrics to plot: 'lle', 'r', 'br'
%                   (default: {'lle', 'r'})
%
% Output:
%   psa         - The loaded ParamSpaceAnalysis object
%   fig_handles - Array of figure handles (one per metric)
%
% See also: ParamSpaceAnalysis2, unit_histogram_patch, load_and_plot_lle_by_stim_period

arguments
    results_dir (1,:) char
    options.NormalizeMode (1,:) char {mustBeMember(options.NormalizeMode, {'count', 'probability'})} = 'count'
    options.Metrics (1,:) cell = {'lle', 'r'}
end

normalize_mode = options.NormalizeMode;
metrics_to_plot = lower(options.Metrics);

if ~exist(results_dir, 'dir')
    error('load_and_make_unit_histograms:NotFound', ...
        'Directory not found: %s', results_dir);
end

fprintf('=== Loading Parameter Space Analysis for Unit Histograms ===\n');
fprintf('Directory: %s\n\n', results_dir);

%% Load PSA object (same as load_and_plot_param_space_analysis)
psa_file = fullfile(results_dir, 'psa_object.mat');
if exist(psa_file, 'file')
    fprintf('Loading PSA object from: %s\n', psa_file);
    loaded = load(psa_file);
    psa = loaded.psa;
    fprintf('Loaded PSA object successfully.\n');
else
    fprintf('No psa_object.mat found. Creating new object and loading results...\n');
    psa = ParamSpaceAnalysis2();
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
    fprintf('Loading results from per-condition MAT files...\n');
    psa.load_results(results_dir);
end

%% Get condition info
condition_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
num_conditions = length(condition_names);

condition_titles = containers.Map(...
    {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});

%% Build metric configuration based on options
% Map short names to internal metric names and display properties
metric_config = struct();
metric_config.lle = struct('field', 'LLE', 'label', '\lambda_1', 'range', [-2.3, 1.5], 'inf_both', true);
metric_config.r = struct('field', 'mean_rate', 'label', 'Mean Firing Rate', 'range', [0, 1], 'inf_both', false);
metric_config.br = struct('field', 'mean_synaptic_output', 'label', 'Mean Synaptic Output', 'range', [0, 1], 'inf_both', false);

% Filter to requested metrics
metrics = {};
metric_labels = {};
metric_ranges = {};
metric_fields = {};
metric_inf_both = {};
for i = 1:length(metrics_to_plot)
    key = metrics_to_plot{i};
    if isfield(metric_config, key)
        cfg = metric_config.(key);
        metrics{end+1} = cfg.field; %#ok<AGROW>
        metric_labels{end+1} = cfg.label; %#ok<AGROW>
        metric_ranges{end+1} = cfg.range; %#ok<AGROW>
        metric_fields{end+1} = cfg.field; %#ok<AGROW>
        metric_inf_both{end+1} = cfg.inf_both; %#ok<AGROW>
    else
        warning('load_and_make_unit_histograms:UnknownMetric', ...
            'Unknown metric: %s. Valid options: lle, r, br', key);
    end
end

if isempty(metrics)
    error('load_and_make_unit_histograms:NoValidMetrics', ...
        'No valid metrics specified.');
end

n_bins = 25;  % Match psa.plot

% Precompute bin edges for each metric (consistent across conditions)
metric_edges = cell(1, length(metrics));
for m_idx = 1:length(metrics)
    hist_range = metric_ranges{m_idx};
    % Compute edges using linspace, add inf bins for overflow
    edges = linspace(hist_range(1), hist_range(2), n_bins + 1);
    if metric_inf_both{m_idx}
        edges = [-inf, edges, inf];  % Both sides for LLE-like metrics
    else
        edges = [edges, inf];  % Just upper overflow for rate-like metrics
    end
    metric_edges{m_idx} = edges;
end

% First pass: collect all f values across all conditions for global normalization
all_f_combined = [];
for c_idx = 1:num_conditions
    cond_name = condition_names{c_idx};
    if isfield(psa.results, cond_name)
        results_cell = psa.results.(cond_name);
        for k = 1:length(results_cell)
            res = results_cell{k};
            if isstruct(res) && isfield(res, 'success') && res.success
                if isfield(res, 'config') && isfield(res.config, 'f')
                    all_f_combined(end+1) = res.config.f;
                elseif isfield(psa.model_defaults, 'f')
                    all_f_combined(end+1) = psa.model_defaults.f;
                end
            end
        end
    end
end

% Compute global f limits
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

% Get colormap (same as beeswarm)
cmap_f = blue_gray_red_colormap(256);

%% Create figures for each metric
fig_handles = gobjects(1, length(metrics));  % Pre-allocate figure handle array

for m_idx = 1:length(metrics)
    metric = metrics{m_idx};
    metric_label = metric_labels{m_idx};
    metric_range = metric_ranges{m_idx};

    fprintf('\nPlotting %s distribution...\n', metric);

    % Create figure with matching dimensions
    fig = figure('Name', sprintf('%s Unit Histogram', metric), ...
        'Position', [100, 100, 300 * num_conditions, 300]);

    for c_idx = 1:num_conditions
        cond_name = condition_names{c_idx};
        ax = subplot(1, num_conditions, c_idx);

        % Extract metric values and f values
        values = [];
        f_values = [];

        if isfield(psa.results, cond_name)
            results_cell = psa.results.(cond_name);
            for k = 1:length(results_cell)
                res = results_cell{k};
                if isstruct(res) && isfield(res, 'success') && res.success
                    if isfield(res, metric) && ~isnan(res.(metric))
                        values(end+1) = res.(metric);

                        % Get f value
                        if isfield(res, 'config') && isfield(res.config, 'f')
                            f_values(end+1) = res.config.f;
                        elseif isfield(psa.model_defaults, 'f')
                            f_values(end+1) = psa.model_defaults.f;
                        else
                            f_values(end+1) = 0.5;  % Default
                        end
                    end
                end
            end
        end

        % Plot unit histogram
        if ~isempty(values)
            unit_histogram_patch(values(:), f_values(:), ...
                'BinEdges', metric_edges{m_idx}, ...
                'SortMode', 'sorted', ...
                'Axes', ax, ...
                'Colormap', cmap_f, ...
                'CLim', [f_min, f_max], ...
                'Normalize', normalize_mode, ...
                'EdgeColor', 'none');

            % Add stability reference line for LLE
            if strcmpi(metric, 'LLE')
                hold(ax, 'on');
                xline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 2);
                hold(ax, 'off');
            end
        end

        % Labels
        if condition_titles.isKey(cond_name)
            title(ax, condition_titles(cond_name), 'FontWeight', 'normal');
        else
            title(ax, strrep(cond_name, '_', ' '), 'FontWeight', 'normal');
        end

        if c_idx == 1
            if strcmpi(normalize_mode, 'probability')
                ylabel(ax, 'Probability');
            else
                ylabel(ax, 'Count');
            end
        end
        xlabel(ax, metric_label);

        % Set x limits based on metric
        xlim(ax, 1.05.*metric_range);
    end

    % Link y-axes
    drawnow;
    ax_handles = findobj(fig, 'Type', 'Axes');
    linkaxes(ax_handles, 'y');

    % Store figure handle
    fig_handles(m_idx) = fig;
end

%% Create colorbar figure for f values (same as beeswarm)
if has_f_variation
    fig_cb = figure('Name', 'f Value Colorbar', ...
        'Position', [500, 200, 300, 300]);

    ax_cb = axes(fig_cb);

    % Create gradient image (low values at bottom, high at top)
    n_colors = 256;
    gradient_img = repmat(linspace(0, 1, n_colors)', 1, 2);
    imagesc(ax_cb, [0 1], [f_min f_max], gradient_img);
    colormap(ax_cb, cmap_f);

    % Configure axes
    ax_cb.XTick = [];
    ax_cb.YDir = 'normal';
    ax_cb.XColor = 'none';  % Hide x-axis completely
    ylabel(ax_cb, 'fraction excitatory', 'FontSize', 12);
    box(ax_cb, 'off');

    % Set aspect ratio
    pbaspect(ax_cb, [0.1 1 1]);


end

fprintf('\nDone!\n');

end

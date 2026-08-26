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
%   ColorBy       - SRNNModel2 parameter to colour points by (default: 'f')
%   ColorFcn      - @(psa, res) -> scalar, used INSTEAD of ColorBy. For a colour
%                   axis that is not a model property -- e.g. the realized E:I
%                   weight balance, which has to be computed from a rebuilt W.
%                   ColorBy still names the axis for messages when this is set.
%   CLim          - [lo hi] fixed colour limits. Default [] derives them from the
%                   data, which makes the colorbar's extent depend on which
%                   networks happened to be sampled; pass an explicit range for a
%                   bar that means the same thing in every run.
%   ColorLabel    - y-label on the colorbar figure (default 'fraction excitatory')
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
    options.LLERange (1,2) double = [-2.3, 1.5]   % LLE histogram range (bin edges)
    options.ColorBy (1,:) char = 'f'
    options.ColorFcn = []
    options.CLim double = []
    options.ColorLabel (1,:) char = 'fraction excitatory'
end

normalize_mode = options.NormalizeMode;
metrics_to_plot = lower(options.Metrics);
color_by = options.ColorBy;

% One accessor for the colour value, so the two passes below cannot disagree
% about what they are colouring by.
if isempty(options.ColorFcn)
    color_of = @(psa, res) psa.effective_param(res, color_by);
else
    color_of = options.ColorFcn;
end

if ~isempty(options.CLim)
    assert(numel(options.CLim) == 2 && options.CLim(2) > options.CLim(1), ...
        'load_and_make_unit_histograms:BadCLim', ...
        'CLim must be [lo hi] with hi > lo; got %s.', mat2str(options.CLim));
end

if ~exist(results_dir, 'dir')
    error('load_and_make_unit_histograms:NotFound', ...
        'Directory not found: %s', results_dir);
end

fprintf('=== Loading Parameter Space Analysis for Unit Histograms ===\n');
fprintf('Directory: %s\n\n', results_dir);

%% Load the run
% from_dir handles what this used to do by hand: read psa_object.mat, fall back
% to the per-condition result files, and consolidate temp_batches/ for a run that
% was interrupted.
psa = ParamSpaceAnalysis2.from_dir(results_dir);

%% Get condition info
condition_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
num_conditions = length(condition_names);

% The shared map, not a local copy. This file held a SIXTH copy of the same
% four-entry table that srnn_condition_titles was created to replace, and it was
% missed when the other five were consolidated -- so the three regimes added
% since (sfa_only_oneTS, std_only_oneTS, sfa3_std1) fell through to their raw
% snake_case names here while reading properly in every other figure. Still
% isKey-guarded below: a saved run directory can name any condition it likes.
condition_titles = srnn_condition_titles();

%% Build metric configuration based on options
% Map short names to internal metric names and display properties
metric_config = struct();
metric_config.lle = struct('field', 'LLE', 'label', '\lambda_1', 'range', options.LLERange, 'inf_both', true);
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

% First pass: collect all colour-parameter values for global normalization
all_f_combined = [];
for c_idx = 1:num_conditions
    cond_name = condition_names{c_idx};
    if isfield(psa.results, cond_name)
        results_cell = psa.results.(cond_name);
        for k = 1:length(results_cell)
            res = results_cell{k};
            if isstruct(res) && isfield(res, 'success') && res.success
                all_f_combined(end+1) = color_of(psa, res); %#ok<AGROW>
            end
        end
    end
end

% Colour limits. An explicit CLim wins outright: derived limits make the bar's
% extent depend on which networks the sweep happened to sample, so a run that
% missed the extremes silently loses its end ticks and the bar stops being
% comparable between runs.
if ~isempty(options.CLim)
    f_min = options.CLim(1);
    f_max = options.CLim(2);
    has_f_variation = true;
    n_out = sum(all_f_combined < f_min | all_f_combined > f_max);
    if n_out > 0
        % Saying so matters: unit_histogram_patch clamps to CLim, so these
        % networks are drawn in the end colour and read as merely extreme
        % rather than off-scale.
        warning('load_and_make_unit_histograms:ColorOutsideCLim', ...
            ['%d of %d networks fall outside CLim [%g %g] (data range ' ...
             '[%.3f %.3f]) and are clamped to the end colours.'], ...
            n_out, numel(all_f_combined), f_min, f_max, ...
            min(all_f_combined), max(all_f_combined));
    end
elseif ~isempty(all_f_combined) && max(all_f_combined) > min(all_f_combined)
    % A constant colour parameter would give a degenerate CLim, so fall back to
    % the unit range as if none were found.
    f_min = min(all_f_combined);
    f_max = max(all_f_combined);
    has_f_variation = true;
else
    f_min = 0;
    f_max = 1;
    has_f_variation = false;
end

if has_f_variation
    fprintf('Coloring by %s: CLim [%.3f, %.3f], data [%.3f, %.3f]\n', ...
        color_by, f_min, f_max, min(all_f_combined), max(all_f_combined));
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
                        values(end+1) = res.(metric); %#ok<AGROW>
                        f_values(end+1) = color_of(psa, res); %#ok<AGROW>
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
    ylabel(ax_cb, options.ColorLabel, 'FontSize', 12);
    box(ax_cb, 'off');

    % Set aspect ratio
    pbaspect(ax_cb, [0.1 1 1]);


end

fprintf('\nDone!\n');

end

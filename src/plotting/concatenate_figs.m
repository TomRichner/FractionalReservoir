function fig_combined = concatenate_figs(fig_handles, direction, options)
% CONCATENATE_FIGS Combine multiple figures into a single tiled figure
%
% Usage:
%   fig_combined = concatenate_figs([fig1, fig2, fig3], 'vertical')
%   fig_combined = concatenate_figs([fig1, fig2], 'horizontal')
%   fig_combined = concatenate_figs(..., 'HideTitlesAfterFirstRow', true)
%
% This function copies the axes from multiple source figures into a new
% combined figure using tiledlayout. All source figures must have the
% same number of axes (subplots).
%
% Input:
%   fig_handles - Array of figure handles to combine
%   direction   - 'vertical' (stack rows) or 'horizontal' (stack columns)
%
% Options:
%   HideTitlesAfterFirstRow - If true, hide titles for rows 2+ (default: false)
%
% Output:
%   fig_combined - Handle to the new combined figure
%
% Example:
%   % Combine 3 figures each with 4 subplots into a 3x4 grid
%   fig_combined = concatenate_figs([fig1, fig2, fig3], 'vertical');
%
% See also: tiledlayout, copyobj

arguments
    fig_handles (1,:) matlab.ui.Figure
    direction (1,:) char {mustBeMember(direction, {'vertical', 'horizontal'})}
    options.HideTitlesAfterFirstRow (1,1) logical = false
    options.CloseSourceFigs (1,1) logical = false
end

%% Validate inputs
n_figs = length(fig_handles);
if n_figs < 2
    error('concatenate_figs:TooFewFigures', ...
        'At least 2 figures are required for concatenation.');
end

% Get axes from all figures (exclude colorbars, legends, etc.)
axes_per_fig = cell(n_figs, 1);
for i = 1:n_figs
    all_children = fig_handles(i).Children;
    % Filter to only Axes objects, excluding legends/colorbars
    axes_per_fig{i} = findobj(all_children, 'Type', 'Axes', '-not', 'Tag', 'Colorbar');
end

% Validate all figures have the same number of axes
n_axes = cellfun(@length, axes_per_fig);
if ~all(n_axes == n_axes(1))
    error('concatenate_figs:AxesMismatch', ...
        'All figures must have the same number of axes. Found: [%s]', ...
        strjoin(string(n_axes), ', '));
end
n_subplots = n_axes(1);

%% Determine grid dimensions
if strcmpi(direction, 'vertical')
    n_rows = n_figs;
    n_cols = n_subplots;
else  % horizontal
    n_rows = n_subplots;
    n_cols = n_figs;
end

%% Create combined figure
% Compute figure size based on original figure dimensions
first_fig_pos = fig_handles(1).Position;
if strcmpi(direction, 'vertical')
    combined_width = first_fig_pos(3);
    combined_height = first_fig_pos(4) * n_figs;
else
    combined_width = first_fig_pos(3) * n_figs;
    combined_height = first_fig_pos(4);
end

fig_combined = figure('Name', 'Combined Figure', ...
    'Position', [100, 100, combined_width, combined_height]);

% Create tiledlayout with no padding between tiles
t = tiledlayout(fig_combined, n_rows, n_cols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% Copy axes into combined figure
for fig_idx = 1:n_figs
    % Get axes in correct order (MATLAB returns children in reverse order)
    src_axes = flip(axes_per_fig{fig_idx});

    for ax_idx = 1:n_subplots
        src_ax = src_axes(ax_idx);

        % Calculate tile index based on direction
        if strcmpi(direction, 'vertical')
            tile_idx = (fig_idx - 1) * n_cols + ax_idx;
        else
            tile_idx = (ax_idx - 1) * n_cols + fig_idx;
        end

        % Create new axes in the tile
        new_ax = nexttile(t, tile_idx);

        % Copy all children from source axes (use allchild to include
        % objects with HandleVisibility='off', like significance brackets)
        children = allchild(src_ax);
        if ~isempty(children)
            copyobj(children, new_ax);
        end

        % Copy axis properties
        new_ax.XLim = src_ax.XLim;
        new_ax.YLim = src_ax.YLim;
        new_ax.XTick = src_ax.XTick;
        new_ax.YTick = src_ax.YTick;
        new_ax.XTickLabel = src_ax.XTickLabel;
        new_ax.YTickLabel = src_ax.YTickLabel;
        new_ax.XLabel.String = src_ax.XLabel.String;
        new_ax.YLabel.String = src_ax.YLabel.String;

        % Copy or hide title based on row position
        if options.HideTitlesAfterFirstRow && fig_idx > 1
            new_ax.Title.String = '';
        else
            new_ax.Title.String = src_ax.Title.String;
            new_ax.Title.FontWeight = src_ax.Title.FontWeight;
        end
        new_ax.Box = src_ax.Box;
        new_ax.XDir = src_ax.XDir;
        new_ax.YDir = src_ax.YDir;
        new_ax.ColorOrder = src_ax.ColorOrder;
        new_ax.XAxis.Visible = src_ax.XAxis.Visible;
        new_ax.YAxis.Visible = src_ax.YAxis.Visible;

        % Copy tick angle if set
        if isprop(src_ax, 'XTickLabelRotation')
            new_ax.XTickLabelRotation = src_ax.XTickLabelRotation;
        end
    end
end

%% Optionally close original figures
if options.CloseSourceFigs
    for i = 1:n_figs
        close(fig_handles(i));
    end
end

fprintf('Created combined figure: %d rows x %d cols\n', n_rows, n_cols);

end

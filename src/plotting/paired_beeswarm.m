function [hs, hl] = paired_beeswarm(data, varargin)
% PAIRED_BEESWARM Plot beeswarm chart with paired lines connecting conditions
%
%   [hs, hl] = paired_beeswarm(data)
%   [hs, hl] = paired_beeswarm(data, 'Name', Value, ...)
%
%   INPUTS:
%       data        - N x C matrix where N = number of subjects/observations
%                     and C = number of conditions. Each row represents a
%                     paired observation across conditions.
%
%   NAME-VALUE PAIRS:
%       'Axes'          - Handle to axes for plotting (default: gca)
%       'Colors'        - N x 3 RGB matrix or 1 x 3 RGB vector for coloring
%                         dots and lines. If N x 3, each row colors that
%                         subject's dot/line. If 1 x 3, all use same color.
%                         (default: uses 'lines' colormap per condition)
%       'Filled'        - Logical, true for filled markers, false for open
%                         (default: true)
%       'MarkerSize'    - Relative size of markers (default: 1)
%       'LineWidth'     - Width of connecting lines (default: 0.5)
%       'Labels'        - Cell array of condition labels (default: {'1','2',...})
%       'Alpha'         - Transparency for markers and lines (default: 0.5)
%       'ShowYAxis'     - Logical, true to show y-axis (default: true)
%       'SortStyle'     - Beeswarm sort style: 'nosort', 'up', 'down', 'fan',
%                         'rand', 'square', 'hex' (default: 'nosort')
%       'CorralStyle'   - How to handle points outside channel: 'none',
%                         'gutter', 'omit', 'random' (default: 'none')
%
%   OUTPUTS:
%       hs - Cell array of scatter handles (one per condition)
%       hl - Array of line handles for paired connections
%
%   EXAMPLE:
%       N = 50;
%       data = [randn(N,1), randn(N,1) + 0.5];
%       colors = parula(N);
%       paired_beeswarm(data, 'Colors', colors);
%
%   See also BEESWARM, LINE

% Parse inputs
p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'Axes', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
addParameter(p, 'Colors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2) == 3));
addParameter(p, 'Filled', true, @islogical);
addParameter(p, 'MarkerSize', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'LineWidth', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Labels', {}, @iscell);
addParameter(p, 'Alpha', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(p, 'ShowYAxis', true, @islogical);
addParameter(p, 'SortStyle', 'nosort', @ischar);
addParameter(p, 'CorralStyle', 'none', @ischar);
parse(p, data, varargin{:});

opts = p.Results;

% Get dimensions
[N, C] = size(data);

% Setup axes
if isempty(opts.Axes)
    ax = gca;
else
    ax = opts.Axes;
end
hold(ax, 'on');

% Generate default labels if not provided
if isempty(opts.Labels)
    opts.Labels = arrayfun(@num2str, 1:C, 'UniformOutput', false);
end

% Determine color mode and setup colors
usePerSubjectColors = false;
if isempty(opts.Colors)
    % Default: use lines colormap per condition (beeswarm default behavior)
    conditionColors = lines(C);
elseif size(opts.Colors, 1) == 1
    % Single color for all
    conditionColors = repmat(opts.Colors, C, 1);
elseif size(opts.Colors, 1) == N
    % Per-subject colors
    usePerSubjectColors = true;
    subjectColors = opts.Colors;
    conditionColors = lines(C);  % Still need for beeswarm's colormap
else
    error('Colors must be 1x3 (single color), Nx3 (per subject), or empty');
end

% Prepare data for beeswarm
% Create x position vector (1, 2, ..., C) and y vector
xPos = repmat(1:C, N, 1);
xPos = xPos(:);
yVals = data(:);
pairID = repmat((1:N)', C, 1);

% Call beeswarm to get the jittered x positions
% beeswarm returns the jittered x coordinates
% Note: use_current_axes=false allows beeswarm to set up proper axis limits,
% which is required for 'hex' and 'square' sort styles to work correctly
xJittered = beeswarm(xPos, yVals, ...
    'sort_style', opts.SortStyle, ...
    'corral_style', opts.CorralStyle, ...
    'dot_size', opts.MarkerSize, ...
    'use_current_axes', false, ...
    'colormap', conditionColors, ...
    'MarkerFaceAlpha', opts.Alpha, ...
    'MarkerEdgeColor', 'none');

% Now we need to replot if using per-subject colors, since beeswarm plots by condition
% First, get the scatter handles that beeswarm created
hs_beeswarm = findobj(ax, 'Type', 'Scatter');
hs = cell(C, 1);
for i = 1:min(length(hs_beeswarm), C)
    hs{i} = hs_beeswarm(C - i + 1);  % They're in reverse order
end

% If using per-subject colors, we need to delete beeswarm's plots and replot
if usePerSubjectColors
    delete(hs_beeswarm);
    hs = cell(N, 1);

    % Plot each point individually with its own color
    for k = 1:N
        for c = 1:C
            idx = (c-1)*N + k;
            if ~isnan(xJittered(idx)) && ~isnan(yVals(idx))
                hs{k} = scatter(ax, xJittered(idx), yVals(idx), ...
                    opts.MarkerSize * 36, subjectColors(k,:), 'filled', ...
                    'MarkerFaceAlpha', opts.Alpha, 'MarkerEdgeColor', 'none');
            end
        end
    end
end

% Draw paired lines
hl = gobjects(N, 1);
for k = 1:N
    idx = (pairID == k) & ~isnan(yVals) & ~isnan(xJittered);
    if nnz(idx) < 2
        continue;
    end
    xs = xJittered(idx);
    ys = yVals(idx);

    % Sort by x position to ensure correct line drawing
    [xs, ord] = sort(xs);
    ys = ys(ord);

    if usePerSubjectColors
        lineColor = subjectColors(k, :);
    elseif ~isempty(opts.Colors) && size(opts.Colors, 1) == 1
        lineColor = opts.Colors;
    else
        lineColor = [0.5 0.5 0.5];  % Gray for default case
    end

    if opts.Alpha < 1
        lineColor = [lineColor, opts.Alpha];  % RGBA
    end
    hl(k) = line(ax, xs, ys, 'LineWidth', opts.LineWidth, 'Color', lineColor);
end

% Put points on top of lines
children = ax.Children;
scatterHandles = findobj(children, 'Type', 'Scatter');
for i = 1:length(scatterHandles)
    uistack(scatterHandles(i), 'top');
end

% Set axis labels
ax.XTick = 1:C;
ax.XTickLabel = opts.Labels;
xlim(ax, [0.5, C + 0.5]);
box(ax, 'on');

% Optionally hide y-axis
if ~opts.ShowYAxis
    ax.YAxis.Visible = 'off';
end

hold(ax, 'off');

end

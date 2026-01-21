function [hs, hl] = paired_swarm(data, varargin)
% PAIRED_SWARM Plot swarm chart with paired lines connecting conditions
%
%   [hs, hl] = paired_swarm(data)
%   [hs, hl] = paired_swarm(data, 'Name', Value, ...)
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
%                         (default: [0 0 0] black)
%       'Filled'        - Logical, true for filled markers, false for open
%                         (default: true)
%       'MarkerSize'    - Size of swarm markers in points^2 (default: 36)
%       'LineWidth'     - Width of connecting lines (default: 0.5)
%       'XJitterWidth'  - Width of x-axis jitter (default: 0.3)
%       'Labels'        - Cell array of condition labels (default: {'1','2',...})
%       'Alpha'         - Transparency for markers and lines (default: 1)
%       'ShowYAxis'     - Logical, true to show y-axis (default: true)
%
%   OUTPUTS:
%       hs - Handle to swarmchart object
%       hl - Array of line handles for paired connections
%
%   EXAMPLE:
%       N = 50;
%       data = [randn(N,1), randn(N,1) + 0.5];
%       colors = parula(N);
%       paired_swarm(data, 'Colors', colors, 'Filled', true);
%
%   See also SWARMCHART, LINE

% Parse inputs
p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'Axes', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
addParameter(p, 'Colors', [0 0 0], @(x) isnumeric(x) && size(x,2) == 3);
addParameter(p, 'Filled', true, @islogical);
addParameter(p, 'MarkerSize', 36, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'LineWidth', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'XJitterWidth', 0.3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Labels', {}, @iscell);
addParameter(p, 'Alpha', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(p, 'ShowYAxis', true, @islogical);
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

% Expand colors if single color provided
if size(opts.Colors, 1) == 1
    colors = repmat(opts.Colors, N, 1);
else
    colors = opts.Colors;
end

% Prepare data for swarmchart
% Reshape data to column vectors with x positions
xPos = repmat(1:C, N, 1);
xPos = xPos(:);
yVals = data(:);
pairID = repmat((1:N)', C, 1);

% Expand colors for all data points
colorVals = repmat(colors, C, 1);

% Create swarmchart
if opts.Filled
    hs = swarmchart(ax, xPos, yVals, opts.MarkerSize, colorVals, 'filled');
else
    hs = swarmchart(ax, xPos, yVals, opts.MarkerSize, colorVals);
end
hs.XJitter = 'density';
hs.XJitterWidth = opts.XJitterWidth;

% Set marker alpha if specified
if opts.Alpha < 1
    hs.MarkerFaceAlpha = opts.Alpha;
    hs.MarkerEdgeAlpha = opts.Alpha;
end

drawnow;  % Ensure jitter is computed

% Extract jittered coordinates
xyz = struct(hs).XYZJittered;  % undocumented, but works
xj = xyz(:,1);
yj = xyz(:,2);

% Draw paired lines
hl = gobjects(N, 1);
for k = 1:N
    idx = (pairID == k) & ~isnan(yj) & ~isnan(xj);
    if nnz(idx) < 2
        continue;
    end
    [xs, ord] = sort(xj(idx));
    ys = yj(idx);
    ys = ys(ord);

    lineColor = colors(k, :);
    if opts.Alpha < 1
        lineColor = [lineColor, opts.Alpha];  % RGBA
    end
    hl(k) = line(ax, xs, ys, 'LineWidth', opts.LineWidth, 'Color', lineColor);
end

% Put points on top of lines
uistack(hs, 'top');

% Set axis labels
ax.XTick = 1:C;
ax.XTickLabel = opts.Labels;
xlim(ax, [0.5, C + 0.5]);
box(ax, 'on');

% Optionally hide y-axis
if ~opts.ShowYAxis
    ax.YAxis.Visible = 'off';
end

end

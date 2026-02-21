function unit_histogram_patch(data, color_values, options)
% UNIT_HISTOGRAM_PATCH Creates a stacked histogram where every element is a patch.
%
% Usage:
%   unit_histogram_patch(data, color_values)
%   unit_histogram_patch(data, color_values, 'BinEdges', edges)
%   unit_histogram_patch(data, color_values, 'Axes', ax, 'Colormap', cmap)
%
% Inputs:
%   data         : Vector of x-axis distribution
%   color_values : Vector of values to determine color
%   options      : Name-value pairs:
%       'BinEdges'  - Explicit bin edges vector (required for consistent bins)
%       'NumBins'   - Number of bins (default 20, ignored if BinEdges provided)
%       'SortMode'  - 'sorted' (gradient) or 'random' (default 'sorted')
%       'Axes'      - Axes handle to plot on (default: new figure)
%       'Colormap'  - Colormap matrix (default: parula)
%       'CLim'      - [min max] color limits (default: [min max] of color_values)
%       'Normalize' - 'count' or 'probability' (default 'count')
%       'EdgeColor' - Edge color for patches (default 'none')

arguments
    data (:,1) double
    color_values (:,1) double
    options.BinEdges (1,:) double = []
    options.NumBins (1,1) double = 20
    options.SortMode (1,:) char {mustBeMember(options.SortMode, {'sorted', 'random'})} = 'sorted'
    options.Axes = []
    options.Colormap (:,3) double = parula(256)
    options.CLim (1,2) double = [NaN NaN]
    options.Normalize (1,:) char {mustBeMember(options.Normalize, {'count', 'probability'})} = 'count'
    options.EdgeColor = 'none'
end

% 1. Discretization - use BinEdges if provided, otherwise NumBins
if ~isempty(options.BinEdges)
    edges = options.BinEdges;
    [counts, ~] = histcounts(data, edges);
else
    [counts, edges] = histcounts(data, options.NumBins);
end

% Handle infinite edges for geometry (same approach as psa.plot)
% Find first non-infinite bin width to use as standard
finite_mask = ~isinf(edges);
finite_edges_only = edges(finite_mask);
if length(finite_edges_only) >= 2
    bin_width = finite_edges_only(2) - finite_edges_only(1);
else
    bin_width = 1;  % Fallback
end

% Create finite edges for plotting geometry
plot_edges = edges;
if isinf(plot_edges(1))
    plot_edges(1) = plot_edges(2) - bin_width;
end
if isinf(plot_edges(end))
    plot_edges(end) = plot_edges(end-1) + bin_width;
end

bin_indices = discretize(data, edges);

% 2. Sorting Logic
T = table(data, color_values, bin_indices);
T = T(~isnan(T.bin_indices), :); % Remove NaNs

if strcmpi(options.SortMode, 'random')
    T.rnd = rand(height(T), 1);
    T = sortrows(T, {'bin_indices', 'rnd'});
else
    % Default: Gradient effect
    T = sortrows(T, {'bin_indices', 'color_values'});
end

% 3. Calculate Stack Heights
[~, ~, unique_idx] = unique(T.bin_indices);
group_counts = accumarray(unique_idx, 1);

T.y_pos = zeros(height(T), 1);
current_idx = 1;
for i = 1:length(group_counts)
    n = group_counts(i);
    T.y_pos(current_idx : current_idx + n - 1) = (0 : n - 1)';
    current_idx = current_idx + n;
end

% 4. Apply probability normalization if requested
if strcmpi(options.Normalize, 'probability')
    total_count = height(T);
    T.y_pos = T.y_pos / total_count;
    y_height = 1 / total_count;
else
    y_height = 1;
end

% 5. Construct Geometry (using plot_edges instead of edges for positions)
x_lefts = plot_edges(T.bin_indices);
x_lefts = x_lefts(:)';

y_bottoms = T.y_pos(:)';

% Create 4xN matrices
X_verts = [x_lefts; x_lefts + bin_width; x_lefts + bin_width; x_lefts];
Y_verts = [y_bottoms; y_bottoms; y_bottoms + y_height; y_bottoms + y_height];

% 6. Map color values to colormap indices
if isnan(options.CLim(1))
    c_min = min(T.color_values);
    c_max = max(T.color_values);
else
    c_min = options.CLim(1);
    c_max = options.CLim(2);
end

% Normalize color values to [0, 1]
if c_max > c_min
    c_normalized = (T.color_values - c_min) / (c_max - c_min);
else
    c_normalized = ones(height(T), 1) * 0.5;
end
c_normalized = max(0, min(1, c_normalized));  % Clamp

% Map to colormap indices
n_colors = size(options.Colormap, 1);
color_indices = round(c_normalized * (n_colors - 1)) + 1;
face_colors = options.Colormap(color_indices, :);

% 7. Rendering
if isempty(options.Axes)
    figure;
    ax = gca;
else
    ax = options.Axes;
end

% Plot each patch with its color
hold(ax, 'on');
for i = 1:size(X_verts, 2)
    patch(ax, X_verts(:,i), Y_verts(:,i), face_colors(i,:), ...
        'EdgeColor', options.EdgeColor);
end
hold(ax, 'off');

% Set axis properties
axis(ax, 'tight');
if strcmpi(options.Normalize, 'probability')
    ylim(ax, [0, max(counts)/height(T) * 1.05]);
else
    ylim(ax, [0, max(counts) * 1.05]);
end
box(ax, 'off');

end
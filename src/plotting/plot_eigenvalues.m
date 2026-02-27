function ax = plot_eigenvalues(eigenvalues, ax, time_value, x_lim, y_lim, circle_params)
% plot_eigenvalues - Plot eigenvalue distribution on complex plane
%
% Syntax:
%   ax = plot_eigenvalues(eigenvalues, ax, time_value)
%   ax = plot_eigenvalues(eigenvalues, ax, time_value, x_lim, y_lim)
%   ax = plot_eigenvalues(eigenvalues, ax, time_value, x_lim, y_lim, circle_params)
%
% Description:
%   Plots eigenvalue distribution on the complex plane with polished styling
%   including LaTeX formatted labels and manual reference axes. Based on the
%   style from RMT_figure.m for a professional appearance.
%
% Inputs:
%   eigenvalues - Complex vector of eigenvalues to plot
%   ax          - Axes handle to plot on (e.g., from subplot)
%   time_value  - Time value to display in title (in seconds)
%   x_lim       - (Optional) [xmin xmax] limits for x-axis
%   y_lim       - (Optional) [ymin ymax] limits for y-axis
%   circle_params - (Optional) struct with fields 'center', 'radius', and
%                   optionally 'outlier_threshold' (default 1.04)
%
% Outputs:
%   ax - Axes handle (returned for convenience)
%
% Example:
%   figure;
%   ax = subplot(2, 2, 1);
%   evals = eig(randn(100));
%   ax = plot_eigenvalues(evals, ax, 5.0);
%
% See also: plot_lyapunov, plot_SRNN_tseries

% Handle optional arguments
if nargin < 6
    circle_params = [];
end
if nargin < 5
    y_lim = [];
end
if nargin < 4
    x_lim = [];
end

% Make the specified axes current
axes(ax);

% Classify and plot eigenvalues using three-tier outlier coloring similar to RMT.plot_spectrum
mSize = 4;
hold on;

has_circle = ~isempty(circle_params) && isfield(circle_params, 'center') && isfield(circle_params, 'radius');

if has_circle
    R = circle_params.radius;
    xc = real(circle_params.center);
    yc = imag(circle_params.center);

    % Read outlier_threshold with default for backward compatibility
    if isfield(circle_params, 'outlier_threshold')
        outlier_threshold = circle_params.outlier_threshold;
    else
        outlier_threshold = 1.04;
    end

    % Compute distances from center for all eigenvalues
    distances = abs(eigenvalues - circle_params.center);

    % Interior eigenvalues (within R): black unfilled circles
    interior_mask = distances <= R;
    interior_eigs = eigenvalues(interior_mask);
    plot(real(interior_eigs), imag(interior_eigs), 'ko', 'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);

    % Near outlier eigenvalues (between R and outlier_threshold*R): black Xs
    near_outlier_mask = (distances > R) & (distances <= outlier_threshold * R);
    near_outlier_eigs = eigenvalues(near_outlier_mask);
    if ~isempty(near_outlier_eigs)
        plot(real(near_outlier_eigs), imag(near_outlier_eigs), 'kx', 'MarkerSize', mSize, 'LineWidth', 0.5);
    end

    % Far outlier eigenvalues (beyond outlier_threshold*R): green filled circles
    far_outlier_mask = distances > outlier_threshold * R;
    far_outlier_eigs = eigenvalues(far_outlier_mask);
    if ~isempty(far_outlier_eigs)
        plot(real(far_outlier_eigs), imag(far_outlier_eigs), 'o', 'MarkerSize', mSize, 'MarkerFaceColor', [0 .7 0], 'MarkerEdgeColor', [0 .7 0]);
    end

    % Draw theoretical radius as solid black circle
    theta = linspace(0, 2*pi, 100);
    plot(xc + R*cos(theta), yc + R*sin(theta), 'k-', 'LineWidth', 2);
else
    % No circle params: plot all eigenvalues as black unfilled circles
    plot(real(eigenvalues), imag(eigenvalues), 'ko', 'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
end

% Get axis limits (use auto-scaled if not provided)
if isempty(x_lim)
    x_lim = xlim;
end

% Ensure right xlim always includes 0
if x_lim(2) < 0
    x_lim(2) = 0.05;
end

if isempty(y_lim)
    y_lim = ylim;
end

% Turn axis off for cleaner appearance
axis off;

% Add manual reference lines at Re=0 and Im=0
hold on;
h_x = plot(x_lim, [0, 0], 'k', 'LineWidth', 1.25);
h_y = plot([0, 0], y_lim, 'k', 'LineWidth', 1.25);

% Move reference lines to bottom layer
uistack([h_x, h_y], 'bottom');

% Add LaTeX formatted axis labels
% Position Re label near right edge of x-axis
text(1.02*x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', ...
    'VerticalAlignment', 'middle');

% Position Im label near top of y-axis
text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

% Restore axis limits (in case they changed)
xlim(x_lim);
ylim(y_lim);

hold off;

% Set equal aspect ratio for proper circular representation
axis equal;
end


function [fig_handle, ax_handles] = plot_SRNN_tseries_minimal(t_out, u, x, r, a, b, params, lya_results, Lya_method)
% PLOT_SRNN_TSERIES_MINIMAL Create minimal time series plots for SRNN simulation
%
% Syntax:
%   [fig_handle, ax_handles] = plot_SRNN_tseries_minimal(t_out, u, x, r, a, b, params, lya_results, Lya_method)
%
% Inputs:
%   t_out       - Time vector from ODE solver
%   u           - External input structure with fields u.E and u.I
%   x           - Dendritic states (n x nt)
%   r           - Firing rates (n x nt) - NOT PLOTTED
%   a           - Adaptation variables (n_a x nt) - NOT PLOTTED
%   b           - STD variables (n_b x nt) - NOT PLOTTED
%   params      - Parameters structure containing network configuration
%   lya_results - Results from Lyapunov exponent computation
%   Lya_method  - String indicating Lyapunov method ('benettin', 'qr', or 'none')
%
% Outputs:
%   fig_handle  - Handle to the created figure
%   ax_handles  - Array of axes handles for all subplots
%
% Description:
%   Creates a tiled layout figure with minimal subplots:
%   1. Dendritic states
%   2. Lyapunov exponent(s) (if enabled)
%   All time series plots are linked along the x-axis for synchronized zooming.

% Determine which subplots are needed
has_lyapunov = ~strcmpi(Lya_method, 'none');

% Calculate total number of subplots
n_plots = 1;  % Always: Dendritic states
if has_lyapunov
    n_plots = n_plots + 1;
end

% Create figure with tiled layout
fig_handle = figure('Position', [200         573        1252         326]); % Slightly shorter height since fewer plots
tiledlayout(n_plots, 1);

% Initialize array to store axes handles
ax_handles = [];

% Always create: Dendritic states
ax_handles(end+1) = nexttile;
plot_mean = false;
if isfield(params, 'plot_mean_dendrite')
    plot_mean = params.plot_mean_dendrite;
end
plot_dendritic_state(t_out, x, plot_mean);
set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

% Conditionally create: Lyapunov exponent(s)
if has_lyapunov
    ax_handles(end+1) = nexttile;
    if strcmpi(Lya_method, 'benettin')
        plot_lyapunov(lya_results, Lya_method, {'filtered', 'EOC', 'value'});
    else
        plot_lyapunov(lya_results, Lya_method);
    end
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Link x-axes of all time series plots
linkaxes(ax_handles,'x');

% Add time scale bar overlay in lower right of last subplot
axes(ax_handles(end));  % Make last subplot current
hold on;

% Calculate scale bar length (round 1/10 of total time range)
scale_bar_length = round(0.1 * (t_out(end) - t_out(1)));

% Get current axis limits
xlims = xlim;
ylims = ylim;

% Position scale bar in lower right corner
% X position: end at 95% of x-axis width
x_end = xlims(1) + 0.95 * (xlims(2) - xlims(1));
x_start = x_end - scale_bar_length;
% Y position: at 10% height from bottom
y_pos = ylims(1) + 0.10 * (ylims(2) - ylims(1));

% Draw scale bar
plot([x_start, x_end], [y_pos, y_pos], 'k-', 'LineWidth', 4);

% Add text label below scale bar
text_x = (x_start + x_end) / 2;  % Center of scale bar
text_y = ylims(1) + 0.05 * (ylims(2) - ylims(1));  % Below scale bar
text(text_x, text_y, sprintf('%d seconds', scale_bar_length), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;

end


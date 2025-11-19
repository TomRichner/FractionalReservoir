function [fig_handle, ax_handles] = plot_SRNN_paired_pulse_tseries(t_out, u, x, r, a, b, params, lya_results, Lya_method)
% PLOT_SRNN_PAIRED_PULSE_TSERIES Create time series plots for paired pulse SRNN simulation
%
% Syntax:
%   [fig_handle, ax_handles] = plot_SRNN_paired_pulse_tseries(t_out, u, x, r, a, b, params, lya_results, Lya_method)
%
% Inputs:
%   t_out       - Time vector from ODE solver
%   u           - External input structure with fields u.E and u.I
%   x           - Dendritic states (n x nt)
%   r           - Firing rates (n x nt)
%   a           - Adaptation variables (n_a x nt)
%   b           - STD variables (n_b x nt)
%   params      - Parameters structure containing network configuration
%   lya_results - Results from Lyapunov exponent computation
%   Lya_method  - String indicating Lyapunov method ('benettin', 'qr', or 'none')
%
% Outputs:
%   fig_handle  - Handle to the created figure
%   ax_handles  - Array of axes handles for all subplots
%
% Description:
%   Creates a tiled layout figure optimized for paired pulse experiments.
%   Uses thinner line widths (1.0) for stimulus and gives dendritic states
%   double the vertical space (2 subplot rows).
%   All time series plots are linked along the x-axis for synchronized zooming.

% Determine which subplots are needed
has_adaptation = params.n_a_E > 0 || params.n_a_I > 0;
has_std = params.n_b_E > 0 || params.n_b_I > 0;
has_lyapunov = ~strcmpi(Lya_method, 'none');

% Calculate total number of subplot rows
% External input: 1 row, Dendritic states: 2 rows, Firing rates: 1 row
n_plots = 4;  % Base: stim (1) + dendrite (2) + firing (1) = 4
if has_adaptation
    n_plots = n_plots + 1;
end
if has_std
    n_plots = n_plots + 1;
end
if has_lyapunov
    n_plots = n_plots + 1;
end

% Create figure with tiled layout
fig_handle = figure('Position', [200 300 1200 800]);
tiledlayout(n_plots, 1);

% Initialize array to store axes handles
ax_handles = [];

% Always create: External input (1 row, thin lines for paired pulse)
ax_handles(end+1) = nexttile;
plot_external_input_paired_pulse(t_out, u);
set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

% Always create: Dendritic states (2 rows)
ax_handles(end+1) = nexttile([2 1]);  % Span 2 rows
plot_mean = false;
if isfield(params, 'plot_mean_dendrite')
    plot_mean = params.plot_mean_dendrite;
end
plot_dendritic_state(t_out, x, plot_mean);
set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

% Always create: Firing rates (1 row)
ax_handles(end+1) = nexttile;
plot_firing_rate(t_out, r);
set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

% Conditionally create: Adaptation variables
if has_adaptation
    ax_handles(end+1) = nexttile;
    plot_adaptation(t_out, a, params);
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Conditionally create: STD variables (b)
if has_std
    ax_handles(end+1) = nexttile;
    plot_std_variable(t_out, b, params);
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

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


function plot_external_input_paired_pulse(t, u)
% plot_external_input_paired_pulse - Plot external input with thin lines for paired pulse
%
% Syntax:
%   plot_external_input_paired_pulse(t, u)
%
% Description:
%   Plots external input time series for excitatory and inhibitory neurons
%   with LineWidth of 1.0, optimized for paired pulse visualization.
%
% Inputs:
%   t - Time vector (nt x 1)
%   u - Struct with fields:
%       .E - External input for E neurons (n_E x nt)
%       .I - External input for I neurons (n_I x nt)

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    % Plot inhibitory neurons first (background layer) with LineWidth 1.0
    plot_lines_with_colormap(t, u.I, cmap_I, 'LineWidth', 1.0);
    
    % Plot excitatory neurons on top with LineWidth 1.0
    hold on;
    plot_lines_with_colormap(t, u.E, cmap_E, 'LineWidth', 1.0);
    hold off;
    ylabel('stim');
    
    % Set yticks to match ylim
    yl = ylim;
    yticks(yl);
end



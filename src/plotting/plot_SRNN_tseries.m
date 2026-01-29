function [fig_handle, ax_handles] = plot_SRNN_tseries(t_out, u, x, r, a, b, br, params, lya_results, Lya_method, T_plot)
% PLOT_SRNN_TSERIES Create comprehensive time series plots for SRNN simulation
%
% Syntax:
%   [fig_handle, ax_handles] = plot_SRNN_tseries(t_out, u, x, r, a, b, br, params, lya_results, Lya_method, T_plot)
%
% Inputs:
%   t_out       - Time vector from ODE solver
%   u           - External input structure with fields u.E and u.I
%   x           - Dendritic states (n x nt)
%   r           - Firing rates (n x nt)
%   a           - Adaptation variables (n_a x nt)
%   b           - STD variables (n_b x nt)
%   br          - Synaptic output (n x nt)
%   params      - Parameters structure containing network configuration
%   lya_results - Results from Lyapunov exponent computation
%   Lya_method  - String indicating Lyapunov method ('benettin', 'qr', or 'none')
%   T_plot      - Optional [start, end] time range for plotting
%
% Outputs:
%   fig_handle  - Handle to the created figure
%   ax_handles  - Array of axes handles for all subplots
%
% Description:
%   Creates a tiled layout figure with conditional subplots based on which
%   features are active in the simulation (adaptation, STD, Lyapunov).
%   Always includes: External input, Dendritic states, and Firing rates.
%   All time series plots are linked along the x-axis for synchronized zooming.

if nargin < 11
    T_plot = [];
end

% Determine which subplots are needed
has_adaptation = params.n_a_E > 0 || params.n_a_I > 0;
has_std_vars = (params.n_b_E > 0 || params.n_b_I > 0) && ~isempty(b); % Checks if 'b' (STD variables) should be plotted
has_synaptic_output = ~isempty(br); % Checks if 'br' (synaptic output) should be plotted
has_lyapunov = ~strcmpi(Lya_method, 'none');
has_firing_rate = ~isempty(r);

% Calculate total number of subplots
n_plots = 2;  % Always: External input, Dendritic states
if has_firing_rate
    n_plots = n_plots + 1;
end
if has_synaptic_output
    n_plots = n_plots + 1; % Synaptic output (br)
end
if has_adaptation
    n_plots = n_plots + 1;
end
if has_std_vars
    n_plots = n_plots + 1; % STD variables (b)
end
if has_lyapunov
    n_plots = n_plots + 1;
end

% Create figure with tiled layout
fig_handle = figure('Position', [413   623   880   640]);
tiledlayout(n_plots, 1);

% Initialize array to store axes handles
ax_handles = [];

% Always create: External input
ax_handles(end+1) = nexttile;
plot_external_input(t_out, u);
set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

% Always create: Dendritic states
ax_handles(end+1) = nexttile;
plot_mean = false;
if isfield(params, 'plot_mean_dendrite')
    plot_mean = params.plot_mean_dendrite;
end
plot_dendritic_state(t_out, x, plot_mean);
set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

% Conditionally create: Firing rates
if has_firing_rate
    ax_handles(end+1) = nexttile;
    plot_firing_rate(t_out, r);
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Conditionally create: Synaptic output (br)
if has_synaptic_output
    ax_handles(end+1) = nexttile;
    plot_synaptic_output(t_out, br);
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Conditionally create: Adaptation variables
if has_adaptation
    ax_handles(end+1) = nexttile;
    plot_adaptation(t_out, a, params);
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Conditionally create: STD variables (b)
if has_std_vars
    ax_handles(end+1) = nexttile;
    plot_std_variable(t_out, b, params);
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Conditionally create: Lyapunov exponent(s)
if has_lyapunov
    ax_handles(end+1) = nexttile;
    if strcmpi(Lya_method, 'benettin')
        plot_lyapunov(lya_results, Lya_method, {'local', 'EOC', 'value'});
    else
        plot_lyapunov(lya_results, Lya_method);
    end
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Link x-axes of all time series plots
linkaxes(ax_handles,'x');

% Apply T_plot limit if provided
if ~isempty(T_plot)
    xlim(ax_handles(end), T_plot);
end

% Add time scale bar overlay in lower right of last subplot
axes(ax_handles(end));  % Make last subplot current
hold on;

% Get current axis limits
xlims = xlim;
ylims = ylim;

% Calculate scale bar length (round 1/10 of visible time range)
scale_bar_length = round(0.1 * (xlims(2) - xlims(1)));
if scale_bar_length < 1
    % For very short durations, don't round to integer
    scale_bar_length = 0.1 * (xlims(2) - xlims(1));
end

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
text(text_x, text_y, sprintf('%g seconds', scale_bar_length), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;

end


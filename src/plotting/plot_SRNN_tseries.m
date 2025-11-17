function [fig_handle, ax_handles] = plot_SRNN_tseries(t_out, u, x, r, a, b, params, lya_results, Lya_method)
% PLOT_SRNN_TSERIES Create comprehensive time series plots for SRNN simulation
%
% Syntax:
%   [fig_handle, ax_handles] = plot_SRNN_tseries(t_out, u, x, r, a, b, params, lya_results, Lya_method)
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
%   Creates a tiled layout figure with conditional subplots based on which
%   features are active in the simulation (adaptation, STD, Lyapunov).
%   Always includes: External input, Dendritic states, and Firing rates.
%   All time series plots are linked along the x-axis for synchronized zooming.

% Determine which subplots are needed
has_adaptation = params.n_a_E > 0 || params.n_a_I > 0;
has_std = params.n_b_E > 0 || params.n_b_I > 0;
has_lyapunov = ~strcmpi(Lya_method, 'none');

% Calculate total number of subplots
n_plots = 3;  % Always: External input, Dendritic states, Firing rates
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

% Always create: External input
ax_handles(end+1) = nexttile;
plot_external_input(t_out, u);

% Always create: Dendritic states
ax_handles(end+1) = nexttile;
plot_dendritic_state(t_out, x);

% Always create: Firing rates
ax_handles(end+1) = nexttile;
plot_firing_rate(t_out, r);

% Conditionally create: Adaptation variables
if has_adaptation
    ax_handles(end+1) = nexttile;
    plot_adaptation(t_out, a, params);
end

% Conditionally create: STD variables (b)
if has_std
    ax_handles(end+1) = nexttile;
    plot_std_variable(t_out, b, params);
end

% Conditionally create: Lyapunov exponent(s)
if has_lyapunov
    ax_handles(end+1) = nexttile;
    plot_lyapunov(lya_results, Lya_method);
end

% Link x-axes of all time series plots
linkaxes(ax_handles,'x');

% Add xlabel only to the last subplot
xlabel(ax_handles(end), 'Time (s)');

end


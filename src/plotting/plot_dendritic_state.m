function plot_dendritic_state(t, x, plot_mean)
% plot_dendritic_state - Plot dendritic states for E and I neurons
%
% Syntax:
%   plot_dendritic_state(t, x, plot_mean)
%
% Description:
%   Plots dendritic state (x) time series for excitatory and inhibitory neurons
%   on the current axes. Uses custom colormaps: reds/magentas for inhibitory
%   and blues/greens for excitatory neurons. Inhibitory neurons are plotted
%   first (background), then excitatory on top.
%
% Inputs:
%   t - Time vector (nt x 1)
%   x - Struct with fields:
%       .E - Dendritic states for E neurons (n_E x nt)
%       .I - Dendritic states for I neurons (n_I x nt)
%   plot_mean - (Optional) Boolean, if true plots mean potential in black (default: false)
%
% Example:
%   subplot(5, 1, 2);
%   plot_dendritic_state(t_out, x, true);

    if nargin < 3
        plot_mean = false;
    end

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    % Plot inhibitory neurons first (background layer)
    plot_lines_with_colormap(t, x.I, cmap_I);
    
    % Plot excitatory neurons on top
    hold on;
    plot_lines_with_colormap(t, x.E, cmap_E);
    
    if plot_mean
        % Calculate mean across excitatory neurons only
        mean_x = mean(x.E, 1);
        plot(t, mean_x, 'k', 'LineWidth', 3);
    end
    
    hold off;
    ylabel('dendrite');
end


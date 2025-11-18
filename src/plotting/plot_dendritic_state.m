function plot_dendritic_state(t, x)
% plot_dendritic_state - Plot dendritic states for E and I neurons
%
% Syntax:
%   plot_dendritic_state(t, x)
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
%
% Example:
%   subplot(5, 1, 2);
%   plot_dendritic_state(t_out, x);

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    % Plot inhibitory neurons first (background layer)
    plot_lines_with_colormap(t, x.I, cmap_I);
    
    % Plot excitatory neurons on top
    hold on;
    plot_lines_with_colormap(t, x.E, cmap_E);
    hold off;
    ylabel('dendrite');
end


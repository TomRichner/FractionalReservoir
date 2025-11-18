function plot_firing_rate(t, r)
% plot_firing_rate - Plot firing rates for E and I neurons
%
% Syntax:
%   plot_firing_rate(t, r)
%
% Description:
%   Plots firing rate (r) time series for excitatory and inhibitory neurons
%   on the current axes. Uses custom colormaps: reds/magentas for inhibitory
%   and blues/greens for excitatory neurons. Inhibitory neurons are plotted
%   first (background), then excitatory on top.
%
% Inputs:
%   t - Time vector (nt x 1)
%   r - Struct with fields:
%       .E - Firing rates for E neurons (n_E x nt)
%       .I - Firing rates for I neurons (n_I x nt)
%
% Example:
%   subplot(5, 1, 4);
%   plot_firing_rate(t_out, r);

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    % Plot inhibitory neurons first (background layer)
    plot_lines_with_colormap(t, r.I, cmap_I);
    
    % Plot excitatory neurons on top
    hold on;
    plot_lines_with_colormap(t, r.E, cmap_E);
    hold off;
    ylabel('firing rate');
    yticks([0, 1]);
    ylim([0, 1]);
end


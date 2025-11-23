function plot_synaptic_output(t, br)
% plot_synaptic_output - Plot synaptic output (br = b .* r) for E and I neurons
%
% Syntax:
%   plot_synaptic_output(t, br)
%
% Description:
%   Plots synaptic output time series for excitatory and inhibitory neurons
%   on the current axes. Uses custom colormaps: reds/magentas for inhibitory
%   and blues/greens for excitatory neurons. Inhibitory neurons are plotted
%   first (background), then excitatory on top.
%
% Inputs:
%   t  - Time vector (nt x 1)
%   br - Struct with fields:
%        .E - Synaptic output for E neurons (n_E x nt)
%        .I - Synaptic output for I neurons (n_I x nt)
%
% Example:
%   subplot(5, 1, 5);
%   plot_synaptic_output(t_out, br);

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    % Plot inhibitory neurons first (background layer)
    plot_lines_with_colormap(t, br.I, cmap_I);
    
    % Plot excitatory neurons on top
    hold on;
    plot_lines_with_colormap(t, br.E, cmap_E);
    hold off;
    ylabel('synaptic output');
    yticks([0, 1]);
    ylim([0, 1]);
end


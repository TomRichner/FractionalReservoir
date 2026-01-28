function plot_external_input(t, u)
% plot_external_input - Plot external input for E and I neurons
%
% Syntax:
%   plot_external_input(t, u)
%
% Description:
%   Plots external input time series for excitatory and inhibitory neurons
%   on the current axes. Uses custom colormaps: reds/magentas for inhibitory
%   and blues/greens for excitatory neurons. Inhibitory neurons are plotted
%   first (background), then excitatory on top.
%
% Inputs:
%   t - Time vector (nt x 1)
%   u - Struct with fields:
%       .E - External input for E neurons (n_E x nt)
%       .I - External input for I neurons (n_I x nt)
%
% Example:
%   subplot(5, 1, 1);
%   plot_external_input(t_out, u);

% Get colormaps (8 colors each)
cmap_I = inhibitory_colormap(8);
cmap_E = excitatory_colormap(8);

% Plot inhibitory neurons first (background layer)
plot_lines_with_colormap(t, u.I, cmap_I);

% Plot excitatory neurons on top
hold on;
plot_lines_with_colormap(t, u.E, cmap_E);
hold off;
ylabel('stim');

% Set yticks to match ylim
yl = ylim;
yticks([-0.5 0 0.5]);
end


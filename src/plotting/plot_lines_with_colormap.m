function plot_lines_with_colormap(t, data, cmap, varargin)
% plot_lines_with_colormap - Plot multiple lines with explicit colormap assignment
%
% Syntax:
%   plot_lines_with_colormap(t, data, cmap)
%   plot_lines_with_colormap(t, data, cmap, 'PropertyName', PropertyValue, ...)
%
% Description:
%   Plots each row of data as a separate line with colors explicitly assigned
%   from the provided colormap. This ensures colors are applied correctly even
%   when using hold on with multiple colormaps. Colors cycle if there are more
%   data rows than colors in the colormap.
%
% Inputs:
%   t       - Time vector (nt x 1)
%   data    - Data matrix (n_lines x nt), each row is plotted as a line
%   cmap    - Colormap matrix (n_colors x 3) with RGB values in [0, 1]
%   varargin - Optional name-value pairs for plot properties (e.g., 'LineWidth', 1.5)
%
% Example:
%   % Plot inhibitory neurons with red colormap
%   cmap_I = inhibitory_colormap(8);
%   plot_lines_with_colormap(t, u.I, cmap_I, 'LineWidth', 1);
%   hold on;
%   
%   % Plot excitatory neurons with blue colormap
%   cmap_E = excitatory_colormap(8);
%   plot_lines_with_colormap(t, u.E, cmap_E, 'LineWidth', 1);
%   hold off;

    % Check if data is empty
    if isempty(data)
        return;
    end
    
    % Get number of lines to plot
    n_lines = size(data, 1);
    n_colors = size(cmap, 1);
    
    % Plot each line with explicit color
    hold on;
    for i = 1:n_lines
        % Cycle through colors if needed
        color_idx = mod(i - 1, n_colors) + 1;
        plot(t, data(i, :), 'Color', cmap(color_idx, :), varargin{:});
    end
end


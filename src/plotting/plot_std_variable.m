function plot_std_variable(t, b, params)
% plot_std_variable - Plot short-term depression variables for E and I neurons
%
% Syntax:
%   plot_std_variable(t, b, params)
%
% Description:
%   Plots short-term depression variable (b) time series for excitatory and
%   inhibitory neurons on the current axes. Uses custom colormaps: reds/magentas
%   for inhibitory and blues/greens for excitatory neurons. Inhibitory neurons
%   are plotted first (background), then excitatory on top.
%
% Inputs:
%   t      - Time vector (nt x 1)
%   b      - Struct with fields:
%            .E - STD variable for E neurons (n_E x nt)
%            .I - STD variable for I neurons (n_I x nt)
%   params - Struct with fields:
%            .n_b_E - Number of STD timescales for E neurons (0 or 1)
%            .n_b_I - Number of STD timescales for I neurons (0 or 1)
%
% Example:
%   subplot(5, 1, 5);
%   plot_std_variable(t_out, b, params);

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    has_std = false;
    
    % Plot inhibitory STD first (background layer)
    if params.n_b_I > 0 && ~isempty(b.I) && size(b.I, 1) > 0
        % Check if b.I is not all ones (actual STD dynamics present)
        if ~all(b.I(:) == 1)
            plot_lines_with_colormap(t, b.I, cmap_I);
            has_std = true;
        end
    end
    
    % Plot excitatory STD on top
    if params.n_b_E > 0 && ~isempty(b.E) && size(b.E, 1) > 0
        % Check if b.E is not all ones (actual STD dynamics present)
        if ~all(b.E(:) == 1)
            if has_std
                hold on;
            end
            plot_lines_with_colormap(t, b.E, cmap_E);
            has_std = true;
        end
    end
    
    if has_std
        hold off;
        ylabel('depression');
        ylim([0, 1]);
        yticks([0, 1]);
    else
        % No STD variables to plot
        text(0.5, 0.5, 'No STD variables', 'HorizontalAlignment', 'center');
        axis off;
    end
end


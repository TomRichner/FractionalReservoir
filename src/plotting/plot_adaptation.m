function plot_adaptation(t, a, params)
% plot_adaptation - Plot adaptation variables for E and I neurons
%
% Syntax:
%   plot_adaptation(t, a, params)
%
% Description:
%   Plots adaptation variable (a) time series for excitatory and inhibitory
%   neurons on the current axes. Handles multiple adaptation timescales per
%   neuron. Uses custom colormaps: reds/magentas for inhibitory and
%   blues/greens for excitatory neurons. Inhibitory neurons are plotted
%   first (background), then excitatory on top.
%
% Inputs:
%   t      - Time vector (nt x 1)
%   a      - Struct with fields:
%            .E - Adaptation for E neurons (n_E x n_a_E x nt) or empty
%            .I - Adaptation for I neurons (n_I x n_a_I x nt) or empty
%   params - Struct with fields:
%            .n_E   - Number of excitatory neurons
%            .n_I   - Number of inhibitory neurons
%            .n_a_E - Number of adaptation timescales for E neurons
%            .n_a_I - Number of adaptation timescales for I neurons
%
% Example:
%   subplot(5, 1, 3);
%   plot_adaptation(t_out, a, params);

    % Get colormaps (8 colors each)
    cmap_I = inhibitory_colormap(8);
    cmap_E = excitatory_colormap(8);
    
    has_adaptation = false;
    
    % Plot inhibitory adaptation first (background layer)
    if ~isempty(a.I) && params.n_a_I > 0
        % Reshape from (n_I x n_a_I x nt) to ((n_I * n_a_I) x nt)
        a_I_reshaped = reshape(a.I, params.n_I * params.n_a_I, []);
        plot_lines_with_colormap(t, a_I_reshaped, cmap_I);
        has_adaptation = true;
    end
    
    % Plot excitatory adaptation on top
    if ~isempty(a.E) && params.n_a_E > 0
        if has_adaptation
            hold on;
        end
        % Reshape from (n_E x n_a_E x nt) to ((n_E * n_a_E) x nt)
        a_E_reshaped = reshape(a.E, params.n_E * params.n_a_E, []);
        plot_lines_with_colormap(t, a_E_reshaped, cmap_E);
        has_adaptation = true;
    end
    
    if has_adaptation
        hold off;
        ylabel('Adaptation a');
    else
        % No adaptation variables to plot
        text(0.5, 0.5, 'No adaptation variables', 'HorizontalAlignment', 'center');
        axis off;
    end
end


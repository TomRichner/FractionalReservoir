function [fig_handle, ax_handles] = plot_SRNN_combined_tseries(runs, gap_duration, subplots_to_plot)
% PLOT_SRNN_COMBINED_TSERIES Combined time series plots for multiple SRNN runs
%
% Syntax:
%   [fig_handle, ax_handles] = plot_SRNN_combined_tseries(runs, gap_duration, subplots_to_plot)
%
% Inputs:
%   runs             - Cell array of structs, each containing:
%                      .plot_data (struct with t, u, x, r, a, b, br)
%                      .params    (struct)
%                      .lya_results (struct)
%                      .Lya_method (string)
%   gap_duration     - Time in seconds to insert between runs (default: 10)
%   subplots_to_plot - (Optional) Cell array of strings specifying which subplots to include
%                      and their order. Valid options:
%                      'u_ex', 'x', 'r', 'br', 'a', 'b', 'lya'
%
% Outputs:
%   fig_handle       - Handle to the created figure
%   ax_handles       - Array of axes handles

if nargin < 2
    gap_duration = 10;
end

if nargin < 3
    subplots_to_plot = {};
end

n_runs = length(runs);
if n_runs == 0
    error('No runs provided');
end

% Determine which subplots to show
if isempty(subplots_to_plot)
    % Auto-detect based on data if not specified
    subplots_to_plot = {'u_ex', 'x'}; % Always include input and dendritic states

    % Check features across all runs
    has_adaptation = false;
    has_std_vars = false;
    has_synaptic_output = false;
    has_lyapunov = false;
    has_firing_rate = false;

    for i = 1:n_runs
        p = runs{i}.params;
        d = runs{i}.plot_data;

        if p.n_a_E > 0 || p.n_a_I > 0, has_adaptation = true; end
        if (p.n_b_E > 0 || p.n_b_I > 0), has_std_vars = true; end
        if ~isempty(d.br), has_synaptic_output = true; end
        if ~strcmpi(runs{i}.Lya_method, 'none'), has_lyapunov = true; end
        if ~isempty(d.r), has_firing_rate = true; end
    end

    if has_firing_rate, subplots_to_plot{end+1} = 'r'; end
    if has_synaptic_output, subplots_to_plot{end+1} = 'br'; end
    if has_adaptation, subplots_to_plot{end+1} = 'a'; end
    if has_std_vars, subplots_to_plot{end+1} = 'b'; end
    if has_lyapunov, subplots_to_plot{end+1} = 'lya'; end
end

n_plots = length(subplots_to_plot);

% Create figure
fig_handle = figure('Position', [200         371        1300         726]);
tiledlayout(n_plots, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
ax_handles = gobjects(n_plots, 1);

% Helper to get shifted time
get_t_shifted = @(t, offset) t - t(1) + offset;

% Iterate through requested subplots
for curr_ax_idx = 1:n_plots
    subplot_type = subplots_to_plot{curr_ax_idx};
    ax_handles(curr_ax_idx) = nexttile;
    hold on;
    offset = 0;

    for i = 1:n_runs
        r = runs{i};
        t_shifted = get_t_shifted(r.plot_data.t, offset);
        p = r.params;

        switch subplot_type
            case 'u_ex'
                plot_external_input(t_shifted, r.plot_data.u);

            case 'x'
                plot_mean = false;
                if isfield(r.params, 'plot_mean_dendrite')
                    plot_mean = r.params.plot_mean_dendrite;
                end
                plot_dendritic_state(t_shifted, r.plot_data.x, plot_mean);

            case 'r'
                if ~isempty(r.plot_data.r)
                    plot_firing_rate(t_shifted, r.plot_data.r);
                end

            case 'br'
                if ~isempty(r.plot_data.br)
                    plot_synaptic_output(t_shifted, r.plot_data.br);
                end

            case 'a'
                a = r.plot_data.a;
                has_run_adaptation = (p.n_a_E > 0 || p.n_a_I > 0);

                if has_run_adaptation
                    cmap_I = inhibitory_colormap(8);
                    cmap_E = excitatory_colormap(8);

                    % I neurons
                    if ~isempty(a.I) && p.n_a_I > 0
                        a_I_sum = sum(a.I, 2);
                        a_I_summed = reshape(a_I_sum, p.n_I, []);
                        plot_lines_with_colormap(t_shifted, a_I_summed, cmap_I);
                    elseif p.n_a_I == 0
                        % If partial adaptation (e.g. only E), I is skipped or zeros?
                        % Assuming if n_a_I == 0 we don't plot anything for I unless we want zeros.
                        % Let's stick to previous logic: plot only what exists.
                    end

                    % E neurons
                    if ~isempty(a.E) && p.n_a_E > 0
                        a_E_sum = sum(a.E, 2);
                        a_E_summed = reshape(a_E_sum, p.n_E, []);
                        plot_lines_with_colormap(t_shifted, a_E_summed, cmap_E);
                    end
                else
                    % Plot zeros - Excitatory Colormap by default for visual consistency
                    cmap_E = excitatory_colormap(8);
                    cmap_I = inhibitory_colormap(8);

                    if p.n_I > 0
                        % Plot I first (background)
                        zeros_I = zeros(p.n_I, length(t_shifted));
                        plot_lines_with_colormap(t_shifted, zeros_I, cmap_I);
                    end

                    if p.n_E > 0
                        % Plot E on top
                        zeros_E = zeros(p.n_E, length(t_shifted));
                        plot_lines_with_colormap(t_shifted, zeros_E, cmap_E);
                    end
                end

            case 'b'
                b = r.plot_data.b;
                has_run_std = (p.n_b_E > 0 || p.n_b_I > 0);

                if has_run_std
                    cmap_I = inhibitory_colormap(8);
                    cmap_E = excitatory_colormap(8);

                    if p.n_b_I > 0 && ~isempty(b.I)
                        plot_lines_with_colormap(t_shifted, b.I, cmap_I);
                    elseif p.n_b_I == 0 && p.n_I > 0
                        ones_I = ones(p.n_I, length(t_shifted));
                        plot_lines_with_colormap(t_shifted, ones_I, cmap_I);
                    end

                    if p.n_b_E > 0 && ~isempty(b.E)
                        plot_lines_with_colormap(t_shifted, b.E, cmap_E);
                    elseif p.n_b_E == 0 && p.n_E > 0
                        ones_E = ones(p.n_E, length(t_shifted));
                        plot_lines_with_colormap(t_shifted, ones_E, cmap_E);
                    end
                else
                    % Plot ones
                    cmap_E = excitatory_colormap(8);
                    cmap_I = inhibitory_colormap(8);

                    if p.n_I > 0
                        % Plot I first
                        ones_I = ones(p.n_I, length(t_shifted));
                        plot_lines_with_colormap(t_shifted, ones_I, cmap_I);
                    end

                    if p.n_E > 0
                        % Plot E on top
                        ones_E = ones(p.n_E, length(t_shifted));
                        plot_lines_with_colormap(t_shifted, ones_E, cmap_E);
                    end
                end
                ylabel('depression');
                ylim([0, 1.05]);
                yticks([0, 1]);

            case 'lya'
                if isfield(r.lya_results, 't_lya')
                    lya_res_shifted = r.lya_results;
                    lya_res_shifted.t_lya = get_t_shifted(lya_res_shifted.t_lya, offset);
                    plot_lyapunov(lya_res_shifted, r.Lya_method, {'local', 'EOC'});
                end

            otherwise
                warning('Unknown subplot type: %s', subplot_type);
        end

        hold on; % Ensure subsequent runs plot on top
        if i < n_runs
            offset = t_shifted(end) + gap_duration;
        end
    end

    if strcmp(subplot_type, 'a')
        ylabel('adaptation');
        yl = ylim;
        ylim([-0.05, yl(2)]);
    end

    hold off;
    set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
end

% Link axes
linkaxes(ax_handles, 'x');

% Add scale bar to last subplot
if ~isempty(ax_handles)
    axes(ax_handles(end));
    hold on;
    xlims = xlim;
    ylims = ylim;
    scale_bar_length = round(0.1 * (xlims(2) - xlims(1)) / 5) * 5;
    if scale_bar_length < 5, scale_bar_length = 5; end

    x_end = xlims(1) + 0.85 * (xlims(2) - xlims(1));
    x_start = x_end - scale_bar_length;
    y_pos = ylims(1) + 0.10 * (ylims(2) - ylims(1));

    plot([x_start, x_end], [y_pos, y_pos], 'k-', 'LineWidth', 4);
    text_x = (x_start + x_end) / 2;
    text_y = ylims(1) + 0.05 * (ylims(2) - ylims(1));
    text(text_x, text_y, sprintf('%g seconds', scale_bar_length), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    hold off;
end

end

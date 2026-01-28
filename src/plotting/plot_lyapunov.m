function plot_lyapunov(lya_results, Lya_method, plot_options)
% plot_lyapunov - Plot Lyapunov exponent(s) on current axes
%
% Syntax:
%   plot_lyapunov(lya_results, Lya_method)
%   plot_lyapunov(lya_results, Lya_method, plot_options)
%
% Description:
%   Plots Lyapunov exponent analysis results on the current axes.
%   Handles two methods:
%     - 'benettin': Plots single largest exponent with filtered version
%     - 'qr': Plots full Lyapunov spectrum with legend
%
% Inputs:
%   lya_results  - Struct containing Lyapunov analysis results:
%                  For 'benettin': LLE, local_lya, t_lya, lya_fs
%                  For 'qr': LE_spectrum, local_LE_spectrum_t, t_lya, params.N_sys_eqs
%   Lya_method   - String: 'benettin' or 'qr'
%   plot_options - (Optional) Cell array of strings specifying what to plot for 'benettin':
%                  'local': Raw local Lyapunov exponent time series
%                  'filtered': Filtered version of local exponent
%                  'asym': Asymptotic line at final LLE value
%                  'EOC': Edge of chaos line at zero
%                  'value': Text showing LLE value
%                  Default: {'local', 'filtered', 'asym', 'EOC', 'value'}
%
% Example:
%   subplot(6, 1, 6);
%   plot_lyapunov(lya_results, 'benettin');
%   plot_lyapunov(lya_results, 'benettin', {'filtered', 'EOC', 'value'});

% Parse input arguments
if nargin < 3
    plot_options = {'local', 'filtered', 'asym', 'EOC', 'value'};
end

% Validate plot_options for benettin method
if strcmpi(Lya_method, 'benettin')
    valid_options = {'local', 'filtered', 'asym', 'EOC', 'value'};
    if ~iscell(plot_options)
        error('plot_options must be a cell array of strings');
    end
    for i = 1:length(plot_options)
        if ~any(strcmpi(plot_options{i}, valid_options))
            error('Invalid plot_option: %s. Valid options are: %s', ...
                plot_options{i}, strjoin(valid_options, ', '));
        end
    end
end

if strcmpi(Lya_method, 'benettin')
    % Plot for Benettin's method (single exponent)

    % Check which components to plot
    plot_local = any(strcmpi('local', plot_options));
    plot_filtered = any(strcmpi('filtered', plot_options));
    plot_asym = any(strcmpi('asym', plot_options));
    plot_EOC = any(strcmpi('EOC', plot_options));
    plot_value = any(strcmpi('value', plot_options));

    % Design filter if needed
    if plot_filtered
        filter_order = 2;
        filter_cutoff = 0.25;  % Normalized cutoff frequency
        [bL, aL] = butter(filter_order, filter_cutoff/(lya_results.lya_fs/2), 'low');
    end

    % Track legend entries and handles
    legend_entries = {};
    plot_started = false;

    % Plot raw local Lyapunov exponent
    if plot_local
        plot(lya_results.t_lya, lya_results.local_lya)
        hold on
        plot_started = true;
        legend_entries{end+1} = 'Local LLE';
    end

    % Plot filtered signal
    if plot_filtered
        if ~plot_started
            hold on
            plot_started = true;
        end
        local_lya_filt = filtfilt(bL, aL, lya_results.local_lya);
        colors = lines(1);
        plot(lya_results.t_lya, local_lya_filt, 'Color', colors(1,:))
        legend_entries{end+1} = 'Local LLE';
    end

    % Add asymptotic line at LLE value
    if plot_asym
        if ~plot_started
            hold on
            plot_started = true;
        end
        % yline(lya_results.LLE, '--r', 'LineWidth', 1.5)
        plot([lya_results.t_lya(1), lya_results.t_lya(end)], ...
            [lya_results.LLE, lya_results.LLE], '--r', 'LineWidth', 1.5);
        % legend_entries{end+1} = 'LLE';
    end

    % Add reference line at zero (EOC)
    if plot_EOC
        if ~plot_started
            hold on
            plot_started = true;
        end
        % yline(0, '--k')
        plot([lya_results.t_lya(1), lya_results.t_lya(end)], [0, 0], '--k');
        % legend_entries{end+1} = 'EOC';
    end

    ylabel('\lambda_1')

    % % Add legend if there are entries
    % if ~isempty(legend_entries)
    %     legend(legend_entries, 'Location', 'best')
    % end

    % Add LLE text just above y=0 line, on the right side
    if plot_value
        ylims = ylim;
        xlims = xlim;
        text_y = 0.05 * (ylims(2) - ylims(1));  % Slightly above zero
        text_x = xlims(2);  % Right edge of x-axis
        text(text_x, text_y, ['$\lambda_1 = ' sprintf('%.2f', lya_results.LLE) '$'], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'Interpreter', 'latex');
    end

    hold off
    box off

elseif strcmpi(Lya_method, 'qr')
    % Plot for QR method (full spectrum)

    % Plot spectrum (reverse order for plotting, then reorder handles)
    plot_data = lya_results.local_LE_spectrum_t(:, end:-1:1);
    line_handles = plot(lya_results.t_lya, plot_data);
    line_handles = line_handles(end:-1:1); % reorder handles to match sorted spectrum order

    hold on
    yline(0, '--k')
    ylabel('\lambda_1')

    % Add legend with final values (most positive exponents on top)
    legend_count = min(5, lya_results.params.N_sys_eqs);
    legend_entries = cell(1, legend_count);
    for i = 1:legend_count
        legend_entries{i} = sprintf('\\lambda_{%d} = %.3f', i, lya_results.LE_spectrum(i));
    end
    legend(line_handles(1:legend_count), legend_entries, 'Location', 'best')
    hold off
    box off
end
end


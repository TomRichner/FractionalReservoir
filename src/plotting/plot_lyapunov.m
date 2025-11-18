function plot_lyapunov(lya_results, Lya_method)
% plot_lyapunov - Plot Lyapunov exponent(s) on current axes
%
% Syntax:
%   plot_lyapunov(lya_results, Lya_method)
%
% Description:
%   Plots Lyapunov exponent analysis results on the current axes.
%   Handles two methods:
%     - 'benettin': Plots single largest exponent with filtered version
%     - 'qr': Plots full Lyapunov spectrum with legend
%
% Inputs:
%   lya_results - Struct containing Lyapunov analysis results:
%                 For 'benettin': LLE, local_lya, t_lya, lya_fs
%                 For 'qr': LE_spectrum, local_LE_spectrum_t, t_lya, params.N_sys_eqs
%   Lya_method  - String: 'benettin' or 'qr'
%
% Example:
%   subplot(6, 1, 6);
%   plot_lyapunov(lya_results, 'benettin');

    if strcmpi(Lya_method, 'benettin')
        % Plot for Benettin's method (single exponent)
        
        % Design low-pass Butterworth filter
        filter_order = 3;
        filter_cutoff = 1;  % Normalized cutoff frequency
        [bL, aL] = butter(filter_order, filter_cutoff/(lya_results.lya_fs/2), 'low');
        
        % Only filter data after specified time
        t_filter_start = 0.5; % seconds
        idx_filter_start = find(lya_results.t_lya >= t_filter_start, 1, 'first');
        
        % Plot raw local Lyapunov exponent
        plot(lya_results.t_lya, lya_results.local_lya)
        hold on
        
        % Plot filtered signal if we have enough data
        if ~isempty(idx_filter_start)
            t_lya_filt = lya_results.t_lya(idx_filter_start:end);
            local_lya_filt = filtfilt(bL, aL, lya_results.local_lya(idx_filter_start:end));
            plot(t_lya_filt, local_lya_filt)
        end
        
        % Add reference line at zero
        yline(0, '--k')
        
        ylabel('Lyapunov Exponent')
        title(sprintf('LLE = %.4f', lya_results.LLE))
        legend('Local', 'Filtered', 'EOC', 'Location', 'best')
        hold off
        
    elseif strcmpi(Lya_method, 'qr')
        % Plot for QR method (full spectrum)
        
        % Plot spectrum (reverse order for plotting, then reorder handles)
        plot_data = lya_results.local_LE_spectrum_t(:, end:-1:1);
        line_handles = plot(lya_results.t_lya, plot_data);
        line_handles = line_handles(end:-1:1); % reorder handles to match sorted spectrum order
        
        hold on
        yline(0, '--k')
        ylabel('Lyapunov Exponents')
        
        % Add legend with final values (most positive exponents on top)
        legend_count = min(5, lya_results.params.N_sys_eqs);
        legend_entries = cell(1, legend_count);
        for i = 1:legend_count
            legend_entries{i} = sprintf('\\lambda_{%d} = %.3f', i, lya_results.LE_spectrum(i));
        end
        legend(line_handles(1:legend_count), legend_entries, 'Location', 'best')
        hold off
    end
end


function [MI_corrected, MI_uncorrected] = mutual_info_SISO(data_in, data_out, n_bins_in, n_bins_out)
% computes mutual information with and without Miller-Madow bias correction. 
% This correction helps to prevent overestimation of MI when the number of
% data samples is small relative to the number of bins.
%
% OUTPUTS:
%   MI_corrected   - Mutual information with Miller-Madow bias correction
%   MI_uncorrected - Raw mutual information without bias correction
%
% data_in and data_out shoudl be channels x time

    [n_ch_in, n_t]   = size(data_in);
    [n_ch_out, n_t_out] = size(data_out);
    
    assert(n_t == n_t_out, 'data_in and data_out must have same number of time samples');
    assert(n_ch_in == 1, 'data_in must have a single channel for SISO');
    assert(n_ch_out == 1, 'data_out must have a single channel for SISO');

    min_in = min(data_in,[],'all');
    max_in = max(data_in,[],'all');

    min_out = min(data_out,[],'all');
    max_out = max(data_out,[],'all');

    edges_in = linspace(min_in, max_in, n_bins_in + 1);
    edges_out = linspace(min_out, max_out, n_bins_out + 1);

    P_in = histcounts(data_in, edges_in);
    P_in = P_in ./ sum(P_in, 'all');

    P_out = histcounts(data_out, edges_out);
    P_out = P_out ./ sum(P_out, 'all');

    P_joint = histcounts2(data_in, data_out, edges_in, edges_out);
    P_joint = P_joint ./ sum(P_joint, 'all');

    MI_plugin = sum(P_joint .* log2(P_joint ./ (P_in' * P_out)), 'all', 'omitnan');

    % Miller-Madow bias correction
    % This corrects for the positive bias of the plug-in MI estimator
    bias_correction = ((n_bins_in - 1) * (n_bins_out - 1)) / (2 * n_t * log(2));
    MI_bias_corrected = MI_plugin - bias_correction;
    
    % MI cannot be negative, so floor at 0.
    MI_uncorrected = max(0, MI_plugin);
    MI_corrected = max(0, MI_bias_corrected);

end
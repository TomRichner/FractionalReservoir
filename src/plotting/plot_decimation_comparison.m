function fig = plot_decimation_comparison(t_orig, x_orig, t_deci, x_deci, params, title_str)
% PLOT_DECIMATION_COMPARISON Overlay original and downsampled time series
%
% Creates a figure comparing original data (thick lines) with downsampled
% data (thin red lines) overlaid on the same axes.
%
% Syntax:
%   fig = plot_decimation_comparison(t_orig, x_orig, t_deci, x_deci, params, title_str)
%
% Inputs:
%   t_orig    - Original time vector (n_samples_orig x 1)
%   x_orig    - Original data [n_samples_orig x n_channels]
%   t_deci    - Downsampled time vector (n_samples_deci x 1)
%   x_deci    - Downsampled data [n_samples_deci x n_channels]
%   params    - Struct with optional plotting parameters:
%               .time_window  - [t_start, t_end] in seconds (default: first 10s)
%               .n_E          - Number of excitatory neurons (default: half)
%               .linewidth_orig - Line width for original data (default: 1.5)
%               .linewidth_deci - Line width for decimated data (default: 0.5)
%               .detrend_mode   - 'none', 'dc', or 'linear' (default: 'linear')
%   title_str - Figure title string
%
% Outputs:
%   fig       - Handle to created figure
%
% Example:
%   params.time_window = [0, 5];
%   params.n_E = 50;
%   fig = plot_decimation_comparison(t, x, t_ds, x_ds, params, 'SRNN_none');
%
% See also: plot_dendritic_state, SRNN_load_and_preprocess

% Parse parameters
if ~isfield(params, 'time_window')
    params.time_window = [t_orig(1), min(t_orig(1) + 10, t_orig(end))];
end
if ~isfield(params, 'n_E')
    n_channels = size(x_orig, 2);
    params.n_E = round(n_channels / 2);
end
if ~isfield(params, 'linewidth_orig')
    params.linewidth_orig = 1.5;
end
if ~isfield(params, 'linewidth_deci')
    params.linewidth_deci = 0.5;
end
if ~isfield(params, 'detrend_mode')
    params.detrend_mode = 'linear';
end

n_E = params.n_E;
n_channels = size(x_orig, 2);
n_I = n_channels - n_E;

% Apply detrending to FULL original data (matching what was done to decimated during processing)
switch params.detrend_mode
    case 'none'
        x_orig_detrended = x_orig;
    case 'dc'
        x_orig_detrended = x_orig - mean(x_orig, 1);
    case 'linear'
        x_orig_detrended = detrend(x_orig);
    otherwise
        x_orig_detrended = x_orig;
end

% Get time window indices
t_win = params.time_window;
orig_mask = t_orig >= t_win(1) & t_orig <= t_win(2);
deci_mask = t_deci >= t_win(1) & t_deci <= t_win(2);

% Window the data AFTER detrending original on full chunk
t_orig_win = t_orig(orig_mask);
t_deci_win = t_deci(deci_mask);
x_orig_all = x_orig_detrended(orig_mask, :);
x_deci_all = x_deci(deci_mask, :);  % Use as-is (already detrended during processing)

% Create figure
fig = figure('Position', [100 100 1400 500]);

% Single panel with all neurons
hold on;
% Plot original data with normal lines (MATLAB default colors cycle)
h_orig = plot(t_orig_win, x_orig_all, 'LineWidth', params.linewidth_orig);
% Plot downsampled data in black
h_deci = plot(t_deci_win, x_deci_all, 'k', 'LineWidth', params.linewidth_deci);
hold off;

ylabel('Dendrite (x)');
xlabel('Time (s)');
title(sprintf('Decimation Comparison: %s (n=%d neurons)', title_str, n_channels), 'Interpreter', 'none');
legend([h_orig(1), h_deci(1)], {'Original (lw 1.5)', 'Downsampled (lw 0.5, black)'}, 'Location', 'northeast');
xlim(t_win);

end

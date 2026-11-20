function [T, M] = outlier_table(data_dir, opts)
%OUTLIER_TABLE The outlier measures of every run in a folder: one row per run, printed, and saved as CSV.
%
%   T = outlier_table()                    data/explorations/outliers
%   [T, M] = outlier_table(data_dir, 'print', true)
%
% T is a table of the scalar measures (outlier_measures, field `row`), M a
% cell of the full measure structs in the same order. The CSV is
% outlier_measures.csv in data_dir. The printed block is the short form;
% rates are in 1/s, dwell times in ms.
%
% See also: outlier_measures, outlier_run, fig_outliers

arguments
    data_dir (1,:) char = ''
    opts.print (1,1) logical = true
end
if isempty(data_dir)
    data_dir = fullfile(fileparts(which('setup_paths')), 'data', 'explorations', 'outliers');
end
files = dir(fullfile(data_dir, 'run_*.mat'));
if isempty(files)
    error('outlier_table:NoRuns', 'No run_*.mat in %s.', data_dir);
end

rows = cell(1, numel(files)); M = cell(1, numel(files));
for k = 1:numel(files)
    L = load(fullfile(files(k).folder, files(k).name), 'R');
    M{k} = outlier_measures(L.R);
    rows{k} = M{k}.row;
end
T = struct2table([rows{:}], 'AsArray', true);
if ischar(T.network);   T.network = {T.network};     end   % a single run gives char columns
if ischar(T.condition); T.condition = {T.condition}; end

% Order: network, seed, then the conditions in the paper's order.
cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
[~, ci] = ismember(T.condition, cond_order);
[~, order] = sortrows(table(T.network, T.seed, ci));
T = T(order, :);
M = M(order);
writetable(T, fullfile(data_dir, 'outlier_measures.csv'));

if opts.print
    fprintf('\nOutlier measures (rates 1/s; dwell ms; M in units of W)\n');
    fprintf('%-10s %-14s %4s | %7s | %6s %6s | %7s %7s %7s | %7s %6s | %5s %5s %5s | %7s | %5s | %6s %6s %6s | %5s %5s | %7s %7s | %7s %5s\n', ...
        'network', 'condition', 'seed', 'lam_1', 'kappa', 'R_eff', 'edge', 'top', 'top95', ...
        'excess', 'ratio', 'f_det', 'f_uns', 'f_both', 'X', 'n_out', 'dw_blk', 'dw_0', 'dw_x', ...
        'same', 'N_eff', 'K1', 'K_err', 'full_ab', 'f_fu');
    for k = 1:height(T)
        r = T(k, :);
        fprintf('%-10s %-14s %4d | %+7.3f | %6.3f %6.3f | %+7.2f %+7.2f %+7.2f | %+7.2f %6.2f | %5.2f %5.2f %5.2f | %7.3f | %5.2f | %6.0f %6.0f %6.0f | %5.2f %5.1f | %+7.2f %+7.2f | %+7.2f %5.2f\n', ...
            r.network{1}, r.condition{1}, r.seed, r.lambda_1, r.kappa, r.R_eff_med, ...
            r.rate_edge_med, r.rate_top_med, r.rate_top_p95, r.excess_med, r.ratio_med, ...
            r.frac_detached, r.frac_unstable, r.frac_both, r.X, r.n_out_mean, ...
            1e3 * r.dwell_bulk_med, 1e3 * r.dwell_zero_med, 1e3 * r.dwell_x_med, ...
            r.frac_same, r.N_eff, r.K1_med, r.K_err_med, r.full_abscissa_med, r.full_frac_unstable);
    end
    fprintf(['\nedge, top, top95: median rate of the bulk edge, median and 95th percentile of the rightmost rate.\n' ...
        'excess: median rate beyond the edge. ratio: median Re lambda_1 / edge. f_det, f_uns, f_both: fraction of\n' ...
        'states with the rightmost beyond the edge, beyond zero, beyond both. X: mean exposure (1/s).\n' ...
        'n_out: mean number of right outliers. dw_*: median dwell beyond the edge, zero, both.\n' ...
        'same: fraction of consecutive states led by the same mode. K1, K_err: K_eff prediction (units of W) and\n' ...
        'its median error where an outlier is predicted. full_ab, f_fu: full-Jacobian abscissa, fraction > 0.\n']);
end
end

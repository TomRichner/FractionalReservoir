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

% Order: network, the conditions in the paper's order, then seed.
cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
[~, ci] = ismember(T.condition, cond_order);
[~, order] = sortrows(table(T.network, ci, T.seed));
T = T(order, :);
M = M(order);
writetable(T, fullfile(data_dir, 'outlier_measures.csv'));

if opts.print
    fprintf('\nOutlier measures (rates 1/s; dwell s)\n');
    fprintf('%-10s %-14s %4s | %7s | %5s %5s | %5s %5s %6s | %5s %5s | %6s %6s | %5s %6s %6s | %6s %5s %6s %6s | %12s %5s\n', ...
        'network', 'condition', 'seed', 'lam_1', 'kappa', 'agree', 'rho', 'rho_u', 'reined', 'on_E', 'on_I', ...
        'edge', 'top', 'f_uns', 'dw_0', 'dw_max', 'full', 'f_fu', 'dw_fu', 'X_full', 'continuation', 'N_eff');
    for k = 1:height(T)
        r = T(k, :);
        fprintf('%-10s %-14s %4d | %+7.3f | %5.3f %5.3f | %5.2f %5.2f %6.2f | %5.2f %5.2f | %+6.2f %+6.2f | %5.2f %6.2f %6.2f | %+6.2f %5.2f %6.2f %6.3f | %+5.2f%+5.2fi %5.1f\n', ...
            r.network{1}, r.condition{1}, r.seed, r.lambda_1, r.kappa, r.agree_re, ...
            r.ratio_med, r.rho_uniform, r.reined, r.frac_on_E, r.frac_on_I, ...
            r.rate_edge_med, r.rate_top_med, r.frac_unstable, r.dwell_zero_med, r.dwell_zero_max, ...
            r.full_abscissa_med, r.full_frac_unstable, r.full_dwell_zero_med, r.X_full, ...
            r.cont_re_med, r.cont_im_med, r.N_eff);
    end
    fprintf(['\nlam_1: leading Lyapunov exponent. kappa, agree: 95th percentile and median of the null''s rightmost real part over R_eff.\n' ...
        'rho: median Re lambda_1 / edge (1 = at the edge of the bulk). rho_u: the same with uniform (rms) gains: the structure alone.\n' ...
        'reined: fraction of the structural detachment removed. on_E, on_I: fraction of neurons off the floor and ceiling.\n' ...
        'edge, top: median frozen-gain rate of the bulk edge and of the rightmost eigenvalue. f_uns, dw_0, dw_max: fraction of states\n' ...
        'with the frozen rightmost beyond zero, median and longest dwell there. full, f_fu, dw_fu, X_full: the FULL Jacobian''s\n' ...
        'rightmost real part (median), fraction of states beyond zero, median dwell there (0.1 s steps), mean positive part.\n' ...
        'continuation: the full-Jacobian eigenvalue carrying the frozen rightmost pattern (median over the dense states).\n' ...
        'N_eff: number of directions that take turns as the rightmost frozen mode.\n']);
end
end

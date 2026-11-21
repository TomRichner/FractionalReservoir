function [T, M, MM, Rs] = outlier_mode_table(data_dir, opts)
%OUTLIER_MODE_TABLE The intended-mode measures of every run in a folder: one row per run, printed, and saved as CSV.
%
%   T = outlier_mode_table()                 data/explorations/outlier_structures
%   [T, M, MM, Rs] = outlier_mode_table(data_dir, 'print', true, 'keep_runs', false)
%
% T is a table of outlier_mode_measures' rows, in the order of the files
% sorted by network, condition (the paper's order) and seed. M and MM are
% cells of the outlier_measures structs and of the mode series; Rs the
% runs themselves when keep_runs is true (about 60 MB each). The CSV is
% outlier_mode_measures.csv in data_dir.
%
% See also: outlier_mode_measures, outlier_run, fig_outlier_structures

arguments
    data_dir (1,:) char = ''
    opts.print (1,1) logical = true
    opts.keep_runs (1,1) logical = false
end
if isempty(data_dir)
    data_dir = fullfile(fileparts(which('setup_paths')), 'data', 'explorations', 'outlier_structures');
end
files = dir(fullfile(data_dir, 'run_*.mat'));
if isempty(files)
    error('outlier_mode_table:NoRuns', 'No run_*.mat in %s.', data_dir);
end

n = numel(files);
rows = cell(1, n); M = cell(1, n); MM = cell(1, n); Rs = cell(1, n);
for k = 1:n
    L = load(fullfile(files(k).folder, files(k).name), 'R');
    M{k} = outlier_measures(L.R);
    [rows{k}, MM{k}] = outlier_mode_measures(L.R, M{k});
    if opts.keep_runs; Rs{k} = L.R; end
end
T = struct2table([rows{:}], 'AsArray', true);
if ischar(T.network);   T.network = {T.network};     end   % a single run gives char columns
if ischar(T.condition); T.condition = {T.condition}; end
if ischar(T.target);    T.target = {T.target};       end

cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
[~, ci] = ismember(T.condition, cond_order);
[~, order] = sortrows(table(T.network, T.seed, ci));
T = T(order, :); M = M(order); MM = MM(order); Rs = Rs(order);
writetable(T, fullfile(data_dir, 'outlier_mode_measures.csv'));

if opts.print
    fprintf('\nIntended-mode measures (rates 1/s, rotation rad/s)\n');
    fprintf('%-24s %4s %-14s | %11s | %7s | %5s %5s %6s | %4s %4s | %6s %6s %6s | %5s %5s %5s | %5s | %6s %6s | %5s %5s | %5s %5s\n', ...
        'network', 'seed', 'condition', 'lam0 / R', 'lam_1', 'rho', 'rho_u', 'reined', 'ov', 'f_rt', ...
        'r_str', 'r_frz', 'r_full', 'w_str', 'w_frz', 'w_ful', 'g_rat', 'dl_obs', 'dl_prd', 'on_tp', 'on_rs', 'on_E', 'on_I');
    for k = 1:height(T)
        r = T(k, :);
        fprintf('%-24s %4d %-14s | %+5.2f%+5.2fi | %+7.3f | %5.2f %5.2f %6.2f | %4.2f %4.2f | %+6.2f %+6.2f %+6.2f | %5.2f %5.2f %5.2f | %5.2f | %+6.2f %+6.2f | %5.2f %5.2f | %5.2f %5.2f\n', ...
            r.network{1}, r.seed, r.condition{1}, r.lam0_re_R, r.lam0_im_R, r.lambda_1, ...
            r.mode_rho, r.mode_rho_uniform, r.mode_reined, r.mode_ov, r.mode_frac_rightmost, ...
            r.rate_struct, r.rate_frozen, r.rate_full, r.rot_struct, r.rot_frozen, r.rot_full, ...
            r.gain_ratio, r.shift_obs_re, r.shift_pred_re, r.on_top, r.on_rest, r.frac_on_E, r.frac_on_I);
    end
    fprintf(['\nlam0 / R: the intended structural eigenvalue over the structural bulk radius. lam_1: leading Lyapunov exponent.\n' ...
        'rho, rho_u, reined: detachment of the mode (Re / bulk edge), the same with uniform gains, the fraction removed.\n' ...
        'ov: overlap of the followed eigenvector with the structural one. f_rt: fraction of states where the mode is the rightmost.\n' ...
        'r_str, r_frz, r_full: real part of the mode with uniform gains, frozen gains, and in the full Jacobian. w_*: its rotation rate.\n' ...
        'g_rat: participation-weighted gain / mean gain. dl_obs, dl_prd: real part of the relative shift of the eigenvalue from the\n' ...
        'uniform-gain picture, observed and predicted to first order. on_tp, on_rs: responsive fraction of the top carriers and the rest.\n']);
end
end

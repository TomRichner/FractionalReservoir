% EXPLORE_SMALL_NETWORK_SETTLING Does starting the Lyapunov computation after a settling period change the top-K spectrum of the 80-neuron network?
%
% WHY. explore_small_network_spectrum starts cold: the network starts at
% t = 0, the tangent basis aligns over 0 to 5 s and the exponents
% accumulate over 5 to 15 s. The slowest SFA rung has time constants of 5
% to 20 s, so that window sits inside the network's transient. A 30 s
% check at K = 100 moved every exponent up by about 0.018 1/s, but it
% changed two things at once: later states and a longer window. This
% script separates them with three runs on the same network and seed
% (same x0; the noise over the common stretch is NOT the same, see below):
%
%   A  cold, as before:       window [5 15] s,   warm-up 5 s
%   B  settled, same lengths: window [65 75] s,  warm-up 5 s  (60 s of
%                             settling, then the same 5 s warm-up and 10 s
%                             window: isolates the settling)
%   C  settled and long:      window [80 140] s, warm-up 20 s (iteration
%                             starts at 60 s; the better estimate)
%
% The runs are rows of a small settings table (label, T_range,
% lya_T_interval, lya_warmup); add a row to add a run. Every run uses
% T_range = [0 T] with T the end of its window, so every run starts from
% the same state x0. The noise is NOT the same over the common stretch:
% sra1 takes two Gaussians per step, and wiener_increments draws them as
% xi1 = randn(n, n_steps) and THEN xi2 = randn(n, n_steps), so xi1 is
% shared over the common steps but xi2 depends on the run's length. The
% script prints how far the rates over the first 15 s differ between runs
% (2026-11-19: up to 0.015 in r, from the first step on).
%
% THE NETWORK AND THE REST OF THE RUN are those of
% explore_small_network_spectrum: small_network_model, the condition
% sfa3_stdc2 by default, seed [7 8], 500 Hz, sra1, topk with lya_dt 0.05 s,
% K = 200.
%
% MEMORY. During run() the class holds the whole trajectory (nt x N
% doubles; at 140 s, 500 Hz and 640 states about 360 MB). The script
% prints the estimate for the longest run and stops if it exceeds
% mem_limit_gb.
%
% WHAT IT PRINTS
%   - per run: lambda_1, range, number positive, wall time (total, and the
%     Lyapunov part as lyapunov_topk times it)
%   - rank-by-rank differences of the sorted spectra for every pair of
%     runs (B - A, C - A, C - B): median, min and max, and the median over
%     the top n_top ranks against the rest
%   - in the settled run (settle_run): population means of the SFA term
%     (c/K) sum_k a_ik, of the chained pool's b1 and b2 (all routes) and of
%     the rate r, over each of stat_windows (the cold window and the start
%     and end of the settled run), read from the stored output (plot_data
%     at fs / plot_deci = 100 Hz)
%   - in the settled run: a within-run convergence indicator, the running
%     finite-time exponents at the midpoint of the window against the full
%     window (lya_results.finite_LE_spectrum_t)
%   - the largest gaps between consecutive sorted exponents in the settled
%     run, with their ranks
%
% THE FIGURE (one, left open, PNG to figs/explorations/small_network/):
%   (a) the sorted spectra of all runs overlaid, exponent against rank
%   (b) histograms of the first run and the settled run on shared bins,
%       drawn as outlines
% Both panels carry a zero line.
%
% POOL: none; pool_set(0) records that this session holds none.
%
% PATH: run setup_paths once in the session first, with the current folder
% at the repository root (it adds scripts/ with genpath, which includes
% this folder). Takes about three minutes on a laptop.
%
% See also: explore_small_network_spectrum, small_network_model, lyapunov_topk

%% Settings
condition  = 'sfa3_stdc2';
lya_K      = 200;
rng_seeds  = [7 8];
fs         = 500;           % Hz
lya_dt     = 0.05;          % s
plot_deci  = 5;             % stored output at fs / plot_deci = 100 Hz
% The runs: label, T_range, lya_T_interval, lya_warmup (one row each).
runs = cell2struct({ ...
    'A', [0 15],  [5 15],   5; ...
    'B', [0 75],  [65 75],  5; ...
    'C', [0 140], [80 140], 20}, ...
    {'label', 'T_range', 'lya_T_interval', 'lya_warmup'}, 2);
settle_run   = 3;           % the run whose trajectory and convergence are examined
stat_windows = {[5 15], [60 80], [120 140]};   % s, windows for the settling means
n_top        = 60;          % ranks 1..n_top against the rest, for the shifts
n_gaps       = 8;
mem_limit_gb = 2;           % stop if the longest run's stored trajectory would exceed this
save_fig     = true;
out_dir      = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', 'small_network');

%% No pool
pool_set(0, 'label', 'small_network_settling');

%% Memory check: the trajectory run() holds for the longest run
n_runs = numel(runs);
m_probe = small_network_model(condition, 'seed', rng_seeds, 'T_range', [0 1], 'fs', fs);
n_state = m_probe.N_sys_eqs;
T_longest = max(arrayfun(@(r) diff(r.T_range), runs));
mem_gb = 8 * (T_longest * fs + 1) * n_state / 1e9;
fprintf('memory: longest run %g s x %g Hz x %d states = %.2f GB as doubles (limit %.1f GB)\n', ...
    T_longest, fs, n_state, mem_gb, mem_limit_gb);
if mem_gb > mem_limit_gb
    error('explore_small_network_settling:Memory', ...
        'The longest run would store %.2f GB; shorten it or raise mem_limit_gb.', mem_gb);
end
clear m_probe

%% Run
lambda    = cell(1, n_runs);
results   = cell(1, n_runs);
plot_data = cell(1, n_runs);
wall_s    = zeros(1, n_runs);
lya_s     = zeros(1, n_runs);
for k = 1:n_runs
    t_start = tic;
    mk = small_network_model(condition, 'seed', rng_seeds, ...
        'T_range', runs(k).T_range, 'fs', fs, 'lya_method', 'topk', ...
        'lya_K', lya_K, 'lya_dt', lya_dt, ...
        'lya_T_interval', runs(k).lya_T_interval, 'lya_warmup', runs(k).lya_warmup, ...
        'model_args', {'plot_deci', plot_deci});
    mk.run();
    wall_s(k) = toc(t_start);
    results{k} = mk.lya_results;
    lya_s(k) = mk.lya_results.seconds;
    lambda{k} = sort(mk.lya_results.LE_spectrum(:), 'descend');
    plot_data{k} = mk.plot_data;
    if k == settle_run
        params_settle = mk.get_params();
        c_eff = mk.c ./ max(1, mk.n_a);
        types = mk.cell_type_names;
    end
    clear mk
    fprintf('run %s done in %.1f s\n', runs(k).label, wall_s(k));
end

%% Per run
fprintf('\n=== n = 80, %s, seed [%d %d], K = %d ===\n', condition, rng_seeds, lya_K);
for k = 1:n_runs
    lam = lambda{k};
    fprintf(['  %s  window [%g %g] s, warm-up %g s, T_range [%g %g]:  lambda_1 %+.4f, ' ...
        'range %+.4f to %+.4f, positive %d;  wall %.1f s (Lyapunov %.1f s)\n'], ...
        runs(k).label, runs(k).lya_T_interval, runs(k).lya_warmup, runs(k).T_range, ...
        lam(1), lam(1), lam(end), nnz(lam > 0), wall_s(k), lya_s(k));
end

% The trajectories should agree over their common start (same seed).
t_common = 15;
for k = 2:n_runs
    r1 = rates_until(plot_data{1}, t_common);
    rk = rates_until(plot_data{k}, t_common);
    fprintf('  rates over 0-%g s, run %s against run %s: max |difference| %.3g\n', ...
        t_common, runs(k).label, runs(1).label, max(abs(r1(:) - rk(:))));
end

%% Rank-by-rank differences
fprintf('\n  rank-by-rank differences of the sorted spectra (1/s):\n');
for k2 = 2:n_runs
    for k1 = 1:k2-1
        d = lambda{k2} - lambda{k1};
        fprintf(['    %s - %s: median %+.4f, min %+.4f (rank %d), max %+.4f (rank %d); ' ...
            'median over ranks 1-%d %+.4f, over %d-%d %+.4f\n'], ...
            runs(k2).label, runs(k1).label, median(d), min(d), find(d == min(d), 1), ...
            max(d), find(d == max(d), 1), n_top, median(d(1:n_top)), ...
            n_top + 1, numel(d), median(d(n_top+1:end)));
    end
end

%% Settling of the trajectory in the settled run
pd = plot_data{settle_run};
fprintf('\n  run %s: population means over the stored output (%g Hz)\n', ...
    runs(settle_run).label, fs / plot_deci);
fprintf('    %-12s %10s %10s %10s %10s\n', 'window (s)', 'SFA term', 'b1', 'b2', 'rate r');
for w = 1:numel(stat_windows)
    win = stat_windows{w};
    keep = pd.t >= win(1) & pd.t <= win(2);
    [sfa_all, r_all] = deal([]);
    for q = 1:numel(types)
        a_q = pd.a.(types{q});                               % n_q x K_a x nt
        if ~isempty(a_q)
            sfa_all = [sfa_all; c_eff(q) * squeeze(sum(a_q(:, :, keep), 2))]; %#ok<AGROW>
        end
        r_all = [r_all; pd.r.(types{q})(:, keep)]; %#ok<AGROW>
    end
    b_comp = {[], []};
    for pre = 1:numel(types)
        for post = 1:numel(types)
            if isfield(pd.b, types{pre}) && isfield(pd.b.(types{pre}), types{post})
                b_route = pd.b.(types{pre}).(types{post});   % n_pre x n_b x nt
                for j = 1:min(2, size(b_route, 2))
                    b_comp{j} = [b_comp{j}; squeeze(b_route(:, j, keep))];
                end
            end
        end
    end
    fprintf('    %-12s %10.4f %10.4f %10.4f %10.4f\n', sprintf('%g-%g', win), ...
        mean_or_nan(sfa_all), mean_or_nan(b_comp{1}), mean_or_nan(b_comp{2}), mean_or_nan(r_all));
end
if params_settle.N_sys_eqs ~= n_state
    fprintf('    (state dimension changed between the probe and the run)\n');
end

%% Within-run convergence in the settled run
res = results{settle_run};
win = runs(settle_run).lya_T_interval;
t_mid = mean(win);
row_mid = find(res.t_lya >= t_mid, 1);
row_end = find(all(isfinite(res.finite_LE_spectrum_t), 2), 1, 'last');
fin_mid = res.finite_LE_spectrum_t(row_mid, :)';     % columns in the final order
fin_end = res.finite_LE_spectrum_t(row_end, :)';
d_col = fin_end - fin_mid;
d_sorted = sort(fin_end, 'descend') - sort(fin_mid, 'descend');
fprintf(['\n  run %s convergence: running exponents at t = %.2f s (window %g-%g s, ' ...
    'half of it accumulated) against t = %.2f s\n'], ...
    runs(settle_run).label, res.t_lya(row_mid), win, res.t_lya(row_end));
fprintf('    lambda_1 %+.4f -> %+.4f (change %+.4f)\n', max(fin_mid), max(fin_end), max(fin_end) - max(fin_mid));
fprintf('    per column: median change %+.4f, max |change| %.4f\n', median(d_col), max(abs(d_col)));
fprintf('    sorted rank by rank: median change %+.4f, max |change| %.4f\n', median(d_sorted), max(abs(d_sorted)));
fprintf('    (final spectrum against the stored LE_spectrum: max |difference| %.3g)\n', ...
    max(abs(sort(fin_end, 'descend') - lambda{settle_run})));

%% Gaps in the settled run
lam = lambda{settle_run};
steps = lam(1:end-1) - lam(2:end);
[step_sorted, step_rank] = sort(steps, 'descend');
fprintf('\n  run %s, largest steps between consecutive exponents (rank i: between i and i+1):\n', ...
    runs(settle_run).label);
for g = 1:min(n_gaps, numel(steps))
    i = step_rank(g);
    fprintf('    rank %3d: %+.4f -> %+.4f, step %.4f\n', i, lam(i), lam(i + 1), step_sorted(g));
end
fprintf('  median step %.4f\n', median(steps));

%% Figure
fig = figure('Name', sprintf('Small network (n = 80): %s, top-%d, settling', condition, lya_K), ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1000 650]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('n = 80, %s, seed [%d %d], K = %d: cold start against settled', ...
    strrep(condition, '_', '\_'), rng_seeds, lya_K));
colors  = lines(n_runs);
markers = {'o', 's', '^', 'd', 'v', 'x'};

ax1 = nexttile(tl);
hold(ax1, 'on');
for k = 1:n_runs
    plot(ax1, 1:numel(lambda{k}), lambda{k}, markers{mod(k - 1, numel(markers)) + 1}, ...
        'MarkerSize', 3.5, 'Color', colors(k, :), 'MarkerFaceColor', colors(k, :), ...
        'DisplayName', sprintf('%s: window %g-%g s, \\lambda_1 = %+.3f', ...
        runs(k).label, runs(k).lya_T_interval, lambda{k}(1)));
end
yline(ax1, 0, 'k--', 'HandleVisibility', 'off');
hold(ax1, 'off');
xlim(ax1, [0, lya_K + 1]);
xlabel(ax1, 'Rank');
ylabel(ax1, 'Lyapunov exponent (1/s)');
title(ax1, '(a) Sorted exponents');
legend(ax1, 'Location', 'southwest');
grid(ax1, 'on');
box(ax1, 'off');

ax2 = nexttile(tl);
both = [lambda{1}; lambda{settle_run}];
range_both = max(both) - min(both);
bin_width = max(0.005, round(range_both * 8 / numel(lambda{1}) / 0.005) * 0.005);
edges = (floor(min(both) / bin_width) * bin_width):bin_width:(ceil(max(both) / bin_width) * bin_width + bin_width);
hold(ax2, 'on');
for k = [1 settle_run]
    histogram(ax2, lambda{k}, edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, ...
        'EdgeColor', colors(k, :), 'DisplayName', sprintf('%s: window %g-%g s', ...
        runs(k).label, runs(k).lya_T_interval));
end
xline(ax2, 0, 'k--', 'HandleVisibility', 'off');
hold(ax2, 'off');
xlim(ax2, [edges(1), max(edges(end), 0) + bin_width]);
xlabel(ax2, 'Lyapunov exponent (1/s)');
ylabel(ax2, 'Count');
title(ax2, sprintf('(b) Histograms, bins %.3g 1/s', bin_width));
legend(ax2, 'Location', 'northwest');
box(ax2, 'off');
fprintf('\n  histogram bins %.3f 1/s (about 8 exponents per bin over the joint range)\n', bin_width);

if save_fig
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    png = fullfile(out_dir, sprintf('small_network_settling_%s_K%d.png', condition, lya_K));
    exportgraphics(fig, png, 'Resolution', 150);
    fprintf('  saved %s\n', png);
end

%% Local functions
function r = rates_until(pd, t_end)
% RATES_UNTIL All neurons' rates (rows) over 0..t_end of the stored output.
keep = pd.t <= t_end + 1e-9;
names = fieldnames(pd.r);
r = [];
for q = 1:numel(names)
    r = [r; pd.r.(names{q})(:, keep)]; %#ok<AGROW>
end
end

function v = mean_or_nan(x)
% MEAN_OR_NAN The mean of every entry, NaN for an empty array.
if isempty(x)
    v = NaN;
else
    v = mean(x, 'all');
end
end

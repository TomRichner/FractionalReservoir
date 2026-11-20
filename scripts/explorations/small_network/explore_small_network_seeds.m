% EXPLORE_SMALL_NETWORK_SEEDS The settled top-K Lyapunov spectrum of the 80-neuron network on ten network seeds, for three conditions.
%
% WHY. On seeds 1 and 7 the chained multiple-timescale network (sfa3_stdc2,
% "MTS") gave nearly the same spectrum (lambda_1 -0.073 and -0.071), while
% no adaptation (-6.0 and +0.17) and one timescale (sfa1_std1, "1TS";
% +0.0025 and -0.41) did not. Tom's observation: in a small network chained
% MTS seems to have a much lower seed-to-seed variance of lambda_1 than the
% other two. This script looks at ten seeds. It is an exploration: nothing
% is tuned.
%
% THE NETWORK: the preset ..._mu7chain_n80 (the paper's network at n = 80;
% srnn_param_preset). Seed k is rng_seeds [k, k+1], as run_lyapunov_spectrum
% numbers them.
%
% THE RUN: run_lyapunov_spectrum with the arguments of the opt-in stage
% lyapunov_spectrum_n80 (lyapunov_spectrum_n80_stage_plan): K = K_max = 240
% (capped at the state dimension, so no adaptation gets all 80), T 140 s,
% exponents accumulated over [80 140] s after a 20 s warm-up, fs 500 Hz,
% the preset's noise only (no noise-off variant). One call per condition,
% into data/explorations/small_network_seeds/<condition>/.
%   The spectrum runs keep no trajectory, so the regime of each seed comes
% from trajectory-only runs of the same models (build_from_preset, same
% preset, condition, rng_seeds, T_range [0 140], fs 500, lya_method 'none',
% output at 100 Hz), in a parfor, saved beside the spectra as
% regime_data.mat. The trajectory does not depend on the Lyapunov method but
% does depend on T_range (the noise draw depends on the run length), so the
% T_range is the spectrum runs'.
%   recompute = false loads what is on disk and computes only what is
% missing; max_new_conditions (default Inf; set it in the workspace before
% running to override) stops the script after that many conditions were
% computed, so each MATLAB call stays short.
%
% WHAT IT PRINTS, per condition: lambda_1 of every seed, its median, mean,
% standard deviation, min and max; per seed the number of positive
% exponents, the drift of lambda_1 over the last quarter of the window
% (lya_summary's lambda_1_drift), and the regime over [80 140] s: mean rate,
% share of neuron-time at the floor (r = 0), at the ceiling (r = 1) and
% between, neurons silent or saturated throughout, the standard deviation
% of the population-mean rate and its dominant frequency. Then the seed-to-
% seed standard deviation of the exponent at fixed ranks.
%
% THE FIGURE (one, 1300 x 850, left open, PNG to
% figs/explorations/small_network/): three columns (no adaptation, 1TS,
% MTS) and a fourth for the comparison. A seed has one color everywhere.
% Row 1: sorted exponents against rank. Row 2, on shared bins (about 8
% exponents per bin per seed, median over seeds): the pooled histogram as
% count per seed (black) with the per-bin min-max over seeds (blue band),
% or, with hist_style 'outlines', per-seed outlines in the seed colors over
% the pooled histogram (gray). Row 3: lambda_1 per seed with the median. Column 4:
% lambda_1 of all three conditions on one symmetric-log axis (linear within
% +-symlog_c, logarithmic beyond).
%
% PATH: run setup_paths once in the session first, with the current folder
% at the repository root. Computing takes about ten minutes on a laptop
% with a pool of 3; loading and drawing a few seconds.
%
% See also: explore_small_network_spread, run_lyapunov_spectrum,
%           lyapunov_spectrum_n80_stage_plan

%% Settings
preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain_n80';
conditions  = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
cond_labels = {'no adaptation', '1TS (sfa1\_std1)', 'chained MTS (sfa3\_stdc2)'};
n_seeds     = 10;
K           = 240;           % K_max = K: no retry doubling
T           = 140;           % s
lya_window  = [80 140];      % s
lya_warmup  = 20;            % s
fs          = 500;           % Hz
plot_deci   = 5;             % trajectory output at 100 Hz
ranks_sd    = [1 10 40 80 160 240];   % ranks for the seed-to-seed SD
recompute   = false;         % true: recompute everything; false: compute only what is missing
if ~exist('max_new_conditions', 'var'), max_new_conditions = Inf; end
n_pool      = 3;
hist_style  = 'band';        % 'band' (pooled, with the min-max over seeds) or 'outlines' (per-seed stairs);
                             % ten outlines are busy, the band reads better
symlog_c    = 0.01;          % 1/s, linear range of the symmetric-log axis
tau_d       = 0.1;           % s, the preset's membrane time constant: pinned neurons give -1/tau_d
pin_tol     = 0.1;           % 1/s, an exponent within this of -1/tau_d counts as a pinned direction
save_fig    = true;
root        = fileparts(which('setup_paths'));
data_root   = fullfile(root, 'data', 'explorations', 'small_network_seeds');
out_dir     = fullfile(root, 'figs', 'explorations', 'small_network');

% Ten distinguishable colors (the Tableau 10 set), one per seed.
seed_colors = [ ...
    0.122 0.467 0.706; 1.000 0.498 0.055; 0.173 0.627 0.173; 0.839 0.153 0.157; ...
    0.580 0.404 0.741; 0.549 0.337 0.294; 0.890 0.467 0.761; 0.498 0.498 0.498; ...
    0.737 0.741 0.133; 0.090 0.745 0.812];

%% Compute what is missing (one run_lyapunov_spectrum call per condition)
n_cond = numel(conditions);
spec_file = @(c) fullfile(data_root, conditions{c}, 'lyapunov_spectrum_data.mat');
reg_file  = @(c) fullfile(data_root, conditions{c}, 'regime_data.mat');
todo = find(arrayfun(@(c) recompute || ~isfile(spec_file(c)) || ~isfile(reg_file(c)), 1:n_cond));
todo = todo(1:min(numel(todo), max_new_conditions));
if ~isempty(todo)
    pool_set(n_pool, 'label', 'small_network_seeds');
    for c = todo
        cname = conditions{c};
        cdir = fileparts(spec_file(c));
        if ~isfolder(cdir), mkdir(cdir); end
        t0 = tic;
        run_lyapunov_spectrum('preset_name', preset_name, 'run_mode', 'production', ...
            'out_dir', cdir, 'verbose', 'minimal', 'conditions', {cname}, ...
            'n_seeds', n_seeds, 'K', K, 'K_max', K, 'T', T, 'lya_T_interval', lya_window, ...
            'lya_warmup', lya_warmup, 'fs', fs, 'noise_off_too', false);
        t_spec = toc(t0);
        t0 = tic;
        regime = cell(1, n_seeds);
        parfor k = 1:n_seeds
            regime{k} = regime_trial(preset_name, cname, [k, k + 1], T, fs, plot_deci, lya_window);
        end
        regime = [regime{:}];
        t_reg = toc(t0);
        save(reg_file(c), 'regime', 'preset_name', 'cname', 'T', 'fs', 'plot_deci', 'lya_window');
        fprintf('[compute] %-14s spectra %.0f s, trajectories %.0f s\n', cname, t_spec, t_reg);
    end
    pool_release();
    if isempty(gcp('nocreate'))
        fprintf('[pool] released; no pool is up\n');
    else
        fprintf('[pool] WARNING: a pool is still up after pool_release\n');
    end
end
if any(arrayfun(@(c) ~isfile(spec_file(c)) || ~isfile(reg_file(c)), 1:n_cond))
    fprintf('[compute] conditions still missing; run again to compute the rest\n');
    return
end

%% Load
D = struct('lam', {}, 'r', {}, 'reg', {}, 'minutes', {});
for c = 1:n_cond
    S = load(spec_file(c), 'results', 'settings');
    R = load(reg_file(c), 'regime');
    runs = S.results(1).noise_on;
    lam = NaN(n_seeds, max(arrayfun(@(r) numel(r.LE_spectrum), runs)));
    for k = 1:n_seeds
        l = sort(runs(k).LE_spectrum(:), 'descend')';
        lam(k, 1:numel(l)) = l;
    end
    D(c) = struct('lam', lam, 'r', runs, 'reg', R.regime, 'minutes', S.settings.minutes);
end

%% Print, per condition
for c = 1:n_cond
    lam = D(c).lam; r = D(c).r; g = D(c).reg;
    l1 = lam(:, 1);
    fprintf('\n=== %s (state dimension %d, K used %d; spectra %.1f min) ===\n', conditions{c}, ...
        r(1).N, r(1).K_used, D(c).minutes);
    fprintf('  lambda_1: median %+.4f, mean %+.4f, SD %.4f, min %+.4f, max %+.4f\n', ...
        median(l1), mean(l1), std(l1), min(l1), max(l1));
    fprintf(['  seed  lambda_1  lambda_2  n_pos  drift_l1   rate   floor  ceil   betw  ' ...
        'silent sat  popSD   f_pk(Hz)  n_pin\n']);
    fprintf('  (n_pin: exponents within %g 1/s of -1/tau_d = %g; drift_l1: lambda_1 change over the last quarter of the window)\n', pin_tol, -1 / tau_d);
    for k = 1:n_seeds
        drift = NaN;
        if isfield(r(k), 'lambda_1_drift'), drift = r(k).lambda_1_drift; end
        fprintf('  %3d  %+8.4f  %+8.4f  %4d  %+8.4f  %6.3f  %5.3f  %5.3f  %5.3f  %4d  %4d  %6.4f  %6.3f   %3d\n', ...
            k, lam(k, 1), lam(k, 2), nnz(lam(k, :) > 0), drift, g(k).rate, g(k).floor, g(k).ceiling, ...
            g(k).between, g(k).n_silent, g(k).n_saturated, g(k).pop_sd, g(k).f_peak, ...
            nnz(abs(lam(k, :) + 1 / tau_d) < pin_tol));
    end
end

%% Seed-to-seed spread across the spectrum
fprintf('\n=== seed-to-seed SD of the exponent at fixed ranks (1/s) ===\n');
fprintf('  %-14s', 'rank'); fprintf(' %8d', ranks_sd); fprintf('   median over all ranks\n');
for c = 1:n_cond
    lam = D(c).lam;
    fprintf('  %-14s', conditions{c});
    for q = ranks_sd
        if q <= size(lam, 2), fprintf(' %8.4f', std(lam(:, q))); else, fprintf(' %8s', '-'); end
    end
    fprintf('   %.4f\n', median(std(lam, 0, 1)));
end
fprintf('  (IQR of lambda_1 over seeds:');
for c = 1:n_cond, fprintf(' %s %.4f;', conditions{c}, iqr(D(c).lam(:, 1))); end
fprintf(')\n');

%% Figure
fig = figure('Name', 'Small network (n = 80): spectra over ten seeds', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [30 40 1300 850]);
tl = tiledlayout(fig, 3, n_cond + 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['n = 80 (mu7chain\\_n80), seeds [k k+1], k = 1..%d; window %g-%g s, ' ...
    'warm-up %g s, K = min(%d, state dimension)'], n_seeds, lya_window, lya_warmup, K));
h_seed = gobjects(1, n_seeds);
for c = 1:n_cond
    lam = D(c).lam;
    K_c = size(lam, 2);
    lo = min(lam(:)); hi = max([lam(:); 0]);
    pad = 0.04 * (hi - lo);
    lim = [lo - pad, hi + pad];

    % Row 1: ranked spectra
    ax = nexttile(tl, c);
    hold(ax, 'on');
    for k = 1:n_seeds
        h_seed(k) = plot(ax, 1:K_c, lam(k, :), '-', 'Color', seed_colors(k, :), 'LineWidth', 0.9, ...
            'DisplayName', sprintf('seed %d', k));
    end
    yline(ax, 0, 'k--', 'HandleVisibility', 'off');
    hold(ax, 'off');
    ylim(ax, lim); xlim(ax, [0, K_c + 1]);
    xlabel(ax, 'Rank'); ylabel(ax, 'Exponent (1/s)');
    title(ax, cond_labels{c});
    grid(ax, 'on'); box(ax, 'off');

    % Row 2: distributions on shared bins (8 exponents per bin per seed, median over seeds)
    ax = nexttile(tl, (n_cond + 1) + c);
    bw_seed = max(0.005, round((max(lam, [], 2) - min(lam, [], 2)) * 8 / K_c / 0.005) * 0.005);
    bw = median(bw_seed);
    edges = (floor(lo / bw) * bw):bw:(ceil(hi / bw) * bw + bw);
    counts = zeros(n_seeds, numel(edges) - 1);
    for k = 1:n_seeds, counts(k, :) = histcounts(lam(k, :), edges); end
    hold(ax, 'on');
    pooled = mean(counts, 1);
    if strcmp(hist_style, 'band')
        xs = reshape([edges(1:end-1); edges(2:end)], 1, []);
        lo_b = reshape(repmat(min(counts, [], 1), 2, 1), 1, []);
        hi_b = reshape(repmat(max(counts, [], 1), 2, 1), 1, []);
        fill(ax, [xs, fliplr(xs)], [lo_b, fliplr(hi_b)], [0.6 0.75 0.95], 'EdgeColor', 'none', ...
            'FaceAlpha', 0.6);
        stairs(ax, [edges, edges(end)], [pooled, pooled(end), 0], 'k', 'LineWidth', 1.2);
    else
        histogram(ax, 'BinEdges', edges, 'BinCounts', pooled, 'FaceColor', [0.8 0.8 0.8], ...
            'EdgeColor', 'none');
        for k = 1:n_seeds
            histogram(ax, 'BinEdges', edges, 'BinCounts', counts(k, :), 'DisplayStyle', 'stairs', ...
                'EdgeColor', seed_colors(k, :), 'LineWidth', 0.8);
        end
    end
    xline(ax, 0, 'k--');
    hold(ax, 'off');
    xlim(ax, lim);
    xlabel(ax, 'Exponent (1/s)'); ylabel(ax, 'Count per seed');
    title(ax, sprintf('bins %.3g 1/s (~8 per bin per seed)', bw), 'FontWeight', 'normal');
    box(ax, 'off');

    % Row 3: lambda_1 per seed
    ax = nexttile(tl, 2 * (n_cond + 1) + c);
    l1 = lam(:, 1);
    hold(ax, 'on');
    yline(ax, 0, 'k--');
    yline(ax, median(l1), '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    scatter(ax, 1:n_seeds, l1, 36, seed_colors, 'filled');
    hold(ax, 'off');
    span = max(l1) - min(l1); span = max(span, 1e-3);
    ylim(ax, [min([l1; 0]) - 0.1 * span, max([l1; 0]) + 0.1 * span]);
    xlim(ax, [0.5, n_seeds + 0.5]); xticks(ax, 1:n_seeds);
    xlabel(ax, 'Seed'); ylabel(ax, '\lambda_1 (1/s)');
    title(ax, sprintf('\\lambda_1: median %+.3g, SD %.3g, range [%+.3g, %+.3g]', ...
        median(l1), std(l1), min(l1), max(l1)), 'FontWeight', 'normal', 'FontSize', 8);
    grid(ax, 'on'); box(ax, 'off');
end

% Column 4: lambda_1 of the three conditions on one symmetric-log axis
ax = nexttile(tl, n_cond + 1, [3 1]);
symlog = @(y) sign(y) .* log10(1 + abs(y) / symlog_c);
hold(ax, 'on');
yline(ax, 0, 'k--');
offs = linspace(-0.3, 0.3, n_seeds);
for c = 1:n_cond
    l1 = D(c).lam(:, 1);
    plot(ax, c + [-0.38 0.38], symlog(median(l1)) * [1 1], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    scatter(ax, c + offs, symlog(l1), 30, seed_colors, 'filled');
end
hold(ax, 'off');
tick_vals = [-10 -1 -0.1 -0.01 0 0.01 0.1 1];
yticks(ax, symlog(tick_vals));
yticklabels(ax, arrayfun(@(v) sprintf('%g', v), tick_vals, 'UniformOutput', false));
all_l1 = cell2mat(arrayfun(@(d) d.lam(:, 1), D(:), 'UniformOutput', false));
ylim(ax, [symlog(min(all_l1)) - 0.2, symlog(max(all_l1)) + 0.2]);
xlim(ax, [0.4, n_cond + 0.6]); xticks(ax, 1:n_cond);
xticklabels(ax, {'none', '1TS', 'MTS'});
ylabel(ax, sprintf('\\lambda_1 (1/s), symmetric log (linear within \\pm%g)', symlog_c));
title(ax, '\lambda_1, all conditions', 'FontWeight', 'normal');
grid(ax, 'on'); box(ax, 'off');

lgd = legend(h_seed, 'Orientation', 'horizontal', 'NumColumns', n_seeds, 'FontSize', 8);
lgd.Layout.Tile = 'south';

if save_fig
    if ~isfolder(out_dir), mkdir(out_dir); end
    png = fullfile(out_dir, sprintf('small_network_seeds_K%d_%s.png', K, hist_style));
    exportgraphics(fig, png, 'Resolution', 150);
    fprintf('\nsaved %s\n', png);
end

%% Local functions
function g = regime_trial(preset_name, cname, seeds, T, fs, plot_deci, window)
% REGIME_TRIAL The regime of one seed over the window, from a trajectory-only run.
rng(0, 'twister');                   % as run_lyapunov_spectrum: a seed must mean one network
t0 = tic;
m = build_from_preset(preset_name, cname, 'verbose', 'near-none', 'rng_seeds', seeds, ...
    'fs', fs, 'T_range', [0 T], 'lya_method', 'none', 'plot_deci', plot_deci);
m.run();
pd = m.plot_data;
keep = pd.t >= window(1) & pd.t <= window(2);
r_all = [];
types = m.cell_type_names;
for q = 1:numel(types)
    r_all = [r_all; pd.r.(types{q})(:, keep)]; %#ok<AGROW>
end
pop = mean(r_all, 1);
dt = median(diff(pd.t));
x = pop - mean(pop);
nfft = numel(x);
P = abs(fft(x)).^2;
f = (0:nfft - 1) / (nfft * dt);
half = 2:floor(nfft / 2);             % skip DC
[~, i_pk] = max(P(half));
g = struct('seeds', seeds, 'rate', mean(r_all, 'all'), ...
    'floor', mean(r_all == 0, 'all'), 'ceiling', mean(r_all == 1, 'all'), ...
    'between', mean(r_all > 0 & r_all < 1, 'all'), ...
    'n_silent', sum(all(r_all == 0, 2)), 'n_saturated', sum(all(r_all == 1, 2)), ...
    'pop_sd', std(pop), 'f_peak', f(half(i_pk)), ...
    'peak_power_frac', P(half(i_pk)) / sum(P(half)), ...
    'neuron_sd', mean(std(r_all, 0, 2)), 'seconds', toc(t0));
end

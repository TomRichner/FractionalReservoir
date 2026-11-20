% EXPLORE_SMALL_NETWORK_SPREAD The settled Lyapunov spectrum of the 80-neuron network for three conditions, with and without the per-neuron SFA and depression spreads.
%
% WHY. The settled window (iteration from 60 s, warm-up 20 s, exponents
% accumulated over 80 to 140 s; explore_small_network_settling) is the
% better estimate of the spectrum. This script takes it for no adaptation,
% one timescale (sfa1_std1) and chained multiple timescales (sfa3_stdc2),
% each with the per-neuron spreads as the preset has them (tau_a_spread
% [0.25 0.25], std_spread [0.25 0.25]) and with both spreads off. Without
% the spreads every neuron of a type has the same time constants, so the
% spectrum can be read against the bare rates 1/tau.
%
% THE NETWORK: small_network_model (n = 80, alpha = 0.2, F = 0.18634),
% seed [7 8]. The setpoint spread sigma_S_c stays as the preset has it in
% every run; only tau_a_spread and std_spread are switched off.
%
% THE RUN: T_range [0 140] s, 500 Hz, sra1, topk with lya_dt 0.05 s,
% window [80 140] s, warm-up 20 s, K = min(K_max, state dimension) passed
% explicitly as lya_K (no_adaptation has 80 states, so K = 80 = all of
% them). Every run has the same T_range, so all six share x0 and the noise.
% The six runs are rows of a settings table and run in a parfor over a
% pool of 3 opened with pool_set and closed with pool_release. Each worker
% does the whole run and its diagnostics and returns only a summary, so the
% 140 s trajectories (about 360 MB each at 640 states) stay on the workers.
%
% WHAT IT PRINTS, per run: state dimension, K, lambda_1, range, number
% positive, Kaplan-Yorke dimension, wall time (total and Lyapunov part);
% the realized time constants (min / median / max over the neurons, as tau
% and as 1/tau); W, the setpoints and x0 against the other run of the same
% condition; the within-run convergence (running exponents at the middle
% of the window against the full window); the population means of the
% rate, the SFA term and the depression states b1, b2 over settle_windows;
% the share of neuron-time over the window at the floor (r = 0), at the
% ceiling (r = 1) and between, and neurons silent or saturated throughout;
% the largest gaps; and for the no-spread runs the shelves (runs of
% consecutive sorted exponents whose steps are all below shelf_tol, at
% least shelf_min long) beside the bare rates. Per condition: the rank-by-
% rank difference without minus with spread. All exponents of the
% no-spread adapting runs are printed in full.
%
% THE FIGURE (one, 1250 x 700, left open, PNG to
% figs/explorations/small_network/): three columns (conditions), two rows.
% Top: sorted exponents against rank, with and without spread. Bottom:
% histograms of the two on shared bins, as outlines. Each column has its
% own exponent range, shared by its two rows. With show_rates, the bare
% rates -1/tau of the no-spread run's time constants are drawn on the top
% panel as thin dotted lines (only those inside the panel's range).
%
% PATH: run setup_paths once in the session first, with the current folder
% at the repository root. Takes about three minutes on a laptop.
%
% See also: explore_small_network_settling, small_network_model, lyapunov_topk

%% Settings
runs = cell2struct({ ...
    'no_adaptation spread',    'no_adaptation', true; ...
    'no_adaptation no spread', 'no_adaptation', false; ...
    'sfa1_std1 spread',        'sfa1_std1',     true; ...
    'sfa1_std1 no spread',     'sfa1_std1',     false; ...
    'sfa3_stdc2 spread',       'sfa3_stdc2',    true; ...
    'sfa3_stdc2 no spread',    'sfa3_stdc2',    false}, ...
    {'label', 'condition', 'spread'}, 2);
common = struct( ...
    'seed',           [7 8], ...
    'T_range',        [0 140], ...       % s
    'fs',             500, ...           % Hz
    'lya_dt',         0.05, ...          % s
    'lya_T_interval', [80 140], ...      % s
    'lya_warmup',     20, ...            % s
    'K_max',          200, ...           % K = min(K_max, state dimension)
    'plot_deci',      5, ...             % stored output at 100 Hz
    'settle_windows', {{[60 80], [120 140]}}, ...  % s, for the population means
    'no_spread_args', {{'tau_a_spread', [0 0], 'std_spread', [0 0]}});
n_top      = 60;        % top band against the rest, for the pair differences
n_gaps     = 6;
shelf_tol  = [];        % 1/s; [] = a quarter of the run's median step
shelf_min  = 3;         % shortest run of exponents that counts as a shelf
show_rates = true;      % bare rates -1/tau on the no-spread top panels
n_pool     = 3;
save_fig   = true;
out_dir    = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', 'small_network');

%% Pool (the project's pool module; never parpool directly)
max_workers = 0;
try
    pool_set(n_pool, 'label', 'small_network');
    p = gcp('nocreate');
    if isempty(p)
        fprintf('[pool] pool_set returned but no pool is up; running serially\n');
    else
        fprintf('[pool] pool up with %d workers\n', p.NumWorkers);
        max_workers = p.NumWorkers;
    end
catch err
    fprintf('[pool] pool_set FAILED: %s\n  running serially\n', err.message);
end

%% Run
n_runs = numel(runs);
out = cell(1, n_runs);
parfor (k = 1:n_runs, max_workers)
    out{k} = run_one(runs(k), common);
end
pool_release();
if isempty(gcp('nocreate'))
    fprintf('[pool] released; no pool is up\n');
else
    fprintf('[pool] WARNING: a pool is still up after pool_release\n');
end

%% Per run
for k = 1:n_runs
    o = out{k};
    lam = o.lambda;
    fprintf('\n=== %s ===\n', runs(k).label);
    fprintf('  state dimension %d, K = %d, wall %.1f s (Lyapunov %.1f s)\n', ...
        o.N, numel(lam), o.wall_s, o.lya_s);
    fprintf('  lambda_1 %+.4f, range %+.4f to %+.4f, positive %d, D_KY %s\n', ...
        lam(1), lam(1), lam(end), nnz(lam > 0), dky_text(o));
    fprintf('  time constants (s): min / median / max over the neurons   [1/tau]\n');
    for j = 1:numel(o.tau)
        t = o.tau(j).values;
        fprintf('    %-22s %.4g / %.4g / %.4g   [%.4f / %.4f / %.4f]\n', o.tau(j).name, ...
            min(t), median(t), max(t), 1 / max(t), 1 / median(t), 1 / min(t));
    end
    fprintf('  convergence: running exponents at t = %.2f s against %.2f s: lambda_1 change %+.4f, max |change| by sorted rank %.4f\n', ...
        o.conv.t_mid, o.conv.t_end, o.conv.d_lambda1, o.conv.max_abs_sorted);
    fprintf('  population means      %10s %10s %10s %10s\n', 'rate r', 'SFA term', 'b1', 'b2');
    for w = 1:numel(o.settle)
        s = o.settle(w);
        fprintf('    %-19s %10.4f %10.4f %10.4f %10.4f\n', sprintf('%g-%g s', s.window), ...
            s.rate, s.sfa, s.b1, s.b2);
    end
    fprintf('  neuron-time %g-%g s: floor %.3f, ceiling %.3f, between %.3f; silent throughout %d, saturated %d of %d\n', ...
        common.lya_T_interval, o.occ.floor, o.occ.ceiling, o.occ.between, o.occ.n_silent, o.occ.n_saturated, o.n);
    steps = lam(1:end-1) - lam(2:end);
    [step_sorted, step_rank] = sort(steps, 'descend');
    fprintf('  largest steps (rank i: between i and i+1); median step %.5f:\n', median(steps));
    for g = 1:min(n_gaps, numel(steps))
        i = step_rank(g);
        fprintf('    rank %3d: %+.4f -> %+.4f, step %.4f\n', i, lam(i), lam(i + 1), step_sorted(g));
    end
    if ~runs(k).spread
        tol = shelf_tol;
        if isempty(tol), tol = median(steps) / 4; end
        sh = find_shelves(lam, tol, shelf_min);
        fprintf('  shelves (every step < %.5f 1/s, at least %d exponents): %d\n', tol, shelf_min, numel(sh));
        for j = 1:numel(sh)
            fprintf('    ranks %3d-%3d (%3d): level %+.4f (%+.4f to %+.4f)\n', sh(j).first, sh(j).last, ...
                sh(j).last - sh(j).first + 1, mean(lam(sh(j).first:sh(j).last)), lam(sh(j).first), lam(sh(j).last));
        end
        fprintf('  bare rates -1/tau (1/s):');
        for j = 1:numel(o.tau)
            fprintf('  %s %+.4f', o.tau(j).name, -1 / median(o.tau(j).values));
        end
        fprintf('\n');
        if ~strcmp(runs(k).condition, 'no_adaptation')
            fprintf('  all %d exponents (1/s), by rank:\n', numel(lam));
            for i0 = 1:10:numel(lam)
                i1 = min(i0 + 9, numel(lam));
                fprintf('   %3d-%3d:', i0, i1);
                fprintf(' %+.4f', lam(i0:i1));
                fprintf('\n');
            end
        end
    end
end

%% Per condition: without spread minus with spread
conditions = unique({runs.condition}, 'stable');
fprintf('\n=== without minus with spread, rank by rank (1/s) ===\n');
for c = 1:numel(conditions)
    i_with = find(strcmp({runs.condition}, conditions{c}) & [runs.spread], 1);
    i_wo   = find(strcmp({runs.condition}, conditions{c}) & ~[runs.spread], 1);
    if isempty(i_with) || isempty(i_wo), continue; end
    a = out{i_with}; b = out{i_wo};
    d = b.lambda - a.lambda;
    nt = min(n_top, numel(d) - 1);
    fprintf(['  %-14s median %+.4f, min %+.4f (rank %d), max %+.4f (rank %d); ' ...
        'median ranks 1-%d %+.4f, %d-%d %+.4f;  identical spectra %d\n'], conditions{c}, ...
        median(d), min(d), find(d == min(d), 1), max(d), find(d == max(d), 1), ...
        nt, median(d(1:nt)), nt + 1, numel(d), median(d(nt+1:end)), isequal(a.lambda, b.lambda));
    fprintf('  %-14s W identical %d, setpoints identical %d, x0 identical %d, rates identical %d\n', '', ...
        isequal(a.W_sum, b.W_sum) && isequal(a.W_nnz, b.W_nnz), isequal(a.S_c, b.S_c), ...
        isequal(a.S0, b.S0), isequal(a.r_check, b.r_check));
end

%% Figure
fig = figure('Name', 'Small network (n = 80): settled spectra, with and without spread', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [40 60 1250 700]);
tl = tiledlayout(fig, 2, numel(conditions), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('n = 80, seed [%d %d], window %g-%g s (warm-up %g s), K = min(%d, state dimension)', ...
    common.seed, common.lya_T_interval, common.lya_warmup, common.K_max));
col_with = [0.85 0.33 0.10];
col_wo   = [0.00 0.45 0.74];
for c = 1:numel(conditions)
    i_with = find(strcmp({runs.condition}, conditions{c}) & [runs.spread], 1);
    i_wo   = find(strcmp({runs.condition}, conditions{c}) & ~[runs.spread], 1);
    a = out{i_with}; b = out{i_wo};
    both = [a.lambda; b.lambda];
    lo = min(both); hi = max([both; 0]);
    pad = 0.04 * (hi - lo);
    lim = [lo - pad, hi + pad];

    ax = nexttile(tl, c);
    hold(ax, 'on');
    plot(ax, 1:numel(a.lambda), a.lambda, 'o', 'MarkerSize', 3, 'Color', col_with, ...
        'MarkerFaceColor', col_with, 'DisplayName', sprintf('spread, \\lambda_1 = %+.3f', a.lambda(1)));
    plot(ax, 1:numel(b.lambda), b.lambda, 's', 'MarkerSize', 3, 'Color', col_wo, ...
        'MarkerFaceColor', col_wo, 'DisplayName', sprintf('no spread, \\lambda_1 = %+.3f', b.lambda(1)));
    yline(ax, 0, 'k--', 'HandleVisibility', 'off');
    if show_rates
        % One line per distinct bare rate; constants that share a value
        % share a line and a label (in sfa1_std1, tau_a1 = tau_drive).
        rates = arrayfun(@(t) -1 / median(t.values), b.tau);
        [levels, ~, which_level] = unique(round(rates, 9));
        for j = 1:numel(levels)
            y = levels(j);
            if y > lim(1) && y < lim(2)
                names = strjoin(strcat('1/', {b.tau(which_level == j).name}), ' = ');
                yline(ax, y, ':', ['bare ' names], 'Color', [0.3 0.3 0.3], ...
                    'LineWidth', 0.75, 'FontSize', 7, 'Interpreter', 'none', ...
                    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', ...
                    'HandleVisibility', 'off');
            end
        end
    end
    hold(ax, 'off');
    ylim(ax, lim);
    xlim(ax, [0, max(numel(a.lambda), numel(b.lambda)) + 1]);
    xlabel(ax, 'Rank');
    ylabel(ax, 'Lyapunov exponent (1/s)');
    title(ax, strrep(conditions{c}, '_', '\_'));
    legend(ax, 'Location', 'southwest', 'FontSize', 8);
    grid(ax, 'on');
    box(ax, 'off');

    ax = nexttile(tl, numel(conditions) + c);
    range_both = max(both) - min(both);
    bw = max(0.005, round(range_both * 8 / numel(a.lambda) / 0.005) * 0.005);
    edges = (floor(min(both) / bw) * bw):bw:(ceil(max(both) / bw) * bw + bw);
    hold(ax, 'on');
    histogram(ax, a.lambda, edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.3, ...
        'EdgeColor', col_with, 'DisplayName', 'spread');
    histogram(ax, b.lambda, edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.3, ...
        'EdgeColor', col_wo, 'DisplayName', 'no spread');
    xline(ax, 0, 'k--', 'HandleVisibility', 'off');
    hold(ax, 'off');
    xlim(ax, lim);
    xlabel(ax, 'Lyapunov exponent (1/s)');
    ylabel(ax, 'Count');
    title(ax, sprintf('bins %.3g 1/s', bw));
    legend(ax, 'Location', 'northwest', 'FontSize', 8);
    box(ax, 'off');
end
if save_fig
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    png = fullfile(out_dir, sprintf('small_network_spread_settled_K%d.png', common.K_max));
    exportgraphics(fig, png, 'Resolution', 150);
    fprintf('\nsaved %s\n', png);
end

%% Local functions
function o = run_one(row, common)
% RUN_ONE Build and run one row, and return its spectrum and diagnostics (no trajectory).
extra = {'plot_deci', common.plot_deci};
if ~row.spread
    extra = [extra, common.no_spread_args];
end
t_start = tic;
probe = small_network_model(row.condition, 'seed', common.seed, 'T_range', [0 1], ...
    'fs', common.fs, 'model_args', extra);
K = min(common.K_max, probe.N_sys_eqs);
clear probe
m = small_network_model(row.condition, 'seed', common.seed, 'T_range', common.T_range, ...
    'fs', common.fs, 'lya_method', 'topk', 'lya_K', K, 'lya_dt', common.lya_dt, ...
    'lya_T_interval', common.lya_T_interval, 'lya_warmup', common.lya_warmup, ...
    'model_args', extra);
m.run();
o.wall_s = toc(t_start);
res = m.lya_results;
o.lya_s = res.seconds;
o.lambda = sort(res.LE_spectrum(:), 'descend');
o.D_KY = res.D_KY;
o.D_KY_resolved = res.D_KY_resolved;
o.N = m.N_sys_eqs;
o.n = m.n;
o.W_sum = full(sum(m.W(:)));
o.W_nnz = nnz(m.W);
o.S_c = m.S_c_vec;
o.S0 = m.S0;

% Convergence: running exponents at the middle of the window against the end.
t_mid = mean(common.lya_T_interval);
row_mid = find(res.t_lya >= t_mid, 1);
row_end = find(all(isfinite(res.finite_LE_spectrum_t), 2), 1, 'last');
f_mid = sort(res.finite_LE_spectrum_t(row_mid, :), 2, 'descend')';
f_end = sort(res.finite_LE_spectrum_t(row_end, :), 2, 'descend')';
o.conv = struct('t_mid', res.t_lya(row_mid), 't_end', res.t_lya(row_end), ...
    'd_lambda1', f_end(1) - f_mid(1), 'max_abs_sorted', max(abs(f_end - f_mid)));

% The realized time constants.
p = m.get_params();
types = m.cell_type_names;
tau = struct('name', {}, 'values', {});
n_rungs = max(cellfun(@(x) size(x, 2), p.tau_a_matrix));
for j = 1:n_rungs
    v = [];
    for q = 1:numel(types)
        if size(p.tau_a_matrix{q}, 2) >= j
            v = [v; p.tau_a_matrix{q}(:, j)]; %#ok<AGROW>
        end
    end
    tau(end + 1) = struct('name', sprintf('tau_a%d', j), 'values', v); %#ok<AGROW>
end
fields = {'tau_rec1', 'tau_rel', 'tau_rec2', 'tau_rest', 'tau_drive'};
for f = 1:numel(fields)
    v = [];
    for pre = 1:numel(types)
        for post = 1:numel(types)
            P = p.route{pre, post}.P;
            if isfield(P, fields{f}) && ~isempty(P.(fields{f}))
                v = [v; P.(fields{f})]; %#ok<AGROW>
                break;   % a neuron's routes share their draws: one route per presynaptic type
            end
        end
    end
    for j = 1:size(v, 2)
        name = fields{f};
        if size(v, 2) > 1, name = sprintf('%s%d', name, j); end
        tau(end + 1) = struct('name', name, 'values', v(:, j)); %#ok<AGROW>
    end
end
tau(end + 1) = struct('name', 'tau_d', 'values', p.tau_d(:));
o.tau = tau;

% Settling and occupancy, from the stored output.
pd = m.plot_data;
c_eff = m.c ./ max(1, m.n_a);
windows = [common.settle_windows, {common.lya_T_interval}];
settle = struct('window', {}, 'rate', {}, 'sfa', {}, 'b1', {}, 'b2', {});
for w = 1:numel(windows)
    keep = pd.t >= windows{w}(1) & pd.t <= windows{w}(2);
    [r_all, sfa_all] = deal([]);
    for q = 1:numel(types)
        r_all = [r_all; pd.r.(types{q})(:, keep)]; %#ok<AGROW>
        a_q = pd.a.(types{q});
        if ~isempty(a_q) && size(a_q, 2) > 0
            sfa_all = [sfa_all; c_eff(q) * reshape(sum(a_q(:, :, keep), 2), size(a_q, 1), [])]; %#ok<AGROW>
        end
    end
    b_comp = {[], []};
    for pre = 1:numel(types)
        for post = 1:numel(types)
            if isfield(pd.b, types{pre}) && isfield(pd.b.(types{pre}), types{post})
                b_route = pd.b.(types{pre}).(types{post});
                for j = 1:min(2, size(b_route, 2))
                    b_comp{j} = [b_comp{j}; reshape(b_route(:, j, keep), size(b_route, 1), [])];
                end
            end
        end
    end
    if w <= numel(common.settle_windows)
        settle(end + 1) = struct('window', windows{w}, 'rate', mean(r_all, 'all'), ...
            'sfa', mean_or_nan(sfa_all), 'b1', mean_or_nan(b_comp{1}), 'b2', mean_or_nan(b_comp{2})); %#ok<AGROW>
    else
        o.occ = struct('floor', mean(r_all == 0, 'all'), 'ceiling', mean(r_all == 1, 'all'), ...
            'between', mean(r_all > 0 & r_all < 1, 'all'), ...
            'n_silent', sum(all(r_all == 0, 2)), 'n_saturated', sum(all(r_all == 1, 2)));
    end
end
o.settle = settle;
o.r_check = pd.r.(types{1})(:, end);   % last stored E rates, to compare runs
end

function sh = find_shelves(lam, tol, min_len)
% FIND_SHELVES Runs of consecutive sorted exponents whose every step is below tol.
small = (lam(1:end-1) - lam(2:end)) < tol;
sh = struct('first', {}, 'last', {});
i = 1;
while i <= numel(small)
    if small(i)
        j = i;
        while j <= numel(small) && small(j)
            j = j + 1;
        end
        if j - i + 1 >= min_len
            sh(end + 1) = struct('first', i, 'last', j); %#ok<AGROW>
        end
        i = j;
    else
        i = i + 1;
    end
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

function s = dky_text(o)
% DKY_TEXT The Kaplan-Yorke dimension, or 'unresolved'.
if o.D_KY_resolved
    s = sprintf('%.3f', o.D_KY);
else
    s = 'unresolved within K';
end
end

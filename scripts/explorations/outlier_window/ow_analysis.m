% OW_ANALYSIS What 20, 50 and 100 s windows would have shown, emulated on the 200 s runs, and what each window costs.
%
% Reads data/explorations/outlier_window/ow_series.mat (ow_extract). Writes
%   data/explorations/outlier_window/ow_truncation_table.md   item 2's table
%   data/explorations/outlier_window/ow_prefix_table.md       item 4: the
%       window-wide quantities of the paper figures on 20 s against 200 s
%   data/explorations/outlier_window/ow_cost_table.md         item 3: cost
%       split and estimated stage wall times (list scheduling on 12 workers)
%   figs/explorations/outlier_window/ow_window_comparison.png the figure
%
% TRUNCATION. A run's series is cut to its first L s (round(L / dt)
% samples) BEFORE the excursions are found (ospec_excursions), so an
% excursion running past L is censored at L, as a run of length L would
% have it. The 200 s runs are one noise realization of 260 s; a real run of
% 60 + L s would be another realization, with the same statistics.
% Per seed the survival is censored_survival's Kaplan-Meier estimate; the
% median curve across seeds follows fig_outlier_dwell's rule (drawn while
% at least half of the seeds with excursions still contribute).

setup_paths();
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dd = fullfile(root, 'data', 'explorations', 'outlier_window');
fd = fullfile(root, 'figs', 'explorations', 'outlier_window');
if ~isfolder(fd); mkdir(fd); end
S0 = load(fullfile(dd, 'ow_series.mat'));
X = S0.X;
nets = {'reference', 'EEx1p5'}; net_txt = {'Reference', 'E-to-E \times 1.5'};
conds = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
st = manuscript_style();
cond_txt = grouped_condition_label(conds, 'sentence');
Ls = [20 50 100 200];
ser = {'frozen', 'full'}; ser_txt = {'frozen gains', 'full Jacobian'};
xnet = {X.network}; xcond = {X.condition};

%% Item 2: truncated excursions
Q = cell(numel(nets), numel(conds), numel(Ls), 2);
for w = 1:numel(nets)
    for c = 1:numel(conds)
        sel = find(strcmp(xnet, nets{w}) & strcmp(xcond, conds{c}));
        for j = 1:numel(Ls)
            for k = 1:2
                E = cell(1, numel(sel));
                for s = 1:numel(sel)
                    x = X(sel(s));
                    if k == 1; t = x.t; y = x.rate_top; dt = x.dt;
                    else; t = x.full_t; y = x.full_re; dt = x.full_dt; end
                    m = min(numel(t), round(Ls(j) / dt));
                    E{s} = ospec_excursions(t(1:m), y(1:m));
                end
                Q{w, c, j, k} = dwell_summary([E{:}]);
            end
        end
    end
end

fid = fopen(fullfile(dd, 'ow_truncation_table.md'), 'w');
fprintf(fid, ['# Excursions beyond zero on the first L s of the 200 s runs (fastmedium, 10 seeds)\n\n' ...
    '`exc/seed`: excursions per seed, median [min, max]. `cut`: share of all excursions censored by the window. ' ...
    '`q25 / median / q75` (s): Kaplan-Meier quartiles on the seeds pooled; `> x`: not reached within the longest ' ...
    'observed x. `beyond`: fraction of time beyond zero, median [min, max] over seeds. `reach` (s): how far the ' ...
    'median survival curve is drawn (fig_outlier_dwell''s rule), and `S at reach` its value there (0: the curve ' ...
    'reaches the bottom; > 0: ends in a triangle). `seeds exc`: seeds with at least one excursion. Steps: frozen ' ...
    '0.05 s, full 0.1 s.\n\n']);
fprintf(fid, '| network | condition | series | L (s) | seeds exc | exc/seed | cut | q25 | median | q75 | beyond | reach | S at reach |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for w = 1:numel(nets)
    for c = 1:numel(conds)
        for k = 1:2
            for j = 1:numel(Ls)
                q = Q{w, c, j, k};
                qt = arrayfun(@(v) qtxt(v, q.pooled_t_max, q.n_exc), q.pooled_q, 'UniformOutput', false);
                reach = '-'; send = '-';
                if ~isempty(q.x); reach = sprintf('%.3g', q.x(end)); send = sprintf('%.2f', q.med(end)); end
                cut = '-'; if q.n_exc > 0; cut = sprintf('%d%% (%d/%d)', round(100 * q.n_cut / q.n_exc), q.n_cut, q.n_exc); end
                fprintf(fid, '| %s | %s | %s | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s |\n', nets{w}, conds{c}, ...
                    ser_txt{k}, Ls(j), q.n_seeds_exc, mmm(q.n_per_seed, '%g'), cut, qt{:}, mmm(q.fraction, '%.3f'), reach, send);
            end
        end
    end
end
fclose(fid);

%% The comparison figure: full-Jacobian median survival, four window lengths
cols = interp1([0 1], [0.72 0.80 0.92; 0.03 0.19 0.42], linspace(0, 1, numel(Ls)));   % light to dark blue
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 st.canvas_w 6.2], 'Visible', 'off');
tl = tiledlayout(fig, numel(nets), numel(conds), 'TileSpacing', 'compact', 'Padding', 'compact');
h = gobjects(1, numel(Ls));
for w = 1:numel(nets)
    for c = 1:numel(conds)
        ax = nexttile(tl); hold(ax, 'on');
        for j = 1:numel(Ls)
            q = Q{w, c, j, 2};
            if isempty(q.x); continue; end
            [sx, sy] = stair_points(q.x, q.med);
            h(j) = plot(ax, sx, sy, '-', 'Color', cols(j, :), 'LineWidth', 2);
            if q.med(end) > 0
                plot(ax, q.x(end), q.med(end), '>', 'MarkerSize', 7, 'MarkerFaceColor', cols(j, :), 'MarkerEdgeColor', cols(j, :));
            end
        end
        set(ax, 'XScale', 'log', 'XLim', [0.08 320], 'YLim', [0 1.04], 'YTick', 0:0.5:1, 'XTick', 10 .^ (-1:2), ...
            'XMinorTick', 'off', 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        if w == 1; title(ax, cond_txt{c}, 'FontSize', st.title_fs, 'FontWeight', 'normal'); end
        if c == 1; ylabel(ax, {net_txt{w}, 'Fraction longer'}, 'FontSize', st.label_fs); end
        if w == numel(nets) && c == 2; xlabel(ax, 'Time beyond zero (s)', 'FontSize', st.label_fs); end
        if all(cellfun(@(q) q.n_exc == 0, Q(w, c, :, 2)))
            text(ax, 1, 0.5, 'No excursion', 'HorizontalAlignment', 'center', 'FontSize', st.tick_fs, 'Color', [0.4 0.4 0.4]);
        elseif Q{w, c, end, 2}.n_exc <= 3
            q = Q{w, c, end, 2};
            none_at = Ls(cellfun(@(z) z.n_exc == 0, squeeze(Q(w, c, :, 2)))');
            txt = {sprintf('%d excursions of %.2g s', q.n_exc, max(q.pooled_t_max)), sprintf('in %d seed at 200 s,', q.n_seeds_exc)};
            if ~isempty(none_at); txt{end + 1} = ['none at ' strjoin(arrayfun(@(v) sprintf('%d', v), none_at, 'UniformOutput', false), ', ') ' s']; end %#ok<AGROW>
            text(ax, 8, 0.5, txt, 'HorizontalAlignment', 'center', 'FontSize', st.tick_fs, 'Color', [0.4 0.4 0.4]);
        end
    end
end
ok = isgraphics(h);
lg = legend(h(ok), arrayfun(@(L) sprintf('%d s', L), Ls(ok), 'UniformOutput', false), 'Orientation', 'horizontal', ...
    'FontSize', st.tick_fs, 'Box', 'off');
lg.Layout.Tile = 'north';
exportgraphics(fig, fullfile(fd, 'ow_window_comparison.png'), 'Resolution', 200, 'BackgroundColor', 'white');
close(fig);

%% Item 4: window-wide quantities, first 20 s against all 200 s
fid = fopen(fullfile(dd, 'ow_prefix_table.md'), 'w');
qn = {'ratio_med', 'ratio_p05', 'ratio_p95', 'frac_on_E', 'gain_ratio_EI', 'frac_unstable'};
fprintf(fid, ['# Window-wide quantities of the paper figures: first 20 s against the whole 200 s (fastmedium, 10 seeds)\n\n' ...
    'Per network and condition: the median over seeds at 20 s and at 200 s, and the median and largest absolute ' ...
    'per-seed change (200 s minus 20 s) relative to the seeds'' spread at 200 s (IQR). `ratio_*`: detachment ' ...
    'Re lambda_1 / R_eff (median, 5th, 95th percentile over the window; fig_outlier_dose left, fig_outlier_plane); ' ...
    '`frac_on_E`: median fraction of E neurons on the slope (dose middle); `gain_ratio_EI` (dose right); ' ...
    '`frac_unstable`: fraction of time the frozen-gain eigenvalue is beyond zero (dose table).\n\n']);
fprintf(fid, '| network | condition | quantity | 20 s | 200 s | median change | largest change | IQR across seeds at 200 s |\n|---|---|---|---|---|---|---|---|\n');
for w = 1:numel(nets)
    for c = 1:numel(conds)
        sel = find(strcmp(xnet, nets{w}) & strcmp(xcond, conds{c}));
        for i = 1:numel(qn)
            a = arrayfun(@(x) x.pre(1).(qn{i}), X(sel));
            b = arrayfun(@(x) x.pre(end).(qn{i}), X(sel));
            d = b - a;
            [~, im] = max(abs(d));
            fprintf(fid, '| %s | %s | %s | %.3g | %.3g | %+.3g | %+.3g | %.3g |\n', nets{w}, conds{c}, qn{i}, ...
                median(a), median(b), median(d), d(im), iqr_(b));
        end
    end
end
fclose(fid);

%% Item 3: cost split, and wall times by list scheduling on 12 workers
is_long = ismember(xnet, nets);
Sx = [X.seconds];
t_idx = [X.seconds_index];
comp = [[Sx.sim]; [Sx.fast]; [Sx.null]; [Sx.full]; [Sx.dense]];
other = t_idx - sum(comp, 1);                     % load, save, setup
fid = fopen(fullfile(dd, 'ow_cost_table.md'), 'w');
fprintf(fid, '# Cost of the outlier_spectra stage, fastmedium (measured from R.seconds and the index)\n\n');
fprintf(fid, '| runs | condition | n | worker-min per run, median [min, max] | sim | frozen pass | null | full pass | dense | other |\n|---|---|---|---|---|---|---|---|---|---|\n');
for g = [true false]
    for c = 1:numel(conds)
        sel = is_long == g & strcmp(xcond, conds{c});
        fprintf(fid, '| %s | %s | %d | %s | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f |\n', ternary(g, 'long (200 s)', 'short (20 s)'), ...
            conds{c}, nnz(sel), mmm(t_idx(sel) / 60, '%.1f'), median(comp(:, sel) / 60, 2), median(other(sel) / 60));
    end
end
fprintf(fid, '\nSum of run times: long %.0f worker-min (%d runs), short %.0f worker-min (%d runs); share long %.1f%%.\n', ...
    sum(t_idx(is_long)) / 60, nnz(is_long), sum(t_idx(~is_long)) / 60, nnz(~is_long), 100 * sum(t_idx(is_long)) / sum(t_idx));
% Per-state costs measured in the runs
for c = 1:numel(conds)
    sel = find(is_long & strcmp(xcond, conds{c}));
    fs = [X(sel).full_sec];
    fz = [Sx(sel).fast] ./ arrayfun(@(x) numel(x.t), X(sel));
    fprintf(fid, '\n%s (N = %d state variables): frozen-gain state %.3f s (median over runs), full-Jacobian state %.3f s (median over states; mean %.3f s).', ...
        conds{c}, X(sel(1)).N, median(fz), median(fs), mean(fs));
end
fprintf(fid, '\n\n');
% Scenarios: the long networks' window L, frozen step ds, full step df
scen = struct('name', {'200 s, 0.05 / 0.1 (as run)', '100 s, 0.05 / 0.1', '50 s, 0.05 / 0.1', '20 s, 0.05 / 0.1', ...
    '200 s, 0.1 / 0.1', '200 s, 0.1 / 0.2', '20 s, 0.01 / 0.1 (fast mode''s steps)', 'no full networks (all 20 s, no full or dense pass)'}, ...
    'L', {200, 100, 50, 20, 200, 200, 20, 20}, 'ds', {0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.01, 0.01}, ...
    'df', {0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1}, 'nofull', {false, false, false, false, false, false, false, true});
wall_measured = schedule(t_idx, is_long, X, 200, 0.05);
fprintf(fid, ['Stage wall time: logged 387.8 min. List scheduling of the measured run times in the stage''s order ' ...
    '(run_outlier_spectra_stage: by its cost key, longest first; one run at a time per worker) gives %.0f min.\n\n'], wall_measured / 60);
fprintf(fid, ['Estimated wall times (12 workers). Each long run''s measured parts are rescaled: simulation by (60 + L)/260, ' ...
    'frozen pass by its number of states (L/ds), full pass by its number (L/df), null and dense passes unchanged; ' ...
    'short runs as measured. EXTRAPOLATED, except the first row.\n\n']);
fprintf(fid, '| long networks: window, frozen step / full step | worker-min, long runs | wall (min), 12 workers | longest run (min) |\n|---|---|---|---|\n');
for i = 1:numel(scen)
    sc = scen(i);
    tl_ = t_idx;
    li = find(is_long);
    for r = li
        x = X(r); s = x.seconds;
        sim = s.sim * (60 + sc.L) / 260;
        fz = s.fast / numel(x.t) * sc.L / sc.ds;
        fl = s.full / numel(x.full_t) * sc.L / sc.df;
        dn = s.dense; nl = s.null;
        if sc.nofull
            % as a short network: the short runs of the same condition
            sh = ~is_long & strcmp(xcond, x.condition);
            tl_(r) = median(t_idx(sh));
        else
            tl_(r) = sim + fz + nl + fl + dn + other(r);
        end
    end
    wl = schedule(tl_, is_long, X, sc.L, sc.ds);
    fprintf(fid, '| %s | %.0f | %.0f | %.0f |\n', sc.name, sum(tl_(is_long)) / 60, wl / 60, max(tl_) / 60);
end
fclose(fid);
fprintf('wrote the tables and %s\n', fullfile(fd, 'ow_window_comparison.png'));

%% ------------------------------------------------------------------------
function T = schedule(t_run, is_long, X, L, ds)
% Wall time of the runs on 12 workers, handed out one at a time in the
% stage's order: its cost key (window / step + 50 window for the full
% networks), descending, stable on the index order.
cost = zeros(1, numel(X));
for r = 1:numel(X)
    if is_long(r); cost(r) = L / ds + 50 * L; else; cost(r) = 20 / 0.01; end
end
[~, order] = sort(cost, 'descend');
free = zeros(1, 12);
for r = order
    [t0, w] = min(free);
    free(w) = t0 + t_run(r);
end
T = max(free);
end

function q = dwell_summary(E)
% fig_outlier_dwell's dwell_curve, copied (an exploration may not change the figure).
n_seeds = numel(E);
K = cell(1, n_seeds);
for s = 1:n_seeds
    K{s} = censored_survival(E(s).duration, E(s).censored);
end
has = cellfun(@(k) k.n > 0, K);
q = struct();
q.n_seeds = n_seeds; q.n_seeds_exc = nnz(has);
q.fraction = [E.fraction];
q.n_per_seed = cellfun(@(k) k.n, K);
q.n_exc = sum(q.n_per_seed);
q.n_cut = sum(arrayfun(@(e) nnz(e.censored), E));
pooled = censored_survival([E.duration], [E.censored]);
q.pooled_q = [pooled.quantile(0.25), pooled.quantile(0.5), pooled.quantile(0.75)];
q.pooled_t_max = pooled.t_max;
q.x = []; q.med = [];
if ~any(has); return; end
x = unique([E.duration]);
S = nan(q.n_seeds_exc, numel(x));
Kh = K(has);
for s = 1:q.n_seeds_exc
    S(s, :) = Kh{s}.at(x);
end
valid = sum(~isnan(S), 1) >= max(1, ceil(q.n_seeds_exc / 2));
last = find(valid, 1, 'last');
x = x(1:last); S = S(:, 1:last);
q.x = [max([E.dt]), x];
q.med = [1, median(S, 1, 'omitnan')];
end

function [sx, sy] = stair_points(x, y)
sx = reshape([x; x], 1, []); sx = sx(2:end);
sy = reshape([y; y], 1, []); sy = sy(1:end - 1);
end

function s = qtxt(v, t_max, n_exc)
if n_exc == 0; s = '-'; elseif isnan(v); s = sprintf('> %.3g', t_max); else; s = sprintf('%.3g', v); end
end

function s = mmm(v, fmt)
v = v(~isnan(v));
if isempty(v); s = '-';
elseif isscalar(v) || max(v) == min(v); s = sprintf(fmt, median(v));
else; s = sprintf([fmt ' [' fmt ', ' fmt ']'], median(v), min(v), max(v));
end
end

function v = iqr_(x)
s = sort(x(~isnan(x))); n = numel(s);
if n < 2; v = 0; return; end
p = @(pp) interp1(((1:n) - 0.5) / n, s, min(max(pp, 0.5 / n), 1 - 0.5 / n));
v = p(0.75) - p(0.25);
end

function out = ternary(c, a, b)
if c; out = a; else; out = b; end
end

function out = fig_outlier_dwell(cfg)
% FIG_OUTLIER_DWELL How long the rightmost eigenvalue stays beyond zero: censored survival curves, median across seeds.
%
%   out = FIG_OUTLIER_DWELL('run_dir', d)
%   out = FIG_OUTLIER_DWELL('run_dir', d, 'band', 'range')
%
% Draws what the opt-in stage outlier_spectra saved (FD-040), through
% outlier_fig_data. Reads one run directory, or the index file given as
% data_file; errors when neither holds the stage; simulates nothing.
% FD-039 (Tom, 2026-11-20: "could become a good figure. It would need more
% seeds and then we could plot the median for each. And we would need to
% run it for much longer"). Supplementary candidate.
%
% ONE PANEL PER NETWORK (the reference and E-to-E x 1.5). An EXCURSION is
% a run of consecutive sampled states at which the rightmost eigenvalue
% has positive real part (ospec_excursions); the curves are the fraction
% of excursions that last longer than the time on the horizontal axis
% (log). Per condition, in its color:
%   dashed, thin    frozen gains: the rightmost eigenvalue of the
%                   frozen-gain effective connectivity (slow variables held)
%   solid, thick    the full Jacobian (slow variables free): its rightmost
%                   network mode; an eigenvalue beyond zero is always a
%                   network mode, so the series is safe for this (FD-040)
% Each curve is the MEDIAN ACROSS SEEDS of the per-seed survival estimate,
% with the spread as a translucent band (cfg.band: 'iqr', the 25th to 75th
% percentile across seeds, or 'range'; with two seeds they are the same).
%
% CUT-OFF EXCURSIONS ARE CENSORED, NOT COUNTED AS COMPLETE. An excursion
% that contains the first or the last sample of the window lasted at least
% its observed length. Per seed the survival is the Kaplan-Meier estimate
% (censored_survival): it drops only where a COMPLETE excursion ended, and
% a cut one counts as still going at its observed length. Chosen because
% the two simpler treatments are both wrong in the direction that matters:
% counting cut excursions as complete draws a condition that is beyond
% zero throughout as a step at the window's length, and dropping them
% removes the longest excursions and shortens every curve. Beyond a seed's
% longest observed excursion, when that one was cut, its estimate says
% nothing and the seed leaves the median; the median is drawn while at
% least half of the seeds that have excursions still contribute.
%
% WHERE A CURVE ENDS ABOVE ZERO it ends in a right-pointing triangle: the
% excursions still going there were cut off by the window, and nothing is
% known beyond. A condition that is beyond zero for the whole window in
% every seed is therefore a level line at 1 ending in a triangle at the
% window's length, not a step down. A condition with no excursion in any
% seed has no curve (the table says so).
%
% Resolution: an excursion's length is its number of samples times the
% sampling step (the series' step for the frozen gains, the full pass's
% coarser step for the full Jacobian), so no curve begins left of its
% step, a shorter excursion can be missed and two separated by less are
% merged.
%
% Nothing on the panels but the network names and the axis names.
% Fig_Outlier_Dwell_table.md has a caption paragraph and, per network,
% condition and series: the fraction of time beyond zero, the number of
% excursions and how many were cut (at the start, at the end, spanning
% the whole window), and the quartiles of the dwell time, across seeds.
%
% See also: outlier_fig_data, censored_survival, ospec_excursions,
%           fig_outlier_plane, fig_outlier_slow_variables

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.networks                  = {'reference', 'EEx1p5'}
    cfg.conditions                = {}
    cfg.seed double               = []    % unused: every seed enters the medians
    cfg.band (1,:) char {mustBeMember(cfg.band, {'iqr', 'range', 'none'})} = 'iqr'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = outlier_fig_data(cfg, 'fig_outlier_dwell', 'dwell', 'networks', cfg.networks);
st = F.st;
[n_net, n_cond, ~] = size(F.D);
series = {'frozen', 'full'};
series_txt = {'frozen gains', 'full Jacobian'};
style = {'--', '-'}; width = [1.25, 2.25];

%% Per network, condition and series: the seeds' estimates and their median
C = cell(n_net, n_cond, 2);
dt_min = Inf; window_max = 0;
for w = 1:n_net
    for c = 1:n_cond
        runs = squeeze(F.D(w, c, :)); runs = runs(~cellfun(@isempty, runs));
        if isempty(runs); continue; end
        runs = [runs{:}];
        for k = 1:2
            E = {runs.(series{k})};
            E = E(~cellfun(@isempty, E));
            if isempty(E); continue; end
            C{w, c, k} = dwell_curve([E{:}], cfg.band);
            dt_min = min(dt_min, C{w, c, k}.dt);
            window_max = max(window_max, C{w, c, k}.window_s);
        end
    end
end
if ~isfinite(dt_min)
    error('fig_outlier_dwell:NoRuns', 'No finished run among the panels.');
end

%% Draw
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.6 * n_net + 0.6, 4.6], 'Visible', vis);
tl = tiledlayout(fig, 1, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');
x_lim = [dt_min / 1.25, 1.6 * window_max];
for w = 1:n_net
    ax = nexttile(tl); hold(ax, 'on');
    for pass = 1:2                                       % every band under every line
        for c = 1:n_cond
            col = F.colors{c};
            for k = 1:2
                q = C{w, c, k};
                if isempty(q) || isempty(q.x); continue; end
                if pass == 1
                    if ~strcmp(cfg.band, 'none') && q.n_seeds_exc > 1
                        [bx, blo] = stair_points(q.x, q.lo); [~, bhi] = stair_points(q.x, q.hi);
                        patch(ax, [bx, fliplr(bx)], [blo, fliplr(bhi)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                    end
                else
                    [sx, sy] = stair_points(q.x, q.med);
                    plot(ax, sx, sy, style{k}, 'Color', col, 'LineWidth', width(k));
                    if q.med(end) > 0                    % still going where the window cut it: nothing known beyond
                        plot(ax, q.x(end), q.med(end), '>', 'MarkerSize', 7, 'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                    end
                end
            end
        end
    end
    set(ax, 'XScale', 'log', 'XLim', x_lim, 'YLim', [0 1.04], 'YTick', 0:0.5:1, 'XTick', 10 .^ (-3:3), ...
        'XMinorTick', 'off', 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    title(ax, F.net_labels{w}, 'FontSize', st.title_fs, 'FontWeight', 'normal');
    xlabel(ax, 'Time beyond zero (s)', 'FontSize', st.label_fs);
    if w == 1; ylabel(ax, 'Fraction of excursions longer', 'FontSize', st.label_fs); end
end

fig_tag = 'Fig_Outlier_Dwell';
out = struct('figs', fig, 'files', {{}}, 'source', F.file, 'curves', {C});
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        write_table(fid, fig_tag, F, C, series_txt, cfg);
        clear closer
        out.files{end + 1} = table_name;
    end
end
for w = 1:n_net
    for c = 1:n_cond
        for k = 1:2
            q = C{w, c, k};
            if isempty(q); continue; end
            vprintf(cfg.verbose, 'verbose', '%s | %s | %s: %d excursions, %d cut, time beyond zero %s\n', F.networks{w}, ...
                F.cond_names{c}, series_txt{k}, q.n_exc, q.n_cut, mmm(q.fraction, '%.3f'));
        end
    end
end
end

%% ------------------------------------------------------------------------
function q = dwell_curve(E, band)
% From the excursions of every seed (a struct array, one ospec_excursions
% output per seed) to the median survival curve and the table's numbers.
n_seeds = numel(E);
K = cell(1, n_seeds);
for s = 1:n_seeds
    K{s} = censored_survival(E(s).duration, E(s).censored);
end
has = cellfun(@(k) k.n > 0, K);
q = struct();
q.n_seeds = n_seeds; q.n_seeds_exc = nnz(has);
q.dt = max([E.dt]); q.window_s = max([E.window_s]);
q.fraction = [E.fraction];
q.n_per_seed = cellfun(@(k) k.n, K);
q.n_exc = sum(q.n_per_seed);
q.n_cut_start = sum(arrayfun(@(e) nnz(e.cut_start & ~e.cut_end), E));
q.n_cut_end = sum(arrayfun(@(e) nnz(e.cut_end & ~e.cut_start), E));
q.n_cut_both = sum(arrayfun(@(e) nnz(e.cut_start & e.cut_end), E));
q.n_cut = q.n_cut_start + q.n_cut_end + q.n_cut_both;
q.seed_median = cellfun(@(k) k.quantile(0.5), K(has));
q.longest = max([E.duration, NaN]);
pooled = censored_survival([E.duration], [E.censored]);
q.pooled_q = [pooled.quantile(0.25), pooled.quantile(0.5), pooled.quantile(0.75)];
q.pooled_t_max = pooled.t_max;
q.x = []; q.med = []; q.lo = []; q.hi = [];
if ~any(has); return; end
% The grid: every observed length of any seed. Each seed's estimate there;
% NaN where it says nothing (beyond its longest excursion, when that was cut).
x = unique([E.duration]);
S = nan(q.n_seeds_exc, numel(x));
Kh = K(has);
for s = 1:q.n_seeds_exc
    S(s, :) = Kh{s}.at(x);
end
valid = sum(~isnan(S), 1) >= max(1, ceil(q.n_seeds_exc / 2));
last = find(valid, 1, 'last');
valid(:) = false; valid(1:last) = true;                 % one stretch, from the shortest length on
x = x(valid); S = S(:, valid);
q.x = [q.dt, x];                                        % no excursion is shorter than one sample
q.med = [1, median(S, 1, 'omitnan')];
if strcmp(band, 'range')
    q.lo = [1, min(S, [], 1)]; q.hi = [1, max(S, [], 1)];
else
    q.lo = [1, quantile_rows(S, 0.25)]; q.hi = [1, quantile_rows(S, 0.75)];
end
end

function v = quantile_rows(S, p)
% The p-quantile down each column, ignoring NaN, by linear interpolation
% between the order statistics at (i - 0.5)/n (the minimum and the maximum
% outside them), without the Statistics Toolbox.
v = nan(1, size(S, 2));
for j = 1:size(S, 2)
    s = sort(S(~isnan(S(:, j)), j));
    n = numel(s);
    if n == 0; continue; end
    pos = min(max(p * n + 0.5, 1), n);
    lo = floor(pos); hi = ceil(pos);
    v(j) = s(lo) + (pos - lo) * (s(hi) - s(lo));
end
end

function [sx, sy] = stair_points(x, y)
% The vertices of a right-continuous step function: level y(i) from x(i)
% to x(i + 1), then the drop at x(i + 1).
sx = reshape([x; x], 1, []); sx = sx(2:end);
sy = reshape([y; y], 1, []); sy = sy(1:end - 1);
end

function write_table(fid, fig_tag, F, C, series_txt, cfg)
[n_net, n_cond, ~] = size(C);
any_curve = C(~cellfun(@isempty, C));
windows = unique(cellfun(@(q) q.window_s, any_curve));
n_seed = size(F.seeds, 1);
switch cfg.band
    case 'iqr';   band_txt = 'Bands: the 25th to 75th percentile across seeds of the per-seed estimates.';
    case 'range'; band_txt = 'Bands: the range across seeds of the per-seed estimates.';
    otherwise;    band_txt = 'The spread across seeds is not drawn.';
end
whole = {};
for w = 1:n_net
    for c = 1:n_cond
        for k = 1:2
            q = C{w, c, k};
            if ~isempty(q) && q.n_exc > 0 && q.n_cut_both == q.n_exc && q.n_seeds_exc == q.n_seeds
                whole{end + 1} = sprintf('%s, %s, %s', F.networks{w}, F.cond_names{c}, series_txt{k}); %#ok<AGROW>
            end
        end
    end
end
if isempty(whole)
    whole_txt = 'No condition is beyond zero for the whole window in every seed.';
else
    whole_txt = ['BEYOND ZERO FOR THE WHOLE WINDOW IN EVERY SEED (a level line at 1 ending in a triangle at the window''s ' ...
        'length; its dwell time is longer than the window, and no more is known): ' strjoin(whole, '; ') '.'];
end
fprintf(fid, '# %s\n\n', fig_tag);
fprintf(fid, ['One panel per network (%s). An excursion is a run of consecutive sampled states at which the rightmost ' ...
    'eigenvalue has positive real part. The curves give the fraction of excursions that last longer than the time on ' ...
    'the horizontal axis (log). Colors are the conditions (%s). DASHED, thin: the rightmost eigenvalue of the ' ...
    'frozen-gain effective connectivity (slow variables held). SOLID, thick: the rightmost network mode of the full ' ...
    'Jacobian (slow variables free). Each curve is the median across seeds of the per-seed Kaplan-Meier estimate. %s ' ...
    'An excursion that contains the first or the last sample of the window was cut off by it and is CENSORED: it ' ...
    'counts as still going at its observed length and is not counted as having ended there. Where a seed''s longest ' ...
    'excursion was cut, its estimate ends there; the median is drawn while at least half of the seeds that have ' ...
    'excursions contribute. A curve that ends above zero ends in a right-pointing TRIANGLE: the excursions still going ' ...
    'there were cut off by the window. A condition without a curve has no excursion in any seed. %s\n\n'], ...
    strjoin(F.networks, ', '), strjoin(F.titles, ', '), band_txt, whole_txt);
fprintf(fid, ['Length of the series %s s (samples times step, after settling); %d seed(s): %s. Lengths are numbers of samples times the sampling ' ...
    'step (`step`), so the resolution is one step. Data: `%s`.\n\n'], strjoin(arrayfun(@(v) sprintf('%g', v), windows, ...
    'UniformOutput', false), ' and '), n_seed, ...
    strjoin(arrayfun(@(k) sprintf('[%d %d]', F.seeds(k, :)), 1:n_seed, 'UniformOutput', false), ', '), F.file);
fprintf(fid, ['`time beyond zero`: fraction of the sampled states with the eigenvalue beyond zero, median over seeds ' ...
    '[min, max]. `excursions`: the total over seeds (and per seed, median [min, max]). `cut`: of those, how many the ' ...
    'window cut off, as start / end / both (both: one excursion spanning the whole window). `dwell q25, median, q75` ' ...
    '(s): quartiles of the dwell time from the Kaplan-Meier estimate on the excursions of all seeds pooled, censoring ' ...
    'respected; `> x` means the estimate does not fall that far within the longest observed excursion x, so the ' ...
    'quartile is longer than that. `median per seed`: the same median from each seed''s own estimate, median over the ' ...
    'seeds where it is reached [min, max], and the number of seeds where it is. `longest`: the longest observed ' ...
    'length, cut or not.\n\n']);
fprintf(fid, ['| network | condition | series | step (s) | seeds | seeds with excursions | time beyond zero | excursions | ' ...
    'cut (start / end / both) | dwell q25 | dwell median | dwell q75 | median per seed | longest (s) |\n']);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for w = 1:n_net
    for c = 1:n_cond
        for k = 1:2
            q = C{w, c, k};
            if isempty(q); continue; end
            qt = arrayfun(@(v) quartile_txt(v, q.pooled_t_max, q.n_exc), q.pooled_q, 'UniformOutput', false);
            if isempty(q.seed_median) || all(isnan(q.seed_median))
                per_seed = '-';
            else
                per_seed = sprintf('%s (%d)', mmm(q.seed_median, '%.3g'), nnz(~isnan(q.seed_median)));
            end
            longest = '-'; if q.n_exc > 0; longest = sprintf('%.3g', q.longest); end
            fprintf(fid, '| %s | %s | %s | %g | %d | %d | %s | %d (%s) | %d (%d / %d / %d) | %s | %s | %s | %s | %s |\n', ...
                F.networks{w}, F.cond_names{c}, series_txt{k}, q.dt, q.n_seeds, q.n_seeds_exc, mmm(q.fraction, '%.3f'), ...
                q.n_exc, mmm(q.n_per_seed, '%g'), q.n_cut, q.n_cut_start, q.n_cut_end, q.n_cut_both, qt{:}, per_seed, longest);
        end
    end
end
no_full = {};
for w = 1:n_net
    for c = 1:n_cond
        if ~isempty(C{w, c, 1}) && isempty(C{w, c, 2})
            no_full{end + 1} = sprintf('%s, %s', F.networks{w}, F.cond_names{c}); %#ok<AGROW>
        end
    end
end
if ~isempty(no_full)
    fprintf(fid, '\nNo full-Jacobian series was computed for: %s (the stage computes it for its full networks only).\n', ...
        strjoin(no_full, '; '));
end
if F.n_failed > 0
    fprintf(fid, '\n%d run(s) of these networks did not finish and are not in the figure.\n', F.n_failed);
end
end

function s = quartile_txt(v, t_max, n_exc)
if n_exc == 0
    s = '-';
elseif isnan(v)
    s = sprintf('> %.3g', t_max);
else
    s = sprintf('%.3g', v);
end
end

function s = mmm(v, fmt)
% median [min, max] over seeds, collapsed when they coincide; NaN ignored.
v = v(~isnan(v));
if isempty(v)
    s = '-';
elseif isscalar(v) || max(v) == min(v)
    s = sprintf(fmt, median(v));
else
    s = sprintf([fmt ' [' fmt ', ' fmt ']'], median(v), min(v), max(v));
end
end

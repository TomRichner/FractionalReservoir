function out = fig_lyapunov_ranked_spectrum(cfg)
% FIG_LYAPUNOV_RANKED_SPECTRUM The ranked top-K Lyapunov spectrum of the 80-neuron network and its histogram, per condition.
%
%   out = FIG_LYAPUNOV_RANKED_SPECTRUM('run_dir', d)
%   out = FIG_LYAPUNOV_RANKED_SPECTRUM('data_file', f, 'out_dir', o)
%
% Draws what the opt-in stage lyapunov_spectrum_n80 saved (FD-032; the top
% 240 exponents per condition on the 80-neuron version of the mu7chain
% network, run_lyapunov_spectrum with the preset's noise). It reads one run
% directory, <run_dir>/lyapunov_spectrum_n80/lyapunov_spectrum_data.mat, or
% the file given, errors when neither holds it, and simulates nothing.
% Supplementary.
%
% ONE COLUMN PER CONDITION, in the condition's color and titled by its name
% (manuscript_style, titles_by_name):
%   top     the exponents sorted in descending order against their rank,
%           one marker series per network seed (a marker shape per seed),
%           and the zero line
%   bottom  the histogram of the same exponents, every seed pooled, drawn
%           as COUNT PER SEED (pooled counts divided by the number of
%           seeds), so it reads the same at one seed or three; the zero
%           line
% Each column has its own exponent limits (no adaptation spans a far wider
% range than the adapting conditions); within a column the exponent axis of
% the two rows (top y, bottom x) is the same range, the column's exponents
% padded by 4% and always including zero.
%
% THE BIN WIDTH, per column, by the Freedman-Diaconis rule on the column's
% pooled exponents with the number of exponents PER SEED as the sample size,
% floored at a fiftieth of the column's range:
%   w = max(2 IQR / m^(1/3), range / 50),   m = the mean number per seed
% (with m the per-seed count, adding seeds does not narrow the bins, which
% is what keeps a count-per-seed histogram comparable across run modes; the
% floor keeps a column whose exponents mostly coincide, no adaptation with
% 70 of 80 at -1/tau_d, from drawing invisible bins). Bin edges are aligned
% to zero. The width of each column is written in the table.
%
% NO STATISTICS ON THE PANELS. Per condition, Fig_Lyapunov_Ranked_Spectrum_table.md
% gives the state dimension N, K used, the seeds, lambda_1, the number of
% positive exponents, the range of the exponents, the convergence drift of
% lambda_1 (|finite-time lambda_1 at the window's end minus at three
% quarters of it|, lya_summary's lambda_1_drift), and the bin width; then
% the run's settings
% (preset, n, T, window, warm-up, fs, K requested).
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions, names drawn as columns in the
% order given; empty draws every condition the .mat holds
% (select_conditions).
%
% See also: run_lyapunov_spectrum, lyapunov_spectrum_n80_stage_plan,
%           fig_lyapunov_spectrum, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.conditions                = {}    % names to draw, in order; empty = all the run holds (select_conditions)
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'lyapunov_spectrum_n80')}, ...
    'lyapunov_spectrum_data.mat', ...
    'Run the stage lyapunov_spectrum_n80 (run_lyapunov_spectrum_n80_stage) first');
D = load(data_file);
[~, idx] = select_conditions({D.results.name}, cfg.conditions);
R = D.results(idx);
n_cond = numel(R);
titles = titles_by_name({R.name}, {R.title});
variant = D.settings.variants{1};    % the stage runs one: the preset's noise
markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Position', [80 80 420 * n_cond, 620], 'Visible', vis);
tl = tiledlayout(fig, 2, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = cell(1, n_cond);
for i = 1:n_cond
    runs = R(i).(variant);
    col = st.condition_color(R(i).name);
    lams = arrayfun(@(r) sort(r.LE_spectrum(:)', 'descend'), runs, 'UniformOutput', false);
    pooled = [lams{:}];
    n_seeds = numel(runs);

    lo = min(min(pooled), 0); hi = max(max(pooled), 0);
    pad = 0.04 * max(hi - lo, eps);
    lim = [lo - pad, hi + pad];
    w = bin_width(pooled, mean(cellfun(@numel, lams)));
    edges = (floor(lim(1) / w):ceil(lim(2) / w)) * w;

    % Top: sorted exponents against rank, one marker series per seed.
    ax = nexttile(tl, i); hold(ax, 'on');
    for s = 1:n_seeds
        lam = lams{s};
        plot(ax, 1:numel(lam), lam, markers{1 + mod(s - 1, numel(markers))}, ...
            'LineStyle', 'none', 'MarkerSize', 3, 'Color', col, 'MarkerFaceColor', 'none', ...
            'LineWidth', 0.75);
    end
    yline(ax, 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw);
    hold(ax, 'off');
    ylim(ax, lim);
    xlim(ax, [0, max(cellfun(@numel, lams)) + 1]);
    title(ax, titles{i}, 'FontWeight', 'normal', 'FontSize', st.title_fs);
    xlabel(ax, 'rank i', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, '\lambda_i (1/s)', 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw); box(ax, 'off');

    % Bottom: histogram, every seed pooled, count per seed.
    ax2 = nexttile(tl, n_cond + i); hold(ax2, 'on');
    counts = histcounts(pooled, edges) / n_seeds;
    histogram(ax2, 'BinEdges', edges, 'BinCounts', counts, 'FaceColor', col, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    xline(ax2, 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw);
    hold(ax2, 'off');
    xlim(ax2, lim);
    xlabel(ax2, '\lambda (1/s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax2, 'count per seed', 'FontSize', st.label_fs); end
    set(ax2, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw); box(ax2, 'off');

    rows{i} = table_row(titles{i}, runs, pooled, w);
end

fig_tag = 'Fig_Lyapunov_Ranked_Spectrum';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        S = D.settings;
        fprintf(fid, '# Ranked Lyapunov spectrum, %d-neuron network\n\n', S.n);
        fprintf(fid, ['Preset `%s`; T = [0, %g] s; exponents accumulated over [%g, %g] s after a %g s ' ...
            'warm-up; fs %g Hz, lya_dt %g s; K requested %d, capped at the state dimension N; ' ...
            'noise: %s; %d seed(s). Data: `%s`.\n\n'], S.preset_name, S.T, S.lya_T_interval, ...
            S.lya_warmup, S.fs, S.lya_dt, S.K, strrep(variant, '_', ' '), S.n_seeds, data_file);
        fprintf(fid, ['Median [min, max] over seeds where seeds differ. lambda_1 and its drift in 1/s ' ...
            '(drift: |finite-time lambda_1 at the end of the window minus at three quarters of it|); ' ...
            'range is the smallest and largest exponent over all seeds; bin width by the ' ...
            'Freedman-Diaconis rule, 2 IQR / m^(1/3), on the pooled exponents with m the ' ...
            'mean number per seed, and at least range / 50; the histogram is counts ' ...
            'per seed.\n\n']);
        fprintf(fid, '| Condition | N | K used | seeds | lambda_1 | n positive | range | lambda_1 drift | bin width |\n');
        fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
end

%% ------------------------------------------------------------------------
function w = bin_width(x, m)
% Freedman-Diaconis on the pooled exponents with m, the number per seed, as
% the sample size, and never narrower than range / 50. The floor matters
% where most exponents coincide: without adaptation 70 of 80 sit at
% -1/tau_d, the IQR is nearly zero, and the bare rule gave bins of 0.0007
% that drew as nothing.
q = quantile(x, [0.25 0.75]);
w = max(2 * (q(2) - q(1)) / m^(1/3), (max(x) - min(x)) / 50);
if ~(w > 0); w = 1e-3; end
end

function s = table_row(ttl, runs, pooled, w)
s = sprintf('| %s | %s | %s | %d | %s | %s | [%.4f, %.4f] | %s | %.4g |', ttl, ...
    mmm([runs.N], '%d'), mmm([runs.K_used], '%d'), numel(runs), mmm([runs.LLE], '%+.4f'), ...
    mmm([runs.n_positive], '%d'), min(pooled), max(pooled), ...
    mmm(abs([runs.lambda_1_drift]), '%.3g'), w);
end

function s = mmm(v, fmt)
% median [min, max] over seeds, collapsed when they coincide.
if numel(v) == 1 || max(v) == min(v)
    s = sprintf(fmt, median(v));
else
    s = sprintf([fmt ' [' fmt ', ' fmt ']'], median(v), min(v), max(v));
end
end

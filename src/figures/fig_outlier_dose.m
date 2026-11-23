function out = fig_outlier_dose(cfg)
% FIG_OUTLIER_DOSE Detachment, responsive E fraction and E/I gain against the structural outlier, across E-to-E levels and seeds.
%
%   out = FIG_OUTLIER_DOSE('run_dir', d)
%   out = FIG_OUTLIER_DOSE('run_dir', d, 'spread', 'band')
%
% Draws what the opt-in stage outlier_spectra saved (FD-040), through
% outlier_fig_data. Reads one run directory, or the index file given as
% data_file; errors when neither holds the stage; simulates nothing.
% FD-037 (Tom, 2026-11-20, on the exploration's outliers_dose: "a good
% figure that should go into the supplement. We will need it with more
% seeds"). Supplementary.
%
% THE NETWORKS are every one in the index that is the reference or an
% E-to-E level (EEx<a>p<b>). HORIZONTAL, all three panels: what the level
% does to the weights, the structural outlier (largest real part of the
% eigenvalues of the realized block-mean matrix) over the structural bulk
% radius; below 1 the outlier is inside the bulk. It differs a little from
% seed to seed, so each level's median point sits at the MEDIAN ACROSS
% SEEDS of that ratio.
%
%   left    detachment of the rightmost frozen-gain eigenvalue: its real
%           part over the THEORETICAL radius R_eff of its own state, the
%           median over the window (log axis; gray line at 1, the radius).
%           Open diamonds on a black line: the same with every gain equal
%           (the rms gain), the weights alone; median over seeds and
%           conditions (it is nearly the same in every condition).
%   middle  the fraction of excitatory neurons off the floor and the
%           ceiling of the nonlinearity (nonzero slope), median over the
%           window
%   right   the mean gain of the excitatory neurons over that of the
%           inhibitory ones, median over the window
%
% Each condition, in its color: the MEDIAN ACROSS SEEDS (thick line, filled
% markers) and the spread.
%
% THE SPREAD (cfg.spread). 'lines' (default): one thin line per seed, each
% point at THAT seed's own structural ratio. Chosen over a band because the
% levels of one seed are the same random network with one block mean
% changed, so a seed's line is a within-network dose curve and shows
% whether every network takes the same course (a band keeps only the
% envelope and would hide, for instance, seeds that cross); because with
% the two seeds of an `essential` run a band is just the two lines filled
% in; and because near the radius the three conditions' bands lie on top of
% each other on the log axis. 'band': the range across seeds (min to max)
% as a translucent band at the median horizontal positions. 'none': the
% medians alone.
%
% The left panel's upper limit leaves headroom above the top diamond; a
% detachment below detachment_floor is drawn at the bottom of the axis.
% The figure is manuscript_style canvas_w wide (or fig_w, which a composite
% sets to its block) and fig_h tall; the vertical names are on two lines and
% the horizontal one, the same for all three panels, is written once, under
% the middle panel. Nothing but axis names on the panels; the conditions are told by color
% (the caption paragraph of Fig_Outlier_Dose_table.md names them), and the
% table gives, per level and condition, median [min, max] over seeds of
% everything drawn, then the per-seed values.
%
% See also: outlier_fig_data, ospec_measures, fig_outlier_plane

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
    cfg.seed double               = []    % unused: every seed is drawn
    cfg.spread (1,:) char {mustBeMember(cfg.spread, {'lines', 'band', 'none'})} = 'lines'
    cfg.detachment_floor (1,1) double {mustBePositive} = 0.25   % the lowest value the left panel's log axis shows
    cfg.fig_w double = []                   % inches; [] = manuscript_style canvas_w (a composite passes its block's width)
    cfg.fig_h (1,1) double {mustBePositive} = 3.1   % inches
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = outlier_fig_data(cfg, 'fig_outlier_dose', 'dose', 'networks', 'dose');
st = F.st;
[n_net, n_cond, n_seed] = size(F.D);
gray = [0.55 0.55 0.55];

%% Per level, condition and seed (NaN where there is no run)
fields = {'structural_ratio', 'ratio_theory_med', 'rho_uniform_theory', 'frac_on_E', 'frac_on_I', 'gain_ratio_EI', ...
    'reined_theory', 'kappa', 'lambda_1'};
V = struct();
for k = 1:numel(fields)
    V.(fields{k}) = nan(n_net, n_cond, n_seed);
end
for i = 1:numel(F.D)
    if isempty(F.D{i}); continue; end
    for k = 1:numel(fields)
        V.(fields{k})(i) = F.D{i}.(fields{k});
    end
end
x_seed = squeeze(median(V.structural_ratio, 2, 'omitnan'));          % n_net x n_seed: the same weights in every condition
x_seed = reshape(x_seed, n_net, n_seed);
x_med = median(x_seed, 2, 'omitnan');                                % where each level is placed
[x_med, order] = sort(x_med);
x_seed = x_seed(order, :);
for k = 1:numel(fields); V.(fields{k}) = V.(fields{k})(order, :, :); end
nets = F.networks(order);
diamond = median(reshape(V.rho_uniform_theory, n_net, []), 2, 'omitnan');   % over conditions and seeds

%% Draw
ykey = {'ratio_theory_med', 'frac_on_E', 'gain_ratio_EI'};
% Two lines each, so that three panels and their names fit across the
% common canvas; the horizontal name is the same for all three and is
% written once, under the middle panel.
ylab = {{'Detachment,', 'Re \lambda_1 / bulk radius'}, {'Responsive fraction', 'of E neurons'}, 'Mean gain, E / I'};
fig_w = cfg.fig_w; if isempty(fig_w); fig_w = st.canvas_w; end
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w cfg.fig_h], 'Visible', vis);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
x_pad = 0.06 * max(max(x_seed(:)) - min([x_seed(:); 0]), eps);
x_lim = [min([0; x_seed(:)]) - x_pad * (min(x_seed(:)) < 0), max(x_seed(:)) + x_pad];
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    Y = V.(ykey{p});                                                 % n_net x n_cond x n_seed
    if p == 1
        % The log axis and its floor. A detachment at or below the floor (a
        % population pinned so hard that the bulk has collapsed onto its
        % center) is drawn AT the floor; the table has the value.
        y_top = max([diamond(:); Y(:); 1], [], 'omitnan');
        y_bot = min([Y(Y > 0); 1], [], 'omitnan');
        y_lim = [max(cfg.detachment_floor, min(0.4, 0.85 * y_bot)), max(2, 1.3 * y_top)];   % headroom above the top diamond
        n_floor = nnz(Y < y_lim(1));
        Y(Y < y_lim(1)) = y_lim(1);                                  % NaN (no run) stays NaN
    end
    if p == 1
        yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
        plot(ax, x_med, diamond, '-d', 'Color', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.25);
    end
    for c = 1:n_cond
        col = F.colors{c};
        Yc = reshape(Y(:, c, :), n_net, n_seed);
        switch cfg.spread
            case 'lines'
                for s = 1:n_seed
                    ok = ~isnan(Yc(:, s));
                    plot(ax, x_seed(ok, s), Yc(ok, s), '-', 'Color', [col 0.45], 'LineWidth', 0.75);
                end
            case 'band'
                lo = min(Yc, [], 2); hi = max(Yc, [], 2);
                ok = ~isnan(lo);
                if nnz(ok) > 1
                    patch(ax, [x_med(ok); flipud(x_med(ok))], [lo(ok); flipud(hi(ok))], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
        end
    end
    for c = 1:n_cond                                                 % the medians on top of every seed's line
        col = F.colors{c};
        ym = median(reshape(Y(:, c, :), n_net, n_seed), 2, 'omitnan');
        ok = ~isnan(ym);
        plot(ax, x_med(ok), ym(ok), '-o', 'Color', col, 'LineWidth', st.line_lw, 'MarkerSize', 8, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'w');
    end
    if p == 1
        ticks = 2 .^ (ceil(log2(y_lim(1))):floor(log2(y_lim(2))));
        set(ax, 'YScale', 'log', 'YLim', y_lim, 'YTick', ticks, 'YMinorTick', 'off');
    else
        ylim(ax, [0, max(1.1, 1.08 * max(Y(:), [], 'omitnan'))]);
    end
    xlim(ax, x_lim);
    set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    if p == 2; xlabel(ax, 'Structural outlier / bulk radius', 'FontSize', st.label_fs); end
    ylabel(ax, ylab{p}, 'FontSize', st.label_fs);
end

fig_tag = 'Fig_Outlier_Dose';
out = struct('figs', fig, 'files', {{}}, 'source', F.file, 'x', x_med, 'networks', {nets}, 'values', V, 'diamond', diamond);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        write_table(fid, fig_tag, F, nets, V, x_seed, x_med, diamond, cfg, y_lim, n_floor);
        clear closer
        out.files{end + 1} = table_name;
    end
end
for w = 1:n_net
    vprintf(cfg.verbose, 'verbose', '%s: structural ratio %.3g, uniform-gain detachment %.3g\n', nets{w}, x_med(w), diamond(w));
end
end

%% ------------------------------------------------------------------------
function write_table(fid, fig_tag, F, nets, V, x_seed, x_med, diamond, cfg, y_lim, n_floor)
[n_net, n_cond, n_seed] = size(V.structural_ratio);
switch cfg.spread
    case 'lines'
        spread_txt = ['Thin lines: each seed, its points at that seed''s own structural ratio (the levels of one seed ' ...
            'are the same random network with the E-to-E mean changed, so a thin line is a dose curve within one network).'];
    case 'band'
        spread_txt = 'Bands: the range across seeds (minimum to maximum), at the median horizontal positions.';
    otherwise
        spread_txt = 'The spread across seeds is not drawn; it is in the table.';
end
first = F.D{find(~cellfun(@isempty, F.D), 1)};
fprintf(fid, '# %s\n\n', fig_tag);
fprintf(fid, ['The E-to-E mean at %d levels (%s; `reference` is 1), each placed on the horizontal axis by what it does ' ...
    'to the weights: the structural outlier (largest real part of the eigenvalues of the realized block-mean matrix) ' ...
    'over the structural bulk radius, the median across seeds. Below 1 the outlier is inside the bulk. Colors are the ' ...
    'conditions (%s). Thick lines with filled markers: the median across the %d seed(s). %s ' ...
    'LEFT: detachment of the rightmost frozen-gain eigenvalue, its real part over the THEORETICAL radius of the bulk ' ...
    'R_eff of the same state, with no finite-size enlargement (median over the %g s window of the per-state ratio; ' ...
    'log axis; gray line: at the theoretical radius). Open diamonds: the same with uniform gains, the weights alone ' ...
    '(median over seeds and conditions; against the theoretical radius it is nearly the same in every condition). ' ...
    'MIDDLE: the fraction of excitatory neurons off the floor and ceiling of the nonlinearity (nonzero slope). ' ...
    'RIGHT: mean gain of the excitatory neurons over that of the inhibitory ones. The left panel shows %.3g to %.3g; ' ...
    '%d run(s) with a detachment below that are drawn at the bottom of the axis (their values are in the table). ' ...
    '%d neurons. Data: `%s`.\n\n'], ...
    n_net, strjoin(nets, ', '), strjoin(F.titles, ', '), n_seed, spread_txt, first.window_s, y_lim, n_floor, first.n, F.file);
fprintf(fid, ['Median over seeds, with [min, max] where they differ. `structural ratio` is the horizontal position; ' ...
    '`uniform` the diamond''s quantity for that condition; `reined` the fraction of the structural detachment the ' ...
    'gains removed, 1 - (detachment - 1)/(uniform - 1), given where the structure has an outlier by the null''s ' ...
    'measure; `lambda_1` the largest Lyapunov exponent (1/s).\n\n']);
fprintf(fid, ['| network | condition | seeds | structural ratio | detachment | uniform | reined | responsive E | responsive I | ' ...
    'gain E/I | kappa | lambda_1 |\n']);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|\n');
g = @(key, w, c) reshape(V.(key)(w, c, :), 1, []);
for w = 1:n_net
    for c = 1:n_cond
        n_here = nnz(~isnan(g('ratio_theory_med', w, c)));
        if n_here == 0; continue; end
        fprintf(fid, '| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', nets{w}, F.cond_names{c}, n_here, ...
            mmm(g('structural_ratio', w, c), '%.3g'), mmm(g('ratio_theory_med', w, c), '%.3f'), ...
            mmm(g('rho_uniform_theory', w, c), '%.3g'), mmm(g('reined_theory', w, c), '%.3f'), ...
            mmm(g('frac_on_E', w, c), '%.3f'), mmm(g('frac_on_I', w, c), '%.3f'), mmm(g('gain_ratio_EI', w, c), '%.3f'), ...
            mmm(g('kappa', w, c), '%.3f'), mmm(g('lambda_1', w, c), '%+.3f'));
    end
end
fprintf(fid, '\nWhere each level is drawn, and the diamonds:\n\n');
fprintf(fid, '| network | structural ratio, median over seeds | per seed | uniform-gain detachment (diamond) |\n|---|---|---|---|\n');
for w = 1:n_net
    fprintf(fid, '| %s | %.4g | %s | %.4g |\n', nets{w}, x_med(w), ...
        strjoin(arrayfun(@(v) sprintf('%.4g', v), x_seed(w, :), 'UniformOutput', false), ', '), diamond(w));
end
fprintf(fid, '\nPer seed (the thin lines):\n\n');
fprintf(fid, '| network | condition | seed | structural ratio | detachment | uniform | responsive E | gain E/I |\n|---|---|---|---|---|---|---|---|\n');
for w = 1:n_net
    for c = 1:n_cond
        for s = 1:n_seed
            if isnan(V.ratio_theory_med(w, c, s)); continue; end
            fprintf(fid, '| %s | %s | %d-%d | %.4g | %.4f | %.4g | %.4f | %.4f |\n', nets{w}, F.cond_names{c}, F.seeds(s, :), ...
                V.structural_ratio(w, c, s), V.ratio_theory_med(w, c, s), V.rho_uniform_theory(w, c, s), ...
                V.frac_on_E(w, c, s), V.gain_ratio_EI(w, c, s));
        end
    end
end
if F.n_failed > 0
    fprintf(fid, '\n%d run(s) of these networks did not finish and are not in the figure.\n', F.n_failed);
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

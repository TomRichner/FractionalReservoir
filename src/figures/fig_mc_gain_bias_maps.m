function out = fig_mc_gain_bias_maps(cfg)
% FIG_MC_GAIN_BIAS_MAPS Memory capacity and lambda_1 on the gain-by-bias grid, one column per condition (block C).
%
%   out = FIG_MC_GAIN_BIAS_MAPS('run_dir', d)
%   out = FIG_MC_GAIN_BIAS_MAPS('data_file', f)
%   out = FIG_MC_GAIN_BIAS_MAPS('data', D, 'save', false)   % D from read_mc_gain_sweep, or a table
%
% Rows (Plan C, section 5):
%   1. memory capacity, median over seeds, one colour scale shared by every
%      condition, white (none) to dark red, with a colour bar;
%   2. lambda_1 with the memory-task input on, median over seeds, a diverging
%      scale centred on zero, shared, with a colour bar;
%   3. (cfg.show_iqr, default true) the interquartile range of memory
%      capacity over seeds, shared scale.
% x is the gain (level_of_chaos) on a log axis, y the bias S_c. The
% lambda_1 = 0 contour (of the median map of the same condition) is drawn
% on every row, green dashed as everywhere in the manuscript. The preset's
% own point (D.preset_gain, D.preset_bias: gain 1, bias 0.2) is a white
% star. A cell with no successful seed is left grey (the axes colour).
%
% SATURATED AND SILENT CELLS are marked lightly so that a low memory
% capacity there is not read as chaos: a cell whose median over seeds of
% the fraction of neurons pinned above 0.9 is >= cfg.flag_fraction (0.5)
% carries a small grey '+', one whose median fraction silent below 0.01 is
% >= cfg.flag_fraction a small grey 'o'. Flags are per condition.
%
% WHICH CONDITIONS: cfg.conditions through select_conditions (empty = all the
% result holds, in its order); titles and colours by name from
% manuscript_style. Works for three and four columns.
%
% The table is read by read_mc_gain_sweep (the one reader), located by
% resolve_mc_gain_sweep.
%
% See also: fig_mc_gain_curves, fig_mc_robustness, mc_gain_bias_statistics

arguments
    cfg.verbose       (1,:) char    = 'minimal'
    cfg.run_dir       (1,:) char    = ''
    cfg.data_file     (1,:) char    = ''
    cfg.data                        = []
    cfg.out_dir       (1,:) char    = ''
    cfg.save          (1,1) logical = true
    cfg.visible       (1,1) logical = true
    cfg.preset_name   (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                  = {}
    cfg.show_iqr      (1,1) logical = true
    cfg.flag_fraction (1,1) double  = 0.5
    cfg.lle_clim                    = []   % [lo hi]; empty = symmetric about 0 from the data
    cfg.fig_tag       (1,:) char    = 'Fig_MC_gain_bias_maps'
end

setup_paths();
D = load_data(cfg);
[conds, ~] = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();

%% The grids
nb = numel(D.biases); ng = numel(D.gains);
MC = nan(nb, ng, nc); L = MC; IQ = MC; PIN = MC; SIL = MC;
for c = 1:nc
    MC(:, :, c)  = mc_gain_grid(D, conds{c}, 'mc_total');
    L(:, :, c)   = mc_gain_grid(D, conds{c}, 'lle_input');
    IQ(:, :, c)  = mc_gain_grid(D, conds{c}, 'mc_total', 'iqr');
    PIN(:, :, c) = mc_gain_grid(D, conds{c}, 'frac_pinned');
    SIL(:, :, c) = mc_gain_grid(D, conds{c}, 'frac_silent');
end
mc_clim = [0, max([MC(:); eps])];
if isempty(cfg.lle_clim)
    m = max([abs(L(:)); eps]);
    lle_clim = [-m m];
else
    lle_clim = cfg.lle_clim;
end
iq_clim = [0, max([IQ(:); eps])];

rows = {struct('V', MC, 'clim', mc_clim, 'cmap', seq_cmap(), ...
               'label', 'Memory capacity', 'cb', 'MC (median)')
        struct('V', L,  'clim', lle_clim, 'cmap', div_cmap(), ...
               'label', '\lambda_1, task input on', 'cb', '\lambda_1 (1/s)')};
if cfg.show_iqr
    rows{end+1} = struct('V', IQ, 'clim', iq_clim, 'cmap', gray_cmap(), ...
                         'label', 'Spread over seeds', 'cb', 'MC IQR');
end
nr = numel(rows);

%% Draw
fig = figure('Color', 'w', 'Position', [60 60 260 * nc + 140, 250 * nr + 60]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
outer = tiledlayout(fig, nr, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
x = log10(D.gains); y = D.biases;
xe = edges(x); ye = edges(y);
[xt, xtl] = gain_ticks(D.gains);
ax_all = gobjects(nr, nc);
for r = 1:nr
    R = rows{r};
    inner = tiledlayout(outer, 1, nc, 'TileSpacing', 'compact', 'Padding', 'tight');
    inner.Layout.Tile = r;
    for c = 1:nc
        ax = nexttile(inner); hold(ax, 'on');
        ax_all(r, c) = ax;
        draw_cells(ax, x, y, R.V(:, :, c));
        colormap(ax, R.cmap); clim(ax, R.clim);
        set(ax, 'Color', [0.85 0.85 0.85], 'Layer', 'top', 'TickDir', 'out', ...
            'FontSize', st.tick_fs - 3, 'XTick', xt, 'XTickLabel', xtl, 'Box', 'on');
        xlim(ax, xe([1 end])); ylim(ax, ye([1 end]));
        mark_flags(ax, x, y, PIN(:, :, c) >= cfg.flag_fraction, SIL(:, :, c) >= cfg.flag_fraction);
        zero_contour(ax, x, y, L(:, :, c), st);
        plot(ax, log10(D.preset_gain), D.preset_bias, 'p', 'MarkerSize', 11, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        if r == 1
            % A qualifier in parentheses goes to a second line, so that four
            % long titles ("... Adaptation (product)") do not run together.
            t = regexprep(cond_title(st, conds{c}), '\s+(\(.*\))$', '\n$1');
            title(ax, strsplit(t, newline), 'FontSize', st.title_fs - 2, ...
                'FontWeight', 'normal', 'Color', cond_color(st, conds{c}));
        end
        if c == 1
            ylabel(ax, {R.label, 'Bias S_c'}, 'FontSize', st.label_fs - 3);
        else
            set(ax, 'YTickLabel', []);
        end
        if r == nr
            xlabel(ax, 'Gain', 'FontSize', st.label_fs - 3);
        else
            set(ax, 'XTickLabel', []);
        end
    end
    cb = colorbar(ax_all(r, nc));
    cb.Layout.Tile = 'east';
    cb.Label.String = R.cb; cb.FontSize = st.tick_fs - 4;
end
title(outer, sprintf(['Median over %d seeds. Green dashed: \\lambda_1 = 0. Star: preset (gain %.3g, bias %.3g). ' ...
    '+ saturated, o silent (>= %.0f%% of neurons). Grey: no successful job.'], ...
    numel(D.seeds), D.preset_gain, D.preset_bias, 100 * cfg.flag_fraction), 'FontSize', st.tick_fs - 4);

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', D.source);
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
end
end

%% ------------------------------------------------------------------------
function D = load_data(cfg)
if ~isempty(cfg.data)
    if isstruct(cfg.data) && isfield(cfg.data, 'gains'); D = cfg.data; else; D = read_mc_gain_sweep(cfg.data); end
else
    D = resolve_mc_gain_sweep(cfg.data_file, cfg.run_dir);
end
end

function e = edges(v)
% Cell edges: midpoints, extended by half a step at each end.
v = v(:)';
if isscalar(v); e = v + [-0.5 0.5]; return; end
m = (v(1:end-1) + v(2:end)) / 2;
e = [v(1) - (m(1) - v(1)), m, v(end) + (v(end) - m(end))];
end

function draw_cells(ax, x, y, V)
xe = edges(x); ye = edges(y);
C = nan(numel(ye), numel(xe));
C(1:end-1, 1:end-1) = V;
[X, Y] = meshgrid(xe, ye);
surface(ax, X, Y, zeros(size(X)), C, 'FaceColor', 'flat', 'EdgeColor', 'none');
end

function mark_flags(ax, x, y, pinned, silent)
[X, Y] = meshgrid(x, y);
g = [0.55 0.55 0.55];
if any(pinned(:))
    plot(ax, X(pinned), Y(pinned), '+', 'Color', g, 'MarkerSize', 5, 'LineWidth', 0.8);
end
if any(silent(:))
    plot(ax, X(silent), Y(silent), 'o', 'Color', g, 'MarkerSize', 4, 'LineWidth', 0.8);
end
end

function zero_contour(ax, x, y, Lc, st)
if numel(x) < 2 || numel(y) < 2 || nnz(isfinite(Lc)) < 4; return; end
if ~(any(Lc(:) > 0) && any(Lc(:) < 0)); return; end
contour(ax, x, y, Lc, [0 0], '--', 'LineColor', st.zeroline_color, 'LineWidth', 1.8);
end

function [xt, xtl] = gain_ticks(gains)
cand = [0.01 0.03 0.1 0.3 1 3 10 30 100];
xt_v = cand(cand >= min(gains) * 0.999 & cand <= max(gains) * 1.001);
if numel(xt_v) < 2; xt_v = gains([1 end]); end
xt = log10(xt_v);
xtl = arrayfun(@(v) sprintf('%g', v), xt_v, 'UniformOutput', false);
end

function c = seq_cmap()
k = [1 1 1; 0.99 0.88 0.55; 0.93 0.55 0.25; 0.70 0.18 0.18; 0.30 0.02 0.20];
c = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

function c = div_cmap()
k = [0.02 0.19 0.38; 0.26 0.58 0.76; 1 1 1; 0.84 0.38 0.30; 0.40 0.00 0.05];
c = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

function c = gray_cmap()
c = repmat(linspace(1, 0.15, 256)', 1, 3);
end

function t = cond_title(st, nm)
if isKey(st.condition_title, nm); t = st.condition_title(nm); else; t = nm; end
end

function c = cond_color(st, nm)
if isKey(st.condition_color, nm); c = st.condition_color(nm); else; c = [0.3 0.3 0.3]; end
end

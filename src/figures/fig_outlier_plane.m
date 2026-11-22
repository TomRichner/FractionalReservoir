function out = fig_outlier_plane(cfg)
% FIG_OUTLIER_PLANE Frozen-gain spectra in the complex plane, networks by conditions, with the neurons on the response curve.
%
%   out = FIG_OUTLIER_PLANE('run_dir', d)
%   out = FIG_OUTLIER_PLANE('run_dir', d, 'seed', [1 2], 'outlier_factor', 1.05)
%
% Draws what the opt-in stage outlier_spectra saved (FD-040), through
% outlier_fig_data. Reads one run directory, or the index file given as
% data_file; errors when neither holds the stage; simulates nothing.
% FD-036 (Tom, 2026-11-20, on the exploration's
% outliers_complex_plane_v3_alt: "looks awesome ... should go in the paper,
% or at least the supplement").
%
% ROWS are the networks (cfg.networks; the reference and E-to-E x 1.5),
% COLUMNS the conditions. Each panel overlays the calibration states of ONE
% seed (cfg.seed; the first seed of the index by default):
%   dots      every eigenvalue of the frozen-gain effective connectivity
%             M = W D at each state, as a rate (lambda - 1)/tau_d in 1/s,
%             pale and translucent, in the condition's color
%   marked    an eigenvalue farther from the center -1/tau_d than
%             outlier_factor (1.03) times the theoretical radius OF ITS OWN
%             STATE: opaque, a little larger, same color, no edge. (The
%             exploration judged against the drawn circle, the median
%             radius, which counted eigenvalues of the states with a
%             smaller bulk as inside and of those with a larger one as
%             outliers: 15.7 against 5.8 per state for single-timescale
%             adaptation in the reference network. Lead's decision.)
%   solid gray circle   the theoretical bulk: radius R_eff/tau_d from the
%             gains, the MEDIAN over the drawn states, no finite-size factor
%   thin dashed circle  that radius times the run's kappa: how far the
%             rightmost eigenvalue of a structureless matrix with the same
%             gains reaches (95th percentile of the null draws)
%   diamond   where the rightmost eigenvalue would be with every gain equal
%             (the rms gain): the outlier the weights alone make. Median
%             over the drawn states.
%   No ringed markers for the rightmost eigenvalue (Tom, 2026-11-20).
%
% AXES (complex_plane_axes; docs/Figure_style_guide.md): through the
% origin, 1:1, IDENTICAL in every panel, ticks only outside the largest
% dashed circle. The imaginary axis is the stability boundary. The window
% is fitted once to the six drawn panels (or given: re_lim, im_lim,
% far_lim): `margin` bulk radii around the center, widened on the right to
% the rightmost eigenvalue and to any diamond within half a radius of it;
% diamonds farther out are shown beyond a BREAK in the real axis. There is
% no break on the left (lead's decision): eigenvalues left of the window
% (with E-to-E x 1.5 and adaptation, one real eigenvalue per state) are not
% drawn; the table counts them, over every seed, and the caption paragraph
% names them.
%
% THE INSET, lower right of each panel (inset = false leaves it out): the
% firing-rate nonlinearity, the piecewise sigmoid centered on 0, and one
% dot per neuron per drawn state at the neuron's operating point, from the
% stage's saved states (nothing is re-simulated). Horizontal: u, the
% argument of the nonlinearity relative to the neuron's own setpoint.
% Vertical: the rate plus a one-sided random offset, excitatory neurons
% from the curve UP and inhibitory ones from the curve DOWN (inset_jitter),
% in the manuscript's one E and one I color (manuscript_style E_color,
% I_color; until 2026-11-22 the dark rose and teal of
% excitatory_colormap(1) / inhibitory_colormap(1)). Neurons beyond |u| = inset_u_lim are drawn in a
% short strip past a gap in the floor or ceiling line: "further out", not
% how far. The only text on it: the letters E and I (inset_marks), each in
% its cell-type color on its side of the curve.
%
% Text on the figure: the condition names above the columns, in their
% colors, the network names beside the rows, the axis names and ticks.
% Fig_Outlier_Plane_table.md has a caption paragraph and, per panel, as
% the median and range over ALL seeds of the index: the radius, kappa, the
% detachment, the diamond, the marked eigenvalues per state (all, right of
% the center, right of zero) with the null's count, the eigenvalues left
% of the window, and the fractions of E and I neurons on the sloped part,
% the floor and the ceiling.
%
% See also: outlier_fig_data, complex_plane_axes, fig_outlier_dose,
%           fig_outlier_slow_variables, fig_outlier_dwell, ospec_load_run

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''    % the stage's ospec_index.mat, instead of run_dir
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the stage's files
    cfg.networks                  = {'reference', 'EEx1p5'}
    cfg.conditions                = {}    % names to draw, in order; empty = all the stage ran
    cfg.seed double               = []    % the seed drawn; [] = the first of the index
    cfg.outlier_factor (1,1) double {mustBePositive} = 1.03
    cfg.re_lim double             = []    % [lo hi] 1/s; [] = fitted to the drawn panels
    cfg.im_lim double             = []
    cfg.far_lim double            = []    % [f1 f2] beyond a break; [] = fitted (a break only if a diamond is far out)
    cfg.margin (1,1) double {mustBePositive} = 1.2     % the fitted window, in largest dashed radii from the center
    cfg.tick_step double          = []
    cfg.dot_size (1,1) double     = 5
    cfg.dot_alpha (1,1) double    = 0.25
    cfg.outlier_size (1,1) double = 10
    cfg.inset (1,1) logical       = true
    cfg.inset_u_lim (1,1) double {mustBePositive} = 1.5
    cfg.inset_jitter (1,1) double = 0.12
    cfg.inset_size (1,1) double   = 1.5
    cfg.inset_alpha (1,1) double  = 0.2
    cfg.inset_seed (1,1) double   = 1
    cfg.inset_height (1,1) double {mustBePositive} = 0.62   % in largest dashed radii
    cfg.inset_marks (1,1) logical = true    % the letters E and I on the inset
    cfg.inset_mark_fs (1,1) double = 11
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = outlier_fig_data(cfg, 'fig_outlier_plane', 'plane', 'networks', cfg.networks);
st = F.st;
[n_net, n_cond, n_seed] = size(F.D);
f = cfg.outlier_factor;
gray = [0.55 0.55 0.55];

%% What each panel draws (the drawn seed), in rates
P = cell(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        d = F.D{w, c, F.i_seed};
        if isempty(d); continue; end
        p = struct('d', d, 'center', d.center, 'rad', median(d.rad), 'rad_cal', d.kappa * median(d.rad), ...
            'diamond', median(d.z_uniform), 'is_out', abs(d.z - d.center) > f .* d.rad);
        P{w, c} = p;
    end
end
drawn = [P{~cellfun(@isempty, P)}];
if isempty(drawn)
    error('fig_outlier_plane:NoSeed', 'No finished run of seed [%d %d] among the panels.', F.seeds(F.i_seed, 1), F.seeds(F.i_seed, 2));
end
if any(abs([drawn.center] - drawn(1).center) > 1e-12)
    error('fig_outlier_plane:Centers', 'The runs have different tau_d; one window does not fit them.');
end
c0 = drawn(1).center;
rad_max = max([drawn.rad_cal]);

%% One window for every panel
W = fit_window(cfg, drawn, c0, rad_max);
axes_args = {'re_lim', W.re_lim, 'im_lim', W.im_lim, 'far_lim', W.far_lim, 'bulk_center', c0, 'bulk_radius', rad_max, ...
    'tick_step', W.tick_step, 'font_size', st.tick_fs, 'line_width', st.axis_lw};
if ~isempty(W.far_lim); axes_args = [axes_args, {'gap', W.gap}]; end
g = complex_plane_axes([], axes_args{:});

%% Panels placed by hand (inches): identical sizes, 1:1
fig_w = 13; left = 0.6; right = 0.12; col_gap = 0.3; top_m = 0.5; row_gap = 0.35; bottom = 0.15;
pw = (fig_w - left - right - (n_cond - 1) * col_gap) / n_cond;
ph = pw * diff(g.y_lim) / diff(g.x_lim);
% A tick whose label a dot or a diamond of ANY panel would lie on is left out
% of every panel (the panels share their ticks); the sidecar names it.
[hidden_re, hidden_im] = covered_ticks(g, P, cfg, st.tick_fs, diff(g.x_lim) / pw);
g.re_ticks = setdiff(g.re_ticks, hidden_re); g.im_ticks = setdiff(g.im_ticks, hidden_im);
g.hidden_re_ticks = hidden_re; g.hidden_im_ticks = hidden_im;
fig_h = top_m + n_net * ph + (n_net - 1) * row_gap + bottom;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', vis);
ei_col = [st.E_color; st.I_color];                % the manuscript's one E and one I color
rs = RandStream('mt19937ar', 'Seed', cfg.inset_seed);   % private: the jitter does not touch the global stream
th = linspace(0, 2 * pi, 361);
n_not_drawn = zeros(n_net, n_cond); n_marked_drawn = zeros(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        ax = axes(fig, 'Units', 'inches', 'Position', ...
            [left + (c - 1) * (pw + col_gap), bottom + (n_net - w) * (ph + row_gap), pw, ph]);
        names = {'re_label', '', 'im_label', ''};       % the names once per row and per column
        if w == n_net; names{2} = 'Re (s^{-1})'; end
        if c == 1;     names{4} = 'Im (s^{-1})'; end
        cp = complex_plane_axes(ax, axes_args{:}, names{:});
        remove_ticks(cp, hidden_re, hidden_im);
        p = P{w, c};
        if ~isempty(p)
            col = F.colors{c};
            plot(ax, cp.x(p.center + p.rad_cal * cos(th)), p.rad_cal * sin(th), '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.75);
            plot(ax, cp.x(p.center + p.rad * cos(th)), p.rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
            z = p.d.z(:); is_out = p.is_out(:);
            shown = real(z) >= g.re_lim(1) & real(z) <= last_re(g) & abs(imag(z)) <= g.y_lim(2);
            n_not_drawn(w, c) = nnz(~shown);
            n_marked_drawn(w, c) = nnz(shown & is_out);
            zs = z(shown & ~is_out); zo = z(shown & is_out);
            scatter(ax, cp.x(real(zs)), imag(zs), cfg.dot_size, col, 'filled', 'MarkerFaceAlpha', cfg.dot_alpha, 'MarkerEdgeColor', 'none');
            scatter(ax, cp.x(real(zo)), imag(zo), cfg.outlier_size, col, 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
            plot(ax, cp.x(p.diamond), 0, 'd', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
        end
        set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);    % plotting must not have moved them
        if cfg.inset && ~isempty(p)
            per = ax.Position(3) / diff(cp.x_lim);      % inches per unit (1:1, so the same for Im)
            rect = [ax.Position(1) + (W.inset_x - cp.x_lim(1)) * per, ax.Position(2), ...
                (cp.x_lim(2) - W.inset_x) * per, cfg.inset_height * rad_max * per];
            draw_inset(fig, rect, p.d, ei_col, rs, cfg, st.axis_lw);
        end
        if w == 1
            title(ax, F.titles{c}, 'FontSize', st.title_fs, 'FontWeight', 'normal', 'Color', F.colors{c});
        end
        if c == 1
            set(ax.YLabel, 'String', F.net_labels{w}, 'FontSize', st.label_fs, 'Visible', 'on', 'Color', 'k');
        end
    end
end

%% The numbers, per panel over every seed
S = cell(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        S{w, c} = panel_numbers(squeeze(F.D(w, c, :)), f, g.re_lim(1), cfg.inset_u_lim);
    end
end

fig_tag = 'Fig_Outlier_Plane';
out = struct('figs', fig, 'files', {{}}, 'source', F.file, 'window', g, 'numbers', {S}, 'seed', F.seeds(F.i_seed, :));
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        write_table(fid, fig_tag, F, P, S, g, W, cfg, rad_max, c0, n_not_drawn, n_marked_drawn);
        clear closer
        out.files{end + 1} = table_name;
    end
end
for w = 1:n_net
    for c = 1:n_cond
        if isempty(S{w, c}); continue; end
        vprintf(cfg.verbose, 'verbose', '%s | %s: marked per state %s (null %s), left of the window %s\n', ...
            F.networks{w}, F.cond_names{c}, mmm(S{w, c}.n_marked, '%.1f'), mmm(S{w, c}.n_marked_null, '%.1f'), ...
            mmm(S{w, c}.n_left, '%.2f'));
    end
end
if n_seed == 1
    vprintf(cfg.verbose, 'verbose', 'One seed in the index: the table''s ranges are that seed''s values.\n');
end
end

%% ------------------------------------------------------------------------
function W = fit_window(cfg, drawn, c0, rad_max)
% The window of every panel: given, or fitted to the drawn panels. Also the
% left edge of the inset (plotted x) and the room the inset needs.
m = cfg.margin;
im_lim = [-m, m] * rad_max;
re_lim = [c0 - m * rad_max, max(c0 + m * rad_max, 0.1 * rad_max)];
if ~isempty(cfg.im_lim); im_lim = cfg.im_lim(:)'; end
if ~isempty(cfg.re_lim); re_lim = cfg.re_lim(:)'; end
far_lim = cfg.far_lim(:)';
pad = 0.08 * rad_max;
diamonds = [drawn.diamond];
if isempty(cfg.re_lim)
    z = [];
    for k = 1:numel(drawn); z = [z; drawn(k).d.z(:)]; end %#ok<AGROW>
    z = z(abs(imag(z)) <= im_lim(2));
    re_lim(2) = max(re_lim(2), max(real(z)) + pad);
    if isempty(far_lim)
        near = diamonds <= re_lim(2) + 0.5 * rad_max;       % a diamond this close widens the window; farther, a break
        while any(near & diamonds + 2 * pad > re_lim(2))
            re_lim(2) = max(diamonds(near)) + 2 * pad;
            near = diamonds <= re_lim(2) + 0.5 * rad_max;
        end
        if any(~near)
            far_lim = [min(diamonds(~near)) - 0.22 * rad_max, max(diamonds(~near)) + 0.22 * rad_max];
        end
    end
end
if ~isempty(far_lim) && far_lim(1) <= re_lim(2)
    error('fig_outlier_plane:Window', 'far_lim must begin right of the main window, which ends at %.3g 1/s.', re_lim(2));
end
visible = diamonds >= re_lim(1) & diamonds <= re_lim(2);
if ~isempty(far_lim); visible = visible | (diamonds >= far_lim(1) & diamonds <= far_lim(2)); end
if any(~visible)
    error('fig_outlier_plane:MarkOutside', 'A uniform-gain diamond (%.3g 1/s) is outside the window: widen re_lim or move far_lim.', ...
        diamonds(find(~visible, 1)));
end
W.gap = 0.16 * rad_max;
% The inset begins just right of the largest dashed circle and runs to the
% right edge; it gets at least 1.1 radii of width (fitted windows only).
W.inset_x = min(c0 + 1.03 * rad_max, re_lim(2));
need = 1.1 * rad_max;
if cfg.inset && isempty(cfg.re_lim) && isempty(cfg.far_lim)
    width = re_lim(2) - W.inset_x;
    if ~isempty(far_lim); width = width + W.gap + diff(far_lim); end
    if width < need
        if isempty(far_lim); re_lim(2) = re_lim(2) + need - width; else; far_lim(2) = far_lim(2) + need - width; end
    end
end
W.re_lim = re_lim; W.im_lim = im_lim; W.far_lim = far_lim;
W.tick_step = cfg.tick_step;
if isempty(W.tick_step); W.tick_step = nice_step(rad_max / 1.8); end
end

function hi = last_re(g)
if isempty(g.far_lim); hi = g.re_lim(2); else; hi = g.far_lim(2); end
end

function draw_inset(fig, rect, d, ei_col, rs, cfg, lw)
% The nonlinearity with every neuron's operating point at every drawn
% state. Horizontal: u in [-L, L], then a gap, then a strip for the neurons
% beyond, on each side. E above the curve, I below it.
L = cfg.inset_u_lim; gap = 0.2 * L; strip = 0.4 * L; jit = cfg.inset_jitter;
u = d.u(:); r = d.rate(:);
ty = repmat(d.type, 1, size(d.u, 2)); ty = ty(:);
ax = axes(fig, 'Units', 'inches', 'Position', rect, 'Visible', 'off', 'Tag', 'fig_outlier_plane_inset');
hold(ax, 'on');
set(ax, 'XLim', [-1 1] * (L + gap + strip + 0.05), 'YLim', [-jit - 0.06, 1 + jit + 0.06], 'Clipping', 'off');
x = u;
lo = u < -L; hi = u > L;
x(lo) = -L - gap - strip * rs.rand(nnz(lo), 1);
x(hi) =  L + gap + strip * rs.rand(nnz(hi), 1);
y = r + jit * rs.rand(numel(r), 1) .* (3 - 2 * ty);     % type 1 (E) up, type 2 (I) down
order = rs.randperm(numel(x));                          % neither type is drawn on top of the other
scatter(ax, x(order), y(order), cfg.inset_size, ei_col(ty(order), :), 'filled', ...
    'MarkerFaceAlpha', cfg.inset_alpha, 'MarkerEdgeColor', 'none');
uu = linspace(-L, L, 400);
plot(ax, uu, SRNNCellTypePairs.piecewiseSigmoid(uu, d.S_a, 0), 'k-', 'LineWidth', 0.75 * lw);
plot(ax, [-L - gap - strip, -L - gap], [0 0], 'k-', 'LineWidth', 0.75 * lw);   % the floor, beyond the gap
plot(ax, [L + gap, L + gap + strip], [1 1], 'k-', 'LineWidth', 0.75 * lw);     % the ceiling, beyond the gap
if cfg.inset_marks
    % Which row of dots is which (Tom, 2026-11-21: the one text allowed on
    % the inset). "E" on the side of the curve its dots are jittered to
    % (above), "I" below, each in its cell-type color, in the two empty
    % corners of the inset.
    mark = {'FontSize', cfg.inset_mark_fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Tag', 'fig_outlier_plane_inset_mark'};
    text(ax, -0.7 * L, 0.5 + 2 * jit, 'E', 'Color', ei_col(1, :), mark{:});
    text(ax,  0.7 * L, 0.5 - 2 * jit, 'I', 'Color', ei_col(2, :), mark{:});
end
end

function s = panel_numbers(runs, f, re_left, u_lim)
% One value per seed of each number of a panel (rows of the table give
% their median and range). Counts are means per state.
runs = runs(~cellfun(@isempty, runs));
if isempty(runs); s = []; return; end
n = numel(runs);
fields = {'rad', 'kappa', 'detachment', 'diamond', 'n_marked', 'n_marked_right', 'n_marked_unstable', 'n_beyond', ...
    'n_marked_null', 'n_left', 'left_lo', 'left_hi', 'E_slope', 'E_floor', 'E_ceil', 'I_slope', 'I_floor', 'I_ceil', 'u_beyond'};
s = cell2struct(repmat({nan(1, n)}, numel(fields), 1), fields, 1);
s.seeds = zeros(n, 2);
for k = 1:n
    d = runs{k};
    s.seeds(k, :) = d.seed;
    dist = abs(d.z - d.center);
    out = dist > f .* d.rad;
    s.rad(k) = median(d.rad);
    s.kappa(k) = d.kappa;
    s.detachment(k) = d.ratio_theory_med;
    s.diamond(k) = median(d.z_uniform);
    s.n_marked(k) = mean(sum(out, 1));
    s.n_marked_right(k) = mean(sum(out & real(d.z) > d.center, 1));
    s.n_marked_unstable(k) = mean(sum(out & real(d.z) > 0, 1));
    s.n_beyond(k) = mean(sum(dist > d.rad, 1));
    out_null = abs(d.z_null - d.center) > f .* reshape(d.rad, 1, 1, []);
    s.n_marked_null(k) = mean(sum(out_null, 1), 'all');
    is_left = real(d.z) < re_left;
    s.n_left(k) = mean(sum(is_left, 1));
    if any(is_left(:))
        s.left_lo(k) = min(real(d.z(is_left))); s.left_hi(k) = max(real(d.z(is_left)));
    end
    names = 'EI';
    for q = 1:2
        is = d.type == q;
        s.([names(q) '_slope'])(k) = mean(d.slope(is, :) > 0, 'all');
        s.([names(q) '_floor'])(k) = mean(d.slope(is, :) == 0 & d.rate(is, :) < 0.5, 'all');
        s.([names(q) '_ceil'])(k)  = mean(d.slope(is, :) == 0 & d.rate(is, :) >= 0.5, 'all');
    end
    s.u_beyond(k) = mean(abs(d.u) > u_lim, 'all');
end
end

function write_table(fid, fig_tag, F, P, S, g, W, cfg, rad_max, c0, n_not_drawn, n_marked_drawn)
f = cfg.outlier_factor;
[n_net, n_cond, n_seed] = size(F.D);
first = F.D{find(~cellfun(@isempty, F.D), 1)};
seed = F.seeds(F.i_seed, :);
if isempty(g.far_lim)
    break_txt = 'There is no break in the real axis.';
else
    break_txt = sprintf(['THE REAL AXIS IS BROKEN: the stretch from %.3g to %.3g 1/s is left out (the gap, with a slanted ' ...
        'stroke on each side), and the segment to its right runs from %.3g to %.3g 1/s at the same scale; it is drawn in ' ...
        'every panel.'], g.re_lim(2), g.far_lim(1), g.far_lim);
end
% The eigenvalues left of the window, named: which panels have them, and where.
left_txt = {};
for w = 1:n_net
    for c = 1:n_cond
        s = S{w, c};
        if isempty(s) || ~any(s.n_left > 0); continue; end
        left_txt{end + 1} = sprintf('%s, %s: %s per state, real parts from %.3g to %.3g 1/s', F.networks{w}, ...
            F.cond_names{c}, mmm(s.n_left, '%.2g'), min(s.left_lo), max(s.left_hi)); %#ok<AGROW>
    end
end
if isempty(left_txt)
    left_txt = 'No eigenvalue of any seed lies left of the window.';
else
    left_txt = ['EIGENVALUES LEFT OF THE WINDOW ARE NOT DRAWN (there is no break on the left). Over all seeds: ' ...
        strjoin(left_txt, '; ') '.'];
end
fprintf(fid, '# %s\n\n', fig_tag);
fprintf(fid, ['Rows: the networks (%s). Columns: the conditions. Each panel overlays the %d calibration states of one ' ...
    'run, network seed [%d %d] (%d neurons; the states are spread evenly over the %g s measured window). DOTS: every ' ...
    'eigenvalue of the frozen-gain effective connectivity (each neuron''s outgoing weights times its gain at that ' ...
    'state, slow variables held), drawn as a rate (lambda - 1)/tau_d in 1/s, so the imaginary axis is the stability ' ...
    'boundary. SOLID GRAY CIRCLE: the theoretical bulk, centered on -1/tau_d = %.3g 1/s, radius R_eff/tau_d predicted ' ...
    'from the gains, the median over the drawn states, with no finite-size factor. THIN DASHED CIRCLE: that radius ' ...
    'times kappa, how far the rightmost eigenvalue of a structureless matrix with the same entries and gains reaches ' ...
    '(95th percentile of the null draws). OPAQUE, LARGER DOTS: eigenvalues farther from the center than %.3g times the ' ...
    'theoretical radius of THEIR OWN state, in any direction (so a marked dot can lie inside the drawn median circle ' ...
    'and an unmarked one outside it). OPEN DIAMOND: where the rightmost eigenvalue would be if every neuron had the ' ...
    'same gain (the rms of the observed gains): the outlier the weights alone make (median over the drawn states). ' ...
    'The rightmost eigenvalue of each state is not ringed.\n\n'], strjoin(F.networks, ', '), first.n_cal, seed, first.n, ...
    first.window_s, c0, f);
fprintf(fid, ['AXES: through the origin (`complex_plane_axes`), one unit of Re as long as one unit of Im, identical in ' ...
    'every panel: real axis [%.3g, %.3g] 1/s, imaginary axis [%.3g, %.3g] 1/s; ticks only outside the largest dashed ' ...
    'circle (radius %.3g 1/s). %s %s %s\n\n'], g.re_lim, g.y_lim, rad_max, break_txt, left_txt, hidden_ticks_text(g));
if cfg.inset
    a = first.S_a / 2;
    fprintf(fid, ['INSET (lower right of each panel; the letters E and I mark the two cell types). The black curve is the firing-rate nonlinearity, ' ...
        'the piecewise sigmoid (S_a = %g): rate 0 (floor) below u = %g, quadratic up to %g, slope 1 up to %g, quadratic ' ...
        'up to %g, rate 1 (ceiling) above. Each dot is one neuron at one of the drawn states, at its operating point; ' ...
        'excitatory neurons in orange-red ABOVE the curve, inhibitory ones in blue BELOW it (the one E and one I color ' ...
        'of `manuscript_style`, E_color and I_color; they are not the condition colors). Horizontal: ' ...
        'u = x - (c/K) sum_k a_k - S_c, the argument of the nonlinearity relative to the neuron''s own setpoint, so ' ...
        'that every neuron lies on the one curve; shown from %g to %g. Vertical: the rate, plus a uniform random ' ...
        'offset of 0 to %g away from the curve (up for E, down for I; fixed seed %d): the dots are ON the curve, and ' ...
        'the band at a flat part is one rate. A neuron beyond the shown range of u is drawn in the short strip past ' ...
        'the gap in the floor line (left) or the ceiling line (right), at a random place in it. The operating points ' ...
        'are the stage''s saved states; nothing is re-simulated.\n\n'], first.S_a, a - 1, -a, a, 1 - a, ...
        -cfg.inset_u_lim, cfg.inset_u_lim, cfg.inset_jitter, cfg.inset_seed);
end
fprintf(fid, ['Numbers per panel: the median over the %d seed(s) of the index, with [min, max] where the seeds differ ' ...
    '(seeds: %s). Counts of eigenvalues are means per state over a run''s calibration states. `radius` is the median ' ...
    'over those states (the solid circle of that seed); `detachment` is Re lambda_1 over the state''s own R_eff, ' ...
    'median over the whole window; `diamond` the uniform-gain rightmost real part as a rate; `marked` the ' ...
    'eigenvalues beyond %.3g own radii, of which `right of center` (Re > -1/tau_d) and `right of zero` (locally ' ...
    'unstable with the slow variables held); `beyond 1.00` those beyond the own radius itself; `null` the count ' ...
    'beyond %.3g own radii for the structureless matrices saved with each state (what chance alone marks); ' ...
    '`left of window` the eigenvalues with real part below %.3g 1/s. Data: `%s`.\n\n'], n_seed, ...
    strjoin(arrayfun(@(k) sprintf('[%d %d]', F.seeds(k, :)), 1:n_seed, 'UniformOutput', false), ', '), f, f, g.re_lim(1), F.file);
fprintf(fid, ['| network | condition | seeds | radius (1/s) | kappa | detachment | diamond (1/s) | marked | right of center | ' ...
    'right of zero | beyond 1.00 | null | left of window |\n']);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for w = 1:n_net
    for c = 1:n_cond
        s = S{w, c};
        if isempty(s); continue; end
        fprintf(fid, '| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', F.networks{w}, F.cond_names{c}, ...
            size(s.seeds, 1), mmm(s.rad, '%.3g'), mmm(s.kappa, '%.3f'), mmm(s.detachment, '%.3f'), mmm(s.diamond, '%.3g'), ...
            mmm(s.n_marked, '%.1f'), mmm(s.n_marked_right, '%.1f'), mmm(s.n_marked_unstable, '%.1f'), ...
            mmm(s.n_beyond, '%.1f'), mmm(s.n_marked_null, '%.1f'), mmm(s.n_left, '%.2f'));
    end
end
fprintf(fid, ['\nThe neurons on the response curve: fractions of neuron-states (all neurons, all calibration states) ' ...
    'on the sloped part (slope > 0), the floor (rate 0) and the ceiling (rate 1), by cell type, and the fraction ' ...
    'beyond the inset''s range of u (drawn in the strips). Median over seeds [min, max].\n\n']);
fprintf(fid, '| network | condition | E sloped | E floor | E ceiling | I sloped | I floor | I ceiling | beyond u range |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
for w = 1:n_net
    for c = 1:n_cond
        s = S{w, c};
        if isempty(s); continue; end
        fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', F.networks{w}, F.cond_names{c}, ...
            mmm(s.E_slope, '%.3f'), mmm(s.E_floor, '%.3f'), mmm(s.E_ceil, '%.3f'), mmm(s.I_slope, '%.3f'), ...
            mmm(s.I_floor, '%.3f'), mmm(s.I_ceil, '%.3f'), mmm(s.u_beyond, '%.3f'));
    end
end
fprintf(fid, ['\nThe drawn panels (seed [%d %d]): the radii of the two circles, the marked dots in view, and the ' ...
    'eigenvalue dots outside the window on any side (not drawn), totals over the %d states.\n\n'], seed, first.n_cal);
fprintf(fid, '| network | condition | solid radius (1/s) | dashed radius (1/s) | diamond (1/s) | marked, in view | not drawn |\n');
fprintf(fid, '|---|---|---|---|---|---|---|\n');
for w = 1:n_net
    for c = 1:n_cond
        p = P{w, c};
        if isempty(p); continue; end
        fprintf(fid, '| %s | %s | %.3g | %.3g | %.3g | %d | %d |\n', F.networks{w}, F.cond_names{c}, p.rad, p.rad_cal, ...
            p.diamond, n_marked_drawn(w, c), n_not_drawn(w, c));
    end
end
if F.n_failed > 0
    fprintf(fid, '\n%d run(s) of these networks did not finish and are not in the numbers.\n', F.n_failed);
end
if isempty(W.far_lim) && isempty(cfg.far_lim)
    fprintf(fid, '\nNo diamond lies more than half a bulk radius beyond the main window, so no break was needed.\n');
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

function s = nice_step(x)
% 1, 2 or 5 times a power of ten, the largest not above x.
e = 10^floor(log10(x));
m = [1 2 5 10];
s = e * m(find(m * e <= x, 1, 'last'));
end

function [hidden_re, hidden_im] = covered_ticks(g, P, cfg, fs, units_per_inch)
% The ticks whose LABEL a drawn dot or a diamond of any panel would lie on
% (polish list of 2026-11-22: marked eigenvalues at the left end of the
% real axis sat on the "-30"). The label boxes are those of
% complex_plane_axes: under the real axis, centered on the tick; left of
% the imaginary axis, centered on the tick's height. Sizes from the font
% (a digit is about 0.56 em wide), a little generous, plus the marker's
% radius.
pt = units_per_inch / 72;                            % one point, in axis units
tl = 0.018 * diff(g.y_lim);                          % complex_plane_axes' tick length
x = []; y = []; r = [];
for k = 1:numel(P)
    p = P{k};
    if isempty(p); continue; end
    z = p.d.z(:); is_out = p.is_out(:);
    keep = real(z) >= g.re_lim(1) & real(z) <= last_re(g) & abs(imag(z)) <= g.y_lim(2);
    z = z(keep); is_out = is_out(keep);
    x = [x; g.x(real(z)); g.x(p.diamond)]; %#ok<AGROW>
    y = [y; imag(z); 0]; %#ok<AGROW>
    r = [r; (sqrt(cfg.dot_size) / 2) * pt * ~is_out + (sqrt(cfg.outlier_size) / 2) * pt * is_out; 4.5 * pt]; %#ok<AGROW>
end
label_w = @(t) 0.56 * fs * pt * (numel(sprintf('%g', abs(t))) + (t < 0));
label_h = 0.85 * fs * pt;
top = -1.5 * tl - 0.1 * fs * pt;                     % the glyphs begin a little under the top of the text box
hidden_re = [];
for t = g.re_ticks
    xc = g.x(t); w = label_w(t);
    if any(abs(x - xc) <= w / 2 + r & y <= top + r & y >= top - label_h - r)
        hidden_re(end + 1) = t; %#ok<AGROW>
    end
end
hidden_im = [];
for t = g.im_ticks
    w = label_w(t);
    if any(x <= -1.6 * tl + r & x >= -1.6 * tl - w - r & abs(y - t) <= label_h / 2 + r)
        hidden_im(end + 1) = t; %#ok<AGROW>
    end
end
end

function remove_ticks(cp, hidden_re, hidden_im)
% Delete the tick marks and labels complex_plane_axes drew at the hidden
% ticks: each label is a text at the tick's place, its mark the handle
% drawn just before it.
h = cp.handles;
for k = 2:numel(h)
    if ~isgraphics(h(k), 'text') || iscell(h(k).String) || ~isgraphics(h(k - 1), 'line'); continue; end
    pos = h(k).Position;
    is_re = strcmp(h(k).HorizontalAlignment, 'center') && any(abs(cp.x(hidden_re) - pos(1)) < 1e-9);
    is_im = strcmp(h(k).HorizontalAlignment, 'right') && any(abs(hidden_im - pos(2)) < 1e-9);
    if is_re || is_im; delete(h([k - 1, k])); end
end
end

function s = hidden_ticks_text(g)
% The sidecar's sentence on the ticks left out.
if isempty(g.hidden_re_ticks) && isempty(g.hidden_im_ticks)
    s = 'No tick was left out.';
    return;
end
parts = {};
if ~isempty(g.hidden_re_ticks); parts{end + 1} = ['on the real axis at ' strjoin(compose('%g', g.hidden_re_ticks), ', ') ' 1/s']; end
if ~isempty(g.hidden_im_ticks); parts{end + 1} = ['on the imaginary axis at ' strjoin(compose('%g', g.hidden_im_ticks), ', ') ' 1/s']; end
s = ['TICKS LEFT OUT, in every panel, because an eigenvalue dot or a diamond of some panel lies where the tick''s ' ...
    'label would be: ' strjoin(parts, ' and ') '. The window is unchanged.'];
end

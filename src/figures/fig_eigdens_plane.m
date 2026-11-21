function out = fig_eigdens_plane(cfg)
% FIG_EIGDENS_PLANE One frozen-gain spectrum per condition in the complex plane of rates, with the theoretical bulk.
%
%   out = FIG_EIGDENS_PLANE('run_dir', d)
%   out = FIG_EIGDENS_PLANE('run_dir', d, 'network', 'EEx1p5')
%
% Draws what the opt-in stage eigenvalue_density saved (FD-042). Reads one
% run directory or the file given; simulates nothing. Supplementary.
%
% WHAT IS DRAWN. Per condition, one panel: the spectrum of the frozen-gain
% effective connectivity M = W D at ONE state, the middle calibration state
% of the run with the lowest network seed (eigdens_run_measures'
% ev_example), as rates of the dendritic block, z = (lambda - 1) / tau_d
% (1/s). The solid gray circle is the THEORETICAL bulk, centered on
% -1/tau_d with radius R_eff / tau_d from that state's own gains. The
% eigenvalues of the zero-gain neurons sit exactly at the center, on top of
% each other. Eigenvalues farther than 1.03 times the radius from the
% center are drawn opaque and a little larger (the standing choice for the
% complex-plane figures, docs/Plan_spectral_supplement.md, section 1).
%
% AXES (complex_plane_axes; docs/Figure_style_guide.md): through the
% origin, 1:1, ONE scale and identical limits in every panel, ticks only
% outside the largest of the drawn bulk circles. The imaginary axis is the
% stability boundary. The window is the largest bulk with a margin on the
% left, top and bottom (1.35 radii from the center, so that ticks fit
% outside the bulk); on the right it reaches the rightmost eigenvalue,
% across a break in the real axis when that eigenvalue lies more than half
% a bulk radius beyond the window. Eigenvalues outside the window on the
% left, top or bottom are not drawn (with E-to-E x 1.5 and adaptation one
% real eigenvalue per state lies far to the left; no left-hand break, by
% the lead's decision): the table names each panel's clipped ones.
%
% One title line per panel, the condition's name in its color; no other
% text on the panels. Fig_Eigdens_Plane<suffix>_table.md: which run and
% state each panel is, the radius, the counts beyond 1.03 radii, the
% clipped eigenvalues and the window. Options as fig_eigdens_radial.
%
% See also: complex_plane_axes, fig_eigdens_radial, eigdens_run_measures

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.network     (1,:) char    = 'reference'
    cfg.conditions                = {}
end

OUT_FACTOR = 1.03;
setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = eigdens_figure_data(cfg, 'fig_eigdens_plane');
st = F.st;
n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];

%% The drawn state of each condition, in rates
P = cell(1, n_c);
for c = 1:n_c
    rs = F.by_cond{c};
    seeds = cat(1, rs.seed);
    [~, k] = sortrows(seeds);
    r = rs(k(1));
    P{c} = struct('z', (double(r.ev_example(:)) - 1) / r.tau_d, 'center', -1 / r.tau_d, ...
        'rad', r.R_eff_example / r.tau_d, 'seed', r.seed, 'k', r.k_example, 't', r.t(r.k_example));
end
P = [P{:}];
if any(abs([P.center] - P(1).center) > 1e-12)
    error('fig_eigdens_plane:Centers', 'The conditions have different tau_d; one window does not fit them.');
end
c0 = P(1).center;
rad_max = max([P.rad]);
z_all = vertcat(P.z);

%% One window for every panel
MARGIN = 1.35;                          % the window around the largest bulk, in radii
re_lim = [c0 - MARGIN * rad_max, max(c0 + MARGIN * rad_max, 0.15 * rad_max)];
% Room for the name of the real axis right of the largest bulk: the name sits
% under the right end of the axis, and with the window ending 0.35 radii
% past the bulk it was drawn on the bulk's edge (seen 2026-11-21). A quarter
% of the window's width is kept clear of the largest circle.
re_lim(2) = max(re_lim(2), (c0 + 1.05 * rad_max - 0.25 * re_lim(1)) / 0.75);
im_lim = [-MARGIN, MARGIN] * rad_max;
tick_step = nice_step(rad_max / 3);     % so a tick or two falls between the bulk and the window's edge
far_lim = [];
re_right = real(z_all(real(z_all) > re_lim(2) & abs(imag(z_all)) <= im_lim(2)));
if ~isempty(re_right)
    pad = 0.08 * rad_max;
    if min(re_right) - re_lim(2) > 0.5 * rad_max
        far_lim = [min(re_right) - pad, max(re_right) + pad];
        inside_gap = real(z_all) > re_lim(2) & real(z_all) < far_lim(1);
        if any(inside_gap)       % a point between the window and the far segment: no break after all
            far_lim = [];
        end
    end
    if isempty(far_lim)
        re_lim(2) = max(re_right) + pad;
    end
end
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 'far_lim', far_lim, 'bulk_center', c0, 'bulk_radius', rad_max, ...
    'tick_step', tick_step, 'font_size', st.tick_fs, 'line_width', st.axis_lw};
g = complex_plane_axes([], axes_args{:});

%% Panels placed by hand (inches): identical sizes, 1:1
fig_w = 12; left = 0.45; right = 0.15; col_gap = 0.3; top_m = 0.5; bottom = 0.2;
pw = (fig_w - left - right - (n_c - 1) * col_gap) / n_c;
ph = pw * diff(g.y_lim) / diff(g.x_lim);
fig_h = top_m + ph + bottom;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', vis);
th = linspace(0, 2 * pi, 361);
clip = cell(1, n_c); n_out = zeros(1, n_c);
for c = 1:n_c
    ax = axes(fig, 'Units', 'inches', 'Position', [left + (c - 1) * (pw + col_gap), bottom, pw, ph]);
    names = {'re_label', '', 'im_label', ''};
    if c == n_c; names{2} = 'Re (s^{-1})'; end
    if c == 1;   names{4} = 'Im (s^{-1})'; end
    cp = complex_plane_axes(ax, axes_args{:}, names{:});
    p = P(c);
    plot(ax, cp.x(p.center + p.rad * cos(th)), p.rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
    shown = real(p.z) >= re_lim(1) & abs(imag(p.z)) <= im_lim(2);
    is_out = abs(p.z - p.center) > OUT_FACTOR * p.rad;
    n_out(c) = nnz(is_out);
    zs = p.z(shown & ~is_out); zo = p.z(shown & is_out);
    scatter(ax, cp.x(real(zs)), imag(zs), 8, F.colors{c}, 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
    scatter(ax, cp.x(real(zo)), imag(zo), 16, F.colors{c}, 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
    clip{c} = p.z(~shown);
    set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);
    text(ax, mean(cp.x_lim), cp.y_lim(2), F.titles{c}, 'Color', F.colors{c}, 'FontSize', st.title_fs, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
end

fig_tag = ['Fig_Eigdens_Plane' F.suffix];
rows = cell(1, n_c);
for c = 1:n_c
    p = P(c);
    cl = 'none';
    if ~isempty(clip{c})
        cl = strjoin(arrayfun(@(z) sprintf('%.3g%+.3gi', real(z), imag(z)), clip{c}(:)', 'UniformOutput', false), ', ');
    end
    rows{c} = sprintf('| %s | %d-%d | %d (t = %.3g s) | %.3g | %d | %s |', F.titles{c}, p.seed(1), p.seed(2), ...
        p.k, p.t, p.rad, n_out(c), cl);
end
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        far_txt = 'no break';
        if ~isempty(g.far_lim); far_txt = sprintf('a break, then [%.3g, %.3g]', g.far_lim); end
        fprintf(fid, '# %s\n\n', fig_tag);
        fprintf(fid, ['One spectrum of the frozen-gain effective connectivity M = W D per condition; network `%s`; the ' ...
            'middle calibration state of the run with the lowest network seed. Drawn as rates of the dendritic block, ' ...
            '(lambda - 1)/tau_d in 1/s, so the imaginary axis is the stability boundary. Solid gray circle: the ' ...
            'theoretical bulk at that state, centered on -1/tau_d = %.3g 1/s, radius R_eff/tau_d from the state''s own ' ...
            'gains. The zero-gain neurons'' eigenvalues sit exactly at the center. Opaque, larger dots: farther than ' ...
            '%.2f radii from the center. Axes through the origin (`complex_plane_axes`), 1:1, identical in every panel: ' ...
            'real axis [%.3g, %.3g] (%s), imaginary axis [%.3g, %.3g]; ticks only outside the largest circle (radius ' ...
            '%.3g). Eigenvalues outside the window on the left, top or bottom are not drawn; they are listed below. ' ...
            'Data: `%s`.\n\n'], F.network, c0, OUT_FACTOR, g.re_lim, far_txt, g.y_lim, rad_max, F.file);
        fprintf(fid, '| Condition | seeds | state | radius (1/s) | beyond %.2f radii | not drawn (outside the window) |\n', OUT_FACTOR);
        fprintf(fid, '|---|---|---|---|---|---|\n');
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
end

%% ------------------------------------------------------------------------
function s = nice_step(x)
% 1, 2 or 5 times a power of ten, the largest not above x.
e = 10^floor(log10(x));
m = [1 2 5 10];
s = e * m(find(m * e <= x, 1, 'last'));
end

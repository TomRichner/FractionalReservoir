function out = fig_supp3_connectivity_gallery(cfg)
% FIG_SUPP3_CONNECTIVITY_GALLERY Supplementary S3: a gallery of connectivity structures, each a weight image above its equal-gain spectrum.
%
%   out = FIG_SUPP3_CONNECTIVITY_GALLERY('run_dir', d)
%   out = FIG_SUPP3_CONNECTIVITY_GALLERY('run_dir', d, 'n_cols', 4)
%
% S3 of the figure plan (future registry name fig_supp3_grouped; Plan 2,
% package B): extends Figure 1A and names the structures S7 uses. Reads the
% gallery the opt-in stage concept_examples saved (cx_load; structure
% only, nothing was simulated for it); errors when the run directory does
% not hold the stage or the stage has no gallery; simulates nothing.
%
% One tile per structure, in the stage's order, n_cols per row:
%   top     a representative 30 x 30 piece of the network's weights
%           (weight_image_sample, draw_weight_image): E columns red, I
%           columns blue, absent connections white, one color range for
%           every tile (w_clim; default the 95th percentile of the nonzero
%           |w| of all pieces)
%   bottom  the eigenvalues of J = (-I + W)/tau_d with every gain 1, in
%           1/s (complex_plane_axes: axes through the origin, 1:1, ONE
%           window and scale for every tile, the imaginary axis the
%           stability boundary, ticks outside the bulk only), the bulk
%           circle in gray, the eigenvalues clear of the bulk (farther from
%           -1/tau_d than R/tau_d + 0.2 of the largest radius) drawn
%           larger. The window holds every eigenvalue on the left (no left
%           break); eigenvalues far to the right (beyond c0 + far_factor
%           radii) are drawn beyond a break, the far segment placed for the
%           tile that has them and left empty in the others.
% The tile's title names the structure (cfg.titles, or by key).
%
% Writes Fig_Supp3_Connectivity_Gallery.{png,svg,fig} and its _table.md.
%
% See also: concept_examples_stage_plan, cx_load, fig_main1_concept_networks,
%           draw_weight_image, complex_plane_axes

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the stage's file
    cfg.titles                    = {}
    cfg.n_cols (1,1) double       = 4
    cfg.w_clim double             = []
    cfg.far_factor (1,1) double   = 2.4   % eigenvalues right of center + far_factor * largest radius go beyond the break
    cfg.tick_step double          = 100
    cfg.re_tick_candidates double = [-150 100 250]   % per tile, those outside its bulk and inside its window are drawn
    cfg.im_tick_candidates double = [-75 75]
    cfg.dot_size (1,1) double     = 5
    cfg.outlier_size (1,1) double = 30
    cfg.fig_width double          = []    % inches; [] (the default): manuscript_style's canvas_w
    cfg.image_size (1,1) double   = 1.2   % inches
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
C = cx_load(cfg.run_dir, cfg.data_file, 'fig_supp3_connectivity_gallery');
X = C.gallery;
m = numel(X);
if m < 1
    error('fig_supp3_connectivity_gallery:NoGallery', 'The stage file %s holds no gallery.', C.file);
end
st = manuscript_style();
fs = st.tick_fs;
red = st.condition_color('no_adaptation');
gray = [0.55 0.55 0.55];
titles = tile_titles(cfg.titles, X);

%% One window for every tile
c0 = X(1).center;
if any(abs([X.center] - c0) > 1e-12)
    error('fig_supp3_connectivity_gallery:Centers', 'The structures have different tau_d; one window does not fit them.');
end
rad_max = max([X.R_rate]);
far_cut = c0 + cfg.far_factor * rad_max;
z_all = vertcat(X.ev);
z_main = z_all(real(z_all) <= far_cut);
pad = 0.1 * rad_max;
re_lim = [min(min(real(z_main)), c0 - rad_max) - pad, max(max(real(z_main)), c0 + rad_max) + pad];
re_lim(2) = max(re_lim(2), 0.15 * rad_max);
if ~isempty(cfg.re_tick_candidates); re_lim(2) = max(re_lim(2), max(cfg.re_tick_candidates(cfg.re_tick_candidates <= far_cut)) + 0.05 * rad_max); end
im_lim = [-1 1] * (max(max(abs(imag(z_all))), rad_max) + pad);
far = cell(1, m);
z_far = z_all(real(z_all) > far_cut);
far_w = 0.7 * rad_max;
if ~isempty(z_far) && max(real(z_far)) - min(real(z_far)) > 0.8 * far_w
    far_w = max(real(z_far)) - min(real(z_far)) + 0.3 * rad_max;
end
for k = 1:m
    z = X(k).ev(real(X(k).ev) > far_cut);
    if isempty(z); continue; end
    mid = mean([min(real(z)), max(real(z))]);
    f1 = max(re_lim(2) + 1e-6 * rad_max, mid - far_w / 2);
    far{k} = [f1, f1 + far_w];
end
has_far = ~cellfun(@isempty, far);
proto_far = [];
if any(has_far); proto_far = far{find(has_far, 1)}; end
geom_args = {'re_lim', re_lim, 'im_lim', im_lim, 'tick_step', cfg.tick_step, 'font_size', fs, 'line_width', st.axis_lw};
if ~isempty(proto_far); geom_args = [geom_args, {'gap', 0.12 * rad_max}]; end
g = complex_plane_axes([], geom_args{:}, 'far_lim', proto_far);

%% Layout (inches)
nc = min(cfg.n_cols, m);
nr = ceil(m / nc);
W_in = cfg.fig_width;
if isempty(W_in); W_in = st.canvas_w; end
left = 0.05; right = 0; col_gap = 0.2;     % the drawing fills the canvas (the axes boxes are wider than their axis lines)
pw = (W_in - left - right - (nc - 1) * col_gap) / nc;
ph = pw * diff(g.y_lim) / diff(g.x_lim);
img = min(cfg.image_size, pw);
title_h = 0.35; gap_is = 0.12; row_gap = 0.45; top = 0.15; bottom = 0.45;
tile_h = title_h + img + gap_is + ph;
H_in = top + nr * tile_h + (nr - 1) * row_gap + bottom;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 W_in H_in], 'Visible', vis, ...
    'PaperPositionMode', 'auto');

w_clim = cfg.w_clim;
if isempty(w_clim)
    v = [];
    for k = 1:m; s = X(k).sample.W; v = [v; abs(s(s ~= 0))]; end %#ok<AGROW>
    w_clim = prctile(v, 95);
end
th = linspace(0, 2 * pi, 361);
n_out = zeros(1, m);
ax_img = gobjects(1, m); ax_spec = gobjects(1, m);
for k = 1:m
    r = ceil(k / nc); c = k - (r - 1) * nc;
    x0 = left + (c - 1) * (pw + col_gap);
    y_top = H_in - top - (r - 1) * (tile_h + row_gap) - title_h;
    ax_img(k) = axes(fig, 'Units', 'inches', 'Position', [x0 + (pw - img) / 2, y_top - img, img, img]);
    draw_weight_image(ax_img(k), X(k).sample.W, X(k).sample.type, w_clim, 'divider', true, 'line_width', st.axis_lw);
    title(ax_img(k), titles{k}, 'FontSize', fs, 'FontWeight', 'normal', 'Color', 'k');

    ax = axes(fig, 'Units', 'inches', 'Position', [x0, y_top - img - gap_is - ph, pw, ph]);
    ax_spec(k) = ax;
    args = [geom_args, {'bulk_center', X(k).center, 'bulk_radius', 1.15 * X(k).R_rate}];
    if ~isempty(proto_far)
        if has_far(k)
            args = [args, {'far_lim', far{k}}]; %#ok<AGROW>
        else
            args = [args, {'far_lim', proto_far, 'draw_far', false}]; %#ok<AGROW>
        end
    end
    names = {'re_label', '', 'im_label', ''};
    if k == 1; names{2} = 'Re (s^{-1})'; end    % the names once, on a tile with nothing at its right end
    if k == 1; names{4} = 'Im (s^{-1})'; end
    args = [args, tick_args(cfg, X(k), re_lim, far{k}, im_lim)]; %#ok<AGROW>
    cp = complex_plane_axes(ax, args{:}, names{:});
    plot(ax, cp.x(X(k).center + X(k).R_rate * cos(th)), X(k).R_rate * sin(th), '-', 'Color', gray, 'LineWidth', 1.25);
    z = X(k).ev;
    is_out = abs(z - X(k).center) > X(k).R_rate + 0.2 * rad_max;
    n_out(k) = nnz(is_out);
    if any(cp.in_break(real(z)))
        error('fig_supp3_connectivity_gallery:InBreak', 'Tile %d (%s) has an eigenvalue in the removed stretch.', k, X(k).key);
    end
    scatter(ax, cp.x(real(z(~is_out))), imag(z(~is_out)), cfg.dot_size, red, 'filled', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.75);
    scatter(ax, cp.x(real(z(is_out))), imag(z(is_out)), cfg.outlier_size, red, 'filled', 'MarkerEdgeColor', 'none');
    set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);
end

fig_tag = 'Fig_Supp3_Connectivity_Gallery';
out = struct('figs', fig, 'files', {{}}, 'source', C.file, 'window', g, 'far', {far}, 'w_clim', w_clim, ...
    'n_out', n_out, 'titles', {titles}, 'axes', struct('image', ax_img, 'spectrum', ax_spec));
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        write_table(fid, fig_tag, C, X, titles, g, w_clim, rad_max);
        clear closer
        out.files{end + 1} = table_name;
    end
end
end

%% ------------------------------------------------------------------------
function t = tile_titles(given, X)
% Given, or by key (cx_gallery_names, which S7 panel B uses too).
m = numel(X);
if ~isempty(given)
    if numel(given) ~= m
        error('fig_supp3_connectivity_gallery:Titles', 'titles has %d entries for %d structures.', numel(given), m);
    end
    t = cellstr(given);
    return
end
t = cx_gallery_names(X);
end

function write_table(fid, fig_tag, C, X, titles, g, w_clim, rad_max)
fprintf(fid, '# %s\n\n', fig_tag);
fprintf(fid, 'Source: `%s` (stage concept_examples, gallery; structure only). Preset `%s`, condition `%s`.\n\n', ...
    C.file, C.settings.preset_name, C.settings.condition);
fprintf(fid, ['Top of each tile: a %d x %d piece of the weights (submatrix on half E, half I neurons, image seed %d), ' ...
    'color range +-%.4g. Bottom: the eigenvalues of J = (-I + W)/tau_d with every gain 1 (1/s); gray circle radius ' ...
    'R/tau_d about -1/tau_d; larger dots: farther from the center than R/tau_d + %.3g 1/s. Window Re [%g, %g], ' ...
    'Im [%g, %g] 1/s in every tile.\n\n'], C.settings.n_image, C.settings.n_image, C.settings.image_seed, w_clim, ...
    0.2 * rad_max, g.re_lim, g.im_lim);
fprintf(fid, '| Structure | Network | Level of chaos | Weight seed | R/tau_d (1/s) | Eigenvalues clear of the bulk (1/s) | Added mode (1/s) | Block means x n_pre: K (EE, EI, IE, II) |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(X)
    x = X(k);
    z = x.ev(abs(x.ev - x.center) > x.R_rate + 0.2 * rad_max);
    z = z(imag(z) >= 0);
    zs = strjoin(arrayfun(@(v) sprintf('%.1f%+.1fi', real(v), imag(v)), z(:)', 'UniformOutput', false), ', ');
    if isempty(zs); zs = 'none'; end
    lr = '';
    if ~isnan(x.lowrank_rate); lr = sprintf('%.1f%+.1fi', real(x.lowrank_rate), imag(x.lowrank_rate)); end
    fprintf(fid, '| %s | `%s` (%s) | %g | [%d %d] | %.2f | %s | %s | %.2f, %.2f, %.2f, %.2f |\n', titles{k}, x.name, ...
        x.label, x.level_of_chaos, x.seed, x.R_rate, zs, lr, x.K_struct(1, 1), x.K_struct(1, 2), x.K_struct(2, 1), x.K_struct(2, 2));
end
fprintf(fid, ['\nOnly the upper half-plane member of a conjugate pair is listed. K is (post, pre), the structural mean ' ...
    'matrix of ospec_weight_stats (realized block mean, zeros included, times the presynaptic count).\n']);
end

function a = tick_args(cfg, x, re_lim, far_k, im_lim)
% The candidate ticks that fall outside this tile's bulk (1.15 radii, as
% its clear interval) and inside its window; [] leaves complex_plane_axes'
% multiples. A tick whose LABEL an eigenvalue clear of the bulk would sit on
% is not drawn in that tile. The label of a tick on the imaginary axis is
% left of the axis and the label of a tick on the real axis below it, so the
% spot is a box on that side of the tick, in units of the tile's radius
% (the fonts are fixed and the window is about six radii wide). Until
% 2026-11-22 only eigenvalues within 0.35 radii of the imaginary axis
% counted, on either side: a conjugate pair at -21 +- 81i 1/s (0.44 radii
% left, where the label is) sat on the '75' and '-75' of its tile.
a = {};
R = x.R_rate;
clear_re = x.center + [-1 1] * 1.15 * R;
clear_im = [-1 1] * 1.15 * R;
marked = x.ev(abs(x.ev - x.center) > 1.05 * R);
if ~isempty(cfg.re_tick_candidates)
    t = cfg.re_tick_candidates;
    ok = (t <= clear_re(1) | t >= clear_re(2)) & t >= re_lim(1) & t <= re_lim(2) & t ~= 0;
    if ~isempty(far_k); ok = (ok & t <= re_lim(2) - 0.4 * R) | (t >= far_k(1) & t <= far_k(2)); end   % no main tick against the break
    for j = 1:numel(t)       % the label: below the axis, about a radius wide
        if any(abs(real(marked) - t(j)) < 0.6 * R & imag(marked) < -0.04 * R & imag(marked) > -0.75 * R); ok(j) = false; end
    end
    a = [a, {'re_ticks', t(ok)}];
end
if ~isempty(cfg.im_tick_candidates)
    t = cfg.im_tick_candidates;
    ok = (t <= clear_im(1) | t >= clear_im(2)) & t >= im_lim(1) & t <= im_lim(2) & t ~= 0;
    for j = 1:numel(t)       % the label: left of the axis, about a radius wide
        if any(real(marked) > -1.1 * R & real(marked) < 0.3 * R & abs(imag(marked) - t(j)) < 0.3 * R); ok(j) = false; end
    end
    a = [a, {'im_ticks', t(ok)}];
end
end

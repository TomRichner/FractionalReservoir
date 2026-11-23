function out = fig_main1_concept_networks(cfg)
% FIG_MAIN1_CONCEPT_NETWORKS Figure 1: four example networks of the paper's ensemble (weights, equal-gain spectrum, activity), and the adaptation diagrams.
%
%   out = FIG_MAIN1_CONCEPT_NETWORKS('run_dir', d)
%   out = FIG_MAIN1_CONCEPT_NETWORKS('run_dir', d, 'preset_name', p, 't_window', [0 6])
%
% The new Figure 1 of the figure plan (registry name fig_main1_grouped;
% Plan 2, package B). Reads what the opt-in stage concept_examples saved
% (cx_load); errors when the run directory does not hold it; simulates
% nothing.
%
% PANEL A, four columns (the stage's examples, in order) by three rows:
%   row 1  a representative 30 x 30 piece of the network's weights
%          (weight_image_sample, draw_weight_image): E columns red, I
%          columns blue, absent connections white, ONE color range for all
%          four (w_clim; default the 95th percentile of the nonzero |w| of
%          the four pieces), a thin gray line between the E and I columns
%   row 2  the eigenvalues of the Jacobian with every gain equal, as rates
%          (lambda(W) - 1)/tau_d in 1/s, through complex_plane_axes: axes
%          through the origin, 1:1, ONE scale in all four panels, the
%          imaginary axis the stability boundary, the bulk circle (radius
%          R/tau_d about -1/tau_d; gray), ticks outside the bulk only. An
%          eigenvalue right of the main window is drawn beyond a break in
%          the real axis; the far segment has the same width in every
%          panel and is placed per column around that column's far
%          eigenvalues (columns without one leave the room empty).
%          Eigenvalues LEFT of the window are not drawn (with the E-to-E
%          mean raised by half, one real eigenvalue near -141 1/s); the
%          table lists them, in a sentence of their own. Dots in the
%          no-adaptation color (these are unadapted networks): an eigenvalue
%          RIGHT of the boundary (Re > 0) filled at full strength, one left
%          of it or on it in a lighter tint (left_tint), drawn first, so
%          which eigenvalues have crossed is seen at a glance (polish list,
%          2026-11-22; until then all were one red). An eigenvalue clear of
%          the bulk (farther from -1/tau_d than R/tau_d + 0.2 of the largest
%          radius) is drawn larger.
%   row 3  the dendritic potential x of n_E E and n_I I neurons of the
%          unadapted run over t_window, in the E and I colors; one y range
%          for the columns whose traces fit within 1.5 times the median span,
%          a range of its own for a column that does not (the runaway); a
%          time scale bar under the first column
% Column headings (cfg.titles) above the first row.
%
% PANEL B, the adaptation diagrams of the present Figure 1C, unchanged
% (draw_adaptation_concepts): SFA shifts the response curve, depression
% scales it; their product is the neuron's gain, which scales its column of
% the weights.
%
% Writes Fig_Main1_Concept_Networks.{png,svg,fig} and
% Fig_Main1_Concept_Networks_table.md (the numbers behind every column).
%
% See also: concept_examples_stage_plan, run_concept_examples_stage, cx_load,
%           fig_supp3_connectivity_gallery, draw_adaptation_concepts,
%           complex_plane_axes, draw_weight_image

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''    % the stage's concept_examples.mat, instead of run_dir
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % panel B's depression endpoint; '' = the stage's preset
    cfg.titles                    = {}    % one per column; {} = by key (below), else the network's label
    cfg.t_window double           = [0 6]
    cfg.n_E (1,1) double          = 4
    cfg.n_I (1,1) double          = 4
    cfg.w_clim double             = []
    cfg.re_lim double             = []    % [lo hi] 1/s; [] = fitted to the bulks
    cfg.im_lim double             = []
    cfg.far_width double          = []    % width of the far segment, 1/s; [] = 0.6 of the largest bulk radius
    cfg.tick_step double          = 50
    cfg.dot_size (1,1) double     = 6
    cfg.outlier_size (1,1) double = 36
    cfg.left_tint (1,1) double {mustBeGreaterThanOrEqual(cfg.left_tint, 0), mustBeLessThanOrEqual(cfg.left_tint, 1)} = 0.6   % share of white in the dots left of the boundary; 0 = one color
    cfg.scale_bar_s (1,1) double  = 2
    cfg.fig_size double           = []   % [width height] inches; [] (the default): manuscript_style's canvas_w by 9.6
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
C = cx_load(cfg.run_dir, cfg.data_file, 'fig_main1_concept_networks');
X = C.examples;
m = numel(X);
if m < 1
    error('fig_main1_concept_networks:NoExamples', 'The stage file %s holds no examples.', C.file);
end
preset_name = cfg.preset_name;
if isempty(preset_name); preset_name = C.settings.preset_name; end
st = manuscript_style();
fs = st.tick_fs;
red = st.condition_color('no_adaptation');
gray = [0.55 0.55 0.55];
pale = red + cfg.left_tint * (1 - red);     % eigenvalues left of the boundary
titles = column_titles(cfg.titles, X);

%% The spectrum window: one geometry for every column
c0 = X(1).center;
if any(abs([X.center] - c0) > 1e-12)
    error('fig_main1_concept_networks:Centers', 'The examples have different tau_d; one window does not fit them.');
end
rad_max = max([X.R_rate]);
re_lim = cfg.re_lim;
if isempty(re_lim)
    hi = c0 + 1.08 * rad_max;
    for k = 1:m                              % eigenvalues just outside the bulk widen the window
        z = X(k).ev(real(X(k).ev) > hi & real(X(k).ev) <= c0 + 1.6 * rad_max);
        if ~isempty(z); hi = max(hi, max(real(z)) + 0.08 * rad_max); end
    end
    re_lim = [c0 - 1.15 * rad_max, max(hi, 0.15 * rad_max)];
end
im_lim = cfg.im_lim;
if isempty(im_lim)
    im_lim = [-1 1] * 1.12 * rad_max;
    for k = 1:m
        z = X(k).ev(real(X(k).ev) >= re_lim(1));
        im_lim = [-1 1] * max(im_lim(2), max(abs(imag(z))) * 1.06);
    end
end
far_w = cfg.far_width;
if isempty(far_w); far_w = 0.6 * rad_max; end
gap = 0.12 * rad_max;
far = cell(1, m);
for k = 1:m
    z = X(k).ev(real(X(k).ev) > re_lim(2));
    if isempty(z); continue; end
    span = [min(real(z)), max(real(z))];
    if diff(span) > 0.8 * far_w
        error('fig_main1_concept_networks:FarSpread', ...
            'Column %d has eigenvalues from %.3g to %.3g 1/s right of the window; the far segment is %.3g wide (far_width).', ...
            k, span(1), span(2), far_w);
    end
    mid = mean(span);
    f1 = max(re_lim(2) + 1e-6 * rad_max, mid - far_w / 2);
    far{k} = [f1, f1 + far_w];
end
has_far = ~cellfun(@isempty, far);
if any(has_far)
    proto_far = far{find(has_far, 1)};
else
    proto_far = [];
end
geom_args = {'re_lim', re_lim, 'im_lim', im_lim, 'tick_step', cfg.tick_step, 'font_size', fs, 'line_width', st.axis_lw};
if ~isempty(proto_far); geom_args = [geom_args, {'gap', gap}]; end
g = complex_plane_axes([], geom_args{:}, 'far_lim', proto_far);

%% Layout (inches)
fig_size = cfg.fig_size;
if isempty(fig_size); fig_size = [st.canvas_w, 9.6]; end
W_in = fig_size(1); H_in = fig_size(2);
left = 0.75; right = 0.2; col_gap = 0.3;
pw = (W_in - left - right - (m - 1) * col_gap) / m;
ph_spec = pw * diff(g.y_lim) / diff(g.x_lim);
img = min(1.75, pw);
ph_act = 1.45;
top = 0.55;
y_img = H_in - top - img;
y_spec = y_img - 0.28 - ph_spec;
y_act = y_spec - 0.42 - ph_act;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 W_in H_in], 'Visible', vis, ...
    'PaperPositionMode', 'auto');
x0 = @(k) left + (k - 1) * (pw + col_gap);

%% Row 1: the weight images
w_clim = cfg.w_clim;
if isempty(w_clim)
    v = [];
    for k = 1:m; s = X(k).sample.W; v = [v; abs(s(s ~= 0))]; end %#ok<AGROW>
    w_clim = prctile(v, 95);
end
ax_img = gobjects(1, m);
for k = 1:m
    ax_img(k) = axes(fig, 'Units', 'inches', 'Position', [x0(k) + (pw - img) / 2, y_img, img, img]);
    draw_weight_image(ax_img(k), X(k).sample.W, X(k).sample.type, w_clim, 'divider', true, 'line_width', st.axis_lw);
    title(ax_img(k), titles{k}, 'FontSize', fs, 'FontWeight', 'normal', 'Color', 'k');
end

%% Row 2: the equal-gain spectra
th = linspace(0, 2 * pi, 361);
n_left = zeros(1, m); n_far = zeros(1, m); n_crossed = zeros(1, m);
ax_spec = gobjects(1, m);
for k = 1:m
    ax = axes(fig, 'Units', 'inches', 'Position', [x0(k), y_spec, pw, ph_spec]);
    ax_spec(k) = ax;
    args = [geom_args, {'bulk_center', X(k).center, 'bulk_radius', 1.15 * X(k).R_rate}];   % tick labels clear of the circle
    if ~isempty(proto_far)
        if has_far(k)
            args = [args, {'far_lim', far{k}}]; %#ok<AGROW>
        else
            args = [args, {'far_lim', proto_far, 'draw_far', false}]; %#ok<AGROW>
        end
    end
    names = {'re_label', '', 'im_label', ''};
    if k == m; names{2} = 'Re (s^{-1})'; end
    if k == 1; names{4} = 'Im (s^{-1})'; end
    args = [args, {'re_ticks', column_re_ticks(X(k), re_lim, far{k}, cfg.tick_step, gap)}]; %#ok<AGROW>
    cp = complex_plane_axes(ax, args{:}, names{:});
    plot(ax, cp.x(X(k).center + X(k).R_rate * cos(th)), X(k).R_rate * sin(th), '-', 'Color', gray, 'LineWidth', 1.25);
    z = X(k).ev;
    shown = real(z) >= re_lim(1) & abs(imag(z)) <= im_lim(2) & ~cp.in_break(real(z));
    if ~isempty(proto_far) && has_far(k)
        shown = shown & real(z) <= far{k}(2);
    elseif ~isempty(proto_far)
        shown = shown & real(z) <= re_lim(2);
    end
    n_left(k) = nnz(real(z) < re_lim(1));
    n_far(k) = nnz(real(z) > re_lim(2) & shown);
    if nnz(~shown) > n_left(k)
        error('fig_main1_concept_networks:NotShown', ...
            'Column %d: %d eigenvalue(s) fall outside the window other than on its left; widen im_lim or far_width.', ...
            k, nnz(~shown) - n_left(k));
    end
    is_out = abs(z - X(k).center) > X(k).R_rate + 0.2 * rad_max;   % clear of the bulk, not its finite-size edge
    % Left of the boundary first and pale, right of it on top at full strength.
    crossed = real(z) > 0;
    n_crossed(k) = nnz(crossed);
    sides = {~crossed, pale, 'eig_left'; crossed, red, 'eig_crossed'};
    for q = 1:2
        zs = z(shown & ~is_out & sides{q, 1}); zo = z(shown & is_out & sides{q, 1});
        scatter(ax, cp.x(real(zs)), imag(zs), cfg.dot_size, sides{q, 2}, 'filled', 'MarkerEdgeColor', 'none', ...
            'Tag', sides{q, 3});
        scatter(ax, cp.x(real(zo)), imag(zo), cfg.outlier_size, sides{q, 2}, 'filled', 'MarkerEdgeColor', 'none', ...
            'Tag', [sides{q, 3} '_outlier']);
    end
    set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);
end

%% Row 3: activity
nE = min(cfg.n_E, nnz([X(1).trace_type] == 1)); nI = min(cfg.n_I, nnz([X(1).trace_type] == 2));
colE = excitatory_colormap(max(nE, 1)); colI = inhibitory_colormap(max(nI, 1));
spans = zeros(m, 2);
sel = cell(1, m);
for k = 1:m
    tt = X(k).t;
    in_w = tt >= cfg.t_window(1) & tt <= cfg.t_window(2);
    if ~any(in_w)
        error('fig_main1_concept_networks:Window', 't_window [%g %g] s lies outside the saved run (0 to %g s).', ...
            cfg.t_window, tt(end));
    end
    iE = find(X(k).trace_type == 1, nE); iI = find(X(k).trace_type == 2, nI);
    sel{k} = struct('in_w', in_w, 'iE', iE, 'iI', iI);
    v = X(k).x_trace([iE; iI], in_w);
    spans(k, :) = [min(v(:)), max(v(:))];
end
width = diff(spans, 1, 2);
own = width > 1.5 * median(width);
shared = [min(spans(~own, 1)), max(spans(~own, 2))];
ax_act = gobjects(1, m);
y_lims = zeros(m, 2);
for k = 1:m
    ax = axes(fig, 'Units', 'inches', 'Position', [x0(k) + 0.25, y_act, pw - 0.25, ph_act]);
    ax_act(k) = ax;
    hold(ax, 'on');
    s = sel{k};
    tt = X(k).t(s.in_w);
    for j = 1:numel(s.iE)
        plot(ax, tt, X(k).x_trace(s.iE(j), s.in_w), '-', 'Color', colE(j, :), 'LineWidth', 0.8);
    end
    for j = 1:numel(s.iI)
        plot(ax, tt, X(k).x_trace(s.iI(j), s.in_w), '-', 'Color', colI(j, :), 'LineWidth', 0.8);
    end
    if own(k); yl = spans(k, :); else; yl = shared; end
    yl = nice_lim(yl);
    y_lims(k, :) = yl;
    set(ax, 'XLim', cfg.t_window, 'YLim', yl, 'FontSize', fs, 'LineWidth', st.axis_lw, 'Box', 'off', ...
        'TickDir', 'out', 'YTick', nice_ticks(yl));
    ax.XAxis.Visible = 'off';
    if k == 1
        ylabel(ax, 'x_i', 'FontSize', st.label_fs, 'Interpreter', 'tex');
        % the time scale bar
        xb = cfg.t_window(1) + [0 cfg.scale_bar_s];
        yb = yl(1) - 0.08 * diff(yl);
        plot(ax, xb, [yb yb], 'k-', 'LineWidth', 2.5, 'Clipping', 'off');
        text(ax, mean(xb), yb - 0.04 * diff(yl), sprintf('%g s', cfg.scale_bar_s), 'FontSize', fs, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
end

%% Panel B: the adaptation diagrams (unchanged drawing)
yB = 0.62; hB = 1.5;
in2n = @(x, y, w, h) [x / W_in, y / H_in, w / W_in, h / H_in];
% Left edges and widths as shares of the canvas width (laid out on 13 in).
xB = [0.75 3.75 7.0 10.1] / 13 * W_in; wB = [2.0 2.3 2.0 2.3] / 13 * W_in;
posB = [in2n(xB(1), yB, wB(1), hB); in2n(xB(2), yB, wB(2), hB); in2n(xB(3), yB, wB(3), hB); in2n(xB(4), yB, wB(4), hB)];
A = draw_adaptation_concepts(fig, preset_name, posB);
set(findall(A.axes, '-property', 'FontSize'), 'FontSize', fs);
for k = 1:numel(A.axes); set([A.axes(k).XLabel, A.axes(k).YLabel], 'FontSize', st.label_fs); end   % axis names at label_fs

%% Panel letters
lab = axes(fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off', 'Tag', 'panel_letters');
text(lab, 0.08 / W_in, 1 - 0.12 / H_in, '(A)', 'FontSize', st.letter_fs, 'VerticalAlignment', 'top', 'Tag', 'letter_A');
text(lab, 0.08 / W_in, (yB + hB + 0.55) / H_in, '(B)', 'FontSize', st.letter_fs, 'VerticalAlignment', 'top', 'Tag', 'letter_B');

%% Save, and the numbers
fig_tag = 'Fig_Main1_Concept_Networks';
out = struct('figs', fig, 'files', {{}}, 'source', C.file, 'window', g, 'far', {far}, 'w_clim', w_clim, ...
    'y_lims', y_lims, 'n_left', n_left, 'n_far', n_far, 'n_crossed', n_crossed, 'titles', {titles}, ...
    'axes', struct('image', ax_img, 'spectrum', ax_spec, 'activity', ax_act, 'concepts', A.axes));
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        write_table(fid, fig_tag, C, X, titles, g, far, w_clim, n_left, n_far, n_crossed, cfg, A);
        clear closer
        out.files{end + 1} = table_name;
    end
end
end

%% ------------------------------------------------------------------------
function t = column_titles(given, X)
DEFAULT = struct('decay', 'Balanced, weak', 'irregular', 'Balanced, strong', ...
    'runaway', ['E-to-E mean ' char(215) '1.5'], 'oscillation', 'Added rotating mode');   % the caption says it is an added structure
m = numel(X);
if ~isempty(given)
    if numel(given) ~= m
        error('fig_main1_concept_networks:Titles', 'titles has %d entries for %d columns.', numel(given), m);
    end
    t = cellstr(given);
    return
end
t = cell(1, m);
for k = 1:m
    if isfield(DEFAULT, X(k).key); t{k} = DEFAULT.(X(k).key); else; t{k} = X(k).label; end
end
end

function yl = nice_lim(yl)
step = nice_step(diff(yl) / 4);
yl = [floor(yl(1) / step) * step, ceil(yl(2) / step) * step];
if diff(yl) == 0; yl = yl + [-1 1]; end
end

function t = nice_ticks(yl)
step = nice_step(diff(yl) / 3);
t = (ceil(yl(1) / step):floor(yl(2) / step)) * step;
end

function s = nice_step(raw)
p = 10 ^ floor(log10(raw));
q = raw / p;
if q >= 5; s = 5 * p; elseif q >= 2; s = 2 * p; else; s = p; end
end

function write_table(fid, fig_tag, C, X, titles, g, far, w_clim, n_left, n_far, n_crossed, cfg, A)
fprintf(fid, '# %s\n\n', fig_tag);
fprintf(fid, 'Source: `%s` (stage concept_examples). Preset `%s`, condition `%s` (no adaptation).\n\n', ...
    C.file, C.settings.preset_name, C.settings.condition);
fprintf(fid, ['Panel A. Row 1: a %d x %d piece of each network''s weights (the submatrix on %d E and %d I ' ...
    'neurons drawn with image seed %d), color range +-%.4g (95th percentile of the nonzero |w| of the four ' ...
    'pieces unless given). Row 2: the eigenvalues of J = (-I + W)/tau_d with every gain 1, in 1/s; gray circle ' ...
    'radius R/tau_d about -1/tau_d, R the structural bulk radius of the drawn matrix. Window Re [%g, %g], ' ...
    'Im [%g, %g] 1/s, one scale in all columns. An eigenvalue right of the stability boundary (Re > 0) is filled ' ...
    'at full strength, one left of it in a lighter tint. Row 3: dendritic potential x of %d E and %d I neurons of the ' ...
    'unadapted run (noise of the preset), t in [%g, %g] s, from the class''s initial state.\n\n'], ...
    C.settings.n_image, C.settings.n_image, C.settings.n_image / 2, C.settings.n_image / 2, C.settings.image_seed, ...
    w_clim, g.re_lim, g.im_lim, cfg.n_E, cfg.n_I, cfg.t_window);
fprintf(fid, ['| Column | Network | Level of chaos | Seed | n | R/tau_d (1/s) | Rightmost (1/s) | Added mode (1/s) | ' ...
    'Beyond the break | Left of the window | Density EE, EI, IE, II | Nonzero mean EE, EI, IE, II |\n']);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(X)
    x = X(k);
    lr = '';
    if ~isnan(x.lowrank_rate); lr = sprintf('%.1f%+.1fi', real(x.lowrank_rate), imag(x.lowrank_rate)); end
    far_s = '';
    if ~isempty(far{k}); far_s = sprintf('%d in [%.0f, %.0f]', n_far(k), far{k}); end
    lefts = '0';
    if n_left(k) > 0
        zl = x.ev(real(x.ev) < g.re_lim(1));
        lefts = sprintf('%d (%s)', n_left(k), strjoin(arrayfun(@(z) sprintf('%.1f%+.1fi', real(z), imag(z)), zl(:)', ...
            'UniformOutput', false), ', '));
    end
    fprintf(fid, '| %s | `%s` (%s) | %g | [%d %d] | %d | %.2f | %.1f%+.1fi | %s | %s | %s | %.3f, %.3f, %.3f, %.3f | %.4f, %.4f, %.4f, %.4f |\n', ...
        strrep(titles{k}, '\times', 'x'), x.name, x.label, x.level_of_chaos, x.seed, x.n, x.R_rate, real(x.rightmost), ...
        abs(imag(x.rightmost)), lr, far_s, lefts, x.density(1, 1), x.density(1, 2), x.density(2, 1), x.density(2, 2), ...
        x.nz_mean(1, 1), x.nz_mean(1, 2), x.nz_mean(2, 1), x.nz_mean(2, 2));
end
fprintf(fid, '\nDensity and nonzero means are (post, pre): EI is I-to-E. Of the whole network, not of the piece shown.\n\n');
fprintf(fid, 'Eigenvalues right of the boundary (Re > 0, drawn filled): %s.\n\n', ...
    strjoin(arrayfun(@(k) sprintf('%s %d of %d', strrep(titles{k}, '\times', 'x'), n_crossed(k), numel(X(k).ev)), ...
    1:numel(X), 'UniformOutput', false), '; '));
if any(n_left > 0)
    fprintf(fid, '**Not drawn (left of the window, Re < %.1f 1/s):**', g.re_lim(1));
    for k = find(n_left > 0)
        zl = X(k).ev(real(X(k).ev) < g.re_lim(1));
        fprintf(fid, ' column "%s": %d eigenvalue(s), at %s 1/s.', strrep(titles{k}, '\times', 'x'), n_left(k), ...
            strjoin(arrayfun(@(z) sprintf('%.1f%+.1fi', real(z), imag(z)), zl(:)', 'UniformOutput', false), ', '));
    end
    fprintf(fid, [' The window keeps one scale for all the columns, and widening it to hold them would shrink ' ...
        'every bulk. They lie far left of the boundary.\n\n']);
else
    fprintf(fid, 'Every eigenvalue is drawn (none lies left of the window).\n\n');
end
fprintf(fid, '## The unadapted runs\n\n');
fprintf(fid, '| Column | T (s) | Mean rate E | Mean rate I | At the floor | At the ceiling | SD of x | lambda_1 (1/s) | over (s) |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(X)
    x = X(k);
    fprintf(fid, '| %s | %g | %.3f | %.3f | %.2f | %.2f | %.3f | %+.2f | [%g, %g] |\n', strrep(titles{k}, '\times', 'x'), ...
        x.t(end), x.mean_rate_E, x.mean_rate_I, x.frac_floor, x.frac_ceiling, x.x_sd, x.LLE, x.lya_interval);
end
fprintf(fid, ['\nMeans, floor and ceiling fractions (rate below 1e-6 or above 1 - 1e-6, over neurons and samples) ' ...
    'and the SD of x (per neuron over time, averaged) are over the second half of the run. lambda_1: the ' ...
    'class''s Benettin estimate.\n\n']);
fprintf(fid, '## The weight pieces\n\n| Column | Density of the piece | Block means of the piece (EE, EI, IE, II) |\n|---|---|---|\n');
for k = 1:numel(X)
    s = X(k).sample;
    fprintf(fid, '| %s | %.3f | %.4f, %.4f, %.4f, %.4f |\n', strrep(titles{k}, '\times', 'x'), s.density_all, ...
        s.block_mean(1, 1), s.block_mean(1, 2), s.block_mean(2, 1), s.block_mean(2, 2));
end
fprintf(fid, '\nPanel B: the adaptation diagrams of the present Figure 1C (draw_adaptation_concepts). %s\n', A.std_note);
end

function t = column_re_ticks(x, re_lim, far_k, step, gap)
% Multiples of step outside the column's bulk (1.15 radii) in the main
% window, and in its far segment; a main tick too close to the first far
% tick in plotted distance gives way to it, so the far value is labeled.
clear_re = x.center + [-1 1] * 1.15 * x.R_rate;
t = (ceil(re_lim(1) / step):floor(re_lim(2) / step)) * step;
t = t(t ~= 0 & (t <= clear_re(1) | t >= clear_re(2)));
if isempty(far_k); return; end
tf = (ceil(far_k(1) / step):floor(far_k(2) / step)) * step;
if isempty(tf)
    tf = round(mean(far_k) / 10) * 10;      % a round value inside a segment narrower than one step
end
x_far = re_lim(2) + gap + (tf(1) - far_k(1));
t = t(x_far - t >= 0.6 * step);
t = [t, tf];
end

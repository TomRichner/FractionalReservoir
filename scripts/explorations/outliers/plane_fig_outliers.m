function file = plane_fig_outliers(data_dir, fig_dir, opts)
%PLANE_FIG_OUTLIERS The complex-plane outlier comparison on shared axes through the origin, with a break in the real axis.
%
%   file = plane_fig_outliers()
%   file = plane_fig_outliers(data_dir, fig_dir, 'seed', 7)
%
% The six panels of fig_outliers' outliers_complex_plane (rows: networks,
% columns: conditions), redrawn in the style Tom asked for on 2026-11-20
% with complex_plane_axes: every panel has the SAME limits and the same
% length per unit of Re and Im; the real and imaginary axes go through the
% origin (the imaginary axis is the stability boundary, so there is no
% dashed zero line); ticks only outside the largest bulk; and the real axis
% is broken so that the uniform-gain diamond of the E-to-E times 1.5
% network (+32 to +39 1/s) is shown beyond a gap instead of setting the
% scale of its row. Writes outliers_complex_plane_v2.png and its sidecar
% .md; fig_outliers and its PNG are left as they are, for comparison.
%
% Same content as fig_outliers, section 1: the spectrum dots, the gray bulk
% circle, the ringed rightmost eigenvalues with their conjugates, the
% uniform-gain diamond. cond_color, net_name and with_conjugates are COPIES
% of fig_outliers' local functions (that file is not edited from here; no
% adaptation is dark red, FR-004).
%
% THE WINDOW is fixed by the options, not fitted to the data, so the
% figure cannot rescale itself between runs: main window re_lim, im_lim,
% and far_lim beyond the break. Anything in the removed stretch of the
% real axis is an ERROR (complex_plane_axes:InBreak): widen the window or
% move the break. Eigenvalues outside the window on the other sides are
% clipped like on any plot, and COUNTED in the sidecar. The ticks are the
% same in all six panels: they clear the largest of the six bulk circles.
%
% This reads the runs and calls outlier_measures on the six it draws; it
% does not call outlier_table, which rewrites outlier_measures.csv. Nothing
% in data_dir is written. The figure is made invisible and closed.
%
% See also: complex_plane_axes, fig_outliers, outlier_measures

arguments
    data_dir (1,:) char = ''
    fig_dir  (1,:) char = ''
    opts.seed double = []
    opts.re_lim (1,2) double = [-32 12]
    opts.im_lim (1,2) double = [-22 22]
    opts.far_lim (1,2) double = [28 42]
    opts.gap (1,1) double = 3
    opts.tick_step (1,1) double = 10
    opts.resolution (1,1) double = 200
    opts.name (1,:) char = 'outliers_complex_plane_v2'
end

root = fileparts(which('setup_paths'));
if isempty(data_dir); data_dir = fullfile(root, 'data', 'explorations', 'outliers'); end
if isempty(fig_dir);  fig_dir  = fullfile(root, 'figs', 'explorations', 'outliers'); end
if ~isfolder(fig_dir); mkdir(fig_dir); end

%% The runs: networks (reference first), the paper's condition order, one seed
files = dir(fullfile(data_dir, 'run_*.mat'));
tok = regexp({files.name}, '^run_(.+)_(no_adaptation|sfa1_std1|sfa3_stdc2)_s(\d+)-(\d+)\.mat$', 'tokens', 'once');
tok = tok(~cellfun(@isempty, tok));
if isempty(tok); error('plane_fig_outliers:NoRuns', 'No run_*.mat in %s.', data_dir); end
net_of = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
cond_of = cellfun(@(t) t{2}, tok, 'UniformOutput', false);
seed_of = cellfun(@(t) str2double(t{3}), tok);
cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
conds = cond_order(ismember(cond_order, cond_of));
nets = unique(net_of, 'stable'); nets = nets(:);
if any(strcmp(nets, 'reference'))
    nets = [{'reference'}; setdiff(nets, {'reference'}, 'stable')];
end
if isempty(opts.seed); seed1 = max(seed_of); else; seed1 = opts.seed(1); end
n_net = numel(nets); n_cond = numel(conds);

P = cell(n_net, n_cond);                                % what each panel draws, in data units
rows = {};
for w = 1:n_net
    for c = 1:n_cond
        f = fullfile(data_dir, outlier_run_name(nets{w}, conds{c}, [seed1 seed1 + 1], ''));
        if ~isfile(f); continue; end
        L = load(f, 'R'); R = L.R;
        m = outlier_measures(R);
        tau = R.tau_d;
        p.ev = ([R.cal.ev_M] - 1) ./ tau;               % n x n_cal, rates
        p.center = -1 / tau;
        p.rad = m.kappa * median(R.R_eff) / tau;
        top = nan(1, size(p.ev, 2)) + 0i;
        for q = 1:size(p.ev, 2)
            [~, i1] = max(real(p.ev(:, q))); top(q) = p.ev(i1, q);
        end
        p.top = with_conjugates(top);
        p.u = (median([R.cal.uniform_re]) - 1) / tau;
        P{w, c} = p;
        rows{end + 1} = m.row; %#ok<AGROW>
    end
end
drawn = P(~cellfun(@isempty, P));
if isempty(drawn); error('plane_fig_outliers:NoSeed', 'No runs of seed %d in %s.', seed1, data_dir); end
rad_max = max(cellfun(@(p) p.rad, drawn));
centers = unique(cellfun(@(p) p.center, drawn));
if numel(centers) > 1
    error('plane_fig_outliers:Centers', 'The bulk centers differ between runs; one shared clear interval does not describe them.');
end

%% Layout: identical panels, placed by hand (inches)
st = manuscript_style();
fs = st.tick_fs - 2;
gray = [0.55 0.55 0.55];
axes_args = {'re_lim', opts.re_lim, 'im_lim', opts.im_lim, 'far_lim', opts.far_lim, 'gap', opts.gap, ...
    'tick_step', opts.tick_step, 'bulk_center', centers, 'bulk_radius', rad_max, ...
    'font_size', fs, 'line_width', st.axis_lw};
g = complex_plane_axes([], axes_args{:});               % the geometry, before any axes exist
fig_w = 12; left = 0.55; right = 0.10; col_gap = 0.25; top_m = 0.45; row_gap = 0.30; bottom = 0.15;
pw = (fig_w - left - right - (n_cond - 1) * col_gap) / n_cond;
ph = pw * diff(g.y_lim) / diff(g.x_lim);                % 1:1
fig_h = top_m + n_net * ph + (n_net - 1) * row_gap + bottom;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', 'off', 'Tag', 'plane_fig_outliers');
closer = onCleanup(@() close(fig));

n_clip = zeros(n_net, n_cond);
clip_re = cell(n_net, n_cond);                          % [min max] of the real parts of the clipped dots
for w = 1:n_net
    for c = 1:n_cond
        ax = axes(fig, 'Units', 'inches', 'Position', ...
            [left + (c - 1) * (pw + col_gap), bottom + (n_net - w) * (ph + row_gap), pw, ph]);
        names = {'re_label', '', 'im_label', ''};         % the names once per row and per column, like the row labels
        if w == n_net; names{2} = 'Re (s^{-1})'; end
        if c == 1;     names{4} = 'Im (s^{-1})'; end
        cp = complex_plane_axes(ax, axes_args{:}, names{:});
        p = P{w, c};
        if ~isempty(p)
            col = cond_color(st, conds{c});
            th = linspace(0, 2 * pi, 200);
            plot(ax, cp.x(p.center + p.rad * cos(th)), p.rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
            scatter(ax, cp.x(real(p.ev(:))), imag(p.ev(:)), 5, col, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
            plot(ax, cp.x(real(p.top)), imag(p.top), 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
            plot(ax, cp.x(p.u), 0, 'd', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
            k = cp.check(real(p.ev(:)));
            n_clip(w, c) = k.n_left + k.n_right + nnz(imag(p.ev(:)) < g.y_lim(1) | imag(p.ev(:)) > g.y_lim(2));
            out = real(p.ev(:)) < g.re_lim(1) | real(p.ev(:)) > g.far_lim(2) | abs(imag(p.ev(:))) > g.y_lim(2);
            if any(out); clip_re{w, c} = [min(real(p.ev(out))), max(real(p.ev(out)))]; end
            marks = [p.top, p.u];                       % the marks that carry the point must all be in view
            km = cp.check(real(marks));
            if km.n_left + km.n_right > 0 || any(abs(imag(marks)) > g.y_lim(2))
                error('plane_fig_outliers:MarkOutside', 'A ringed eigenvalue or a diamond of %s, %s is outside the window.', nets{w}, conds{c});
            end
        end
        set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);    % plotting must not have moved them
        if w == 1
            title(ax, st.condition_title(conds{c}), 'FontSize', st.title_fs, 'FontWeight', 'normal');
        end
        if c == 1
            set(ax.YLabel, 'String', net_name(nets{w}), 'FontSize', st.label_fs - 2, 'Visible', 'on', 'Color', 'k');
        end
    end
end

file = fullfile(fig_dir, [opts.name '.png']);
exportgraphics(fig, file, 'Resolution', opts.resolution);
fprintf('[plane_fig_outliers] wrote %s\n', file);

%% Sidecar
T = struct2table([rows{:}], 'AsArray', true);
if ischar(T.network);   T.network = {T.network};     end
if ischar(T.condition); T.condition = {T.condition}; end
cols = {'lambda_1', 'kappa', 'agree_re', 'R_eff_med', 'rate_edge_med', 'rate_top_med', 'ratio_med', 'rho_uniform', ...
    'reined', 'frac_detached', 'K1_med', 'K_err_med', 'coh1_med', 'frac_on_E', 'frac_on_I', 'gain_ratio_EI', ...
    'n_unstable_mean', 'shuffle_pct_med'};
fid = fopen(fullfile(fig_dir, [opts.name '.md']), 'w');
fid_closer = onCleanup(@() fclose(fid));
fprintf(fid, '# %s\n\n', opts.name);
fprintf(fid, ['The same six panels as `outliers_complex_plane` (rows: networks, columns: conditions), on axes that are ' ...
    'IDENTICAL in every panel and drawn through the origin (`complex_plane_axes`). One unit of Re is as long as one unit ' ...
    'of Im everywhere, and as long in one panel as in another. The vertical axis is the imaginary axis, Re = 0: the ' ...
    'stability boundary (there is no separate dashed line). THE REAL AXIS IS BROKEN: the stretch from %g to %g 1/s is ' ...
    'left out (the gap with a slanted stroke on each side), and the segment to its right runs from %g to %g 1/s at the ' ...
    'same scale. The break is drawn in every panel; only the E-to-E times 1.5 row has anything beyond it. Ticks are ' ...
    'placed only outside the largest of the six bulk circles (radius %.3g 1/s, the reference network without ' ...
    'adaptation), so they are the same in every panel: Re at %s 1/s, Im at %s 1/s.\n\n'], ...
    g.re_lim(2), g.far_lim(1), g.far_lim(1), g.far_lim(2), rad_max, list(g.re_ticks), list(g.im_ticks));
fprintf(fid, ['Dots: every eigenvalue of the frozen-gain effective connectivity (each neuron''s outgoing weights times ' ...
    'its current gain, slow variables held), drawn as a rate (lambda - 1)/tau_d, at 12 states spread over 20 s after ' ...
    '60 s of settling, overlaid. Gray circle: the edge of the bulk, centered on -1/tau_d = %g 1/s: the radius ' ...
    'predicted from the gains (median over the window) times the null-calibrated factor kappa. Ringed dots: the ' ...
    'rightmost eigenvalue of each state; when it is complex it is one of a conjugate pair and both members are ringed. ' ...
    'Open diamond: where the rightmost eigenvalue would be if every neuron had the same gain (the rms of the observed ' ...
    'gains): the outlier the weights alone make. In the reference network the diamond is inside the circle and the ' ...
    'ringed dots are the edge of the bulk. With E-to-E times 1.5 the diamond is beyond the break (at %s 1/s in the ' ...
    'three conditions) and the actual rightmost eigenvalue is back at the bulk. Seed %d.\n\n'], ...
    centers, list(cellfun(@(p) p.u, P(end, ~cellfun(@isempty, P(end, :))))), seed1);
clip_txt = cell(n_net, n_cond);
for k = 1:numel(clip_txt)
    if isempty(clip_re{k}); clip_txt{k} = '0'; else; clip_txt{k} = sprintf('%d (Re from %.3g to %.3g)', n_clip(k), clip_re{k}); end
end
fprintf(fid, ['Eigenvalue dots outside the window (clipped, not drawn), per panel in reading order: %s. No ringed ' ...
    'eigenvalue and no diamond is outside it, and nothing lies in the removed stretch (either would be an error).\n\n'], ...
    strjoin(reshape(clip_txt.', 1, []), '; '));
fprintf(fid, ['Rates are in 1/s, eigenvalues of the effective connectivity in units of W. ' ...
    'Definitions: `scripts/explorations/outliers/outlier_measures.m`.\n\n']);
head = [{'network', 'condition', 'seed'}, cols];
fprintf(fid, '| %s |\n', strjoin(head, ' | '));
fprintf(fid, '|%s\n', repmat('---|', 1, numel(head)));
for k = 1:height(T)
    cells = cellfun(@(q) sprintf('%.3g', T.(q)(k)), cols, 'UniformOutput', false);
    fprintf(fid, '| %s | %s | %d | %s |\n', T.network{k}, T.condition{k}, T.seed(k), strjoin(cells, ' | '));
end
end

%% ------------------------------------------------------------------------
function s = list(v)
s = strjoin(arrayfun(@(x) sprintf('%.3g', x), v, 'UniformOutput', false), ', ');
end

function col = cond_color(st, name)
% Copy of fig_outliers' cond_color (FR-004: no adaptation in dark red).
if strcmp(name, 'no_adaptation')
    col = [0.75 0.10 0.15];
else
    col = st.condition_color(name);
end
end

function name = net_name(net)
% Copy of fig_outliers' net_name.
tok = regexp(net, '^EEx(\d+)p(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    name = sprintf('E-to-E \\times %s.%s', tok{1}, tok{2});
elseif strcmp(net, 'EI3to2')
    name = 'E:I 3:2';
else
    name = net;
end
end

function z = with_conjugates(z)
% Copy of fig_outliers' with_conjugates: the conjugate of each complex
% eigenvalue appended; a real one is kept once, on the real axis.
z = z(:).';
is_real = abs(imag(z)) <= 1e-9 * max(abs(z), 1);
z(is_real) = real(z(is_real));
z = [z, conj(z(~is_real))];
end

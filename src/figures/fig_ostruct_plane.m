function out = fig_ostruct_plane(cfg)
% FIG_OSTRUCT_PLANE The intended mode of each network of one group in the complex plane of rates: rows networks, columns conditions.
%
%   out = FIG_OSTRUCT_PLANE('run_dir', d, 'group', 'rotation')     'rotation' | 'local' | 'blockmean'
%
% Draws what the opt-in stage outlier_structures saved (FD-043). Reads one
% run directory or the file given; simulates nothing. Supplementary.
%
% WHAT IS DRAWN, for ONE SEED (the first realization of every network: the
% first seed row, or a fixed-weight network's own noise), per network
% (row) and condition (column), as rates of the dendritic block,
% z = (lambda - 1)/tau_d in 1/s, so the imaginary axis is the stability
% boundary:
%   * the whole spectrum of the frozen-gain effective connectivity M at the
%     calibration states, overlaid: small faint dots in the condition color;
%   * those of its eigenvalues farther than 1.03 times their OWN state's
%     theoretical radius from the center: the same dots, less faint
%     (decided after looking at the tiny run: drawn in full they compete
%     with the mode; left out, a panel hides that the bulk has other
%     outliers);
%   * the solid gray circle: the theoretical bulk, centered on -1/tau_d,
%     radius the median over the calibration states of R_eff/tau_d;
%   * the INTENDED mode at each calibration state (the eigenvalue of M that
%     carries the structural mode's eigenvector), with its conjugate when
%     complex: opaque dots in the condition color, slightly larger, no ring;
%   * the open black diamond(s): the same mode with every gain set to the
%     rms gain (median over the window), what the weights alone give.
%
% AXES (complex_plane_axes; docs/Figure_style_guide.md): through the
% origin, 1:1, ONE scale and identical limits in every panel of the
% figure, ticks only outside the largest drawn circle. The window covers
% the largest bulk with a margin; on the right it reaches every drawn
% point, across a break in the real axis when the points beyond the bulk
% leave a stretch wider than half a bulk radius empty (typically a
% uniform-gain diamond far out). NOTHING IS DRAWN LEFT OF THE WINDOW: a
% population mode whose carriers are pinned goes far to the left (with
% E-to-E x 1.5 and adaptation, -59 to -42 1/s at 500 neurons); there is no
% left-hand break (lead's decision, docs/Plan_spectral_supplement.md), and
% the sidecar names and counts such modes and eigenvalues per panel.
%
% One title line per column (the condition, in its color) and one name per
% row (the network, at the left); nothing else on the panels.
% Fig_Ostruct_Plane_<group>_table.md: the caption, the window and the
% break, the counts per panel (mode points left of the window, other
% eigenvalues outside it, outliers beyond 1.03 radii), and the numbers over
% ALL seeds (median and range) with the structural check.
%
% Options: group (required), run_dir, data_file, out_dir, save, visible,
% verbose, preset_name (unused), conditions.
%
% See also: complex_plane_axes, ostruct_figure_data, fig_ostruct_summary

arguments
    cfg.group       (1,:) char {mustBeMember(cfg.group, {'rotation', 'local', 'blockmean'})} = 'rotation'
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
end

OUT_FACTOR = 1.03;
MARGIN = 1.35;                           % the window around the largest bulk, in radii
setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = ostruct_figure_data(cfg, 'fig_ostruct_plane');
st = F.st;
g_idx = find(strcmp(F.group_names, cfg.group));
w_idx = find(F.group == g_idx);
if isempty(w_idx)
    error('fig_ostruct_plane:NoNetwork', 'The stage in %s has no finished run of a %s network.', F.folder, cfg.group);
end
n_r = numel(w_idx); n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];

%% The drawn seed: the first realization of each network
first = min(F.T.realization(ismember(F.T.key, F.keys(w_idx))));
D = ostruct_load(F.folder, 'keep_runs', true, 'measures', false, 'keys', F.keys(w_idx), 'realization', first);
P = cell(n_r, n_c);
for i = 1:n_r
    for c = 1:n_c
        k = find(strcmp(D.runs.key, F.keys(w_idx(i))) & strcmp(D.runs.condition, F.cond_names{c}) & ...
            D.runs.ok & D.runs.realization == first);
        if isempty(k); continue; end
        P{i, c} = panel_data(D.Rs{D.idx == k(1)}, F.cond_names{c}, OUT_FACTOR);
    end
end
have = ~cellfun(@isempty, P);
PP = [P{have}];
c0 = PP(1).center;
if any(abs([PP.center] - c0) > 1e-12)
    error('fig_ostruct_plane:Centers', 'The runs have different tau_d; one window does not fit them.');
end
rad_max = max([PP.rad]);

%% One window for every panel
re_lim = [c0 - MARGIN * rad_max, max(c0 + MARGIN * rad_max, 0.15 * rad_max)];
z_marks = [vertcat(PP.mode); vertcat(PP.diamond)];
z_all = [vertcat(PP.z); z_marks];
im_lim = [-1 1] * max(MARGIN * rad_max, 1.08 * max(abs(imag(z_marks(real(z_marks) >= re_lim(1))))));
% Ticks come only outside the bulk: widen the window to the left and in
% height so that at least one multiple of the step lies outside it.
% A step near two thirds of the bulk radius: with half of it (5 1/s at 500
% neurons) the labels 10, 15, 20, 25 ran together and under the diamonds.
tick_step = nice_step(rad_max / 1.5);
im_lim = [-1 1] * max(im_lim(2), (floor(rad_max / tick_step) + 1) * tick_step * 1.02);
re_lim(1) = min(re_lim(1), floor((c0 - rad_max) / tick_step) * tick_step - 0.02 * tick_step);
pad = 0.08 * rad_max;
right = sort(real(z_all(real(z_all) > re_lim(2) & abs(imag(z_all)) <= im_lim(2))));
far_lim = [];
if ~isempty(right)
    edges = [re_lim(2); right(:)];
    [gap, at] = max(diff(edges));
    if gap > 0.5 * rad_max
        re_lim(2) = max(re_lim(2), edges(at) + pad);
        % The far segment: room for its marks and a tick or two, at least
        % 0.4 bulk radii wide, the points clear of the break marks.
        far_lim = [edges(at + 1) - 0.15 * rad_max, right(end) + pad];
        far_lim(2) = max(far_lim(2), far_lim(1) + 0.4 * rad_max);
        t_far = ceil(far_lim(1) / tick_step) * tick_step;   % a tick in the far segment says where it is
        if t_far > far_lim(2); far_lim(2) = t_far + 0.05 * tick_step; end
        if far_lim(1) <= re_lim(2) + pad          % too close to make a break worth it
            re_lim(2) = far_lim(2); far_lim = [];
        end
    else
        re_lim(2) = right(end) + pad;
    end
end
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 'far_lim', far_lim, 'bulk_center', c0, 'bulk_radius', rad_max, ...
    'tick_step', tick_step, 'font_size', st.tick_fs - 2, 'line_width', st.axis_lw};
g = complex_plane_axes([], axes_args{:});

%% Panels placed by hand (inches): identical sizes, 1:1
fig_w = 12; left = 0.75; right_m = 0.15; col_gap = 0.3; row_gap = 0.35; top_m = 0.45; bottom = 0.3;
pw = (fig_w - left - right_m - (n_c - 1) * col_gap) / n_c;
ph = pw * diff(g.y_lim) / diff(g.x_lim);
fig_h = top_m + n_r * ph + (n_r - 1) * row_gap + bottom;
vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', vis);
th = linspace(0, 2 * pi, 361);
counts = cell(n_r, n_c);
for i = 1:n_r
    for c = 1:n_c
        y0 = fig_h - top_m - i * ph - (i - 1) * row_gap;
        ax = axes(fig, 'Units', 'inches', 'Position', [left + (c - 1) * (pw + col_gap), y0, pw, ph]);
        names = {'re_label', '', 'im_label', ''};
        if c == n_c && i == n_r; names{2} = 'Re (s^{-1})'; end
        if c == 1 && i == 1;     names{4} = 'Im (s^{-1})'; end
        cp = complex_plane_axes(ax, axes_args{:}, names{:});
        set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);
        if c == 1
            text(ax, cp.x_lim(1) - 0.04 * diff(cp.x_lim), 0, F.short{w_idx(i)}, 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', st.tick_fs - 2, 'Interpreter', 'tex');
        end
        if i == 1
            text(ax, mean(cp.x_lim), cp.y_lim(2), F.titles{c}, 'Color', F.colors{c}, 'FontSize', st.title_fs, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
        end
        p = P{i, c};
        if isempty(p); continue; end
        col = F.colors{c};
        in_win = @(z) real(z) >= re_lim(1) & real(z) <= last_re(g) & abs(imag(z)) <= im_lim(2);
        plot(ax, cp.x(p.center + p.rad * cos(th)), p.rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
        zb = p.z(~p.is_out & in_win(p.z)); zo = p.z(p.is_out & in_win(p.z));
        scatter(ax, cp.x(real(zb)), imag(zb), 5, col, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
        scatter(ax, cp.x(real(zo)), imag(zo), 9, col, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
        zm = p.mode(in_win(p.mode));
        plot(ax, cp.x(real(zm)), imag(zm), 'o', 'MarkerSize', 5.5, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
        zd = p.diamond(in_win(p.diamond));
        plot(ax, cp.x(real(zd)), imag(zd), 'd', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
        counts{i, c} = struct('mode_left', nnz(real(p.mode_states) < re_lim(1)), 'n_states', numel(p.mode_states), ...
            'mode_left_re', real(p.mode_states(real(p.mode_states) < re_lim(1))), ...
            'mode_out', nnz(~in_win(p.mode)), 'eig_out', nnz(~in_win(p.z)), 'n_eig', numel(p.z), ...
            'n_beyond', nnz(p.is_out), 'diamond_out', nnz(~in_win(p.diamond)), 'seed', p.seed, 'rad', p.rad);
    end
end

%% Save, and the sidecar
fig_tag = ['Fig_Ostruct_Plane_' cfg.group];
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    group_text = struct( ...
        'rotation', ['Low-rank structure ADDED to the weights on random patterns of neurons (dense, ignoring Dale''s ' ...
        'law): a real eigenvalue (rank one) and conjugate pairs of increasing imaginary part (rank two) at about 1.4 ' ...
        'structural bulk radii, then a real one and a pair at about 2.5 (row names: the structural eigenvalue in bulk ' ...
        'radii, median over seeds).'], ...
        'local', ['Networks with NO structure put in: weight draws chosen for a local outlier of the balanced sparse ' ...
        'Dale matrix (real, reference network on weight seed 27; complex, indegree 200 on weight seed 25), and the ' ...
        'reference network on new weights, whose intended mode is simply the rightmost eigenvalue of W, a bulk-edge mode.'], ...
        'blockmean', ['Outliers made by the block means: single routes changed (real outliers) and the E-I loop ' ...
        'strengthened together with E-to-E (a complex pair).']);
    far_txt = 'no break';
    if ~isempty(g.far_lim); far_txt = sprintf('a break, then [%.3g, %.3g]', g.far_lim); end
    caption = sprintf(['%s Rows: networks; columns: conditions. One seed drawn (realization %d: the first seed row, ' ...
        'or a fixed-weight network''s own noise); the numbers below are over all seeds. Drawn as rates, ' ...
        '(lambda - 1)/tau_d in 1/s, so the imaginary axis is the stability boundary. Faint dots: every eigenvalue of ' ...
        'the frozen-gain effective connectivity (each neuron''s outgoing weights times its current gain, slow ' ...
        'variables held) at the %d calibration states, overlaid; less faint: those farther than %.2f times their own ' ...
        'state''s theoretical radius from the center. Solid gray circle: the theoretical bulk (median radius over the ' ...
        'states), centered on -1/tau_d = %.3g 1/s. Opaque dots: the INTENDED mode at each calibration state (the ' ...
        'eigenvalue whose eigenvector overlaps most with the structural mode''s; mode_ov is the median overlap, and ' ...
        'where it is low the dot marks only the best match), with its conjugate when complex. Open diamonds: the same ' ...
        'mode with every gain set to the rms gain, what the weights alone give. Axes through the origin ' ...
        '(complex_plane_axes), 1:1, identical in every panel: real axis [%.3g, %.3g] (%s), imaginary axis [%.3g, %.3g]; ' ...
        'ticks only outside the largest circle (radius %.3g 1/s). Nothing is drawn left of the window: see the counts.'], ...
        group_text.(cfg.group), first, PP(1).n_cal, OUT_FACTOR, c0, g.re_lim, far_txt, g.y_lim, rad_max);
    notes = {count_table(counts, F, w_idx, re_lim)};
    cols = {'lam0_re_R', 'lam0_im_R', 'coh0', 'cond0', 'lambda_1', 'mode_rho', 'mode_rho_uniform', 'mode_reined', ...
        'mode_ov', 'mode_frac_rightmost', 'rate_struct', 'rate_frozen', 'rate_full', 'rot_struct', 'rot_frozen', ...
        'rot_full', 'frac_on_E', 'frac_on_I', 'on_top', 'on_rest', 'n_out103_mean'};
    name = ostruct_write_table(out_dir, fig_tag, caption, F, cols, w_idx, notes);
    if ~isempty(name); out.files{end + 1} = name; end
end
vprintf(cfg.verbose, 'verbose', '%s: %d networks x %d conditions\n', fig_tag, n_r, n_c);
end

%% ------------------------------------------------------------------------
function p = panel_data(R, condition, out_factor)
% What one panel draws, in rates.
tau = R.tau_d;
z = (R.cal.ev_M - 1) ./ tau;                       % n x n_cal
rad_k = R.cal.R_eff ./ tau;
center = -1 / tau;
is_out = abs(z - center) > out_factor .* rad_k;
if isempty(R.dense)
    error('fig_ostruct_plane:NoDense', 'Run %s | %s has no dense pass: the intended mode at the calibration states is not saved.', ...
        R.network, condition);
end
m = ([R.dense.mode_lam_M].' - 1) ./ tau;          % Im >= 0, one per state
m_all = [m; conj(m(abs(imag(m)) > 1e-9))];
u = (R.ref.lam0 * median(R.mode.d_rms) - 1) / tau;
u_all = u;
if abs(imag(u)) > 1e-9; u_all = [u; conj(u)]; end
p = struct('z', z(:), 'is_out', is_out(:), 'rad', median(rad_k), 'center', center, ...
    'mode', m_all, 'mode_states', m, 'diamond', u_all, 'seed', R.seed, 'n_cal', numel(R.cal.t));
end

function hi = last_re(g)
if isempty(g.far_lim); hi = g.re_lim(2); else; hi = g.far_lim(2); end
end

function s = count_table(counts, F, w_idx, re_lim)
% The counts per panel of the drawn seed, as a markdown table.
lines = {sprintf(['Counts for the drawn seed, per panel. "Mode left of the window": calibration states at which the ' ...
    'intended mode lies left of Re = %.3g 1/s (a population mode whose carriers are pinned), with its real parts; ' ...
    '"not drawn": eigenvalues and mode points outside the window on any side; "beyond 1.03 R": eigenvalues of M ' ...
    'farther than 1.03 times their state''s theoretical radius from the center (over all calibration states).'], re_lim(1)), ...
    '', '| network | condition | seeds | mode left of the window | its real parts (1/s) | mode points not drawn | eigenvalues not drawn | beyond 1.03 R | diamonds not drawn |', ...
    '|---|---|---|---|---|---|---|---|---|'};
for i = 1:size(counts, 1)
    for c = 1:size(counts, 2)
        k = counts{i, c};
        if isempty(k); continue; end
        re = 'none';
        if ~isempty(k.mode_left_re)
            re = sprintf('%.3g to %.3g', min(k.mode_left_re), max(k.mode_left_re));
        end
        lines{end + 1} = sprintf('| %s | %s | %d-%d | %d of %d | %s | %d | %d of %d | %d | %d |', F.keys{w_idx(i)}, ...
            F.cond_names{c}, k.seed(1), k.seed(2), k.mode_left, k.n_states, re, k.mode_out, k.eig_out, k.n_eig, ...
            k.n_beyond, k.diamond_out); %#ok<AGROW>
    end
end
s = strjoin(lines, newline);
end

function s = nice_step(x)
% 1, 2 or 5 times a power of ten, the largest not above x.
e = 10^floor(log10(x));
m = [1 2 5 10];
s = e * m(find(m * e <= x, 1, 'last'));
end

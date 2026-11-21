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
% scale of its row. Writes outliers_complex_plane_v3.png and its sidecar
% .md by default; outliers_complex_plane_v2 (commit 91ec441) is reproduced
% exactly by plane_fig_outliers('', '', 'seed', 7, 'ring_rightmost', true),
% with the other options below off. fig_outliers and its PNG are left as
% they are, for comparison.
%
% Same content as fig_outliers, section 1: the spectrum dots, the gray bulk
% circle, the ringed rightmost eigenvalues with their conjugates, the
% uniform-gain diamond (the ringed dots only with 'ring_rightmost', true:
% Tom, 2026-11-20, finds them distracting and too like Figure 3's Ritz
% values, so they are OFF by default and the rightmost eigenvalues are then
% ordinary dots). cond_color, net_name and with_conjugates are COPIES
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
% TWO CIRCLES ('theory_circle', true; off by default; Tom, 2026-11-20:
% scaling the circle up by kappa makes close outliers look like bulk). The
% solid gray circle is then the theoretical bulk, radius R_eff from the
% gains with no kappa, and the calibrated edge kappa R_eff, the only circle
% otherwise, is kept as a thin dashed lighter one. The sidecar gains the
% detachment of the rightmost eigenvalue against each (computed here, not
% in outlier_measures). Ticks still clear the largest calibrated circle.
%
% MARKED OUTLIERS ('mark_outliers', true; off by default; Tom, 2026-11-20,
% instead of the rings). Every eigenvalue beyond outlier_factor (1.03:
% "Let's just mark those that are beyond 103% of the theoretical radius")
% times the theoretical circle AS DRAWN (distance from -1/tau_d against the
% median R_eff, any direction; the circle itself stays at 1.00 and none is
% drawn at 1.03) is drawn opaque and a little larger (outlier_size), in the
% condition color, with no edge; the sidecar counts them, against the
% run's own null spectra too.
%
% THE INSET ('inset', true; off by default, and then the figure is v2
% exactly; with either option on, the default name is outliers_complex_plane_v3). Tom,
% 2026-11-20: in the lower right of each panel, the model's firing-rate
% nonlinearity with one dot per neuron per state at the neuron's operating
% point, colored by cell type (the single E and I colors of
% excitatory_colormap(1) and inhibitory_colormap(1), a magenta and a teal,
% as in fig_main2_dynamics' legend), for the same 12 states
% whose spectra the panel overlays, so that one sees how many neurons are
% on the sloped part and how many on the flat parts.
%   horizontal  u = x - (c/K) sum_k a_k - S_c_i: the argument the
%               nonlinearity receives (dendritic state minus the adaptation
%               current), relative to the neuron's own setpoint S_c_i, so
%               that all neurons lie on ONE curve, the piecewise sigmoid
%               centered on 0 (SRNNCellTypePairs.piecewiseSigmoid(u, S_a,
%               0); the class's curve is a pure translation in S_c)
%   vertical    the rate r = phi_0(u), from 0 (floor) to 1 (ceiling), plus
%               a uniform jitter of +/- inset_jitter so that piles can be
%               seen (fixed seed inset_seed, a private stream). With
%               'inset_split' the jitter is one-sided by type instead:
%               excitatory neurons from the curve up to + inset_jitter,
%               inhibitory ones from the curve down to - inset_jitter, so
%               that a pile of both does not blend into one color.
%   window      u in [-inset_u_lim, inset_u_lim], the same in all panels.
%               A neuron beyond it is NOT dropped: it is drawn in a strip
%               past a gap in the floor (left) or ceiling (right) line, at
%               a uniformly random horizontal place in the strip (its
%               distance is not shown, only that it is beyond). Counted in
%               the sidecar.
%   no text. The curve is drawn over the dots; its flat parts are the floor
%               and the ceiling.
% The states come from plane_states (a re-simulation of each saved run,
% verified against it; files in states_dir). This function simulates
% nothing: a missing states file is an error that names the call to make.
% 'inset_snapshots' (indices into the 12 states), 'inset_n_neurons' (a
% random subset, the same neurons in every state, drawn per panel in equal
% proportion from each type), 'inset_size' and 'inset_alpha' give variants;
% the sidecar's fractions are always over all neurons and all 12 states,
% and it says what is drawn.
%
% See also: complex_plane_axes, fig_outliers, outlier_measures, plane_states

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
    opts.name (1,:) char = ''
    opts.inset (1,1) logical = false
    opts.theory_circle (1,1) logical = false
    opts.ring_rightmost (1,1) logical = false
    opts.mark_outliers (1,1) logical = false
    opts.outlier_factor (1,1) double = 1.03   % Tom, 2026-11-20: "Let's just mark those that are beyond 103% of the theoretical radius."
    opts.outlier_size (1,1) double = 10
    opts.states_dir (1,:) char = ''
    opts.inset_u_lim (1,1) double = 1.5
    opts.inset_jitter (1,1) double = 0.12
    opts.inset_split (1,1) logical = false
    opts.inset_seed (1,1) double = 1
    opts.inset_snapshots double = []
    opts.inset_n_neurons (1,1) double = Inf
    opts.inset_size (1,1) double = 1.5
    opts.inset_alpha (1,1) double = 0.2
    opts.inset_re_left (1,1) double = 8.5
    opts.inset_im (1,2) double = [-22 -10.5]
end
if isempty(opts.name)
    if opts.ring_rightmost && ~opts.inset && ~opts.theory_circle && ~opts.mark_outliers
        opts.name = 'outliers_complex_plane_v2';
    else
        opts.name = 'outliers_complex_plane_v3';
    end
end

root = fileparts(which('setup_paths'));
if isempty(data_dir); data_dir = fullfile(root, 'data', 'explorations', 'outliers'); end
if isempty(fig_dir);  fig_dir  = fullfile(root, 'figs', 'explorations', 'outliers'); end
if ~isfolder(fig_dir); mkdir(fig_dir); end
states_dir = opts.states_dir;
if isempty(states_dir); states_dir = fullfile(root, 'data', 'explorations', 'outliers_plane'); end

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
        p.rad_theory = median(R.R_eff) / tau;           % no kappa: the radius the gains predict
        p.det = detachment(R, m, opts.outlier_factor);
        top = nan(1, size(p.ev, 2)) + 0i;
        for q = 1:size(p.ev, 2)
            [~, i1] = max(real(p.ev(:, q))); top(q) = p.ev(i1, q);
        end
        p.top = with_conjugates(top);
        p.u = (median([R.cal.uniform_re]) - 1) / tau;
        if opts.inset
            p.inset = load_states(states_dir, nets{w}, conds{c}, [seed1 seed1 + 1], R, m.row);
        end
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

ei_col = [excitatory_colormap(1); inhibitory_colormap(1)];   % the established E and I colors
rs = RandStream('mt19937ar', 'Seed', opts.inset_seed);  % private: the jitter does not touch the global stream
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
            if opts.theory_circle                       % solid: the theoretical bulk; thin dashed, lighter: the calibrated edge
                plot(ax, cp.x(p.center + p.rad * cos(th)), p.rad * sin(th), '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.75);
                plot(ax, cp.x(p.center + p.rad_theory * cos(th)), p.rad_theory * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
            else
                plot(ax, cp.x(p.center + p.rad * cos(th)), p.rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
            end
            scatter(ax, cp.x(real(p.ev(:))), imag(p.ev(:)), 5, col, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
            if opts.mark_outliers                       % beyond outlier_factor times the solid circle AS DRAWN: opaque, a little larger, no edge
                is_out = abs(p.ev(:) - p.center) > opts.outlier_factor * p.rad_theory;
                scatter(ax, cp.x(real(p.ev(is_out))), imag(p.ev(is_out)), opts.outlier_size, col, 'filled', ...
                    'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
            end
            if opts.ring_rightmost
                plot(ax, cp.x(real(p.top)), imag(p.top), 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
            end
            plot(ax, cp.x(p.u), 0, 'd', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
            k = cp.check(real(p.ev(:)));
            n_clip(w, c) = k.n_left + k.n_right + nnz(imag(p.ev(:)) < g.y_lim(1) | imag(p.ev(:)) > g.y_lim(2));
            out = real(p.ev(:)) < g.re_lim(1) | real(p.ev(:)) > g.far_lim(2) | abs(imag(p.ev(:))) > g.y_lim(2);
            if any(out); clip_re{w, c} = [min(real(p.ev(out))), max(real(p.ev(out)))]; end
            marks = [p.top, p.u];                       % the marks that carry the point must all be in view
            km = cp.check(real(marks));
            if km.n_left + km.n_right > 0 || any(abs(imag(marks)) > g.y_lim(2))
                error('plane_fig_outliers:MarkOutside', 'A rightmost eigenvalue or a diamond of %s, %s is outside the window.', nets{w}, conds{c});
            end
        end
        set(ax, 'XLim', cp.x_lim, 'YLim', cp.y_lim);    % plotting must not have moved them
        if opts.inset && ~isempty(p)
            % The inset's rectangle, given in the panel's data units, in inches.
            per = ax.Position(3) / diff(cp.x_lim);      % inches per unit (1:1, so the same for Im)
            rect = [ax.Position(1) + (cp.x(opts.inset_re_left) - cp.x_lim(1)) * per, ...
                ax.Position(2) + (opts.inset_im(1) - cp.y_lim(1)) * per, ...
                (cp.x_lim(2) - cp.x(opts.inset_re_left)) * per, diff(opts.inset_im) * per];
            draw_inset(fig, rect, p.inset, ei_col, rs, opts, st.axis_lw);
        end
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
if opts.theory_circle
    circle_txt = sprintf(['Two circles (Tom, 2026-11-20), both centered on -1/tau_d = %g 1/s. SOLID gray: the bulk, the ' ...
        'radius R_eff predicted from the gains alone (median over the window), with no calibration factor. THIN DASHED, ' ...
        'lighter: that radius times kappa (the 95th percentile of the rightmost eigenvalue of the null draws, relative ' ...
        'to R_eff): how far the outermost eigenvalue of a structureless matrix with these gains reaches by chance; it ' ...
        'was the only circle of `outliers_complex_plane_v2`. An eigenvalue between the two is outside the theoretical ' ...
        'bulk but within what chance reaches. '], centers);
else
    circle_txt = sprintf(['Gray circle: the edge of the bulk, centered on -1/tau_d = %g 1/s: the radius ' ...
        'predicted from the gains (median over the window) times the null-calibrated factor kappa. '], centers);
end
if opts.ring_rightmost
    ring_txt = ['Ringed dots: the ' ...
        'rightmost eigenvalue of each state; when it is complex it is one of a conjugate pair and both members are ringed. '];
    ring_ref = 'ringed dots are';
else
    ring_txt = ['The rightmost eigenvalue of each state is not marked: it is an ordinary dot like the others (the ' ...
        'option `ring_rightmost` draws it large and ringed, as in `outliers_complex_plane_v2`). '];
    ring_ref = 'rightmost eigenvalues are at';
end
fprintf(fid, ['Dots: every eigenvalue of the frozen-gain effective connectivity (each neuron''s outgoing weights times ' ...
    'its current gain, slow variables held), drawn as a rate (lambda - 1)/tau_d, at 12 states spread over 20 s after ' ...
    '60 s of settling, overlaid. %s%s' ...
    'Open diamond: where the rightmost eigenvalue would be if every neuron had the same gain (the rms of the observed ' ...
    'gains): the outlier the weights alone make. In the reference network the diamond is inside the circle and the ' ...
    '%s the edge of the bulk. With E-to-E times 1.5 the diamond is beyond the break (at %s 1/s in the ' ...
    'three conditions) and the actual rightmost eigenvalue is back at the bulk. Seed %d.\n\n'], ...
    circle_txt, ring_txt, ring_ref, list(cellfun(@(p) p.u, P(end, ~cellfun(@isempty, P(end, :))))), seed1);
clip_txt = cell(n_net, n_cond);
if opts.ring_rightmost; clip_what = 'ringed'; else; clip_what = 'rightmost'; end
for k = 1:numel(clip_txt)
    if isempty(clip_re{k}); clip_txt{k} = '0'; else; clip_txt{k} = sprintf('%d (Re from %.3g to %.3g)', n_clip(k), clip_re{k}); end
end
fprintf(fid, ['Eigenvalue dots outside the window (clipped, not drawn), per panel in reading order: %s. No %s ' ...
    'eigenvalue and no diamond is outside it, and nothing lies in the removed stretch (either would be an error).\n\n'], ...
    strjoin(reshape(clip_txt.', 1, []), '; '), clip_what);
if opts.theory_circle
    write_circle_text(fid, P, nets, conds);
end
if opts.mark_outliers
    write_outlier_text(fid, P, nets, conds, opts);
end
if opts.inset
    write_inset_text(fid, P, nets, conds, opts);
end
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
function d = detachment(R, m, factor)
% The rightmost eigenvalue against the two circles (computed here, not in
% outlier_measures). Units of W, measured from the center of the bulk (0 in
% units of W, -1/tau_d as a rate), so a ratio of 1 is on the circle.
ok = R.R_eff > 1e-6;
d.ratio_cal_med = median(m.lam1(ok) ./ (m.kappa .* R.R_eff(ok)));    % outlier_measures' ratio_med
d.ratio_theory_med = median(m.lam1(ok) ./ R.R_eff(ok));              % the same against R_eff alone
d.frac_window_theory = mean(m.lam1 > R.R_eff);                       % all states of the window, each with its own R_eff
d.frac_window_cal = mean(m.lam1 > m.kappa .* R.R_eff);
% The 12 drawn states, against the circles AS DRAWN (median radius over the
% window, distance from the center), and against each state's own radius.
ev = [R.cal.ev_M];                                                   % n x n_cal
[~, i1] = max(real(ev), [], 1);
top = ev(sub2ind(size(ev), i1, 1:size(ev, 2)));
R_med = median(R.R_eff);
d.n_states = size(ev, 2);
d.frac12_theory_drawn = mean(abs(top) > R_med);
d.frac12_cal_drawn = mean(abs(top) > m.kappa * R_med);
d.frac12_theory_re = mean(real(top) > [R.cal.R_eff]);
d.frac12_cal_re = mean(real(top) > m.kappa .* [R.cal.R_eff]);
d.n_between = mean(sum(abs(ev) > R_med & abs(ev) <= m.kappa * R_med, 1));   % eigenvalues per state between the circles
d.n_beyond = mean(sum(abs(ev) > m.kappa * R_med, 1));                       % and outside the dashed one
d.kappa = m.kappa;
% Every eigenvalue outside the theoretical circle (any direction), mean per
% state: against the circle as drawn (R_med) and against each state's own
% R_eff; right of the center (Re > 0 in units of W) and right of zero rate
% (Re > 1). The same for the null draws saved with the run (ev_null: three
% structureless matrices with the observed gains per state).
R_own = [R.cal.R_eff];
d.out_marked = count_out(ev, factor * R_med);                         % what is marked
d.out_marked_own = count_out(ev, factor * R_own);
d.out_drawn = count_out(ev, R_med);                                   % at the circle itself
d.out_own = count_out(ev, R_own);
nul = cat(3, R.cal.ev_null);                                          % n x n_draws x n_cal
d.n_null_draws = size(nul, 2);
acc = zeros(size(nul, 2), 3); acc_own = acc; acc_f = acc;
for j = 1:size(nul, 2)
    acc_f(j, :) = count_out(squeeze(nul(:, j, :)), factor * R_med);
    acc(j, :) = count_out(squeeze(nul(:, j, :)), R_med);
    acc_own(j, :) = count_out(squeeze(nul(:, j, :)), R_own);
end
d.null_marked = mean(acc_f, 1);
d.null_drawn = mean(acc, 1);
d.null_own = mean(acc_own, 1);
end

function c = count_out(ev, rad)
% [all, right of the center, right of zero rate], mean per state (column).
out = abs(ev) > rad;
c = [mean(sum(out, 1)), mean(sum(out & real(ev) > 0, 1)), mean(sum(out & real(ev) > 1, 1))];
end

function write_outlier_text(fid, P, nets, conds, opts)
drawn = P(~cellfun(@isempty, P));
f = opts.outlier_factor;
fprintf(fid, ['MARKED OUTLIERS (Tom, 2026-11-20). Every eigenvalue further than %g times the unscaled theoretical radius ' ...
    'R_eff from the center -1/tau_d, in any direction, is drawn opaque and a little larger (marker area %g pt^2 ' ...
    'against 5), in the condition color, with no edge. The threshold is judged against the circle AS DRAWN (the median ' ...
    'R_eff over the window); the solid circle itself is at the radius (1.00), and no circle is drawn at %g. The ' ...
    'eigenvalues inside the threshold stay pale and translucent. Counts are means per state over the 12 drawn states: ' ...
    'all beyond the threshold, those of them right of the center (Re > -1/tau_d), and those right of zero (Re > 0, ' ...
    'locally unstable with the slow variables frozen). `own R_eff` judges each state against %g times its own radius ' ...
    'instead of the median. The same three counts follow at 1.00 (everything outside the solid circle). `null` is the ' ...
    'count for the structureless matrices saved with each run (`ev_null`: %d draws per state, the weights centered ' ...
    'and permuted within blocks, the observed gains): what chance alone puts beyond the threshold.\n\n'], ...
    f, opts.outlier_size, f, f, drawn{1}.det.n_null_draws);
fprintf(fid, ['| network | condition | beyond %g R (marked) | right of center | right of zero | beyond %g R (own R_eff) | ' ...
    'null: beyond %g R | null: right of center | null: right of zero | beyond 1.00 R | right of center | right of zero | ' ...
    'null: beyond 1.00 R |\n'], f, f, f);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for w = 1:size(P, 1)
    for c = 1:size(P, 2)
        if isempty(P{w, c}); continue; end
        d = P{w, c}.det;
        fprintf(fid, '| %s | %s | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f |\n', ...
            nets{w}, conds{c}, d.out_marked, d.out_marked_own(1), d.null_marked, d.out_drawn, d.null_drawn(1));
    end
end
fprintf(fid, '\n');
end

function write_circle_text(fid, P, nets, conds)
fprintf(fid, ['The rightmost eigenvalue against each circle. `detachment` is Re lambda_1 over the radius, median over ' ...
    'all 2001 states of the window (1 on the circle; `calibrated` is `ratio_med` of the table below, `theoretical` is ' ...
    'the same against R_eff alone). `outside, window` is the fraction of those states with Re lambda_1 beyond the ' ...
    'state''s own radius. `outside, 12 drawn` is the fraction of the 12 drawn states whose rightmost eigenvalue ' ...
    'lies outside the circle AS DRAWN (distance from the center against the median radius); in ' ...
    'parentheses the same by real part against the state''s own radius. The last two columns count all eigenvalues, ' ...
    'mean per drawn state: between the two circles, and outside the dashed one.\n\n']);
fprintf(fid, ['| network | condition | kappa | detachment, theoretical | detachment, calibrated | outside theoretical, window | ' ...
    'outside calibrated, window | outside theoretical, 12 drawn | outside calibrated, 12 drawn | eigenvalues between the circles | ' ...
    'eigenvalues outside the dashed circle |\n']);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|\n');
for w = 1:size(P, 1)
    for c = 1:size(P, 2)
        if isempty(P{w, c}); continue; end
        d = P{w, c}.det;
        fprintf(fid, '| %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f (%.3f) | %.3f (%.3f) | %.1f | %.1f |\n', ...
            nets{w}, conds{c}, d.kappa, d.ratio_theory_med, d.ratio_cal_med, d.frac_window_theory, d.frac_window_cal, ...
            d.frac12_theory_drawn, d.frac12_theory_re, d.frac12_cal_drawn, d.frac12_cal_re, d.n_between, d.n_beyond);
    end
end
fprintf(fid, '\n');
end

function q = load_states(states_dir, net, cond, seed, R, row)
% The operating points of one run (plane_states), checked against the saved
% run, and the fractions on each part of the curve by cell type.
f = fullfile(states_dir, plane_states_name(net, cond, seed));
if ~isfile(f)
    error('plane_fig_outliers:NoStates', ['No %s. The saved runs do not hold the states; make it with ' ...
        'plane_states(''%s'', ''%s'', [%d %d]) (a re-simulation, about half a minute).'], f, net, cond, seed);
end
L = load(f, 'Q'); Q = L.Q;
worst = max(struct2array(Q.check));
if worst > 1e-6 || ~isequal(size(Q.t), size(R.cal_t)) || any(abs(Q.t - R.cal_t) > 1e-9)
    error('plane_fig_outliers:StatesDiffer', ['%s is not the saved run''s trajectory (largest difference %.3g): ' ...
        'the inset would show other states than the spectra.'], f, worst);
end
q = struct('u', Q.u, 'r', Q.r, 'type', Q.type, 'S_a', Q.S_a, 'check', Q.check, 'worst', worst);
n_types = numel(Q.type_names);
if n_types ~= 2 || ~isequal(Q.type_names(:)', {'E', 'I'})
    error('plane_fig_outliers:Types', 'The inset colors E and I; %s has types %s.', f, strjoin(Q.type_names, ', '));
end
sample_of = round((R.cal_t - R.t(1)) / R.dt) + 1;
for k = 1:n_types
    is = Q.type == k;
    q.frac_slope(k) = mean(Q.rp(is, :) > 0, 'all');
    q.frac_floor(k) = mean(Q.rp(is, :) == 0 & Q.r(is, :) == 0, 'all');
    q.frac_ceil(k)  = mean(Q.rp(is, :) == 0 & Q.r(is, :) == 1, 'all');
    q.n(k) = nnz(is) * size(Q.u, 2);
    q.frac_on_saved12(k) = mean(R.frac_on(k, sample_of));     % the saved run, same 12 states
    q.frac_on_all(k) = mean(R.frac_on(k, :));                 % the saved run, all states of the window
end
q.frac_on_measures = [row.frac_on_E, row.frac_on_I];         % outlier_measures: median over the window
if any(abs(q.frac_slope + q.frac_floor + q.frac_ceil - 1) > 1e-12)
    error('plane_fig_outliers:Parts', 'Sloped, floor and ceiling do not add to 1 in %s.', f);
end
end

function draw_inset(fig, rect, q, ei_col, rs, opts, lw)
% The nonlinearity with every neuron's operating point. Horizontal: u in
% [-L, L], then a gap, then a strip for the neurons beyond, on each side.
L = opts.inset_u_lim; gap = 0.3; strip = 0.6; jit = opts.inset_jitter;
[n, n_cal] = size(q.u);
cols = 1:n_cal;
if ~isempty(opts.inset_snapshots); cols = opts.inset_snapshots(:)'; end
keep = true(n, 1);
if isfinite(opts.inset_n_neurons) && opts.inset_n_neurons < n
    keep(:) = false;
    for k = 1:2                                         % the same share of each type
        ids = find(q.type == k);
        m = round(opts.inset_n_neurons * numel(ids) / n);
        keep(ids(rs.randperm(numel(ids), m))) = true;
    end
end
u = q.u(keep, cols); r = q.r(keep, cols);
ty = repmat(q.type(keep), 1, numel(cols));
u = u(:); r = r(:); ty = ty(:);
x = u;
lo = u < -L; hi = u > L;
x(lo) = -L - gap - strip * rs.rand(nnz(lo), 1);
x(hi) =  L + gap + strip * rs.rand(nnz(hi), 1);
if opts.inset_split                                     % E above the curve, I below it: the piles do not mix
    y = r + jit * rs.rand(numel(r), 1) .* (3 - 2 * ty);
else
    y = r + jit * (2 * rs.rand(numel(r), 1) - 1);
end
order = rs.randperm(numel(x));                          % neither type is drawn on top of the other
ax = axes(fig, 'Units', 'inches', 'Position', rect, 'Visible', 'off', 'Tag', 'plane_fig_outliers_inset');
hold(ax, 'on');
scatter(ax, x(order), y(order), opts.inset_size, ei_col(ty(order), :), 'filled', ...
    'MarkerFaceAlpha', opts.inset_alpha, 'MarkerEdgeColor', 'none');
uu = linspace(-L, L, 400);
plot(ax, uu, SRNNCellTypePairs.piecewiseSigmoid(uu, q.S_a, 0), 'k-', 'LineWidth', 0.75 * lw);
plot(ax, [-L - gap - strip, -L - gap], [0 0], 'k-', 'LineWidth', 0.75 * lw);   % the floor, beyond the gap
plot(ax, [L + gap, L + gap + strip], [1 1], 'k-', 'LineWidth', 0.75 * lw);     % the ceiling, beyond the gap
set(ax, 'XLim', [-1 1] * (L + gap + strip + 0.05), 'YLim', [-jit - 0.06, 1 + jit + 0.06], 'Clipping', 'off');
end

function write_inset_text(fid, P, nets, conds, opts)
drawn = P(~cellfun(@isempty, P));
q1 = drawn{1}.inset;
[n, n_cal] = size(q1.u);
cols = 1:n_cal;
if ~isempty(opts.inset_snapshots); cols = opts.inset_snapshots(:)'; end
if isfinite(opts.inset_n_neurons) && opts.inset_n_neurons < n
    who = sprintf(['a random subset of %d of the %d neurons (the same neurons in every state, in the proportion ' ...
        'of the types; drawn per panel with the jitter''s stream)'], opts.inset_n_neurons, n);
else
    who = sprintf('all %d neurons', n);
end
if numel(cols) == n_cal
    when = sprintf('each of the %d states whose spectra the panel overlays', n_cal);
else
    when = sprintf('state %s of the %d whose spectra the panel overlays', list(cols), n_cal);
end
a = q1.S_a / 2;
if opts.inset_split
    how = sprintf(['a uniform random offset by cell type, excitatory neurons UPWARD by 0 to %g and inhibitory neurons ' ...
        'DOWNWARD by 0 to %g (in units of the rate; fixed seed %d), so that the two types lie on either side of the ' ...
        'curve and a pile of both can be read'], opts.inset_jitter, opts.inset_jitter, opts.inset_seed);
else
    how = sprintf('a uniform random jitter of +/- %g (in units of the rate; fixed seed %d) so that piles can be seen', ...
        opts.inset_jitter, opts.inset_seed);
end
fprintf(fid, ['INSET (lower right of each panel; no text on it). The black curve is the model''s firing-rate ' ...
    'nonlinearity, the piecewise sigmoid (`SRNNCellTypePairs.piecewiseSigmoid`, S_a = %g): rate exactly 0 (floor, ' ...
    'slope 0) below u = %g, quadratic up to %g, slope 1 up to %g, quadratic up to %g, rate exactly 1 (ceiling, slope 0) ' ...
    'above. Each dot is one neuron at one state, at its operating point: %s, at %s. Magenta: excitatory; teal: inhibitory ' ...
    '(the first colors of `excitatory_colormap` and `inhibitory_colormap`; these are the CELL TYPE colors and have ' ...
    'nothing to do with the condition colors of the eigenvalue dots). Horizontal: u = x - (c/K) sum_k a_k - S_c, the ' ...
    'argument the nonlinearity receives (dendritic state x minus the neuron''s adaptation current) relative to the ' ...
    'neuron''s own setpoint S_c (drawn per neuron in this preset), so that every neuron lies on the one curve ' ...
    'centered on u = 0. Vertical: the rate r = phi(u), 0 to 1, plus %s: the dots are NOT off the curve, and the band at a flat ' ...
    'part is one rate. The inset shows u from %g to %g, the same in all panels. A neuron beyond that is drawn in the ' ...
    'short strip past the gap in the floor line (left, u < %g) or in the ceiling line (right, u > %g), at a random ' ...
    'horizontal place in the strip: the strip says "further out", not how far. Marker area %g pt^2, opacity %g, ' ...
    'the two types drawn in random order.\n\n'], ...
    q1.S_a, a - 1, -a, a, 1 - a, who, when, how, ...
    -opts.inset_u_lim, opts.inset_u_lim, -opts.inset_u_lim, opts.inset_u_lim, opts.inset_size, opts.inset_alpha);
fprintf(fid, ['The states are a re-simulation of each saved run (`plane_states`; the runs save no states), ' ...
    'checked against it at the 12 states: the mean and rms gains by route, the rates and the fractions with nonzero ' ...
    'slope are identical (largest difference %.3g), and all 500 eigenvalues of the effective connectivity agree to ' ...
    '%.3g (units of W).\n\n'], ...
    max(cellfun(@(p) max([p.inset.check.dbar, p.inset.check.drms, p.inset.check.rate, p.inset.check.frac_on]), drawn)), ...
    max(cellfun(@(p) p.inset.check.ev_M, drawn)));
fprintf(fid, ['Fractions of neuron-states on each part of the curve, over ALL neurons and all %d states whatever ' ...
    'the inset draws (sloped: slope > 0; floor: rate 0; ceiling: rate 1), the range of u, and the share drawn in the ' ...
    'strips. `saved, 12 states` is the saved run''s `frac_on` averaged over the same states (it must equal ' ...
    '`sloped`); `frac_on (measures)` is `outlier_measures`'' `frac_on_E` / `frac_on_I`, the MEDIAN over all 2001 ' ...
    'states of the 20-s window, so it can differ a little from the 12-state mean.\n\n'], n_cal);
fprintf(fid, '| network | condition | type | neuron-states | sloped | floor | ceiling | saved, 12 states | frac_on (measures) | u min | u max | in left strip | in right strip |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
names = 'EI';
for w = 1:size(P, 1)
    for c = 1:size(P, 2)
        if isempty(P{w, c}); continue; end
        q = P{w, c}.inset;
        for k = 1:2
            uk = q.u(q.type == k, :);
            fprintf(fid, '| %s | %s | %s | %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.3f | %.3f |\n', ...
                nets{w}, conds{c}, names(k), q.n(k), q.frac_slope(k), q.frac_floor(k), q.frac_ceil(k), ...
                q.frac_on_saved12(k), q.frac_on_measures(k), min(uk(:)), max(uk(:)), ...
                mean(uk(:) < -opts.inset_u_lim), mean(uk(:) > opts.inset_u_lim));
        end
    end
end
fprintf(fid, '\n');
end

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

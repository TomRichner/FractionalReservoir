function out = fig_main6_amplification(cfg)
% FIG_MAIN6_AMPLIFICATION Main Figure 6: selective amplification, in pattern (non-normal margin, directional gains) and in timescale (frequency response).
%
%   out = FIG_MAIN6_AMPLIFICATION('run_dir', d)
%   out = FIG_MAIN6_AMPLIFICATION('run_dir', d_spectral, 'transient_run_dir', d_paper)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md), registry name
% fig_main6_grouped. Reads saved stages only; simulates nothing; errors
% through resolve_data_file when a stage is missing.
%
%   (A)  the non-normal margin omega(J_xx) - alpha(J_xx) of the dendritic
%        block per sampled state, per condition: 5 to 95 % (line),
%        interquartile range (box), median (stage eig_heatmap; the margin
%        panel of fig_transient_amplification, the former Figure 5A).
%   (B)  the gain of the active dendritic propagator over time after a
%        perturbation, per condition, for five directions (worst case,
%        noise rms, excite E / inhibit I, excite E / excite I, Lyapunov
%        direction), medians over the regular states (stage transient_gain;
%        the former Figure 5B, drawn as fig_grouped_main's case 5 drew it).
%   (C)  the bulk radius of the effective connectivity M(i omega) against
%        frequency, reference network (fig_freqresp_network, left panel).
%   (D)  the loop gain against frequency, reference network (the same
%        figure's right panel). The middle panel is not in the main figure.
%        The bars at the foot of the stand-alone panel are not drawn here.
%        (C) and (D) together take the figure's width.
%   (E)  the complex-plane view of the reference network, one panel per
%        condition (fig_freqresp_plane, its top row): the full Jacobian's
%        eigenvalues, the frozen-gain edge (dashed) and the edge once the
%        slow variables respond (colored). One seed and state: the stage's
%        first seed (cfg.seed picks another).
%
% (A) and (B) are drawn here from the stages' files (the code of the
% present Figure 5 lives in fig_grouped_main's dispatch, which this file
% does not call, so that the dispatch can be pointed here); (C) to (E) are
% the frequency-response figures' own drawings, copied (composite_figure).
%
% WHERE THE STAGES ARE. One run directory holds all three in a registered
% run. transient_run_dir names another run for (A) and (B): while the
% figure was built, the spectral stages were in one run and the paper's
% transient-gain run in another.
%
% Fig_Main6_Grouped_table.md: a caption paragraph, the numbers of (A) and
% (B), and the two frequency-response sidecars.
%
% Options: run_dir, transient_run_dir, out_dir, save, visible, verbose,
% preset_name (unused), conditions (select_conditions), network (the
% network of C to E, 'reference'), seed, gain_window_s (the stretch of
% the horizon shown in B, 2 s), fig_tag.
%
% See also: fig_transient_amplification, fig_transient_gain,
%           fig_freqresp_network, fig_freqresp_plane, composite_figure

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.run_dir     (1,:) char    = ''
    cfg.transient_run_dir (1,:) char = ''   % the run holding eig_heatmap and transient_gain; '' = run_dir
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}    % names to draw, in order; empty = all a stage holds
    cfg.network     (1,:) char    = 'reference'
    cfg.seed double               = []    % the seed drawn in (E); [] = the first of the stage
    cfg.gain_window_s (1,1) double {mustBePositive} = 2
    cfg.fig_tag     (1,:) char    = 'Fig_Main6_Grouped'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
% A source figure that is already drawn when a later one fails is closed,
% not left behind.
before = findall(groot, 'Type', 'figure');
try
    out = build(cfg, out_dir);
catch err
    delete(setdiff(findall(groot, 'Type', 'figure'), before));
    rethrow(err);
end
end

function out = build(cfg, out_dir)
tr_dir = cfg.transient_run_dir;
if isempty(tr_dir); tr_dir = cfg.run_dir; end

% (C) and (D) take the whole width of the figure, as (A)+(B) and (E) do:
% the network figure is drawn that wide (two panels, their height kept), not
% scaled up, so the row is also shorter than it was at 10.5 in (polish list,
% 2026-11-22: the row left a third of the width empty). The bars at the foot
% of (D) are left out: nothing on the panel says what they are; they stay in
% the stand-alone figure and in the sidecar.
% The common canvas (Plan 2, canvas standard): every row is canvas_w less
% composite_figure's side margins; the type is set by role at the end
% (canvas_type).
% The Re name of (E) sits ABOVE the right end of the last plane's real axis,
% as in Figure 5 (Plan 3 fix X11: under the tick labels it was about 0.4 in
% below its axis), raised clear of what that panel draws under it
% (lift_re_label, below); so the foot needs no room for it.
st = manuscript_style(); MARGIN = [0.15 0.15 0.1 0.1];
ROW_W = st.canvas_w - MARGIN(1) - MARGIN(2); ppi = get(groot, 'ScreenPixelsPerInch');
[f_ab, t_ab, src_ab] = draw_margin_and_gain(cfg, tr_dir, ROW_W);
[f_net, t_net, o_net] = composite_part(@fig_freqresp_network, cfg, 'networks', {cfg.network}, 'panels', [1 3], ...
    'bulk_bars', false, 'fig_size', round([ROW_W * ppi, 2.65 * ppi]));   % 2.65 in tall: the figure fits canvas_max_h
[f_plane, t_plane] = composite_part(@fig_freqresp_plane, cfg, 'networks', {cfg.network}, 'seed', cfg.seed, ...
    're_label_side', 'above');

[fig, B] = composite_figure({struct('fig', f_ab, 'width', ROW_W), struct('fig', f_net, 'width', ROW_W), ...
    struct('fig', f_plane, 'width', ROW_W)}, 'visible', cfg.visible, 'margin', MARGIN, 'row_gap', 0.2, 'width', st.canvas_w);
ab = B{1}(1).ax;                                     % margin, then the gain panels, left to right
composite_letter(fig, B{1}(1), 'A', ab(1));
composite_letter(fig, B{1}(1), 'B', ab(2));
net = B{2}(1).ax;
composite_letter(fig, B{2}(1), 'C', net(1));
composite_letter(fig, B{2}(1), 'D', net(2));
composite_letter(fig, B{3}(1), 'E');
canvas_type(fig);
lift_re_label(fig);

% The dendritic time constant, for the caption's corner frequency (Tom's
% question on panel C, 2026-11-22).
src_net = o_net.source; if iscell(src_net); src_net = src_net{1}; end
R = load(src_net, 'runs');
tau_d = R.runs(1).tau_d;
clear R
caption = sprintf(['Main Figure 6. (A) Non-normal margin of the dendritic block of the Jacobian, numerical minus ' ...
    'spectral abscissa, per sampled state and condition: line 5 to 95 %%, box the interquartile range, bar the median. ' ...
    '(B) Gain of the active dendritic propagator against time after a perturbation, medians over the regular ' ...
    'states, for five directions (legend); the first %g s of the saved horizon; dashed gray: unity. (C) Bulk radius ' ...
    'of the effective connectivity M(i omega) against frequency, network `%s`: median over seeds and calibration ' ...
    'states, quartiles as a band. In plain words, the bulk radius of M(i omega) is the radius of the effective ' ...
    'connectivity as a perturbation at that frequency sees it: the slow variables follow a slow perturbation and ' ...
    'not a fast one. Gray: |1 + i omega tau_d|, the threshold, the modulus an eigenvalue must reach for a mode at that ' ...
    'frequency to be marginal; it is 1 at low frequency and rises above 1/(2 pi tau_d) = %.2g Hz (tau_d = %g s) ' ...
    'because the dendrite filters fast input, so a fast mode needs a larger radius to grow. (D) Loop gain, the bulk radius over that modulus; above 1 (gray line) the bulk holds ' ...
    'growing modes at that frequency; triangles: the frequency of the leftover mode (the bars of the stand-alone ' ...
    'figure, the rotation rates of the frozen bulk, are not drawn here; their numbers are in the table below). (E) The same in the complex plane, network `%s`, one seed and state: gray dots the ' ...
    'full Jacobian''s eigenvalues (upper half plane), dashed the edge of the bulk with the adaptive variables frozen, ' ...
    'colored line the edge once they respond, filled dot the leftover mode. Conditions in their colors throughout. ' ...
    '(A) and (B) are read from `%s`; (C) to (E) from the stage frequency_response of `%s`.'], ...
    cfg.gain_window_s, cfg.network, 1 / (2 * pi * tau_d), tau_d, cfg.network, tr_dir, cfg.run_dir);
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Main Figure 6: selective amplification', caption, ...
    {'Panels A and B (eig_heatmap, transient_gain)', t_ab; 'Panels C and D (fig_freqresp_network, panels 1 and 3)', t_net; ...
    'Panel E (fig_freqresp_plane)', t_plane}, [src_ab, {o_net.source}]);
end

%% ------------------------------------------------------------------------
function [fig, txt, sources] = draw_margin_and_gain(cfg, run_dir, fig_w)
% Panels (A) and (B) on a figure of their own, laid out in inches: one
% margin panel and one gain panel per condition, the directions' legend
% under the gain panels. Returns the numbers as markdown.
st = manuscript_style();
f_eig = resolve_data_file('', run_dir, {fullfile(run_dir, 'eig_heatmap')}, 'eig_heatmap_data.mat', ...
    'Run the stage eig_heatmap (run_eig_heatmap) first');
f_gain = resolve_data_file('', run_dir, {fullfile(run_dir, 'transient_gain')}, 'transient_gain_data.mat', ...
    'Run the stage transient_gain first');
sources = {f_eig, f_gain};
E = load(f_eig, 'cond_names', 'num_abscissa_by_cond', 'spec_abscissa_by_cond');
if ~isfield(E, 'num_abscissa_by_cond')
    error('fig_main6_amplification:OldData', '%s predates the numerical-abscissa fields; rerun run_eig_heatmap.', f_eig);
end
[m_names, m_idx] = select_conditions(E.cond_names, cfg.conditions);
G = load(f_gain, 'results', 'settings');
[g_names, g_idx] = select_conditions({G.results.name}, cfg.conditions);
results = G.results(g_idx);
n_m = numel(m_names); n_g = numel(g_names);

% Geometry (inches), made for the 9 in canvas (Plan 3, 2026-11-22): the
% margin panel narrow (0.4 in per condition), the gain panels as wide as
% the rest allows, one axis name under the middle gain panel (three of them
% ran into each other at this width), the key under (B).
fig_h = 3.6; bottom = 1.3; ph = 1.7;
m_left = 0.95; m_w = 0.4 * n_m;
g_left = m_left + m_w + 0.95; g_gap = 0.25;
g_w = (fig_w - 0.1 - g_left - (n_g - 1) * g_gap) / n_g;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_w fig_h], 'Visible', 'off');
ax_style = {'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1};

% (A) the margin.
ax = axes(fig, 'Units', 'inches', 'Position', [m_left bottom m_w ph], 'Tag', 'main6_margin'); hold(ax, 'on');
ticks = grouped_condition_label(m_names, 'tick');
lines = {'| condition | states | 5 % | 25 % | median | 75 % | 95 % |', '|---|---|---|---|---|---|---|'};
top = 0;
for i = 1:n_m
    col = st.condition_color(m_names{i});
    m = E.num_abscissa_by_cond{m_idx(i)} - E.spec_abscissa_by_cond{m_idx(i)};
    q = prctile(m, [5 25 50 75 95]);
    plot(ax, [i i], q([1 5]), '-', 'Color', col, 'LineWidth', 1.2);
    patch(ax, i + 0.25 * [-1 1 1 -1], q([2 2 4 4]), col, 'FaceAlpha', 0.35, 'EdgeColor', col);
    plot(ax, i + 0.25 * [-1 1], q([3 3]), '-', 'Color', col, 'LineWidth', 2.5);
    ticks{i} = sprintf('\\color[rgb]{%.4f,%.4f,%.4f}%s', col, ticks{i});
    lines{end + 1} = sprintf('| %s | %d | %.3g | %.3g | %.3g | %.3g | %.3g |', m_names{i}, numel(m), q); %#ok<AGROW>
    top = max(top, q(5));
end
y_hi = max(100, 25 * ceil(1.05 * top / 25));        % [0, 100] with ticks at 0 to 75, as the former Figure 5A; wider if the data need it
set(ax, ax_style{:}, 'XLim', [0.4, n_m + 0.6], 'XTick', 1:n_m, 'XTickLabel', ticks, 'TickLabelInterpreter', 'tex', ...
    'XTickLabelRotation', 45, 'YLim', [0 y_hi], 'YTick', 0:25:(y_hi - 25));
ylabel(ax, {'Non-normal amplification', 'margin (s^{-1})'}, 'FontSize', st.tick_fs);

% (B) the directional gains.
keys = {'G_worst', 'G_noise', 'G_ei_diff', 'G_ei_sum', 'G_lyap'};
labels = {'Worst case', 'Noise rms', 'Excite E / inhibit I', 'Excite E / excite I', 'Lyapunov direction'};
colors = [0 0 0; 0 0.45 0.7; 0.85 0.33 0.1; 0 0.6 0.5; 0.55 0.2 0.65];
styles = {'-', '--', '-.', ':', '--'}; marks = {'none', 'none', 'o', 's', '^'};
ia = find(strcmp(G.settings.variants, 'active'));
if ~isscalar(ia)
    error('fig_main6_amplification:NoActive', '%s has no ''active'' propagator variant.', f_gain);
end
gain_ax = gobjects(1, n_g); hh = gobjects(1, 5);
glines = {sprintf('| condition | direction | regular states | largest median gain in the first %g s | median gain at %g s |', ...
    cfg.gain_window_s, cfg.gain_window_s), '|---|---|---|---|---|'};
for i = 1:n_g
    ax = axes(fig, 'Units', 'inches', 'Position', [g_left + (i - 1) * (g_w + g_gap), bottom, g_w, ph], ...
        'Tag', sprintf('main6_gain_%s', g_names{i}));
    gain_ax(i) = ax; hold(ax, 'on');
    smp = results(i).samples(strcmp({results(i).samples.kind}, 'regular'));
    t = smp(1).t;
    shown = t <= cfg.gain_window_s + 1e-9;
    for k = 1:5
        Y = cell2mat(arrayfun(@(s) s.(keys{k})(ia, :), smp(:), 'UniformOutput', false));
        md = median(Y, 1);
        hh(k) = plot(ax, t, md, 'Color', colors(k, :), 'LineStyle', styles{k}, 'LineWidth', 1.8, 'Marker', marks{k}, ...
            'MarkerIndices', unique(round(linspace(2, numel(t), 8))), 'MarkerSize', 4);
        glines{end + 1} = sprintf('| %s | %s | %d | %.3g | %.3g |', g_names{i}, labels{k}, numel(smp), ...
            max(md(shown)), md(find(shown, 1, 'last'))); %#ok<AGROW>
    end
    unity = yline(ax, 1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    set(ax, ax_style{:}, 'YScale', 'log', 'XLim', [0 cfg.gain_window_s], 'XTick', 0:1:cfg.gain_window_s, 'YMinorTick', 'off');
    title(ax, two_lines(g_names{i}), 'FontWeight', 'normal', 'FontSize', st.tick_fs, 'Color', st.condition_color(g_names{i}));
    if i == ceil(n_g / 2)
        xlabel(ax, 'Time after perturbation (s)', 'FontSize', st.tick_fs);
    end
    if i == 1
        ylabel(ax, 'Active dendritic gain', 'FontSize', st.tick_fs);
    end
end
linkaxes(gain_ax, 'y');
% The decades, every other one when there are many: at 1.7 in MATLAB's own
% choice kept only 10^0 (9 in canvas, 2026-11-22).
yl = gain_ax(1).YLim;
dec = ceil(log10(yl(1))):floor(log10(yl(2)));
if numel(dec) > 3; dec = dec(mod(dec, 2) == 0); end
set(gain_ax, 'YTick', 10 .^ dec);
set(gain_ax(2:end), 'YTickLabel', []);
lg = legend(gain_ax(min(2, n_g)), [hh unity], [labels {'Unity gain'}], 'NumColumns', 3, 'Box', 'off', ...
    'FontSize', st.tick_fs, 'Units', 'inches');
drawnow;
p = lg.Position;
x_mid = g_left + (n_g * g_w + (n_g - 1) * g_gap) / 2;
lg.Position = [min(max(x_mid - p(3) / 2, 0.1), fig_w - p(3) - 0.05), 0.02, p(3), p(4)];

txt = sprintf(['# Numbers of panels A and B\n\nPanel A: the non-normal margin omega(J_xx) - alpha(J_xx) in 1/s over the ' ...
    'states the stage eig_heatmap sampled (`%s`).\n\n%s\n\nPanel B: the gain of the active dendritic propagator, ' ...
    'medians over the regular states of the stage transient_gain (`%s`; saved horizon %g s).\n\n%s\n'], ...
    f_eig, strjoin(lines, newline), f_gain, G.settings.horizon_s, strjoin(glines, newline));
end

function t = two_lines(name)
% The condition's sentence-case title, broken before "adaptation" when it
% has a word in front of it (the panels are narrower than the long names).
t = grouped_condition_label({name}, 'sentence');
t = t{1};
k = strfind(t, ' adaptation');
if ~isempty(k) && ~startsWith(t, 'No ')
    t = {t(1:k(1) - 1), t(k(1) + 1:end)};
end
end

function lift_re_label(fig)
% The Re name of (E), above the right end of the last plane's real axis, is
% raised until nothing that plane draws under it (the eigenvalues, the edge
% curves, the frozen circle) touches it: 0.05 in clear of the highest such
% point within its width (Plan 3 fix X11). The axis lines and their tick
% labels (Tag complex_plane_axes) are not counted.
h = findall(fig, 'Type', 'text', 'String', 'Re (s^{-1})');
if ~isscalar(h); return; end
ax = h.Parent;
drawnow;
e = h.Extent;                                        % data units
xs = []; ys = [];
for g = reshape(findall(ax, 'Type', 'line'), 1, [])
    if strcmp(g.Tag, 'complex_plane_axes'); continue; end
    xs = [xs, g.XData(:)']; ys = [ys, g.YData(:)']; %#ok<AGROW>
end
for g = reshape(findall(ax, 'Type', 'contour'), 1, [])
    M = g.ContourMatrix; k = 1;
    while k < size(M, 2)
        n = M(2, k); xs = [xs, M(1, k + 1:k + n)]; ys = [ys, M(2, k + 1:k + n)]; %#ok<AGROW>
        k = k + n + 1;
    end
end
px = getpixelposition(ax); ppi = get(groot, 'ScreenPixelsPerInch');
pad_x = 0.05 * ppi * diff(ax.XLim) / px(3); pad_y = 0.05 * ppi * diff(ax.YLim) / px(4);
% First to the right, by at most 0.4 in past the axis end, so that it
% starts clear of what is drawn beside it at its height (the frozen
% circle's foot); only what is still under it then raises it.
low = ys >= 0 & ys <= e(2) + e(4) + pad_y;
x_end = h.Position(1);
if any(low & xs > e(1) - pad_x)
    right = max(xs(low)) + pad_x + e(3);
    if right > x_end && right - x_end <= 0.4 * ppi * diff(ax.XLim) / px(3)
        h.Position(1) = right; drawnow; e = h.Extent;
    end
end
under = xs >= e(1) & xs <= e(1) + e(3) & ys >= 0;
if ~any(under); return; end
lift = max(ys(under)) + pad_y - e(2);
if lift > 0; h.Position(2) = h.Position(2) + lift; end
end

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

[f_ab, t_ab, src_ab] = draw_margin_and_gain(cfg, tr_dir);
[f_net, t_net, o_net] = composite_part(@fig_freqresp_network, cfg, 'networks', {cfg.network}, 'panels', [1 3]);
[f_plane, t_plane] = composite_part(@fig_freqresp_plane, cfg, 'networks', {cfg.network}, 'seed', cfg.seed);

[fig, B] = composite_figure({struct('fig', f_ab, 'width', []), struct('fig', f_net, 'width', 10.5), ...
    struct('fig', f_plane, 'width', [])}, 'visible', cfg.visible);
ab = B{1}(1).ax;                                     % margin, then the gain panels, left to right
composite_letter(fig, B{1}(1), 'A', ab(1));
composite_letter(fig, B{1}(1), 'B', ab(2));
net = B{2}(1).ax;
composite_letter(fig, B{2}(1), 'C', net(1));
composite_letter(fig, B{2}(1), 'D', net(2));
composite_letter(fig, B{3}(1), 'E');

caption = sprintf(['Main Figure 6. (A) Non-normal margin of the dendritic block of the Jacobian, numerical minus ' ...
    'spectral abscissa, per sampled state and condition: line 5 to 95 %%, box the interquartile range, bar the median. ' ...
    '(B) Gain of the active dendritic propagator against time after a perturbation, medians over the regular ' ...
    'states, for five directions (legend); the first %g s of the saved horizon; dashed gray: unity. (C) Bulk radius ' ...
    'of the effective connectivity M(i omega) against frequency, network `%s`: median over seeds and calibration ' ...
    'states, quartiles as a band; gray: |1 + i omega tau_d|, the modulus an eigenvalue must reach for a mode at that ' ...
    'frequency to be marginal. (D) Loop gain, the bulk radius over that modulus; above 1 (gray line) the bulk holds ' ...
    'growing modes at that frequency; triangles: the frequency of the leftover mode; bars: the range of rotation ' ...
    'rates of the frozen bulk. (E) The same in the complex plane, network `%s`, one seed and state: gray dots the ' ...
    'full Jacobian''s eigenvalues (upper half plane), dashed the edge of the bulk with the adaptive variables frozen, ' ...
    'colored line the edge once they respond, filled dot the leftover mode. Conditions in their colors throughout. ' ...
    '(A) and (B) are read from `%s`; (C) to (E) from the stage frequency_response of `%s`.'], ...
    cfg.gain_window_s, cfg.network, cfg.network, tr_dir, cfg.run_dir);
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Main Figure 6: selective amplification', caption, ...
    {'Panels A and B (eig_heatmap, transient_gain)', t_ab; 'Panels C and D (fig_freqresp_network, panels 1 and 3)', t_net; ...
    'Panel E (fig_freqresp_plane)', t_plane}, [src_ab, {o_net.source}]);
end

%% ------------------------------------------------------------------------
function [fig, txt, sources] = draw_margin_and_gain(cfg, run_dir)
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

% Geometry (inches).
fig_w = 13.3; fig_h = 4.35; bottom = 1.45; ph = 2.25;
m_left = 1.0; m_w = 0.55 + 0.5 * n_m;
g_left = m_left + m_w + 1.25; g_gap = 0.4;
g_w = (fig_w - 0.2 - g_left - (n_g - 1) * g_gap) / n_g;
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
    'XTickLabelRotation', 30, 'YLim', [0 y_hi], 'YTick', 0:25:(y_hi - 25));
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
    xlabel(ax, 'Time after perturbation (s)', 'FontSize', st.tick_fs);
    if i == 1
        ylabel(ax, 'Active dendritic gain', 'FontSize', st.tick_fs);
    end
end
linkaxes(gain_ax, 'y');
set(gain_ax(2:end), 'YTickLabel', []);
lg = legend(gain_ax(min(2, n_g)), [hh unity], [labels {'Unity gain'}], 'NumColumns', 3, 'Box', 'off', ...
    'FontSize', st.tick_fs, 'Units', 'inches');
drawnow;
p = lg.Position;
x_mid = g_left + (n_g * g_w + (n_g - 1) * g_gap) / 2;
lg.Position = [min(max(x_mid - p(3) / 2, 0.1), fig_w - p(3) - 0.05), 0.08, p(3), p(4)];

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

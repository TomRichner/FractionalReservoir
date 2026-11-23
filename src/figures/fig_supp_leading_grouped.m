function out = fig_supp_leading_grouped(cfg)
% FIG_SUPP_LEADING_GROUPED Supplementary figure: what the slowly decaying direction is made of; saved data only.
%
%   out = FIG_SUPP_LEADING_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%
% The three panels on the make-up of the leading Lyapunov direction,
% together again (Tom's review, 2026-11-22: Plan2_figure_feedback.md,
% Figure 7 notes 1 and 5, S9 note 2, Figure 4 note 7; Plan3.md, question
% 10). Every panel is drawn by the function that owns it, invisibly and
% unsaved, and its axes are copied here (composite_tools); nothing is
% simulated. On the 9 in canvas:
%   (A) lambda_1 against the longest SFA time constant: the first figure of
%       fig_sfa_EOC_allStd (the tau sweep), with the display rules of
%       Figure 7's former panel A (linear axis, lower limit -0.2, no mean
%       markers).
%   (B) the shares of the leading direction's squared norm in the
%       spike-frequency adaptation, depression and dendritic coordinates
%       against the same time constant: the second figure of
%       fig_sfa_EOC_allStd (S9's former panel B), beside (A) on the same x
%       axis, its key at its right.
%   (C) the dendritic share of the leading perturbation against lambda_1
%       over the joint sample: fig_leading_direction (Figure 4's former
%       panel C), one panel per condition, no line at the uncoupled rate
%       (Tom, 2026-11-21).
% The axis names are in sentence case (Tom, 2026-11-23, G1); lambda_1 is
% named 'Growth rate, lambda_1 (s^-1)' on every panel (G3).
%
% WHICH RUN. run_dir holds the tau sweep (tau_sensitivity_*) and the joint
% sample (param_space_*): both are block A's, in the mu7chain run and in
% the A + C bundle's run directory.
%
% Options
%   run_dir       above
%   preset_name   the bundle's preset (fig_sfa_EOC_allStd)
%   conditions    names to draw in (C), in order; empty = all the joint
%                 sample holds. (A) and (B) draw the one condition of the
%                 selection the tau sweep holds.
%   out_dir, save, visible, verbose, fig_tag
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, fig_sfa_EOC_allStd, fig_leading_direction,
%           fig_main7_memory, fig_supp9_grouped

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''
    cfg.conditions                = {}
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.fig_tag     (1,:) char    = 'Fig_Supp_Leading_Grouped'
end
setup_paths();
T = composite_tools();
st = manuscript_style();
quiet = {'save', false, 'visible', false};

%% The source figures
S_ld = fig_leading_direction('run_dir', cfg.run_dir, 'conditions', cfg.conditions, 'verbose', cfg.verbose, quiet{:});
names = cellfun(@(s) s, {S_ld.stats.name}, 'UniformOutput', false);
nc = numel(names);
tau_args = [{'run_dir', cfg.run_dir, 'conditions', names}, quiet];
if ~isempty(cfg.preset_name); tau_args = [tau_args, {'preset_name', cfg.preset_name}]; end
S_tau = fig_sfa_EOC_allStd(tau_args{:});
src_figs = [S_ld.figs(:); S_tau.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));
assert(numel(S_tau.figs) >= 2, 'fig_supp_leading_grouped:NoBlocks', ...
    'fig_sfa_EOC_allStd returned no leading-vector figure (the tau sweep saved no block fractions).');

%% The canvas, laid out in pixels on canvas_w itself (fit has nothing to stretch)
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
Y_C = 62; H_C = 190;                       % (C), its two-line titles above
Y_AB = Y_C + H_C + 48 + 24 + 74; H_AB = 190;   % (A), (B): above the titles of (C) and their own two-line x labels
H = Y_AB + H_AB + 30 + 18;                 % the title of (A), and its letter
fig = T.canvas(W, H, cfg.visible);

% (A): the lambda_1 panel of the tau sweep (its first figure).
A = T.copy(S_tau.figs(1), fig);
ax_a = A.axes;
assert(isscalar(ax_a), 'fig_supp_leading_grouped:TauPanels', ...
    ['Expected one lambda_1 panel from the tau sweep, found %d: the sweep holds more than one of the ' ...
     'selected conditions. Select one.'], numel(ax_a));
delete(A.legends);
% (B): the leading-vector shares (the second figure of the tau sweep).
B = T.copy(S_tau.figs(2), fig);
assert(isscalar(B.axes), 'fig_supp_leading_grouped:TauPanels', ...
    'Expected one leading-vector panel from the tau sweep, found %d. Select one of its conditions.', numel(B.axes));
ax_b = B.axes;
% (C): one panel per condition.
C = T.copy(S_ld.axes, fig, false);
[~, order] = sort(arrayfun(@(a) a.Position(1), C.axes)); ax_c = C.axes(order);

%% Geometry (pixels)
MA = 108; WA = 228; GAP = 18; MB = 84; WB = 192;   % the margins hold a y label of two lines beside tick labels like '-0.2'
set(ax_a, 'Units', 'pixels', 'PositionConstraint', 'innerposition', 'Position', [MA, Y_AB, WA, H_AB]);
XB = MA + WA + GAP + MB;
set(ax_b, 'Units', 'pixels', 'PositionConstraint', 'innerposition', 'Position', [XB, Y_AB, WB, H_AB]);
MC = 114; CGAP = 36; RIGHT = 30;      % (C): a y label of three lines; the end labels of neighbours apart;
                                      % RIGHT: half the last tick label '10' (at 14 it touched the edge)
WC = (W - MC - (nc - 1) * CGAP - RIGHT) / nc;
for c = 1:nc
    set(ax_c(c), 'Units', 'pixels', 'PositionConstraint', 'innerposition', ...
        'Position', [MC + (c - 1) * (WC + CGAP), Y_C, WC, H_C]);
end

%% Style
T.style(fig);
lam_name = [st.label_growth ' (s^{-1})'];     % G3 (Tom, 2026-11-23): one name for lambda_1
tau_label = {'Longest SFA time', 'constant (s)'};
% (A), as Figure 7 drew it (style_main6_grouped; Tom, 2026-11-17).
tau_name = tau_condition(S_tau.source, names);
set(ax_a, 'XScale', 'linear');
xlabel(ax_a, tau_label, 'Interpreter', 'none');
ylabel(ax_a, {[extractBefore(st.label_growth, ',') ','], [st.label_lle ' (s^{-1})']}, 'Interpreter', 'tex');   % two lines: one is taller than the panel
title(ax_a, T.condition_title(tau_name), 'FontWeight', 'normal', 'Interpreter', 'tex', ...
    'Color', st.condition_color(tau_name));
zero_ref = findall(ax_a, 'Type', 'constantline');
zero_ref = zero_ref(arrayfun(@(h) h.Value == 0, zero_ref));
% Alpha 1: a constant line is drawn at 0.7 by default, so it came out a
% paler green than the dashed lines of (C) (Plan 3 review).
set(zero_ref, 'Color', st.zeroline_color, 'LineStyle', '--', 'LineWidth', st.zeroline_lw, 'Alpha', 1);
delete(findall(ax_a, 'Type', 'line', 'Marker', 'o', 'LineStyle', 'none'));
ax_a.YLim(1) = -0.2;
% (B), as S9 drew it: the state coordinates' palette and markers of
% style_main6_grouped (their own encoding, apart from the condition
% colours), on the x axis of (A).
set(ax_b, 'XScale', 'linear', 'XLim', ax_a.XLim, 'XTick', ax_a.XTick);
xlabel(ax_b, tau_label, 'Interpreter', 'none');
ylabel(ax_b, {'Leading-vector', 'squared-norm fraction'}, 'Interpreter', 'none');
title(ax_b, '');
ax_b.YTick = [0 0.5 1];
% The key: spike-frequency adaptation written out (rule 7 of "For every
% figure": SFA alone is meant), depression as on S8, the dendritic state.
state_names = {'SFA', 'STD', 'x'};
state_keys = {sprintf('Spike-frequency\nadaptation'), 'Depression', 'Dendritic state'};
state_colors = [213 94 0; 0 158 115; 142 68 173] / 255;
state_markers = {'o', 's', '^'};
h_state = gobjects(0); keys = {};
for k = 1:numel(state_names)
    h = findall(ax_b, 'Type', 'line', 'DisplayName', state_names{k});
    if k == 3 && isempty(h); h = findall(ax_b, 'Type', 'line', 'DisplayName', 'Dendritic State'); end
    if isempty(h); continue; end
    set(h, 'Color', state_colors(k, :), 'MarkerFaceColor', state_colors(k, :), 'Marker', state_markers{k}, ...
        'MarkerSize', 5, 'DisplayName', state_keys{k});
    h_state(end + 1) = h(1); keys{end + 1} = state_keys{k}; %#ok<AGROW>
end
delete(B.legends);
lg = legend(ax_b, h_state, keys, 'Box', 'off', 'FontSize', st.tick_fs, 'Location', 'none', 'Units', 'pixels', ...
    'Interpreter', 'none');
lg.ItemTokenSize = [22 18];
% (C): titles on two lines in the condition colours, the names in sentence
% case and lambda_1 named as in (A).
for c = 1:nc
    title(ax_c(c), T.condition_title(names{c}, true), 'FontWeight', 'normal', 'FontSize', st.title_fs, ...
        'Color', st.condition_color(names{c}), 'Interpreter', 'tex');
    xlabel(ax_c(c), lam_name, 'Interpreter', 'tex');
    % At about 230 px a panel, the source's seven labels ran together
    % ('-0.1 0 0.1'): every tick stays, the ones at +-0.1 unlabelled.
    lbl = cellstr(ax_c(c).XTickLabel);
    lbl(ismember(lbl, {'-0.1', '0.1'})) = {''};
    ax_c(c).XTickLabel = lbl;
    if c > 1; ax_c(c).YTickLabel = {}; ylabel(ax_c(c), ''); end
end
ylabel(ax_c(1), {'Dendritic share', 'of the leading', 'perturbation'}, 'Interpreter', 'none');   % on two lines it was taller than the panel
% The y labels beside their tick labels (the copies kept the sources' places).
set([ax_a.YLabel, ax_c(1).YLabel], 'PositionMode', 'auto');

%% The common canvas (the layout is on canvas_w already: fit moves nothing)
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);
for pass = 1:2               % the key settles its size only after a redraw: beside (B), at its top
    drawnow;
    p = lg.Position;
    lg.Position = [XB + WB + 6, Y_AB + H_AB - p(4), p(3), p(4)];
end

%% Panel letters
T.letter(ax_a, 'A', -MA + 2, 30);
T.letter(ax_b, 'B', -MB + 2, 30);
T.letter(ax_c(1), 'C', -MC + 2, 50);
T.margin(fig, 4);            % the x labels of (C) are the lowest thing drawn

%% Save
panels = struct( ...
    'A', 'fig_sfa_EOC_allStd (first figure): lambda_1 against the longest SFA time constant, tau sweep', ...
    'B', 'fig_sfa_EOC_allStd (second figure): leading-vector squared-norm fractions against the longest SFA time constant', ...
    'C', 'fig_leading_direction: dendritic share of the leading perturbation against lambda_1, joint sample');
sources = [{S_tau.source}, cellstr(string(S_ld.source(~cellfun(@isempty, S_ld.source))))];
parts = cell(0, 2);
if cfg.save
    parts = {['Panel A, ' panels.A], T.axes_table(ax_a, 'The median over repetitions at each level of the sweep'); ...
             ['Panel B, ' panels.B], T.axes_table(ax_b, 'The squared-norm fraction of the leading Lyapunov vector in each block of state variables'); ...
             ['Panel C, ' panels.C], S_ld.table_text};
end
caption = sprintf(['Supplementary figure on the leading direction (drawn by fig_supp_leading_grouped; run `%s`; ' ...
    'conditions %s). (A) The leading Lyapunov exponent lambda_1 against the longest SFA time constant in the tau ' ...
    'sweep (%s): median line, interquartile band and gray binned counts; the lower y limit is -0.2 and the mean ' ...
    'markers of the source figure are not drawn. (B) The share of the squared norm of the leading Lyapunov vector ' ...
    'in the spike-frequency adaptation, depression and dendritic coordinates against the same time constant (the ' ...
    'state coordinates have their own colours and markers). (C) For each condition, the share of the leading ' ...
    'perturbation in the dendritic coordinates against lambda_1, one point per configuration of the joint sample ' ...
    '(symmetric-logarithmic axis, linear within 0.05 1/s of zero); dashed: zero growth; the uncoupled rate is not ' ...
    'drawn (it is in the table). Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, strjoin(names, ', '), tau_name, getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary figure: the leading direction', caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function name = tau_condition(tau_dir, requested)
% The condition the tau sweep drew: the selection intersected with what the
% sweep holds (fig_sfa_EOC_allStd's rule), which the caller asserted is one.
f = fullfile(tau_dir, 'param_space_summary.mat');
if isfile(f)
    S = load(f, 'summary_data'); conds = S.summary_data.conditions;
else
    psa = ParamSpaceAnalysis2.from_dir(tau_dir); conds = psa.conditions;
end
held = cellfun(@(c) c.name, conds, 'UniformOutput', false);
hit = requested(ismember(requested, held));
assert(isscalar(hit), 'fig_supp_leading_grouped:TauCondition', ...
    'The tau sweep holds %d of the selected conditions; panels (A) and (B) draw one.', numel(hit));
name = hit{1};
end

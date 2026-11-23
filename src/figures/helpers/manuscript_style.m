function s = manuscript_style()
% MANUSCRIPT_STYLE Shared look for every Stability_Manuscript figure.
%
%   s = MANUSCRIPT_STYLE();
%   set(ax, 'FontSize', s.tick_fs);
%   xlabel(ax, '\lambda_1', 'FontSize', s.label_fs, 'Interpreter', 'tex');
%
% Returns constants; it does NOT touch global state. That is the whole point.
%
% WHY THIS EXISTS, AND WHY IT RETURNS RATHER THAN SETS
%
% tick_fs = 14 and label_fs = 15.4 were duplicated in nine figure scripts, the
% Okabe-Ito condition palette in three, and the condition display names in about
% five. Any of them could drift and nothing would notice.
%
% Worse, two scripts (plot_memory_capacity_combined and
% Fig_memory_capacity_example) set FOUR GRAPHICS-ROOT DEFAULTS and never
% restored them:
%
%   set(0,'DefaultAxesFontSize',14); set(0,'DefaultAxesLineWidth',1.0);
%   set(0,'DefaultTextInterpreter','none'); set(0,'DefaultLegendInterpreter','none');
%
% Run standalone that is harmless. Run from a master script that renders fifteen
% figures in sequence, every figure AFTER memory capacity silently inherits them
% -- including DefaultTextInterpreter = 'none', which would break the
% \lambda_1 and \mu_{EE} tex labels the sensitivity sheets depend on. The bug
% only appears once the figures are batched, and order-dependent figure output
% is miserable to diagnose.
%
% So: this function hands back values to apply per-object. If you genuinely need
% root defaults (some third-party plotters read them), use
% with_manuscript_defaults below, which restores them via onCleanup.
%
% See also: with_manuscript_defaults, write_figure_manifest, save_figure_stable

s = struct();

% --- Type sizes. tick_fs is the MC figure's DefaultAxesFontSize; label_fs is
% that times MATLAB's 1.1 label multiplier, which is how the two families came
% to match in the first place.
s.tick_fs  = 14;
s.label_fs = 15.4;
s.title_fs = 14;

% --- The common canvas (Tom, 2026-11-22: "The figures should have similar
% canvas widths so that the letters appear proportionally the same"). Every
% main and supplementary composite is laid out on one canvas width, with the
% type sizes above in points and never scaled with a block, so letters and
% labels have the same size relative to the figure in every figure. The
% drawing fills the width (the cropped image is within canvas_tol of it).
% Height is free up to canvas_max_h, a page at the same scale.
%
% 9 in, not the working width of 13 in (Tom, 2026-11-22, during the figure
% review: "in general, we should probably get all figures down to a canvas
% width of 9 inches"; on Figure 5, 13 in read as too big a canvas, with the
% panel letters too small). The type sizes stay in points, so letters and
% labels are about 13/9 = 1.4 times larger relative to the figure, and every
% layout is made again for the narrower canvas, not scaled. canvas_max_h and
% canvas_tol keep their proportion to the width (17 and 0.5 at 13 in).
% Read these numbers from here; never type them into a figure function.
s.canvas_w     = 9;      % inches
s.canvas_tol   = 0.35;   % inches the cropped image may fall short of canvas_w
s.canvas_max_h = 11.8;   % inches
s.letter_fs    = 14;     % panel letters, (A), (B), ...

% --- Line weights
s.axis_lw     = 1.0;    % axis lines + tick marks
s.line_lw     = 2.0;    % data traces
s.median_lw   = 3.0;    % median overlay on the sensitivity sheets
s.zeroline_lw = 2.0;    % the green dashed lambda_1 = 0 line

% --- Marks
% Opaque dark blue, not transparent pure blue: plot_sensitivity draws the median
% at alpha 0.35 so it survives near-black cells, and the colormap no longer ramps
% into black, so the transparency is not needed and only muddies the line.
s.median_color   = [0.00 0.00 0.55];
s.zeroline_color = [0.00 0.60 0.00];
s.divider_color  = [0.80 0.80 0.80];   % vertical gray column dividers
s.default_mark_color = [0.55 0.35 0.35];  % reddish-gray preset-default tick
s.default_mark_frac  = 0.05;              % its height, as a fraction of the y-range
s.default_mark_lw    = 2.0;

% --- Adaptation conditions. Okabe-Ito, colorblind-safe. Reddish-purple rather
% than bluish-green keeps all four hues separated (sky-blue and bluish-green were
% too close) and avoids green, which is reserved for the zero line.
%
% ONE EXCEPTION, BY DECISION (Tom, 2026-11-11, FD-007): the product MTS,
% sfa3_std2, is green [0 0.7 0]. Green was kept out of this palette for
% red-green colour blindness (and because the lambda_1 = 0 line is green).
% Product MTS appears only in the four-condition figure set, which is for
% Tom's own review, not the paper, so that reason does not apply there. The
% paper's three-condition set has no green condition: its MTS is the chain,
% sfa3_stdc2, in the blue sfa3_std2 had until then. A figure regenerated
% from a run made before the chain (which has only sfa3_std2) therefore
% draws its MTS column in green.
%
% Keyed BY CONDITION NAME so a run that declares its conditions in a different
% order cannot silently recolour a figure. Both the sweep pipeline's snake_case
% names and the memory-capacity figures' display names are present, because the
% two halves of the paper name the same four regimes differently.
% THE ONE-TIMESCALE REGIMES SHARE THEIR MECHANISM'S HUE, lightened. A reader
% should see at a glance that sfa_only_oneTS and sfa_only are the same mechanism
% differing in timescale structure, not two unrelated conditions -- which giving
% them fresh Okabe-Ito hues would have implied. sfa3_std1 sits between SFA+STD
% and the STD family for the same reason.
%
% These maps are indexed UNGUARDED by fig_sensitivity_medians and
% fig_memory_capacity_example, so a missing key throws
% MATLAB:Containers:Map:NoKey rather than degrading. Adding a condition means
% adding it here.
% TR, 2026-09-15: the multi-timescale blue is Okabe-Ito BLUE (#0072B2), not the
% sky blue (#56B4E9) used until then, which read too light on white; the
% single-timescale orange is the same hue darkened by 20% (0.8 x #E69F00).
sfa_hue  = 0.8 * [0.902 0.624 0.000];
std_hue  = [0.000 0.447 0.698];
both_hue = [0.800 0.475 0.655];
lighten  = @(c) c + 0.45 * (1 - c);

% Keys are the sfaX_stdY names; the pre-2026-09-03 names are kept below so a
% figure regenerated from an older run still gets its colours instead of
% throwing NoKey.
%
% THE PAPER'S TWO ADAPTING REGIMES GET THE TWO SATURATED HUES (TR, 2026-09-09):
% sfa1_std1 orange, sfa3_std2 sky blue. In the 3-condition figures that is
% black / orange / blue (dark red / orange / blue since FR-004, below) --
% three fully distinct colours for the single-vs-multi comparison, rather
% than one purple at two lightnesses.
%
% THE COST, accepted: hue no longer says which MECHANISM a regime carries. It
% used to -- SFA orange, STD blue, both purple, the one-timescale variant
% lighter -- and the 4- and 7-condition sheets relied on that. In those sheets
% sfa3_std2 now shares sky blue with sfa0_std2, and sfa1_std1 shares orange with
% sfa3_std0. Those sets are exploratory and not in the paper; the paper's
% comparison is what these two rows are tuned for.
%
% sfa3_std1 keeps the old scheme (light purple) and is the one row still
% following it; sfa3_std3 keeps purple.
%
% SINCE 2026-11-12 (FD-007) the blue belongs to the CHAIN, sfa3_stdc2, which
% took over the product's colour, position and title as the paper's MTS; the
% product sfa3_std2 is green (see the exception above). So in the 4- and
% 7-condition product sheets sfa3_std2 is now green, not blue.
%
% NO ADAPTATION IS DARK RED, not black (Tom, 2026-11-20; manuscript FR-004).
% The value is the red the Ritz overlay of fig_eig_heatmap substituted for
% it from FD-033 on, where black could not be read on the gray eigenvalue
% density; it read well beside the orange and blue of the adapting regimes,
% and it makes the three conditions three hues. Red and blue therefore also
% serve the E and I neurons (E reddish, I bluish): Tom accepted the double
% use as long as the two are kept apart on a panel. The display name
% 'Baseline' below follows it. Figures drawn before this change show the
% condition in black.
no_adaptation_red = [0.75 0.10 0.15];
cond_colors = { ...
    'no_adaptation',    no_adaptation_red; ...
    'sfa1_std0',        lighten(sfa_hue); ...
    'sfa3_std0',        sfa_hue; ...
    'sfa0_std1',        lighten(std_hue); ...
    'sfa0_std2',        std_hue; ...
    'sfa1_std1',        sfa_hue; ...          % the paper's single-timescale regime
    'sfa3_std1',        lighten(both_hue); ...
    'sfa3_stdc2',       std_hue; ...          % the paper's multi-timescale regime (chain)
    'sfa3_std2',        [0 0.7 0]; ...        % product MTS, four-condition set only
    'sfa0_stdc2',       std_hue; ...          % the chain alone (single-neuron preset)
    'sfa3_std3',        both_hue; ...
    'sfa0_std1_stf1',   lighten(std_hue); ...
    'sfa1_std1_stf1',   lighten(both_hue); ...
    ... % legacy names, for runs saved before the rename
    'sfa_only_oneTS',   lighten(sfa_hue); ...
    'sfa_only',         sfa_hue; ...
    'std_only_oneTS',   lighten(std_hue); ...
    'std_only',         std_hue; ...
    'sfa_and_std',      both_hue; ...
    ... % display strings, used where a figure keys on the printed label
    'Baseline',         no_adaptation_red; ...
    'SFA',              sfa_hue; ...
    'STD',              std_hue; ...
    'SFA+STD',          both_hue};
s.condition_color = containers.Map(cond_colors(:, 1), cond_colors(:, 2));

% Shared with ParamSpaceAnalysis2's plotters, so a manuscript sheet and a sweep
% sheet title the same regime identically.
s.condition_title = srnn_condition_titles();

cond_short = { ...
    'no_adaptation',    'Baseline'; ...
    'sfa1_std0',        'SFA-1'; ...
    'sfa3_std0',        'SFA'; ...
    'sfa0_std1',        'STD-1'; ...
    'sfa0_std2',        'STD'; ...
    'sfa1_std1',        'SFA-1+STD-1'; ...
    'sfa3_std1',        'SFA+STD-1'; ...
    'sfa3_stdc2',       'SFA+STD'; ...        % Tom's words: the chain is "the" STD
    'sfa3_std2',        'SFA+prodSTD'; ...    % ... and the product is named as such
    'sfa0_stdc2',       'STD'; ...            % the chain alone, by the same rule
    'sfa3_std3',        'SFA+STD-3'; ...
    'sfa0_std1_stf1',   'STD-1+STF'; ...
    'sfa1_std1_stf1',   'SFA-1+STD-1+STF'; ...
    ... % legacy
    'sfa_only_oneTS',   'SFA-1'; ...
    'sfa_only',         'SFA'; ...
    'std_only_oneTS',   'STD-1'; ...
    'std_only',         'STD'; ...
    'sfa_and_std',      'SFA+STD'};
s.condition_short = containers.Map(cond_short(:, 1), cond_short(:, 2));

% --- The two cell types. ONE E color and ONE I color for every manuscript
% figure that marks a neuron, a trace or a statistic by its cell type with a
% single color (polish list of 2026-11-22: Figure 8 took the model's orange
% and blue from SRNNCellTypePairs.type_colors, the insets of Figure 5 the
% dark rose and teal of excitatory_colormap(1) / inhibitory_colormap(1)).
% E reddish, I bluish (docs/Figure_style_guide.md). The values are the
% model's own two anchors, SRNNCellTypePairs.type_colors(2), so the model's
% time-series plots and the manuscript figures agree;
% test_composite_spectral_figures pins the equality. The I blue is close to
% the multiple-timescale condition blue and the E orange-red lies between
% the no-adaptation red and the single-timescale yellow: a figure that uses
% both encodings keeps them apart on the panel (an inset, a separate axes),
% as Tom accepted on 2026-11-20. The many-neuron palettes
% (excitatory_colormap, inhibitory_colormap) are a different thing: shades
% that tell the neurons of one type apart.
s.E_color = [0.8500 0.3250 0.0980];
s.I_color = [0.0000 0.4470 0.7410];

% --- Metric labels, so the axis of a given quantity reads the same everywhere.
s.label_lle  = '\lambda_1';
s.label_rate = 'Mean Pop. Firing Rate';

% --- Saving
s.save_types = {'png', 'svg', 'fig'};
end

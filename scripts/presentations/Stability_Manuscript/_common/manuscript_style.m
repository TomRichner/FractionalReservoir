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
% too close) and avoids green entirely, which is reserved for the zero line.
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
sfa_hue  = [0.902 0.624 0.000];
std_hue  = [0.337 0.706 0.914];
both_hue = [0.800 0.475 0.655];
lighten  = @(c) c + 0.45 * (1 - c);

s.condition_color = containers.Map( ...
    {'no_adaptation', 'sfa_only_oneTS', 'sfa_only', 'std_only_oneTS', ...
     'std_only', 'sfa3_std1', 'sfa_and_std', ...
     'Baseline',      'SFA',    'STD',      'SFA+STD'}, ...
    {[0 0 0], lighten(sfa_hue), sfa_hue, lighten(std_hue), ...
     std_hue, lighten(both_hue), both_hue, ...
     [0 0 0], sfa_hue, std_hue, both_hue});

% Shared with ParamSpaceAnalysis2's plotters, so a manuscript sheet and a sweep
% sheet title the same regime identically.
s.condition_title = srnn_condition_titles();

s.condition_short = containers.Map( ...
    {'no_adaptation', 'sfa_only_oneTS', 'sfa_only', 'std_only_oneTS', ...
     'std_only', 'sfa3_std1', 'sfa_and_std'}, ...
    {'Baseline', 'SFA-1', 'SFA', 'STD-1', 'STD', 'SFA+STD-1', 'SFA+STD'});

% --- Metric labels, so the axis of a given quantity reads the same everywhere.
s.label_lle  = '\lambda_1';
s.label_rate = 'Mean Firing Rate';

% --- Saving
s.save_types = {'png', 'svg', 'fig'};
end

function out = fig_FI_curve(cfg)
% FIG_FI_CURVE How SFA and STD reshape the F-I curve (conceptual).
%
%   out = FIG_FI_CURVE()
%   out = FIG_FI_CURVE('out_dir', d, 'save', false)
%
% From the model, the synaptic output of a neuron is
%       br = b * phi(x - c*a)
% The two adaptation mechanisms enter in structurally different ways, so they
% deform the input-output (F-I) curve differently:
%   * SFA (a: 0 -> 1) subtracts c*a INSIDE phi  -> HORIZONTAL shift (bias/threshold)
%   * STD (b: 1 -> 0) multiplies the whole curve -> VERTICAL rescale (gain/slope)
%
% ANALYTIC: it evaluates the static br(x) equation; no simulation is run. Left
% panel varies a (b fixed = 1); right panel varies b (a fixed = 0).
%
% DELIBERATELY NOT PRESET-DRIVEN. Every other figure was moved onto the paper's
% preset; this one was NOT, by decision. It is a schematic of the FORM of the
% equations rather than a picture of a particular network, so it keeps the
% logistic phi at S_c = 0.4 and an ENLARGED c = 0.6. The model's own c_E =
% 0.5/3 is far too small to see, and the point being made -- shift versus
% rescale -- is independent of the operating point. Do not "fix" this into
% consistency with the preset; the inconsistency is the design.
%
% See also: fig_SFA_steady_state, fig_STD_steady_state, paper_config

arguments
    cfg.out_dir (1,:) char    = ''
    cfg.save    (1,1) logical = true
    cfg.visible (1,1) logical = true
    cfg.run_dir     (1,:) char = ''    % unused; accepted for a uniform call
    cfg.preset_name (1,:) char = ''    % unused; see the note above
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

%% ---- Parameters -----------------------------------------------------------
S_c = 0.4;                                     % activation center (class default)
phi = @(x) SRNNModel2.logisticSigmoid(x, S_c); % exactly the model's activation

x = linspace(-0.6, 1.4, 400);                  % input (spans the active region)

% Conceptual adaptation strength. The model default is c_E = 0.15/3 ~ 0.05, too
% small to see; enlarged here for visibility. The functional form b*phi(x - c*a)
% is exact -- only the magnitude of c is chosen for clarity.
c_concept = 0.6;

a_levels = [0, 0.25, 0.5, 0.75, 1.0];   % SFA panel: b fixed = 1
b_levels = [1.0, 0.75, 0.5, 0.25, 0.0]; % STD panel: a fixed = 0

% Type sizes and weights come from the shared manuscript_style so the nine
% figures that used to hardcode tick_fs = 14 / label_fs = 15.4 cannot drift.
tick_fs   = st.tick_fs;
label_fs  = st.label_fs;
title_fs  = 16;          % this figure alone uses a larger panel title
lw        = st.line_lw;

% Sequential color ramps (dark -> light) indexed by mechanism strength.
sfa_cmap = copper(numel(a_levels) + 1);   % warm ramp for SFA
std_cmap = winter(numel(b_levels) + 1);   % cool ramp for STD

%% ---- Figure ---------------------------------------------------------------
fig = figure('Color', 'white', 'Position', [4429, 565, 623, 322]);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Left panel: SFA -> horizontal bias shift (b = 1) -------------------------
ax1 = nexttile(tl); hold(ax1, 'on');
for k = 1:numel(a_levels)
    a = a_levels(k);
    y = 1.0 * phi(x - c_concept * a);   % b = 1
    plot(ax1, x, y, 'LineWidth', lw, 'Color', sfa_cmap(k, :));
end
hold(ax1, 'off');
box(ax1, 'off');
set(ax1, 'FontSize', tick_fs);
xlabel(ax1, 'x', 'FontSize', label_fs);
ylabel(ax1, 'synaptic output  b\cdot r', 'FontSize', label_fs);
title(ax1, 'SFA: bias shift', 'FontWeight', 'normal', 'FontSize', title_fs);

%% Right panel: STD -> vertical gain rescale (a = 0) ------------------------
ax2 = nexttile(tl); hold(ax2, 'on');
for k = 1:numel(b_levels)
    b = b_levels(k);
    y = b * phi(x);   % a = 0
    plot(ax2, x, y, 'LineWidth', lw, 'Color', std_cmap(k, :));
end
hold(ax2, 'off');
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs);
xlabel(ax2, 'x', 'FontSize', label_fs);
title(ax2, 'STD: gain rescale', 'FontWeight', 'normal', 'FontSize', title_fs);

%% Shared scaling + annotations ---------------------------------------------
linkaxes([ax1, ax2], 'xy');
xlim(ax1, [x(1), x(end)]);
ylim(ax1, [0, 1]);
set([ax1, ax2], 'YTick', [0, 1], 'XTick', [0, 1]);

% SFA: 'Increasing SFA' with a right-pointing arrow BELOW the text (bias shifts right).
annotation(fig, 'textbox', [0.045, 0.58, 0.28, 0.08], ...
    'String', 'Increasing SFA', 'EdgeColor', 'none', 'FitBoxToText', 'off', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 12, 'Color', [0 0 0]);
annotation(fig, 'arrow', [0.115, 0.255], [0.565, 0.565], ...
    'Color', [0.35 0.2 0], 'HeadStyle', 'vback2', 'LineWidth', 1);

% STD: 'Increasing STD' with a downward arrow below the text (gain shrinks, black text).
annotation(fig, 'textarrow', [0.70, 0.70], [0.72, 0.46], ...
    'String', 'Increasing STD', 'FontSize', 12, ...
    'TextColor', [0 0 0], 'HeadStyle', 'vback2', 'Color', [0 0.25 0.5]);

%% ---- Save (stable filenames) ----------------------------------------------
fig_tag = 'Fig_FI_curve';

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% ---- Save ------------------------------------------------------------------
out = struct('figs', fig, 'files', {{}}, 'source', 'analytic (no simulation)');
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);

end
end




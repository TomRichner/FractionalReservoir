close all
clc

% Fig_FI_curve.m
% Stability_Manuscript conceptual figure: how SFA and STD reshape the F-I curve.
%
% From SRNNModel2 the synaptic output of a neuron is
%       br = b * phi(x - c*a)
% with the class-default logistic activation phi(x) = 1/(1+exp(-4*(x - S_c))).
% The two adaptation mechanisms enter in structurally different ways, so they
% deform the input-output (F-I) curve differently:
%   * SFA (a: 0 -> 1) subtracts c*a INSIDE phi  -> HORIZONTAL shift (bias/threshold)
%   * STD (b: 1 -> 0) multiplies the whole curve -> VERTICAL rescale (gain/slope)
%
% This is an analytic figure -- it just evaluates the static br(x) equation, no
% SRNNModel2 simulation is run. Left panel varies a (b fixed = 1); right panel
% varies b (a fixed = 0).

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_FI_curve -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

setup_paths();

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

tick_fs   = 14;    % tick numbers   -- matches other Stability_Manuscript fig_ figures
label_fs  = 15.4;  % axis labels    -- matches other fig_ figures
title_fs  = 16;    % panel titles
lw        = 2;

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
old = dir(fullfile(this_dir, [fig_tag '*']));
for a = 1:numel(old)
    fp = fullfile(old(a).folder, old(a).name);
    if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
        delete(fp);
    end
end
save_some_figs_to_folder_2(this_dir, fig_tag, fig.Number, {'png', 'svg', 'fig'});
num = num2str(fig.Number);
movefile(fullfile(this_dir, [fig_tag '_figure_' num '.png']), fullfile(this_dir, [fig_tag '.png']));
movefile(fullfile(this_dir, [fig_tag '_figure_' num '.svg']), fullfile(this_dir, [fig_tag '.svg']));
movefile(fullfile(this_dir, [fig_tag '_f_' num '.fig']),      fullfile(this_dir, [fig_tag '.fig']));

% Git provenance alongside the figure.
capture_git_provenance(this_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(this_dir, 'README_fig_FI_curve.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: SFA & STD reshape the F-I curve\n');
fprintf(fid, '===========================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'WHAT IT SHOWS\n');
fprintf(fid, '  Conceptual, ANALYTIC figure (no SRNNModel2 simulation): plots the\n');
fprintf(fid, '  synaptic output  br = b * phi(x - c*a)  vs input x, with\n');
fprintf(fid, '  phi(x) = 1/(1+exp(-4*(x - S_c))), S_c = %.2f (class-default logistic).\n\n', S_c);
fprintf(fid, '  Left panel  (SFA): b = 1, a swept 0->1. Subtracting c*a inside phi\n');
fprintf(fid, '    shifts the curve HORIZONTALLY -> raises threshold / changes bias.\n');
fprintf(fid, '  Right panel (STD): a = 0, b swept 1->0. The b prefactor rescales the\n');
fprintf(fid, '    curve VERTICALLY -> reduces gain/slope and the saturation ceiling.\n\n');
fprintf(fid, 'CONCEPTUAL PARAMETER\n');
fprintf(fid, '  c = %.2f is enlarged for visibility (model default c_E = 0.15/3 ~ 0.05\n', c_concept);
fprintf(fid, '  is too small to see). The functional form b*phi(x - c*a) is exact.\n');
fprintf(fid, '  a levels: %s ; b levels: %s\n\n', mat2str(a_levels), mat2str(b_levels));
fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
fprintf(fid, '  Fig_FI_curve.png\n  Fig_FI_curve.svg\n  Fig_FI_curve.fig\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);

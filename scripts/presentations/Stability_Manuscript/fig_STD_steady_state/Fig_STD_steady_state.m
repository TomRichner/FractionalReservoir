close all
clc

% Fig_STD_steady_state.m
% Stability_Manuscript conceptual figure: the steady-state depression of the
% DUAL-TIMESCALE STD preset, and the synaptic output it leaves behind.
%
% Setting db/dt = 0 in
%       db_k/dt = (1 - b_k)/tau_rec_k - b_k*r/tau_rel_k
% gives, for each timescale k,
%       b_k(r) = 1 / (1 + r * tau_rec_k / tau_rel_k)
% and the synapse multiplies them, so it sees prod_k b_k(r).
%
% Analytic figure -- no simulation. The taus are READ OFF THE PRESET rather than
% hardcoded, so this cannot drift away from what is actually being swept.
%
% Left panel:  prod(b) vs r  -- how much gain is left at a given rate.
% Right panel: prod(b)*r vs r -- the synaptic output itself, which is what the
%   recurrent sum actually receives. It is NON-MONOTONIC: depression eventually
%   outruns the rate driving it, so beyond a peak, firing harder delivers LESS
%   to the network. That peak is the interesting quantity here.
%
% TWO VERSIONS are written, differing only in the right panel's y range: the
% square one (y over [0, 1], daspect [1 1 1]) and a zoomed one that fills the
% panel with the curves and drops the 1:1 aspect ratio to do it. See the
% comment above the figure loop.
%
% See also: Fig_FI_curve (the same mechanisms as a deformation of the F-I curve)

this_dir = fileparts(mfilename('fullpath'));
setup_paths();
project_root = fileparts(which('setup_paths'));

%% ---- Parameters -----------------------------------------------------------
preset_name = 'celltype_pairs_Sc0p2_noise0p025_dualStd';

% Pull the depression timescales out of the preset's own STD routes. They live
% on the conditions, not on the model_defaults struct: synapse_config can only
% reach the model through a condition, so that is where a preset puts them.
[~, ~, conditions] = srnn_param_preset(preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
sc = conditions{strcmp(cond_names, 'sfa_and_std')}.synapse_config;
route = sc.E.E.std;                 % uniform across all four routes in this preset
tau_rec = route.tau_rec(:)';
tau_rel = route.tau_rel(:)';
ratio = tau_rec ./ tau_rel;         % the only combination that sets the steady state

r = linspace(0, 1, 400);            % rate, over the full range of the nonlinearity

% b_k(r) for each timescale, then their product. rows = timescale, cols = r.
b_each = 1 ./ (1 + ratio(:) * r);
b_prod = prod(b_each, 1);

tick_fs  = 14;    % tick numbers -- matches the other Stability_Manuscript figures
label_fs = 15.4;  % axis labels
title_fs = 16;    % panel titles
lw       = 2;

prod_color     = [0.85 0.325 0.098];   % warm, matching the E colour used elsewhere
single_color   = [0.5 0.5 0.5];
identity_color = [0.55 0.80 0.55];     % faint green: y = r, the undepressed synapse

%% ---- Derived quantities ---------------------------------------------------
output = b_prod .* r;
output_single = b_each(1, :) .* r;
[peak_val, peak_idx] = max(output);
r_peak = r(peak_idx);

% The zoomed version's ceiling: just above the largest value a SINGLE b_k ever
% delivers, which on this preset is its value at r = 1 (b_k*r rises
% monotonically to the asymptote tau_rel/tau_rec, so the right edge is the max).
% Rounded up to the next hundredth so the tick is a round number rather than
% 0.1111..., and so the curve does not touch the top of the axes.
zoom_ymax = ceil(max(output_single) * 100) / 100;

%% ---- Figures --------------------------------------------------------------
% TWO VERSIONS of the same figure, differing only in the right panel's y range:
%
%   'square' : y over [0, 1] with daspect [1 1 1]. The identity line y = r is a
%              literal 45 degrees and a vertical distance can be compared
%              against a horizontal one by eye. Honest but small: both curves
%              sit in the bottom tenth of the panel.
%   'zoom'   : y over [0, zoom_ymax], filling the panel with the curves
%              themselves. The 1:1 aspect ratio is DELIBERATELY DROPPED here --
%              it cannot survive a y range 9x shorter than x's without squashing
%              the panel to a sliver, so the identity line is no longer 45
%              degrees and no cross-axis angle in this version means anything.
%              What it buys is the turnover, which is the point of the panel and
%              is nearly invisible at full scale.
%
% Size copied from Fig_FI_curve so the two figures match, but the POSITION is
% computed rather than copied: that script hardcodes x = 4429, a second-monitor
% coordinate from the machine it was written on, which lands the window
% off-screen on a single 1920-wide display and looks like the figure never
% opened. Centre it on whatever screen is actually attached.
variants = {'square', 'zoom'};
fig_size = [623, 322];
scr = get(groot, 'ScreenSize');
figs = gobjects(1, numel(variants));

for v = 1:numel(variants)
is_zoom = strcmp(variants{v}, 'zoom');

% Offset each window so the second does not land exactly on the first.
fig_pos = [scr(1:2) + max((scr(3:4) - fig_size)/2, 0) + (v - 1)*[40, -40], fig_size];
fig = figure('Color', 'white', 'Position', fig_pos);
figs(v) = fig;
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Left panel: prod(b) vs r --------------------------------------------------
% Identical in both versions: prod(b) is a gain, already bounded by [0, 1], so
% there is nothing to zoom into.
ax1 = nexttile(tl); hold(ax1, 'on');
% One component b_k, thin and dashed, so the gap to the product reads as the
% cost of the SECOND timescale rather than as depression in general. Both
% timescales share tau_rec/tau_rel here, so the components coincide -- drawing
% one is drawing both (the code plots b_1).
plot(ax1, r, b_each(1, :), '--', 'LineWidth', 1, 'Color', single_color);
plot(ax1, r, b_prod, 'LineWidth', lw, 'Color', prod_color);
hold(ax1, 'off');
box(ax1, 'off');
set(ax1, 'FontSize', tick_fs);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'depression  $\prod_k b_k$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax1, 'Steady-state depression', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax1, [0, 1]); ylim(ax1, [0, 1]);
set(ax1, 'XTick', [0, 1], 'YTick', [0, 1]);
legend(ax1, {'single $b_k$', '$\prod_k b_k$'}, 'Interpreter', 'latex', ...
    'Box', 'off', 'FontSize', 12, 'Location', 'northeast');

%% Right panel: prod(b)*r vs r -----------------------------------------------
ax2 = nexttile(tl); hold(ax2, 'on');
% Identity first, so it sits UNDER the data. y = r is the undepressed synapse
% (b == 1): the vertical gap down to each curve is exactly what depression costs
% at that rate. In the zoomed version it leaves the top of the panel almost
% immediately -- that steep exit IS the depression, seen edge-on.
plot(ax2, [0, 1], [0, 1], '-', 'LineWidth', 1, 'Color', identity_color);
% Same pairing as the left panel: one timescale dashed, the product solid. The
% two curves differ in KIND here, not just in scale -- a single b_k rises
% monotonically toward its asymptote 1/(tau_rec/tau_rel), so more rate always
% delivers more, whereas the product turns over.
plot(ax2, r, output_single, '--', 'LineWidth', 1, 'Color', single_color);
plot(ax2, r, output, 'LineWidth', lw, 'Color', prod_color);
plot(ax2, r_peak, peak_val, 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', prod_color, 'MarkerEdgeColor', 'none');
hold(ax2, 'off');
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs);
xlabel(ax2, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax2, 'synaptic output  $\prod_k b_k \cdot r$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax2, 'Delivered output', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax2, [0, 1]);
if is_zoom
    % No daspect: see the header comment above. The y tick is zoom_ymax rather
    % than 1, which is the one number a reader must notice to avoid reading this
    % panel as if it were the square one.
    ylim(ax2, [0, zoom_ymax]);
    set(ax2, 'XTick', [0, 1], 'YTick', [0, zoom_ymax]);
else
    % Equal pixel scaling on both axes, which forces y to share x's [0, 1]
    % range: the curves then sit low in a square panel, which is the honest
    % picture -- at every rate the depressed synapse delivers a small fraction
    % of what an undepressed one would.
    ylim(ax2, [0, 1]);
    daspect(ax2, [1 1 1]);
    set(ax2, 'XTick', [0, 1], 'YTick', [0, 1]);
end

% Mark where the output peaks -- past this rate the synapse delivers less.
% Above-right of the marker: there is headroom there in the square version
% (the title is far off at y = 1) and in the zoomed one (the peak sits below
% the single-b_k curve it is being contrasted with).
text(ax2, r_peak, peak_val, sprintf('  r = %.2f', r_peak), ...
    'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

end

%% ---- Save (stable filenames) ----------------------------------------------
fig_tag_base = 'Fig_STD_steady_state';
fig_tags = {fig_tag_base, [fig_tag_base '_zoom']};

% One sweep for both versions, BEFORE any saving: the base tag is a prefix of
% the zoom tag, so a per-variant delete of '<base>*' would wipe the zoom files
% written on the previous pass.
old = dir(fullfile(this_dir, [fig_tag_base '*']));
for a = 1:numel(old)
    fp = fullfile(old(a).folder, old(a).name);
    if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
        delete(fp);
    end
end

for v = 1:numel(variants)
    fig_tag = fig_tags{v};
    save_some_figs_to_folder_2(this_dir, fig_tag, figs(v).Number, {'png', 'svg', 'fig'});
    num = num2str(figs(v).Number);
    movefile(fullfile(this_dir, [fig_tag '_figure_' num '.png']), fullfile(this_dir, [fig_tag '.png']));
    movefile(fullfile(this_dir, [fig_tag '_figure_' num '.svg']), fullfile(this_dir, [fig_tag '.svg']));
    movefile(fullfile(this_dir, [fig_tag '_f_' num '.fig']),      fullfile(this_dir, [fig_tag '.fig']));
end

capture_git_provenance(this_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(this_dir, 'README_fig_STD_steady_state.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: dual-timescale STD steady state\n');
fprintf(fid, '===========================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'WHAT IT SHOWS\n');
fprintf(fid, '  Conceptual, ANALYTIC figure (no simulation). Setting db/dt = 0 in\n');
fprintf(fid, '    db_k/dt = (1 - b_k)/tau_rec_k - b_k*r/tau_rel_k\n');
fprintf(fid, '  gives b_k(r) = 1/(1 + r*tau_rec_k/tau_rel_k); the synapse multiplies\n');
fprintf(fid, '  the timescales, so it sees prod_k b_k(r).\n\n');
fprintf(fid, '  Left  : prod(b) vs r, with one component b_k dashed for reference.\n');
fprintf(fid, '          Both timescales share tau_rec/tau_rel here, so b_1 and b_2\n');
fprintf(fid, '          coincide at steady state and the dashed curve is either.\n');
fprintf(fid, '  Right : prod(b)*r vs r -- what the recurrent sum actually receives,\n');
fprintf(fid, '          again with one component b_k*r dashed for reference. The\n');
fprintf(fid, '          faint green identity line y = r is the UNDEPRESSED synapse\n');
fprintf(fid, '          (b == 1); the vertical gap down to a curve is what\n');
fprintf(fid, '          depression costs at that rate.\n\n');

fprintf(fid, 'TWO VERSIONS, differing ONLY in the right panel''s y range\n');
fprintf(fid, '  %s      : y over [0, 1] with equal pixel scaling on both\n', fig_tag_base);
fprintf(fid, '      axes (daspect [1 1 1]), so the identity line is a literal 45\n');
fprintf(fid, '      degrees and a vertical distance can be compared against a\n');
fprintf(fid, '      horizontal one by eye. Honest but small: both curves sit in the\n');
fprintf(fid, '      bottom tenth of the panel.\n');
fprintf(fid, '  %s_zoom : y over [0, %.2f] -- just above the largest value a\n', ...
    fig_tag_base, zoom_ymax);
fprintf(fid, '      SINGLE b_k ever delivers (max b_k*r = %.4f, at r = 1). The 1:1\n', ...
    max(output_single));
fprintf(fid, '      aspect ratio is DELIBERATELY DROPPED: it cannot survive a y range\n');
fprintf(fid, '      9x shorter than x''s without squashing the panel to a sliver, so\n');
fprintf(fid, '      no cross-axis angle in this version means anything and the\n');
fprintf(fid, '      identity line leaves the top of the panel at r = %.2f. What it\n', zoom_ymax);
fprintf(fid, '      buys is the turnover, which is the point of the panel and is\n');
fprintf(fid, '      nearly invisible at full scale.\n');
fprintf(fid, '  The left panel is identical in both: prod(b) is a gain, already\n');
fprintf(fid, '  bounded by [0, 1], so there is nothing to zoom into.\n\n');

fprintf(fid, 'PRESET (taus are read from it, not hardcoded)\n');
fprintf(fid, '  %s\n', preset_name);
fprintf(fid, '  tau_rec = %s\n', mat2str(tau_rec));
fprintf(fid, '  tau_rel = %s\n', mat2str(tau_rel));
fprintf(fid, '  tau_rec/tau_rel = %s  (the only combination that sets the steady state)\n\n', ...
    mat2str(ratio));

fprintf(fid, 'NUMBERS WORTH QUOTING\n');
fprintf(fid, '  prod(b) at r = 0.3 : %.4f\n', prod(1 ./ (1 + 0.3*ratio)));
fprintf(fid, '  prod(b) at r = 1   : %.4f   (a %.0fx reduction in synaptic gain)\n', ...
    prod(1 ./ (1 + ratio)), 1/prod(1 ./ (1 + ratio)));
fprintf(fid, '  output prod(b)*r peaks at r = %.3f, value %.4f\n', r_peak, peak_val);
fprintf(fid, '  Beyond that peak the output DECREASES: depression outruns the rate\n');
fprintf(fid, '  driving it, so firing harder delivers less to the network.\n');
fprintf(fid, '  A SINGLE timescale does not do this -- b_k*r rises monotonically to\n');
fprintf(fid, '  an asymptote of %.4f (= tau_rel/tau_rec), reaching %.4f at r = 1,\n', ...
    1/ratio(1), max(output_single));
fprintf(fid, '  i.e. %.1fx the product''s peak. The turnover is specific to n_b > 1.\n\n', ...
    max(output_single) / peak_val);

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for a = 1:numel(fig_tags)
    fprintf(fid, '  %s.png\n  %s.svg\n  %s.fig\n', ...
        fig_tags{a}, fig_tags{a}, fig_tags{a});
end

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);

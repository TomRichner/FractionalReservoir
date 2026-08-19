close all
clc

% Fig_SFA_steady_state.m
% Stability_Manuscript conceptual figure: the steady-state SFA current, and why
% the NUMBER of adaptation timescales does not change it.
%
% Setting da/dt = 0 in
%       da_k/dt = (-a_k + r) / tau_k
% gives a_k = r for EVERY timescale, whatever tau_k is. So
%       c * sum_k a_k  =  c * n_a * r
% (plotted as the positive adaptation variable; it enters the model as
% phi(x - c*sum a), i.e. subtracted, so a larger value here means MORE
% suppression)
% and with c split as a budget over the timescales (c = 0.5/3 for n_a = 3
% against c = 0.5 for n_a = 1) the two are the same line. That is the whole
% content of the left panel, and it is the point: the 0.15/3 convention used
% throughout this project exists precisely so that adding timescales does not
% silently change the operating point.
%
% This is the OPPOSITE of STD (see fig_STD_steady_state). SFA enters as a SUM,
% so a budget-split c makes the count invisible at steady state. STD enters as a
% PRODUCT, so n_b changes the steady state no matter how the taus are chosen --
% two timescales square the depression.
%
% Where the count does show up is the TRANSIENT, which is the right panel.
%
% Analytic figure -- no simulation.
%
% See also: Fig_STD_steady_state, Fig_FI_curve

this_dir = fileparts(mfilename('fullpath'));
setup_paths();
project_root = fileparts(which('setup_paths'));

%% ---- Parameters -----------------------------------------------------------
% Timescales follow the model's own rule, complete_type_defaults:
%   tau_a = logspace(log10(0.25), log10(10), n_a)
% which for n_a = 1 collapses to the slow end, 10 s.
tau_3 = logspace(log10(0.25), log10(10), 3);
tau_1 = logspace(log10(0.25), log10(10), 1);

% c as a BUDGET: the same total 0.5, split over however many timescales there
% are. This is the 0.15/3 convention of the sweep scripts, scaled up here.
c_3 = 0.5 / 3;
c_1 = 0.5;

r = linspace(0, 1, 400);       % rate, over the full range of the nonlinearity

% Steady state: a_k = r for every k, so sum(a) = n_a * r.
sfa_3 = c_3 * numel(tau_3) * r;
sfa_1 = c_1 * numel(tau_1) * r;

r_step = 1;                     % rate step for the transient panel: the full range
                               % of the nonlinearity, so the y-axes of the two
                               % panels share a scale and the transient lands on
                               % the left panel's endpoint.
t = linspace(0, 20, 800);

tick_fs  = 14;    % tick numbers -- matches the other Stability_Manuscript figures
label_fs = 15.4;  % axis labels
title_fs = 16;    % panel titles
lw       = 2;

three_color = [0.85 0.325 0.098];   % warm, matching the E colour used elsewhere
one_color   = [0.5 0.5 0.5];

%% ---- Figure ---------------------------------------------------------------
% Position computed rather than hardcoded: Fig_FI_curve's [4429 ...] is a
% second-monitor coordinate and lands off-screen on a single display.
fig_size = [623, 322];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', ...
    'Position', [scr(1:2) + max((scr(3:4) - fig_size)/2, 0), fig_size]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Left panel: steady-state SFA current vs rate ------------------------------
ax1 = nexttile(tl); hold(ax1, 'on');
% Drawn thick-then-dashed BECAUSE they coincide exactly. A single line would
% look like one case was forgotten; the dashed overlay is the evidence that both
% were computed and landed on top of each other.
plot(ax1, r, sfa_3, 'LineWidth', lw + 1, 'Color', three_color);
plot(ax1, r, sfa_1, '--', 'LineWidth', 1.2, 'Color', one_color);
hold(ax1, 'off');
box(ax1, 'off');
set(ax1, 'FontSize', tick_fs);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'SFA current  c\Sigmaa', 'FontSize', label_fs);
title(ax1, 'Steady-state SFA', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax1, [0, 1]); ylim(ax1, [0, 1]);
set(ax1, 'XTick', [0, 1], 'YTick', [0, 1]);
legend(ax1, ...
    {sprintf('n_a = 3,  c = %.3g', c_3), sprintf('n_a = 1,  c = %.3g', c_1)}, ...
    'Box', 'off', 'FontSize', 12, 'Location', 'northeast');

%% Right panel: the transient, where the count DOES matter ------------------
% Step r: 0 -> r_step at t = 0, from a = 0. Then a_k(t) = r_step*(1 - e^{-t/tau_k}).
sfa_3_t = c_3 * r_step * sum(1 - exp(-t ./ tau_3(:)), 1);
sfa_1_t = c_1 * r_step * sum(1 - exp(-t ./ tau_1(:)), 1);

ax2 = nexttile(tl); hold(ax2, 'on');
plot(ax2, t, sfa_3_t, 'LineWidth', lw, 'Color', three_color);
plot(ax2, t, sfa_1_t, '--', 'LineWidth', 1.2, 'Color', one_color);
yline(ax2, c_3 * numel(tau_3) * r_step, ':', 'Color', [0.6 0.6 0.6]);
hold(ax2, 'off');
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs);
xlabel(ax2, 'time (s)', 'FontSize', label_fs);
ylabel(ax2, 'SFA current  c\Sigmaa', 'FontSize', label_fs);
title(ax2, sprintf('Transient (step to r = %.2g)', r_step), ...
    'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax2, [0, t(end)]);
% Same [0, 1] range as the left panel, so the two are directly comparable and
% the transient visibly lands on the left panel's endpoint.
ylim(ax2, [0, 1]);
set(ax2, 'XTick', [0, t(end)], 'YTick', [0, 1]);

%% ---- Save (stable filenames) ----------------------------------------------
fig_tag = 'Fig_SFA_steady_state';
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

capture_git_provenance(this_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(this_dir, 'README_fig_SFA_steady_state.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: SFA steady state, 1 vs 3 timescales\n');
fprintf(fid, '===============================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'WHAT IT SHOWS\n');
fprintf(fid, '  Conceptual, ANALYTIC figure (no simulation). Setting da/dt = 0 in\n');
fprintf(fid, '    da_k/dt = (-a_k + r)/tau_k\n');
fprintf(fid, '  gives a_k = r for EVERY timescale, so c*sum(a) = c*n_a*r.\n\n');
fprintf(fid, '  Left  : c*sum(a) vs r, for n_a = 3 (c = %.4g) and n_a = 1 (c = %.4g).\n', c_3, c_1);
fprintf(fid, '  Right : the transient response to a step r: 0 -> %.2g.\n\n', r_step);

fprintf(fid, 'THE RESULT: THE TWO LINES IN THE LEFT PANEL ARE IDENTICAL\n');
fprintf(fid, '  Both are %.3g*r exactly. Because a_k = r independent of tau_k, the\n', c_3*numel(tau_3));
fprintf(fid, '  steady-state SFA current depends only on the PRODUCT c*n_a -- and\n');
fprintf(fid, '  splitting c as a budget (0.5/3 over three timescales vs 0.5 over one)\n');
fprintf(fid, '  holds that product fixed. The dashed line is drawn over the solid one\n');
fprintf(fid, '  to show both were computed; they coincide to machine precision.\n\n');
fprintf(fid, '  This is why the project uses c_E = 0.15/3: adding timescales does not\n');
fprintf(fid, '  silently move the operating point.\n\n');

fprintf(fid, 'CONTRAST WITH STD (see fig_STD_steady_state)\n');
fprintf(fid, '  SFA enters the dynamics as a SUM, so a budget-split c makes the count\n');
fprintf(fid, '  invisible at steady state. STD enters as a PRODUCT, so n_b = 2 squares\n');
fprintf(fid, '  the depression no matter how the taus are chosen -- there is no\n');
fprintf(fid, '  budget-split of tau that would make dual STD match single STD.\n\n');

fprintf(fid, 'WHERE THE COUNT DOES SHOW UP: THE TRANSIENT\n');
fprintf(fid, '  tau_a (n_a=3) = %s\n', mat2str(round(tau_3, 4)));
fprintf(fid, '  tau_a (n_a=1) = %s   (logspace collapses to the slow end at n = 1)\n', mat2str(tau_1));
fprintf(fid, '  Three timescales give a multi-exponential approach -- fast partial\n');
fprintf(fid, '  adaptation followed by a slow tail -- against a single exponential.\n');
fprintf(fid, '  Same destination, different route. Both settle to %.4f.\n\n', ...
    c_3*numel(tau_3)*r_step);

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
fprintf(fid, '  %s.png\n  %s.svg\n  %s.fig\n', fig_tag, fig_tag, fig_tag);

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);

% Report the headline result to the console too.
fprintf('\nmax |difference| between the two steady-state lines: %.3g\n', ...
    max(abs(sfa_3 - sfa_1)));

function fig = so_seeds_figure(Rs, S, out_base)
% SO_SEEDS_FIGURE Per weight seed: time beyond Re = 1 before and during the step, and bursts against crossings.
%
%   fig = so_seeds_figure(Rs, S, out_base)
%
% Rs: cell of so_run_network results (the first is the stimulation
% figure's network, seed 42, drawn in black; the others grey).
%   (A) fraction of time the rightmost eigenvalue of M has Re > 1: baseline
%       (10-40 s), during the step (50-100 s), and during the step with
%       adaptation frozen at its pre-step values
%   (B) mean Re of the rightmost eigenvalue of M during the step, per seed,
%       filled where its eigenvector is population-coherent (coherence
%       > 0.3); Re = 1 dotted
%   (C) baseline: bursts, bursts with a crossing (Re > 1 from 0.6 s before
%       onset to 0.1 s after the end), crossings outside any burst
% Writes out_base.png (300 dpi) and .fig.
%
% See also: run_stim_outlier, so_crossings

fs = 14;
ns = numel(Rs);
seed = cellfun(@(R) R.seeds(1), Rs);
fb = cellfun(@(R) R.intact.meas.before.frac_beyond, Rs);
fd = cellfun(@(R) R.intact.meas.during.frac_beyond, Rs);
ff = cellfun(@(R) R.frozen.meas.during.frac_beyond, Rs);
re_d = cellfun(@(R) R.intact.meas.during.mean_re1, Rs);
coh_d = cellfun(@(R) R.intact.meas.during.mean_coh1, Rs);
nb = cellfun(@(R) R.intact.meas.before.n_bursts, Rs);
nbx = cellfun(@(R) R.intact.meas.before.n_bursts_with_cross, Rs);
nc = cellfun(@(R) R.intact.meas.before.n_cross, Rs);
ncb = cellfun(@(R) R.intact.meas.before.n_cross_in_burst, Rs);

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 1500 470]);
set(fig, 'DefaultAxesFontSize', fs, 'DefaultAxesLineWidth', 1, 'DefaultAxesTickDir', 'out', ...
    'DefaultAxesBox', 'off', 'DefaultAxesLabelFontSizeMultiplier', 1, 'DefaultTextFontSize', fs);

ax1 = axes(fig, 'Position', [0.07 0.18 0.24 0.72]); hold(ax1, 'on');
for i = ns:-1:1
    c = [0.6 0.6 0.6]; lw = 1; if i == 1; c = [0 0 0]; lw = 2; end
    jit = (i - 1) * 0.012;
    plot(ax1, [1 2 3] + jit, [fb(i) fd(i) ff(i)], '-o', 'Color', c, 'LineWidth', lw, ...
        'MarkerFaceColor', c, 'MarkerSize', 6);
end
xlim(ax1, [0.7 3.4]); ylim(ax1, [-0.03 1.03]);
set(ax1, 'XTick', 1:3, 'XTickLabel', {'baseline', 'step', 'step, frozen'});
ylabel(ax1, 'fraction of time Re \lambda_1(M) > 1');

ax2 = axes(fig, 'Position', [0.40 0.18 0.24 0.72]); hold(ax2, 'on');
yline(ax2, 1, ':', 'Color', [0.8 0 0], 'LineWidth', 1.5);
for i = 1:ns
    c = [0.4 0.4 0.4]; if i == 1; c = [0 0 0]; end
    fc = 'w'; if coh_d(i) > 0.3; fc = c; end
    plot(ax2, i, re_d(i), 'o', 'MarkerEdgeColor', c, 'MarkerFaceColor', fc, 'MarkerSize', 9, 'LineWidth', 1.5);
end
xlim(ax2, [0.5 ns + 0.5]); ylim(ax2, [0.6 1.4]);
set(ax2, 'XTick', 1:ns, 'XTickLabel', arrayfun(@num2str, seed, 'UniformOutput', false));
xlabel(ax2, 'seed'); ylabel(ax2, 'Re \lambda_1(M), step');

ax3 = axes(fig, 'Position', [0.72 0.18 0.26 0.72]); hold(ax3, 'on');
b = bar(ax3, 1:ns, [nb; nbx; nc - ncb]', 'grouped', 'EdgeColor', 'none');
b(1).FaceColor = [0.3 0.3 0.3]; b(2).FaceColor = [0.95 0.55 0.10]; b(3).FaceColor = [0.3 0.55 0.85];
legend(ax3, {'bursts', 'bursts with crossing', 'crossings outside bursts'}, 'Box', 'off', ...
    'Location', 'northwest', 'FontSize', fs - 1);
xlim(ax3, [0.4 ns + 0.6]); ylim(ax3, [0 max([nb nc - ncb]) * 1.5 + 1]);
set(ax3, 'XTick', 1:ns, 'XTickLabel', arrayfun(@num2str, seed, 'UniformOutput', false));
xlabel(ax3, 'seed'); ylabel(ax3, 'count, baseline 10-40 s');

lab = {'(A)', ax1; '(B)', ax2; '(C)', ax3};
for i = 1:3
    p = lab{i, 2}.Position;
    annotation(fig, 'textbox', [p(1) - 0.065, p(2) + p(4) + 0.07, 0.05, 0.05], 'String', lab{i, 1}, ...
        'EdgeColor', 'none', 'FontSize', fs, 'VerticalAlignment', 'top');
end

if nargin >= 3 && ~isempty(out_base)
    exportgraphics(fig, [out_base '.png'], 'Resolution', 300);
    set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
    savefig(fig, [out_base '.fig']);
    close(fig);
end
end

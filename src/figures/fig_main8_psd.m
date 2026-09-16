function out = fig_main8_psd(source_fig_root)
% FIG_MAIN8_PSD Model PSD panel from the CURRENT run's native bursting_psd.fig.
%
%   out = FIG_MAIN8_PSD(source_fig_root)
%
% Reads <source_fig_root>/fig_stim_engages_adaptation/bursting_psd.fig -- the
% figure fig_stim_engages_adaptation saved earlier in the same figure pass --
% and copies its axes (the per-level PSD curves, legend labels 'no-stim' /
% 'stim') into panel A of a two-panel canvas; panel B is the intentionally
% empty clinical box. No PSD is re-estimated and nothing is digitized.
%
% Until 2026-09-15 this read a checksummed PNG from one archived run and
% calibrated its raster pixels into native axes; that checksum could never
% match a figure from a new run, so the grouped figure failed on every run
% except the one it was written against (TR: grouped figures must be built
% from the run they accompany). The native .fig has always been saved beside
% the PNG, so the raster path was unnecessary.
%
% See also: fig_stim_engages_adaptation, fig_grouped_main

source = fullfile(source_fig_root, 'fig_stim_engages_adaptation', 'bursting_psd.fig');
assert(isfile(source), 'fig_main8_psd:MissingSource', ...
    'No native PSD figure at %s (fig_stim_engages_adaptation must run in the same pass).', source);
src = openfig(source, 'invisible');
guard = onCleanup(@() close(src));
ax0 = findall(src, 'Type', 'axes');
ax0 = ax0(arrayfun(@(a) ~isempty(findall(a, 'Type', 'line')), ax0));
assert(isscalar(ax0), 'fig_main8_psd:BadSource', 'Expected one PSD axes in %s, found %d.', source, numel(ax0));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 780 456]);
ax  = copyobj(ax0, fig);
set(ax, 'Position', [.105 .16 .38 .74], 'Tag', 'psd_model', 'FontSize', 14, 'LineWidth', 1, 'Box', 'off');
title(ax, '');
xlabel(ax, 'frequency (Hz)', 'FontSize', 14);
ylabel(ax, 'Power spectral density of dendritic potential, x', 'FontSize', 14, 'Interpreter', 'none');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 14, 'Interpreter', 'none');

bx = axes(fig, 'Position', [.58 .16 .38 .74], 'Tag', 'psd_human_empty', ...
    'XTick', [], 'YTick', [], 'XLim', [0 1], 'YLim', [0 1], ...
    'Box', 'on', 'LineWidth', 1, 'FontSize', 14);
for a = [ax bx]
    label = '(A)'; if a == bx, label = '(B)'; end
    text(a, -.10, 1.06, label, 'Units', 'normalized', 'FontSize', 14, ...
        'Clipping', 'off', 'VerticalAlignment', 'bottom', 'FontWeight', 'normal', 'Tag', 'panel_label');
end
notes = {'Model PSD: the native axes of `fig_stim_engages_adaptation/bursting_psd.fig` from this run, copied unchanged; no PSD is re-estimated.', ...
    'Panel B is an intentionally empty box. No patient data are shown. No titles; labels (A)/(B), 14-point fonts and axes linewidth 1.0; legend inside panel A at upper right.'};
out = struct('figs', fig, 'source', {{source}}, 'notes', {notes});
end

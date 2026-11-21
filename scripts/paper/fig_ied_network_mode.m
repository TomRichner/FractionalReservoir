function out = fig_ied_network_mode(opts)
%FIG_IED_NETWORK_MODE Appendix: the full Jacobian's network mode in the modular network, at baseline and under stimulation, and before discharges.
%
%   out = fig_ied_network_mode('data_dir', run_dir, 'out_dir', d)
%
% A figure of the ied_stim bundle (called as its other entries are, by
% ied_stim_figures), drawn from its network_mode stage
% (run_ied_network_mode_stage) under <data_dir>/network_mode; simulates
% nothing. Kept apart from Figure 8 at Tom's word ("make the I and J ...
% keep it separate as I'm not sure we will keep it"; mockup panels I, J).
%
%   A, B  the real part of the rightmost network mode of the full Jacobian
%         (1/s; the boundary is 0) every full_dt s through the first
%         window of the baseline block (A) and of the stimulation block
%         (B) of the stage's first network; discharges as thin vertical
%         lines. Same y scale.
%   C     around discharges (baseline blocks, every network): the median
%         over all discharges of that real part against the time from the
%         discharge peak, and the same at control times; thin lines are
%         the per-network medians.
%   D     share of the sampled times the mode is unstable, baseline and
%         stimulation, one line per network.
% Writes fig_ied_network_mode.{png,fig,svg} and _table.md (the stage's per-
% network numbers and its cost).
%
% See also: run_ied_network_mode_stage, ied_network_mode_stage_plan

arguments
    opts.data_dir (1,:) char {mustBeNonempty}
    opts.out_dir (1,:) char {mustBeNonempty}
    opts.visible (1,1) logical = false
end
idx = resolve_data_file('', opts.data_dir, {fullfile(opts.data_dir, 'network_mode')}, 'ied_netmode_index.mat', ...
    'the network_mode stage of the ied_stim bundle (run_ied_network_mode_stage(ied_stim_config(mode)))');
L = load(idx, 'runs', 'settings'); T = L.runs; S = L.settings;
if ~all(T.ok)
    error('fig_ied_network_mode:FailedRuns', 'The index %s lists failed networks.', idx);
end
folder = fileparts(idx);
Ns = cell(1, height(T));
for k = 1:height(T)
    M = load(fullfile(folder, T.file{k}), 'N'); Ns{k} = M.N;
end
if ~isfolder(opts.out_dir); mkdir(opts.out_dir); end
purple = [0.494 0.247 0.596]; green = [0 0.62 0.45];
fs = 9; fa = 8;
vis = 'off'; if opts.visible; vis = 'on'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 7 4.4], 'Visible', vis);
mk = @(p) axes(fig, 'Units', 'inches', 'Position', p);
N1 = Ns{1};
yl = mode_lim(cellfun(@(N) [real(N.block(1).net); real(N.block(2).net)], Ns, 'UniformOutput', false));

%% A, B
ax = gobjects(1, 2); col = {purple, green};
for b = 1:2
    Bk = N1.block(b);
    ax(b) = mk([0.6 + (b - 1) * 3.3, 2.65, 2.95, 1.45]); hold(ax(b), 'on');
    tt = Bk.t - S.sim.settle;
    ev = Bk.events.time - S.sim.settle;
    ev = ev(ev >= tt(1) & ev <= tt(end));
    for e = ev(:)'
        plot(ax(b), [e e], yl, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.6);
    end
    yline(ax(b), 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6);
    plot(ax(b), tt, real(Bk.net), '-o', 'Color', col{b}, 'MarkerFaceColor', col{b}, 'MarkerSize', 2.5, 'LineWidth', 0.9);
    ylim(ax(b), yl); xlim(ax(b), [tt(1) - 0.5, tt(end) + 0.5]);
    xlabel(ax(b), 'Time in block (s)');
end
ylabel(ax(1), 'Re \lambda_{net} (s^{-1})');
title(ax(1), 'No stimulation', 'Color', purple, 'FontWeight', 'normal');
title(ax(2), sprintf('Stimulation (uniform input %.2g)', S.sim.amplitude), 'Color', green, 'FontWeight', 'normal');

%% C
axC = mk([0.6 0.5 2.95 1.55]); hold(axC, 'on');
lags = N1.event.lags;
EV = cell2mat(cellfun(@(N) real(N.event.ev.net), Ns, 'UniformOutput', false)');
CT = cell2mat(cellfun(@(N) real(N.event.ctrl.net), Ns, 'UniformOutput', false)');
for k = 1:numel(Ns)
    plot(axC, lags, median(real(Ns{k}.event.ev.net), 1), '-', 'Color', [purple 0.35], 'LineWidth', 0.6);
end
yline(axC, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6);
h1 = plot(axC, lags, median(EV, 1), '-o', 'Color', purple, 'MarkerFaceColor', purple, 'MarkerSize', 2.5, 'LineWidth', 1.1);
h2 = plot(axC, lags, median(CT, 1), '--s', 'Color', [0.45 0.45 0.45], 'MarkerSize', 2.5, 'LineWidth', 0.9);
xlim(axC, lags([1 end])); xlabel(axC, 'Time from discharge peak (s)'); ylabel(axC, 'Re \lambda_{net} (s^{-1})');
lg = legend(axC, [h1 h2], {sprintf('%d discharges', size(EV, 1)), sprintf('%d control times', size(CT, 1))}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', fa); lg.ItemTokenSize = [12 8];

%% D
axD = mk([4.75 0.5 1.6 1.55]); hold(axD, 'on');
for k = 1:height(T)
    plot(axD, [1 2], [T.frac_unstable_base(k) T.frac_unstable_stim(k)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7);
end
scatter(axD, ones(height(T), 1), T.frac_unstable_base, 16, purple, 'filled');
scatter(axD, 2 * ones(height(T), 1), T.frac_unstable_stim, 16, green, 'filled');
xlim(axD, [0.6 2.4]); xticks(axD, [1 2]); xticklabels(axD, {'No stim.', 'Stim.'}); ylim(axD, [0 1]);
ylabel(axD, 'Share of time unstable');

all_ax = [ax axC axD];
set(all_ax, 'FontSize', fs, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 0.6, ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1);
lets = {ax(1), 'A', 0.02; ax(2), 'B', 3.4; axC, 'C', 0.02; axD, 'D', 4.1};
for k = 1:4
    p = lets{k, 1}.Position;
    annotation(fig, 'textbox', 'Units', 'inches', 'Position', [lets{k, 3}, p(2) + p(4) + 0.05, 0.4, 0.2], ...
        'String', ['(' lets{k, 2} ')'], 'FontSize', 12, 'LineStyle', 'none', 'VerticalAlignment', 'bottom', ...
        'Margin', 0, 'Tag', 'letter');
end

tag = 'fig_ied_network_mode';
exportgraphics(fig, fullfile(opts.out_dir, [tag '.png']), 'Resolution', 300);
savefig(fig, fullfile(opts.out_dir, [tag '.fig']));
print(fig, fullfile(opts.out_dir, [tag '.svg']), '-dsvg', '-vector');
write_table(fullfile(opts.out_dir, [tag '_table.md']), idx, T, S, size(EV, 1), lags, median(EV, 1), median(CT, 1));
out = struct('fig', fig, 'files', {fullfile(opts.out_dir, strcat(tag, {'.png', '.fig', '.svg', '_table.md'}))}, ...
    'out_dir', opts.out_dir);
if ~opts.visible; close(fig); out.fig = []; end
end

function yl = mode_lim(c)
v = cat(1, c{:}); v = v(isfinite(v));
yl = [floor(min(v)) - 0.5, ceil(max(v)) + 0.5];
end

function write_table(file, idx, T, S, n_ev, lags, med_ev, med_ct)
fid = fopen(file, 'w'); guard = onCleanup(@() fclose(fid));
fprintf(fid, ['# fig_ied_network_mode: the full Jacobian''s network mode in the modular network\n\n' ...
    'Source: `%s`. ied_stim run mode %s; blocks %s of the bundle (baseline, first stimulation block), ' ...
    'uniform input %.2g, fs %g Hz. Full Jacobian (dense eig) every %g s over the first %g s of each analysed ' ...
    'block; around the first %d discharges of the baseline block and as many control times at least %g s ' ...
    'from every discharge, lags %g to %g s by %g s. Network mode: at least %g of the eigenvector''s norm^2 in ' ...
    'the dendritic rows. Stage time %.1f min on %d worker(s).\n\n'], idx, S.run_mode, mat2str(S.blocks), ...
    S.sim.amplitude, S.sim.fs, S.full_dt, S.window_s, S.n_events, S.ctrl_gap, S.lags(1), S.lags(end), ...
    S.lags(2) - S.lags(1), S.x_thr, S.minutes, S.n_workers);
fprintf(fid, ['| Seed | Discharges no stim / stim (whole block) | Share unstable no stim / stim | Median Re (1/s) ' ...
    'no stim / stim | Median coherence (10 pop.) | Discharges with Re > 0 in [-0.3, 0] s | Control times with Re > 0 | ' ...
    'Peak of median (1/s) at lag (s) | Repro. vs bundle block | Seconds (sim / eig) |\n|%s\n'], repmat('---|', 1, 10));
for k = 1:height(T)
    fprintf(fid, '| %d | %d / %d | %.2f / %.2f | %+.2f / %+.2f | %.2f / %.2f | %.2f of %d | %.2f | %+.2f at %+.2f | %.1g | %.0f / %.0f |\n', ...
        T.seed(k), T.events_base(k), T.events_stim(k), T.frac_unstable_base(k), T.frac_unstable_stim(k), ...
        T.median_re_base(k), T.median_re_stim(k), T.median_coh_base(k), T.median_coh_stim(k), ...
        T.events_unstable_pre(k), T.n_events_dense(k), T.ctrl_unstable_pre(k), T.peak_median_re(k), T.peak_lag(k), ...
        T.repro_max_diff(k), T.sim_seconds(k), T.eig_seconds(k));
end
fprintf(fid, '\nPanel C, pooled medians over %d discharges (1/s), by lag:\n\n| Lag (s) |%s\n|---|%s\n| Discharges |%s\n| Control |%s\n', ...
    n_ev, sprintf(' %+.2f |', lags), repmat('---|', 1, numel(lags)), sprintf(' %+.2f |', med_ev), sprintf(' %+.2f |', med_ct));
end

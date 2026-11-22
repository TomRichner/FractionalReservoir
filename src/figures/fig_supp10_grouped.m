function out = fig_supp10_grouped(cfg)
%FIG_SUPP10_GROUPED Supplementary S10, stimulation: adaptation is recruited, and it is needed.
%
%   out = FIG_SUPP10_GROUPED('run_dir', d, 'ied_run_dir', e)
%
% Four panels (figures_main_supp_appen.md, S10):
%   A  the modular (IED) model: change of the mean firing rate, the SFA
%      feedback and the depression 1 - b_1 under stimulation, one point per
%      network from its block means, the median as a line. The numbers of
%      fig_ied_stim_blocks' supp_state panel, computed here from the same
%      summary.mat of an ied_stim run (ied_run_dir) by the same rule.
%   B  the bursting network (burst_mode's first seed) around the step
%      onset, adaptation intact (solid) against frozen at its value before
%      the step (dashed): the mean E rate, the mean SFA term (c/K) sum a of
%      the E neurons and the mean depression factor of the E-to-E route.
%      Replaces block B pilot 2's p2_frozen_bursting figure
%      (scripts/explorations/block_B_pilots/p2_frozen_state), which drove
%      half the E neurons only; here the stimulus is the bursting figures'
%      step to every neuron, the one of Figure 8 A.
%   C  every network of burst_mode: left, bursts before the step and those
%      preceded by a crossing of the full Jacobian's network mode; right,
%      share of time that mode is beyond the boundary before and during the
%      step, gains frozen and in the full Jacobian, one line per network.
%   D  synchrony rho of the E population's dendritic potentials (the
%      bursting figure's chi2 measure, floor-corrected), before and during
%      the step, one line per network.
%
% Reads two run directories and simulates nothing: run_dir holds the
% burst_mode stage (a paper bundle's run), ied_run_dir an ied_stim run
% (its summary.mat). Each errors through resolve_data_file when missing.
% Writes Fig_Supp10_Grouped.{png,svg,fig} and _table.md.
%
% See also: fig_burst_mode, bm_load_runs, fig_ied_stim_blocks

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''      % burst_mode_index.mat, instead of run_dir
    cfg.ied_run_dir (1,:) char    = ''      % an ied_stim run directory (summary.mat)
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''      % unused
    cfg.conditions                = {}      % unused
end
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
B = bm_load_runs(cfg.run_dir, cfg.data_file, 'first');
ied_file = resolve_data_file('', cfg.ied_run_dir, {cfg.ied_run_dir}, 'summary.mat', ...
    'the ied_stim bundle (ied_stim_run), and pass its run directory as ied_run_dir');
L = load(ied_file, 'rows');
st = manuscript_style();
R = B.R{1}; T = B.runs;
fs = 9; fa = 8; H = R.t_step;
ec = SRNNCellTypePairs.type_colors(2);
mode_col = 0.6 * st.condition_color(R.condition);
fg_col = [0.62 0.62 0.62];
purple = [0.494 0.247 0.596]; green = [0 0.62 0.45];

vis = 'on'; if ~cfg.visible; vis = 'off'; end
W = 7.0; Ht = 5.6;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 W Ht], 'Visible', vis);
mk = @(p) axes(fig, 'Units', 'inches', 'Position', p);

%% A: state change in the modular model
S = L.rows; seeds = unique(S.seed, 'stable')'; nS = numel(seeds);
state = zeros(nS, 3, 2);
for i = 1:nS
    for c = 1:2
        Tc = S(S.seed == seeds(i) & S.stim == (c == 2), :);
        state(i, :, c) = [mean(Tc.mean_rate) mean(Tc.sfa_feedback) mean(Tc.depression)];
    end
end
pct = 100 * (state(:, :, 2) ./ state(:, :, 1) - 1);
axA = mk([0.62 3.25 2.3 1.95]); hold(axA, 'on');
rng(7, 'twister');
for c = 1:3
    jit = 0.18 * (rand(nS, 1) - 0.5);
    scatter(axA, c + jit, pct(:, c), 14, green, 'filled', 'MarkerFaceAlpha', 0.8);
    plot(axA, [c - 0.3, c + 0.3], median(pct(:, c)) * [1 1], '-', 'Color', 0.55 * green, 'LineWidth', 1.6);
end
yline(axA, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6);
xlim(axA, [0.4 3.6]); xticks(axA, 1:3); xticklabels(axA, {'Rate', 'SFA', 'Depression'});
ylabel(axA, 'Change with stimulation (%)');

%% B: onset, adaptation intact against frozen
win = [H - 5, H + 15];
I = R.intact; Z = R.frozen;
ki = I.t >= win(1) & I.t <= win(2); kz = Z.t >= H & Z.t <= win(2);
hB = 0.53; gB = 0.11; xB = 3.95; wB = 2.9;
yB = 3.25 + [2 * (hB + gB), hB + gB, 0];
iEE = find(strcmp(I.dep_names, 'E to E'), 1);
series = {I.rate(:, 1), Z.rate(:, 1), 'Rate E'; I.adapt(:, 1), Z.adapt(:, 1), 'SFA E'; ...
    I.dep(:, iEE), Z.dep(:, iEE), 'Depr. E\rightarrowE'};
axB = gobjects(1, 3);
for r = 1:3
    axB(r) = mk([xB, yB(r), wB, hB]); hold(axB(r), 'on');
    xline(axB(r), H, ':', 'Color', [0 0 0], 'LineWidth', 0.6);
    hi = plot(axB(r), I.t(ki), series{r, 1}(ki), '-', 'Color', ec(1, :), 'LineWidth', 0.9);
    hz = plot(axB(r), Z.t(kz), series{r, 2}(kz), '--', 'Color', ec(1, :), 'LineWidth', 0.9);
    xlim(axB(r), win); ylabel(axB(r), series{r, 3});
    if r < 3; axB(r).XTickLabel = {}; end
    if r == 3
        % The depression factor is at most 1. A floor a tenth below the
        % deepest point of the window: until 2026-11-22 the automatic limit
        % sat at that point (about 0.47) and the dip after the step was
        % drawn through the time axis.
        lo = min([series{r, 1}(ki); series{r, 2}(kz)]);
        ylim(axB(r), [0.1 * floor(10 * lo - 0.5), 1]);
    end
    if r == 1
        lg = legend(axB(r), [hi hz], {'Intact', 'Adaptation frozen'}, 'Location', 'northeast', ...
            'Box', 'off', 'FontSize', fa, 'AutoUpdate', 'off'); lg.ItemTokenSize = [12 8];
        ylim(axB(r), [0 1]);
    end
end
xlabel(axB(3), 'Time (s)');

%% C: per network
nN = height(T);
axC1 = mk([0.62 0.55 1.75 1.95]); hold(axC1, 'on');
h1 = scatter(axC1, 1:nN, T.n_bursts_before, 30, [0.3 0.3 0.3], 'LineWidth', 0.8);
h2 = scatter(axC1, 1:nN, T.n_bursts_with_cross_full, 12, mode_col, 'filled');
xlim(axC1, [0.4 nN + 0.6]); xticks(axC1, 1:nN);
xticklabels(axC1, arrayfun(@(s) sprintf('%d', s), T.seed1, 'UniformOutput', false));
ylim(axC1, [0, max(T.n_bursts_before) + 4]);
axC1.XTickLabelRotation = 0;
xlabel(axC1, 'Network (seed)'); ylabel(axC1, 'Bursts before the step');
lg = legend(axC1, [h1 h2], {'Bursts', 'Preceded by a crossing'}, 'Location', 'northwest', 'Box', 'off', ...
    'FontSize', fa); lg.ItemTokenSize = [8 8];
axC2 = mk([2.95 0.55 1.45 1.95]); hold(axC2, 'on');
pairs = {T.frac_beyond_before_frozen_gain, T.frac_beyond_during_frozen_gain, fg_col; ...
    T.frac_beyond_before_full, T.frac_beyond_during_full, mode_col};
hc = gobjects(1, 2);
for p = 1:2
    off = (p - 1.5) * 0.12;
    for i = 1:nN
        hc(p) = plot(axC2, [1 2] + off, [pairs{p, 1}(i) pairs{p, 2}(i)], '-o', 'Color', pairs{p, 3}, ...
            'MarkerFaceColor', pairs{p, 3}, 'MarkerSize', 3, 'LineWidth', 0.8);
    end
end
xlim(axC2, [0.6 2.4]); ylim(axC2, [0 1]); xticks(axC2, [1 2]); xticklabels(axC2, {'Before', 'Step'});
ylabel(axC2, 'Share of time unstable');
lg = legend(axC2, hc, {'Gains frozen', 'Full Jacobian'}, 'Location', 'northoutside', 'Box', 'off', 'FontSize', fa);
lg.ItemTokenSize = [12 8];

%% D: synchrony
axD = mk([5.25 0.55 1.6 1.95]); hold(axD, 'on');
for i = 1:nN
    plot(axD, [1 2], [T.rho_before(i) T.rho_during(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7);
end
scatter(axD, ones(nN, 1), T.rho_before, 16, purple, 'filled');
scatter(axD, 2 * ones(nN, 1), T.rho_during, 16, green, 'filled');
xlim(axD, [0.6 2.4]); xticks(axD, [1 2]); xticklabels(axD, {'Before', 'Step'});
ylim(axD, [0 max([T.rho_before; T.rho_during]) * 1.1]);
ylabel(axD, 'Synchrony \rho');

%% Style and letters
all_ax = [axA axB axC1 axC2 axD];
set(all_ax, 'FontSize', fs, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 0.6, ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1);
letter(fig, axA, 'A', 0.02); letter(fig, axB(1), 'B', xB - 0.55);
letter(fig, axC1, 'C', 0.02); letter(fig, axD, 'D', 5.25 - 0.55);

fig_tag = 'Fig_Supp10_Grouped';
out = struct('figs', fig, 'files', {{}}, 'source', {{B.index_file, ied_file}});
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    tfile = fullfile(out_dir, [fig_tag '_table.md']);
    bm_write_table(tfile, B, [], 'Supplementary S10, panels B to D');
    fid = fopen(tfile, 'a'); guard = onCleanup(@() fclose(fid));
    fprintf(fid, ['\n## Panel A: modular model, change with stimulation (%%), per network from its block means\n\n' ...
        'Source: `%s`.\n\n| Seed | Rate | SFA | Depression 1 - b_1 |\n|---|---|---|---|\n'], ied_file);
    for i = 1:nS
        fprintf(fid, '| %d | %+.1f | %+.1f | %+.1f |\n', seeds(i), pct(i, :));
    end
    fprintf(fid, '| median | %+.1f | %+.1f | %+.1f |\n', median(pct, 1));
    fprintf(fid, ['\nPanel B: network %s, window %s s; intact and frozen from t = %.3f s (the last sample ' ...
        'before the step). Drawn: mean E rate, mean (c/K) sum a of E, mean E-to-E depression factor.\n'], ...
        mat2str(R.seed), mat2str(win), Z.t_freeze);
    out.files{end + 1} = tfile;
end
end

function letter(fig, ax, ch, x_in)
p = ax.Position;
annotation(fig, 'textbox', 'Units', 'inches', 'Position', [x_in, p(2) + p(4) + 0.05, 0.4, 0.2], ...
    'String', ['(' ch ')'], 'FontSize', 12, 'FontName', 'Helvetica', 'LineStyle', 'none', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Margin', 0, 'Tag', 'letter');
end

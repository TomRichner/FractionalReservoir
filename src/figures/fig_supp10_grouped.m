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
H = R.t_step;
ec = [st.E_color; st.I_color];                    % the manuscript's one E and one I color
mode_col = 0.6 * st.condition_color(R.condition);
fg_col = [0.62 0.62 0.62];
purple = [0.494 0.247 0.596]; green = [0 0.62 0.45];
lw = st.line_lw; thin = 0.6 * st.line_lw;         % traces; connecting and reference lines

% THE COMMON CANVAS (Plan 2, canvas standard, 2026-11-22). Until then S10
% was drawn 7.0 x 5.6 in with 9-pt type (letters 12 pt), so it printed its
% letters twice as large relative to the figure as a 13-in composite does.
% It is laid out again at canvas_w with the type by role (tick_fs,
% label_fs, letter_fs; canvas_type at the end) and the lines and marks to
% go with it, not enlarged: two rows, A beside B, then C (two parts) and D.
% Geometry in inches.
vis = 'on'; if ~cfg.visible; vis = 'off'; end
W = st.canvas_w;
y2 = 1.0; h2 = 2.8;                               % row 2 (C, D): axes bottom and height
hC2 = 2.2;                                        % C right: room above it for its key
y1 = y2 + h2 + 1.1; h1 = 3.2;                     % row 1 (A, B): above row 2's letters and B's time axis
Ht = y1 + h1 + 0.55;
xL = 1.05; wL = 3.6;                              % left column (A, C left)
xB = 6.1; wB = W - 0.2 - xB;                      % B, C right
xD = 9.85; wD = W - 0.2 - xD;
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
axA = mk([xL y1 wL h1]); hold(axA, 'on');
rng(7, 'twister');
for c = 1:3
    jit = 0.18 * (rand(nS, 1) - 0.5);
    scatter(axA, c + jit, pct(:, c), 50, green, 'filled', 'MarkerFaceAlpha', 0.8);
    plot(axA, [c - 0.3, c + 0.3], median(pct(:, c)) * [1 1], '-', 'Color', 0.55 * green, 'LineWidth', 1.5 * lw);
end
yline(axA, 0, ':', 'Color', [0 0 0], 'LineWidth', st.axis_lw);
xlim(axA, [0.4 3.6]); xticks(axA, 1:3); xticklabels(axA, {'Rate', 'SFA', 'Depression'});
ylabel(axA, 'Change with stimulation (%)');

%% B: onset, adaptation intact against frozen
win = [H - 5, H + 15];
I = R.intact; Z = R.frozen;
ki = I.t >= win(1) & I.t <= win(2); kz = Z.t >= H & Z.t <= win(2);
gB = 0.22; hB = (h1 - 2 * gB) / 3;          % the gap keeps a row's top tick label off the bottom one above
yB = y1 + [2 * (hB + gB), hB + gB, 0];
iEE = find(strcmp(I.dep_names, 'E to E'), 1);
series = {I.rate(:, 1), Z.rate(:, 1), 'Rate E'; I.adapt(:, 1), Z.adapt(:, 1), 'SFA E'; ...
    I.dep(:, iEE), Z.dep(:, iEE), {'Depression', 'E\rightarrowE'}};
axB = gobjects(1, 3);
lgB = [];
for r = 1:3
    axB(r) = mk([xB, yB(r), wB, hB]); hold(axB(r), 'on');
    xline(axB(r), H, ':', 'Color', [0 0 0], 'LineWidth', st.axis_lw);
    hi = plot(axB(r), I.t(ki), series{r, 1}(ki), '-', 'Color', ec(1, :), 'LineWidth', lw);
    hz = plot(axB(r), Z.t(kz), series{r, 2}(kz), '--', 'Color', ec(1, :), 'LineWidth', lw);
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
        % The key above the panel, in one row: inside, at 14 pt, it covered
        % the frozen arm (canvas standard).
        lgB = legend(axB(r), [hi hz], {'Intact', 'Adaptation frozen'}, 'Orientation', 'horizontal', ...
            'Box', 'off', 'AutoUpdate', 'off');
        lgB.ItemTokenSize = [24 10];
        ylim(axB(r), [0 1]);
    end
end
xlabel(axB(3), 'Time (s)');

%% C: per network
nN = height(T);
axC1 = mk([xL y2 wL h2]); hold(axC1, 'on');
h1s = scatter(axC1, 1:nN, T.n_bursts_before, 110, [0.3 0.3 0.3], 'LineWidth', 1.2);
h2s = scatter(axC1, 1:nN, T.n_bursts_with_cross_full, 45, mode_col, 'filled');
xlim(axC1, [0.4 nN + 0.6]); xticks(axC1, 1:nN);
xticklabels(axC1, arrayfun(@(s) sprintf('%d', s), T.seed1, 'UniformOutput', false));
ylim(axC1, [0, max(T.n_bursts_before) + 8]);
axC1.XTickLabelRotation = 0;
xlabel(axC1, 'Network (seed)'); ylabel(axC1, 'Bursts before the step');
lg = legend(axC1, [h1s h2s], {'Bursts', 'Preceded by a crossing'}, 'Location', 'northwest', 'Box', 'off');
lg.ItemTokenSize = [12 12];
axC2 = mk([xB y2 2.4 hC2]); hold(axC2, 'on');
pairs = {T.frac_beyond_before_frozen_gain, T.frac_beyond_during_frozen_gain, fg_col; ...
    T.frac_beyond_before_full, T.frac_beyond_during_full, mode_col};
hc = gobjects(1, 2);
for p = 1:2
    off = (p - 1.5) * 0.12;
    for i = 1:nN
        hc(p) = plot(axC2, [1 2] + off, [pairs{p, 1}(i) pairs{p, 2}(i)], '-o', 'Color', pairs{p, 3}, ...
            'MarkerFaceColor', pairs{p, 3}, 'MarkerSize', 6, 'LineWidth', thin);
    end
end
xlim(axC2, [0.6 2.4]); ylim(axC2, [0 1]); xticks(axC2, [1 2]); xticklabels(axC2, {'Before', 'Step'});
ylabel(axC2, 'Share of time unstable');
lgC = legend(axC2, hc, {'Gains frozen', 'Full Jacobian'}, 'Box', 'off');
lgC.ItemTokenSize = [24 10];

%% D: synchrony
axD = mk([xD y2 wD h2]); hold(axD, 'on');
for i = 1:nN
    plot(axD, [1 2], [T.rho_before(i) T.rho_during(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', thin);
end
scatter(axD, ones(nN, 1), T.rho_before, 50, purple, 'filled');
scatter(axD, 2 * ones(nN, 1), T.rho_during, 50, green, 'filled');
xlim(axD, [0.6 2.4]); xticks(axD, [1 2]); xticklabels(axD, {'Before', 'Step'});
ylim(axD, [0 max([T.rho_before; T.rho_during]) * 1.1]);
ylabel(axD, 'Synchrony \rho');

%% Style, the keys above B and C, and the letters
all_ax = [axA axB axC1 axC2 axD];
set(all_ax, 'TickDir', 'out', 'Box', 'off', 'LineWidth', st.axis_lw, 'TitleFontSizeMultiplier', 1, ...
    'LabelFontSizeMultiplier', 1);
canvas_type(fig);
drawnow;
% B's three y labels in one column, left of the widest tick labels.
lab = [axB.YLabel];
set(lab, 'Units', 'normalized');
xs = arrayfun(@(t) t.Position(1), lab);
for k = 1:3; lab(k).Position(1) = min(xs); end
set([lgB lgC], 'Units', 'inches', 'Location', 'none');
lgB.Position(1:2) = [xB + wB - lgB.Position(3), y1 + h1 + 0.02];   % above B, at its right end
lgC.Position(1:2) = [xB, y2 + hC2 + 0.05];                         % above C right, at its left edge
letter(fig, 'A', 0.05, y1 + h1 + 0.05); letter(fig, 'B', xB - 0.95, y1 + h1 + 0.05);
letter(fig, 'C', 0.05, y2 + h2 + 0.05); letter(fig, 'D', xD - 0.95, y2 + h2 + 0.05);

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

function letter(fig, ch, x_in, y_in)
% The panel letter, its bottom left at (x_in, y_in) inches, at letter_fs.
st = manuscript_style();
annotation(fig, 'textbox', 'Units', 'inches', 'Position', [x_in, y_in, 0.5, 0.3], ...
    'String', ['(' ch ')'], 'FontSize', st.letter_fs, 'FontName', 'Helvetica', 'LineStyle', 'none', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Margin', 0, 'Tag', 'letter');
end

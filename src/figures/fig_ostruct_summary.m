function out = fig_ostruct_summary(cfg)
% FIG_OSTRUCT_SUMMARY Every network of the outlier_structures stage: detachment of the intended mode, responsive fractions, the mode's rate frozen and full.
%
%   out = FIG_OSTRUCT_SUMMARY('run_dir', d)
%
% Draws what the opt-in stage outlier_structures saved (FD-043). Reads one
% run directory or the file given (the stage's ostruct_index.mat);
% simulates nothing. Supplementary.
%
% WHAT IS DRAWN. One column per network, in three groups separated by
% dotted lines (added low-rank structure by strength and rotation; local
% outliers and the reference; block-mean outliers), and within a network
% the conditions side by side in their colors. Every marker is the MEDIAN
% OVER SEEDS of a run's value, and the bar through it the RANGE over seeds
% (no bar with one seed).
%   top     the detachment of the intended mode: the real part of its
%           eigenvalue in the frozen-gain connectivity over the theoretical
%           bulk radius (median over the window), on a log axis; the gray
%           line is the radius; the open black diamonds the same with every
%           gain made equal (the uniform-gain value: what the weights alone
%           give). Values below the axis's floor are drawn at the floor
%           (the sidecar counts them).
%   middle  the responsive fraction of the E neurons (filled) and of the
%           I neurons (open), side by side: detachment at or below 1 with a low
%           responsive fraction is a pinned network, not a controlled one.
%   bottom  the real part of the mode as a rate, with the slow variables
%           frozen (open) and free (filled: its continuation in the full
%           Jacobian), joined by a thin line; dashed: zero.
% No text on the panels but the axis labels and the network names under
% the bottom row. Fig_Ostruct_Summary_table.md: the caption paragraph and
% the numbers (median and range over seeds per network and condition, the
% structural check, runs not drawn).
%
% Options: run_dir, data_file, out_dir, save, visible, verbose,
% preset_name (unused), conditions (select_conditions).
%
% See also: ostruct_figure_data, ostruct_mode_measures, fig_ostruct_plane

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
end

RHO_FLOOR = 0.2; RHO_CEIL = 8;
SUB = 0.06;                              % offset of the two quantities within a condition's slot
setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = ostruct_figure_data(cfg, 'fig_ostruct_summary');
st = F.st;
n_net = numel(F.keys); n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];
off = linspace(-0.25, 0.25, n_c);
if n_c == 1; off = 0; end

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 max(9, 1.15 * n_net + 2) 10], 'Visible', vis);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
n_floor = 0;
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    for w = 1:n_net
        for c = 1:n_c
            if F.n_seed(w, c) == 0; continue; end
            x = w + off(c); col = F.colors{c};
            xa = x - SUB; xb = x + SUB;                       % the two quantities of a slot side by side
            switch p
                case 1
                    clampv = @(v) min(max(v, RHO_FLOOR), RHO_CEIL);
                    md = F.med.mode_rho(w, c);
                    n_floor = n_floor + (md < RHO_FLOOR);
                    range_bar(ax, xa, clampv(F.lo.mode_rho(w, c)), clampv(F.hi.mode_rho(w, c)), col);
                    plot(ax, xa, clampv(md), 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                    range_bar(ax, xb, clampv(F.lo.mode_rho_uniform(w, c)), clampv(F.hi.mode_rho_uniform(w, c)), [0 0 0]);
                    plot(ax, xb, clampv(F.med.mode_rho_uniform(w, c)), 'd', 'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
                case 2
                    range_bar(ax, xa, F.lo.frac_on_E(w, c), F.hi.frac_on_E(w, c), col);
                    range_bar(ax, xb, F.lo.frac_on_I(w, c), F.hi.frac_on_I(w, c), col);
                    plot(ax, xa, F.med.frac_on_E(w, c), 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                    plot(ax, xb, F.med.frac_on_I(w, c), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col, 'LineWidth', 1.5);
                otherwise
                    plot(ax, [xa xb], [F.med.rate_frozen(w, c), F.med.rate_full(w, c)], '-', 'Color', col, 'LineWidth', 0.75);
                    range_bar(ax, xa, F.lo.rate_frozen(w, c), F.hi.rate_frozen(w, c), col);
                    range_bar(ax, xb, F.lo.rate_full(w, c), F.hi.rate_full(w, c), col);
                    plot(ax, xa, F.med.rate_frozen(w, c), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col, 'LineWidth', 1.5);
                    plot(ax, xb, F.med.rate_full(w, c), 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
            end
        end
    end
    seps = find(diff(F.group) ~= 0) + 0.5;
    if ~isempty(seps); xline(ax, seps, ':', 'Color', gray, 'LineWidth', 1.0); end
    if p == 1
        yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
        set(ax, 'YScale', 'log', 'YLim', [RHO_FLOOR * 0.95, RHO_CEIL], 'YTick', [0.25 0.5 1 2 4 8]);
        ylabel(ax, 'Re \lambda_{mode} / bulk radius', 'FontSize', st.label_fs);
    elseif p == 2
        ylim(ax, [0 1.05]);
        ylabel(ax, 'responsive fraction', 'FontSize', st.label_fs);
    else
        yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
        ylabel(ax, 'rate of the mode (s^{-1})', 'FontSize', st.label_fs);
    end
    set(ax, 'XTick', 1:n_net, 'XTickLabel', {}, 'XLim', [0.4, n_net + 0.6], 'FontSize', st.tick_fs, ...
        'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    if p == 3; set(ax, 'XTickLabel', F.short, 'XTickLabelRotation', 35, 'TickLabelInterpreter', 'tex'); end
end

fig_tag = 'Fig_Ostruct_Summary';
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    caption = sprintf(['Is adaptation''s reining in of outliers limited to outliers made by a mismatch of the E and I ' ...
        'block means? Every network of the outlier_structures stage, in three groups separated by dotted lines: ' ...
        'low-rank structure added to the weights on random patterns of neurons (a real eigenvalue and conjugate pairs ' ...
        'of increasing rotation at about 1.4 bulk radii, then two at about 2.5), local outliers of the balanced sparse ' ...
        'Dale matrix with the reference network on new weights for comparison, and outliers of the block means (single ' ...
        'routes changed; the E-I loop strengthened with E-to-E). Within a network the conditions in their colors. Each ' ...
        'marker is the median over %s of the run''s value, the bar the range over seeds. TOP: detachment of the ' ...
        'INTENDED mode (the eigenvalue of the frozen-gain effective connectivity that carries the structural mode''s ' ...
        'eigenvector), its real part over the theoretical bulk radius of each state, median over the window; log axis, ' ...
        'values below %.2g drawn at %.2g (%d network-condition medians). Gray line: the radius. Open diamonds: the same ' ...
        'with every gain set to the rms gain, what the weights alone give. MIDDLE: the fraction of E (filled) and I ' ...
        '(open) neurons off the floor and the ceiling of the nonlinearity. BOTTOM: the real part of the mode as a ' ...
        'rate, slow variables frozen (open) and free (filled, the continuation in the full Jacobian, median over the ' ...
        'calibration states); dashed: zero.'], F.seeds_text, RHO_FLOOR, RHO_FLOOR, n_floor);
    cols = {'lam0_re_R', 'lam0_im_R', 'lambda_1', 'mode_rho', 'mode_rho_p05', 'mode_rho_p95', 'mode_rho_uniform', ...
        'mode_reined', 'mode_ov', 'mode_frac_rightmost', 'rate_struct', 'rate_frozen', 'rate_full', 'frac_on_E', ...
        'frac_on_I', 'on_top', 'on_rest', 'right_rho', 'right_rho_uniform', 'full_abscissa', 'full_frac_unstable'};
    name = ostruct_write_table(out_dir, fig_tag, caption, F, cols, 1:n_net);
    if ~isempty(name); out.files{end + 1} = name; end
end
vprintf(cfg.verbose, 'verbose', '%s: %d networks, %d conditions\n', fig_tag, n_net, n_c);
end

%% ------------------------------------------------------------------------
function range_bar(ax, x, lo, hi, col)
% The range over seeds; nothing when it is a point.
if isfinite(lo) && isfinite(hi) && hi > lo
    plot(ax, [x x], [lo hi], '-', 'Color', col, 'LineWidth', 2.0);
end
end


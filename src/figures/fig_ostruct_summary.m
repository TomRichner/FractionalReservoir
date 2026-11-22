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
%           bulk radius (median over the window), on a LINEAR axis that
%           holds every value; the shaded band from -1 to 1 is the bulk
%           (the eigenvalues are centered on 0 in these units), the gray
%           line at 1 its right edge; the open black diamonds the same with every
%           gain made equal (the uniform-gain value: what the weights alone
%           give). A mode left of the band (below -1) has been pushed out
%           of the bulk on the far side: its carriers are pinned. (Until
%           2026-11-21 the axis was logarithmic and every value below 0.2,
%           the negative ones included, was drawn at the floor.)
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
% preset_name (unused), conditions (select_conditions), and
%   networks   keys of the stage's networks to draw, left to right in the
%              order given ({} : every network, in the stage's order, as
%              before the option existed). fig_supp7_grouped passes the
%              structures of the gallery S3, so that its panel B summarizes
%              what S3 shows.
%   labels     the names under the bottom row, one per drawn network ({} :
%              the stage's short names)
%   fig_size   [width height] inches ([] : 1.15 in per network plus 2, at
%              least 9, by 10)
%   compact    false. true: smaller markers set closer, for a panel with
%              less than about an inch per network.
%
% See also: ostruct_figure_data, ostruct_mode_measures, fig_ostruct_plane

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
    cfg.networks    (1,:) cell    = {}
    cfg.labels      (1,:) cell    = {}
    cfg.fig_size    double        = []
    cfg.compact     (1,1) logical = false
end

SUB = 0.06;                              % offset of the two quantities within a condition's slot
SPREAD = 0.25;                           % the conditions of a network span +-SPREAD about its tick
MS = [8 7];                              % marker sizes: filled, open
if cfg.compact; SUB = 0.07; SPREAD = 0.27; MS = [6 5.5]; end
setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = ostruct_figure_data(rmfield(cfg, {'networks', 'labels', 'fig_size', 'compact'}), 'fig_ostruct_summary');
st = F.st;
n_c = numel(F.cond_names);
% The networks drawn: indices into F.keys, in drawing order.
sel = 1:numel(F.keys);
if ~isempty(cfg.networks)
    [found, sel] = ismember(cfg.networks, F.keys);
    if ~all(found)
        error('fig_ostruct_summary:NoNetwork', 'The outlier_structures stage in %s has no network %s; it has: %s.', ...
            F.file, strjoin(cfg.networks(~found), ', '), strjoin(F.keys(:)', ', '));
    end
end
n_net = numel(sel);
names = F.short(sel);
if ~isempty(cfg.labels)
    if numel(cfg.labels) ~= n_net
        error('fig_ostruct_summary:Labels', 'labels has %d entries for %d networks drawn.', numel(cfg.labels), n_net);
    end
    names = cfg.labels;
end
gray = [0.55 0.55 0.55];
off = linspace(-SPREAD, SPREAD, n_c);
if n_c == 1; off = 0; end

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig_size = [max(9, 1.15 * n_net + 2), 10];
if ~isempty(cfg.fig_size); fig_size = cfg.fig_size(:)'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 fig_size], 'Visible', vis);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
rho_all = [F.lo.mode_rho(sel, :); F.hi.mode_rho(sel, :); F.lo.mode_rho_uniform(sel, :); F.hi.mode_rho_uniform(sel, :)];
rho_all = rho_all(isfinite(rho_all));
rho_lim = [min([rho_all; -1]), max([rho_all; 1])];
rho_lim = rho_lim + [-1 1] * 0.06 * diff(rho_lim);
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    if p == 1                                 % under the markers: the bulk, -1 to 1
        patch(ax, [0.4, n_net + 0.6, n_net + 0.6, 0.4], [-1 -1 1 1], [0.93 0.93 0.93], 'EdgeColor', 'none');
        yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
    end
    for w = 1:n_net
        wi = sel(w);
        for c = 1:n_c
            if F.n_seed(wi, c) == 0; continue; end
            x = w + off(c); col = F.colors{c};
            xa = x - SUB; xb = x + SUB;                       % the two quantities of a slot side by side
            switch p
                case 1
                    md = F.med.mode_rho(wi, c);
                    range_bar(ax, xa, F.lo.mode_rho(wi, c), F.hi.mode_rho(wi, c), col);
                    plot(ax, xa, md, 'o', 'MarkerSize', MS(1), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                    range_bar(ax, xb, F.lo.mode_rho_uniform(wi, c), F.hi.mode_rho_uniform(wi, c), [0 0 0]);
                    plot(ax, xb, F.med.mode_rho_uniform(wi, c), 'd', 'MarkerSize', MS(2), 'MarkerFaceColor', 'w', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
                case 2
                    range_bar(ax, xa, F.lo.frac_on_E(wi, c), F.hi.frac_on_E(wi, c), col);
                    range_bar(ax, xb, F.lo.frac_on_I(wi, c), F.hi.frac_on_I(wi, c), col);
                    plot(ax, xa, F.med.frac_on_E(wi, c), 'o', 'MarkerSize', MS(1), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                    plot(ax, xb, F.med.frac_on_I(wi, c), 'o', 'MarkerSize', MS(2), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col, 'LineWidth', 1.5);
                otherwise
                    plot(ax, [xa xb], [F.med.rate_frozen(wi, c), F.med.rate_full(wi, c)], '-', 'Color', col, 'LineWidth', 0.75);
                    range_bar(ax, xa, F.lo.rate_frozen(wi, c), F.hi.rate_frozen(wi, c), col);
                    range_bar(ax, xb, F.lo.rate_full(wi, c), F.hi.rate_full(wi, c), col);
                    plot(ax, xa, F.med.rate_frozen(wi, c), 'o', 'MarkerSize', MS(2), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col, 'LineWidth', 1.5);
                    plot(ax, xb, F.med.rate_full(wi, c), 'o', 'MarkerSize', MS(1), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
            end
        end
    end
    seps = find(diff(F.group(sel)) ~= 0) + 0.5;
    if ~isempty(seps); xline(ax, seps, ':', 'Color', gray, 'LineWidth', 1.0); end
    if p == 1
        set(ax, 'YLim', rho_lim, 'Layer', 'top');
        ylabel(ax, {'Re \lambda_{mode}', '/ bulk radius'}, 'FontSize', st.label_fs);
    elseif p == 2
        ylim(ax, [0 1.05]);
        ylabel(ax, 'Responsive fraction', 'FontSize', st.label_fs);
    else
        yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
        ylabel(ax, {'Rate of the mode', '(s^{-1})'}, 'FontSize', st.label_fs);
    end
    set(ax, 'XTick', 1:n_net, 'XTickLabel', {}, 'XLim', [0.4, n_net + 0.6], 'FontSize', st.tick_fs, ...
        'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    if p == 3; set(ax, 'XTickLabel', names, 'XTickLabelRotation', 35, 'TickLabelInterpreter', 'tex'); end
end

fig_tag = 'Fig_Ostruct_Summary';
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    caption = sprintf(['Is adaptation''s reining in of outliers limited to outliers made by a mismatch of the E and I ' ...
        'block means? Every network of the outlier_structures stage, in three groups separated by dotted lines: ' ...
        'low-rank structure added to the weights on random patterns of neurons (a real eigenvalue and conjugate pairs ' ...
        'of increasing rotation at about 1.4 bulk radii, then two at about 2.5), local outliers of the balanced sparse ' ...
        'Dale matrix with the reference network on new weights for comparison, and outliers of the block means (single ' ...
        'routes changed; the E-I loop strengthened with E-to-E). Within a network the conditions in their colors. Each ' ...
        'marker is the median over %s of the run''s value, the bar the range over seeds. TOP: detachment of the ' ...
        'INTENDED mode (the eigenvalue of the frozen-gain effective connectivity that carries the structural mode''s ' ...
        'eigenvector), its real part over the theoretical bulk radius of each state, median over the window; LINEAR axis, ' ...
        'every value at its place (drawn range %.3g to %.3g; %d network-condition medians are below -1). Shaded band: ' ...
        'the bulk, from -1 to 1 in these units (the eigenvalues are centered on 0); gray line: its right edge. A mode ' ...
        'below -1 has left the bulk on the far side, which happens when the neurons that carry it are pinned at the ' ...
        'floor or ceiling. Open diamonds: the same ' ...
        'with every gain set to the rms gain, what the weights alone give. MIDDLE: the fraction of E (filled) and I ' ...
        '(open) neurons off the floor and the ceiling of the nonlinearity. BOTTOM: the real part of the mode as a ' ...
        'rate, slow variables frozen (open) and free (filled, the continuation in the full Jacobian, median over the ' ...
        'calibration states); dashed: zero.'], F.seeds_text, rho_lim, nnz(F.med.mode_rho(sel, :) < -1));
    cols = {'lam0_re_R', 'lam0_im_R', 'lambda_1', 'mode_rho', 'mode_rho_p05', 'mode_rho_p95', 'mode_rho_uniform', ...
        'mode_reined', 'mode_ov', 'mode_frac_rightmost', 'rate_struct', 'rate_frozen', 'rate_full', 'frac_on_E', ...
        'frac_on_I', 'on_top', 'on_rest', 'right_rho', 'right_rho_uniform', 'full_abscissa', 'full_frac_unstable'};
    notes = {};
    if ~isempty(cfg.networks)
        left_out = setdiff(F.keys(:)', F.keys(sel), 'stable');
        notes = {sprintf(['DRAWN HERE: %d of the stage''s %d networks, left to right: %s (under the names %s). Where ' ...
            'the caption says "every network" read these; the dotted lines fall where the group changes in this ' ...
            'order. Not drawn and not tabulated here (the stand-alone fig_ostruct_summary has them): %s.'], ...
            n_net, numel(F.keys), strjoin(F.keys(sel), ', '), strjoin(strrep(names, '|', '/'), ', '), ...
            strjoin([left_out, repmat({'none'}, 1, isempty(left_out))], ', '))};
    end
    name = ostruct_write_table(out_dir, fig_tag, caption, F, cols, sel, notes);
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


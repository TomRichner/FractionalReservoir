function out = fig_eigdens_measures(cfg)
% FIG_EIGDENS_MEASURES The four eigenvalue-density counts by condition: a box over the networks, the reference values beside it.
%
%   out = FIG_EIGDENS_MEASURES('run_dir', d)
%   out = FIG_EIGDENS_MEASURES('run_dir', d, 'network', 'EEx1p5')
%
% Draws what the opt-in stage eigenvalue_density saved (FD-042). Reads one
% run directory or the file given; simulates nothing. Supplementary.
% Definitions of the measures: eigdens_measures. Each network (seed) gives
% one value, the mean over its calibration states; the conditions are on
% the horizontal axis by their names (grouped_condition_label 'tick').
%
% BOX AND WHISKER OVER THE NETWORKS (Tom's notes of 2026-11-22, S6 panel
% C): per condition, the measured values of its networks as a box in the
% condition's color (the box the 25th to 75th percentile, MATLAB's
% prctile; the thick bar the median; whiskers to the most extreme network
% within 1.5 interquartile ranges of the box; a network beyond them a
% dot), and right of the box, apart from it, the REFERENCE values as open
% gray diamonds, one per network, in one column. With few networks the
% box says little: with 2 (the `fast` run) prctile puts the box on the two
% values themselves and the bar on their mean, and the whiskers have no
% length; with 1 the box is a bar. Until 2026-11-22 the networks were
% dots spread within a condition with the diamond on top of each.
%
%   1  rim occupancy: the fraction of the eigenvalues of M = W D beyond
%      0.8 R_eff, over the 0.36 of a uniform disk. Box: measured; diamond:
%      the theory with the saved gains; dotted: uniform.
%   2  the number of eigenvalues of the frozen-gain dendritic block with
%      rate in (-2, 0] 1/s ((lambda - 1)/tau_d, tau_d = 0.1 s). Box:
%      measured; diamond: a uniform disk of the same radius.
%   3  the number with rate above 0. Box: frozen gains; diamond: the full
%      Jacobian at the same states (slow variables free).
%   4  the number of eigenvalues of the FULL Jacobian with rate in
%      (-2, 0] 1/s whose eigenvector lies mostly (at least half its norm^2)
%      in the dendritic variables: the fast modes, without the slow
%      variables' own eigenvalues. No reference: the box alone.
% Panels 3 (diamonds) and 4 are empty for a run without the dense pass.
%
% No text on the panels; Fig_Eigdens_Measures<suffix>_table.md has every
% run (seed) of the network with every measure, and the definitions.
% Options as fig_eigdens_radial, and fig_size (inches [w h]; NaN: 34 x 11
% cm).
%
% See also: eigdens_measures, eigdens_run_measures, fig_eigdens_radial

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.network     (1,:) char    = 'reference'
    cfg.conditions                = {}
    cfg.fig_size    (1,2) double  = [NaN NaN]   % inches [w h]; NaN = 34 x 11 cm (fig_supp6_grouped lays it out at the canvas width)
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = eigdens_figure_data(cfg, 'fig_eigdens_measures');
st = F.st; S = F.D.settings.measures;
n_c = numel(F.cond_names);
uniform_shell = 1 - S.shell_r^2;
gray = [0.45 0.45 0.45];
ax_style = {'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out'};
mrow = @(r, M, name) mean(r.(M)(strcmp(r.measure_names, name), :), 'omitnan');
band_lbl = sprintf('%g to 0 s^{-1}', -S.band / S.tau_d);
% The axis names in two lines: at the canvas width a panel is about 1.2 in
% wide and its axis about 1.5 in tall.
panels = { ...
    @(r) mrow(r, 'm_M', 'shell') / uniform_shell, @(r) mrow(r, 'm_theory', 'shell') / uniform_shell, {'Rim occupancy /', 'uniform disk'}; ...
    @(r) mrow(r, 'm_M', 'n_band'),                @(r) mrow(r, 'm_uniform', 'n_band'),               {'Eigenvalues,', band_lbl}; ...
    @(r) mrow(r, 'm_M', 'n_unstable'),            @(r) mean(r.full.n_unstable),                     {'Eigenvalues', 'above 0'}; ...
    @(r) mean(r.full.n_band_x),                   @(r) NaN,                                         {'Full Jacobian,', band_lbl}};

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 34 11]);
if all(isfinite(cfg.fig_size)); fig.Units = 'inches'; fig.Position(3:4) = cfg.fig_size; end
tl = tiledlayout(fig, 1, size(panels, 1), 'TileSpacing', 'compact', 'Padding', 'compact');
for q = 1:size(panels, 1)
    ax = nexttile(tl); hold(ax, 'on');
    ax.Tag = sprintf('eigdens_measures_%d', q);
    if q == 1; plot(ax, [0.4 n_c + 0.6], [1 1], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8); end
    y_hi = 0;
    for c = 1:n_c
        rs = F.by_cond{c};
        ym = arrayfun(panels{q, 1}, rs); yc = arrayfun(panels{q, 2}, rs);
        has_ref = any(isfinite(yc));
        x_box = c - 0.14 * has_ref;                   % the box left of the reference column, or centered when there is none
        box_whisker(ax, x_box, ym(isfinite(ym)), F.colors{c});
        if has_ref
            plot(ax, repmat(c + 0.24, 1, nnz(isfinite(yc))), yc(isfinite(yc)), 'd', 'Color', gray, ...
                'MarkerFaceColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.2, 'Tag', 'eigdens_reference');
        end
        y_hi = max([y_hi, yc, ym], [], 'omitnan');
    end
    if q == 1; y_hi = max(y_hi, 1); end
    if ~(y_hi > 0); y_hi = 1; end
    set(ax, ax_style{:}, 'XLim', [0.4 n_c + 0.6], 'XTick', 1:n_c, 'XTickLabel', F.ticks, ...
        'YLim', [0 y_hi * 1.08], 'XTickLabelRotation', 45);
    ylabel(ax, panels{q, 3}, 'FontSize', st.label_fs);
end

fig_tag = ['Fig_Eigdens_Measures' F.suffix];
runs = [F.by_cond{:}];
rows = arrayfun(@(r) table_row(r, uniform_shell), runs, 'UniformOutput', false);
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        fprintf(fid, '# %s\n\n', fig_tag);
        fprintf(fid, ['The eigenvalue-density measures by condition; network `%s`; one marker per seed (the mean over its ' ...
            'calibration states). Definitions: `eigdens_measures`. PANEL 1: fraction of the eigenvalues of M = W D beyond ' ...
            '%.1f R_eff over the %.2f of a uniform disk; filled: measured; open diamond: theory with the saved gains; ' ...
            'dotted: uniform. PANEL 2: eigenvalues of the frozen-gain dendritic block with rate in (%g, 0] 1/s; filled: ' ...
            'measured; open diamond: a uniform disk of the same radius. PANEL 3: eigenvalues with rate above 0; filled: ' ...
            'frozen gains; open diamond: the full Jacobian at the same states. PANEL 4: eigenvalues of the full Jacobian ' ...
            'in (%g, 0] 1/s whose eigenvector has at least half its norm^2 in the dendritic variables (the fast modes). ' ...
            'Data: `%s`.\n\n'], F.network, S.shell_r, uniform_shell, -S.band / S.tau_d, -S.band / S.tau_d, F.file);
        fprintf(fid, ['Per run, means over its states. `p0`: zero-gain fraction; `PR`: <d^2>^2/<d^4>; `R_eff`: theoretical ' ...
            'bulk radius (units of W; the boundary is Re = 1); `rim`: measured / null draws / theory / W under a uniform ' ...
            'gain; `r50`: median modulus over R_eff, measured / theory (uniform disk 0.707); `band` and `unst`: measured / ' ...
            'theory / uniform disk; `full`: the full Jacobian, all eigenvalues in the band / those mostly dendritic / ' ...
            'unstable / unstable and mostly dendritic; `lambda_1`: the run''s leading Lyapunov exponent (NaN: not ' ...
            'computed); `gap`: largest |measured - theory| of the cumulative radial fraction inside the bulk.\n\n']);
        fprintf(fid, ['| condition | seeds | states | p0 | PR | R_eff | rim | r50 | band | unst | full | lambda_1 | gap |\n' ...
            '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n']);
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
end

%% ------------------------------------------------------------------------
function s = table_row(r, uniform_shell)
m = @(M, name) mean(r.(M)(strcmp(r.measure_names, name), :), 'omitnan');
ttl = grouped_condition_label({r.condition}, 'sentence');
s = sprintf(['| %s | %d-%d | %d | %.3f | %.3f | %.3f | %.3f / %.3f / %.3f / %.3f | %.3f / %.3f | ' ...
    '%.1f / %.1f / %.1f | %.1f / %.1f / %.1f | %.1f / %.1f / %.1f / %.1f | %.3g | %.3f |'], ...
    ttl{1}, r.seed(1), r.seed(2), numel(r.R_eff), mean(r.p0), mean(r.pr, 'omitnan'), mean(r.R_eff), ...
    m('m_M', 'shell') / uniform_shell, m('m_null', 'shell') / uniform_shell, ...
    m('m_theory', 'shell') / uniform_shell, m('m_unigain', 'shell') / uniform_shell, ...
    m('m_M', 'r50'), m('m_theory', 'r50'), ...
    m('m_M', 'n_band'), m('m_theory', 'n_band'), m('m_uniform', 'n_band'), ...
    m('m_M', 'n_unstable'), m('m_theory', 'n_unstable'), m('m_uniform', 'n_unstable'), ...
    mean(r.full.n_band), mean(r.full.n_band_x), mean(r.full.n_unstable), mean(r.full.n_unstable_x), ...
    r.lambda_1, max(r.theory_gap));
end

function box_whisker(ax, x, v, col)
% A box over the values v at x: 25th to 75th percentile (prctile), the
% median as a thick bar, whiskers to the most extreme value within 1.5
% interquartile ranges of the box, values beyond them as dots.
if isempty(v); return; end
q = prctile(v, [25 50 75]);
iqr_ = q(3) - q(1);
inside = v(v >= q(1) - 1.5 * iqr_ & v <= q(3) + 1.5 * iqr_);
w = 0.13;
plot(ax, [x x], [min(inside) q(1)], '-', 'Color', col, 'LineWidth', 1.2, 'Tag', 'eigdens_whisker');
plot(ax, [x x], [q(3) max(inside)], '-', 'Color', col, 'LineWidth', 1.2, 'Tag', 'eigdens_whisker');
patch(ax, x + w * [-1 1 1 -1], q([1 1 3 3]), col, 'FaceAlpha', 0.35, 'EdgeColor', col, 'LineWidth', 1, 'Tag', 'eigdens_box');
plot(ax, x + w * [-1 1], q([2 2]), '-', 'Color', col, 'LineWidth', 2.5, 'Tag', 'eigdens_median');
out = v(v < min(inside) | v > max(inside));
if ~isempty(out)
    plot(ax, repmat(x, size(out)), out, '.', 'Color', col, 'MarkerSize', 10, 'Tag', 'eigdens_outlier');
end
end

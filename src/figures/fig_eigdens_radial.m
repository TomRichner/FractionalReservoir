function out = fig_eigdens_radial(cfg)
% FIG_EIGDENS_RADIAL How the eigenvalues of the frozen-gain connectivity are spread in radius, by condition, against the theory.
%
%   out = FIG_EIGDENS_RADIAL('run_dir', d)
%   out = FIG_EIGDENS_RADIAL('run_dir', d, 'network', 'EEx1p5')
%   out = FIG_EIGDENS_RADIAL('data_file', f, 'out_dir', o)
%
% Draws what the opt-in stage eigenvalue_density saved (FD-042;
% <run_dir>/eigenvalue_density/eigenvalue_density_data.mat). Reads one run
% directory or the file given, errors when neither holds it, and simulates
% nothing. Supplementary.
%
% THE OBJECT. M = W diag(d), the frozen-gain effective connectivity at a
% calibration state; radii are in units of the theoretical bulk radius
% R_eff of THAT state (from its gains and the block variances of W).
%
% LEFT: the fraction of eigenvalues of M within radius r. RIGHT: the density
% relative to a uniform disk of radius R_eff, in equal-area bins of
% (r / R_eff)^2, the eigenvalues exactly at the origin (the zero-gain
% neurons) left out, log axis. In both, per condition in its color: the
% MEDIAN over every seed and calibration state of the network (thick) with
% the interquartile range shaded; thin dashed black: the median of the
% theory evaluated with the saved gains at the same states
% (eigdens_radial_theory_blocks); dotted gray: a uniform disk ((r/R)^2 on
% the left, 1 on the right), which is what W alone gives in the paper's
% network (Harris et al. 2023, Eq. 20, g = 0).
%
% NO TEXT ON THE PANELS. The colors are the conditions' (manuscript_style);
% the numbers (per condition: seeds, states, the zero-gain fraction, R_eff,
% the largest gap between measured and theory, the rim occupancy) are in
% Fig_Eigdens_Radial<suffix>_table.md.
%
% OPTIONS: network (default 'reference'; any network the stage's file
% holds, e.g. 'EEx1p5'; the tag then ends in _<network>), conditions
% (names to draw, in order; empty: all the network has), run_dir,
% data_file, out_dir, save, visible, verbose, preset_name (unused).
%
% See also: run_eigenvalue_density_stage, eigdens_run_measures,
%           eigdens_figure_data, fig_eigdens_gains, fig_eigdens_measures

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the networks are in the data
    cfg.network     (1,:) char    = 'reference'
    cfg.panels      (1,:) cell    = {'cumulative', 'density'}   % which of the two panels, left to right
    cfg.fig_size    (1,2) double  = [NaN NaN]   % inches [w h]; NaN = 13 cm per panel by 11 cm (fig_supp6_grouped lays it out at its width)
    cfg.conditions                = {}
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = eigdens_figure_data(cfg, 'fig_eigdens_radial');
st = F.st; G = F.D.settings.grids; S = F.D.settings.measures;
u = G.u; u2e = G.u2_edges; u2c = (u2e(1:end - 1) + u2e(2:end)) / 2;
n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];
ax_style = {'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out'};
stack = @(rs, f) cat(1, rs.(f));               % seeds x states, as rows

vis = 'on'; if ~cfg.visible; vis = 'off'; end
panels = cfg.panels(ismember(cfg.panels, {'cumulative', 'density'}));
if isempty(panels) || numel(panels) ~= numel(cfg.panels)
    error('fig_eigdens_radial:Panels', 'panels holds ''cumulative'' and/or ''density''.');
end
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 13 * numel(panels) 11]);
if all(isfinite(cfg.fig_size)); fig.Units = 'inches'; fig.Position(3:4) = cfg.fig_size; end
tl = tiledlayout(fig, 1, numel(panels), 'TileSpacing', 'compact', 'Padding', 'compact');

% Left: cumulative fraction within r.
if ismember('cumulative', panels)
ax = nexttile(tl); hold(ax, 'on');
plot(ax, u, min(u.^2, 1), ':', 'Color', gray, 'LineWidth', 1.8);
for c = 1:n_c
    [m, lo, hi] = med_iqr(stack(F.by_cond{c}, 'cdf_M'));
    band(ax, u, lo, hi, F.colors{c});
    plot(ax, u, m, '-', 'Color', F.colors{c}, 'LineWidth', 3);
end
for c = 1:n_c
    plot(ax, u, med_iqr(stack(F.by_cond{c}, 'cdf_theory')), '--', 'Color', 'k', 'LineWidth', 1.1);
end
set(ax, ax_style{:}, 'XLim', [0 1.2], 'YLim', [0 1], 'XTick', 0:0.5:1, 'YTick', 0:0.2:1, 'XTickLabelRotation', 0);   % three x ticks: a panel of S6 is about 1.3 in wide
xlabel(ax, 'r / R_{eff}', 'FontSize', st.label_fs);
ylabel(ax, {'Fraction of eigenvalues', 'within r'}, 'FontSize', st.label_fs);
end

% Right: density relative to a uniform disk.
y_floor = 0.03;
if ismember('density', panels)
ax = nexttile(tl); hold(ax, 'on');
plot(ax, [0 1.44], [1 1], ':', 'Color', gray, 'LineWidth', 1.8);
inside = u2c < 1;
for c = 1:n_c
    [m, lo, hi] = med_iqr(stack(F.by_cond{c}, 'dens_M'));
    [xs, ylo] = stairs(u2e, [max(lo, y_floor), max(lo(end), y_floor)]);
    [~, yhi] = stairs(u2e, [max(hi, y_floor), max(hi(end), y_floor)]);
    band(ax, xs', ylo', yhi', F.colors{c});
    stairs(ax, u2e, [max(m, y_floor / 2), NaN], '-', 'Color', F.colors{c}, 'LineWidth', 2.5);
end
for c = 1:n_c
    th = med_iqr(stack(F.by_cond{c}, 'dens_theory'));
    plot(ax, u2c(inside), th(inside), '--', 'Color', 'k', 'LineWidth', 1.1);
end
set(ax, ax_style{:}, 'XLim', [0 1.2], 'YScale', 'log', 'YLim', [y_floor 5], 'XTick', 0:0.2:1.2, ...
    'YTick', [0.03 0.1 0.3 1 3], 'YTickLabel', {'0.03', '0.1', '0.3', '1', '3'});
xlabel(ax, '(r / R_{eff})^2', 'FontSize', st.label_fs);
ylabel(ax, 'Density / uniform disk', 'FontSize', st.label_fs);
end

fig_tag = ['Fig_Eigdens_Radial' F.suffix];
rows = cell(1, n_c);
uniform_shell = 1 - S.shell_r^2;
for c = 1:n_c
    rs = F.by_cond{c};
    shell = stack(rs, 'm_M'); shell = shell(1:4:end, :);          % m_M is 4 x K per run: rows 1, 5, ... are 'shell'
    shell_th = stack(rs, 'm_theory'); shell_th = shell_th(1:4:end, :);
    rows{c} = sprintf('| %s | %d | %d | %.3f | %.3f | %.3f | %.3f / %.3f | %.4f |', F.titles{c}, numel(rs), ...
        numel(stack(rs, 'R_eff')), med(stack(rs, 'p0')), med(stack(rs, 'pr')), med(stack(rs, 'R_eff')), ...
        med(shell) / uniform_shell, med(shell_th) / uniform_shell, max(stack(rs, 'theory_gap'), [], 'all'));
end
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        fprintf(fid, '# %s\n\n', fig_tag);
        fprintf(fid, ['How the eigenvalues of the frozen-gain effective connectivity M = W D are spread in radius; ' ...
            'network `%s`. Radii in units of the theoretical bulk radius R_eff of each state (from its saved gains and ' ...
            'the block variances of W). LEFT: fraction of eigenvalues within r. RIGHT: density relative to a uniform disk ' ...
            'of radius R_eff in equal-area bins of (r/R_eff)^2 of width %.2g, the eigenvalues exactly at the origin (the ' ...
            'zero-gain neurons) left out, log axis (values below %.2g drawn at the axis). Thick lines: the median over every ' ...
            'seed and calibration state, per condition in its color; shaded: the interquartile range over the same ' ...
            'states. Thin dashed black: the median of the random-matrix theory with the saved gains at the same states ' ...
            '(`eigdens_radial_theory_blocks`; Ahmadian et al. 2015, Eqs. 2.13 to 2.18, generalized to block variances). ' ...
            'Dotted gray: a uniform disk, which is what W alone gives here (Harris et al. 2023, Eq. 20: g = %.3g). The value ' ...
            'at r = 0 on the left is the fraction of neurons with zero gain. A curve above the dotted one on the left is a ' ...
            'spectrum concentrated toward the center. Beyond 1 are the eigenvalues outside the theoretical bulk (finite ' ...
            'size, and the local outliers of a matrix with E and I means). Data: `%s`.\n\n'], ...
            F.network, diff(u2e(1:2)), y_floor, F.net.harris.g, F.file);
        fprintf(fid, ['Medians over seeds and states. `p0`: fraction of neurons with zero gain. `PR`: <d^2>^2/<d^4>, the ' ...
            'density at the rim relative to uniform for a single-variance W. `rim`: fraction of eigenvalues beyond %.1f R_eff ' ...
            'over the %.2f of a uniform disk, measured / theory. `largest gap`: the largest |measured - theory| of the ' ...
            'cumulative fraction inside the bulk at any single state.\n\n'], S.shell_r, uniform_shell);
        fprintf(fid, '| Condition | runs (seeds) | states | p0 | PR | R_eff | rim, measured / theory | largest gap |\n');
        fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
end

%% ------------------------------------------------------------------------
function [m, lo, hi] = med_iqr(A)
% Median and quartiles down the rows (seeds x states), NaN rows ignored.
q = quantile(A, [0.5 0.25 0.75], 1);
m = q(1, :); lo = q(2, :); hi = q(3, :);
end

function v = med(A)
v = median(A(:), 'omitnan');
end

function band(ax, x, lo, hi, col)
x = x(:)'; lo = lo(:)'; hi = hi(:)';
ok = isfinite(lo) & isfinite(hi);
if ~any(ok); return; end
fill(ax, [x(ok), fliplr(x(ok))], [lo(ok), fliplr(hi(ok))], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end

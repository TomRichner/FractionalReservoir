function out = fig_eigdens_gains(cfg)
% FIG_EIGDENS_GAINS The gain distributions behind the eigenvalue densities, and the participation measure.
%
%   out = FIG_EIGDENS_GAINS('run_dir', d)
%   out = FIG_EIGDENS_GAINS('run_dir', d, 'network', 'EEx1p5')
%
% Draws what the opt-in stage eigenvalue_density saved (FD-042): the
% per-neuron gains d (slope of the nonlinearity times the available vesicle
% fraction) at the calibration states, as outlier_spectra saved them.
% Reads one run directory or the file given; simulates nothing.
% Supplementary.
%
% LEFT: the cumulative distribution of d over every neuron, calibration
% state and seed of the network, per condition in its color; the step at
% d = 0 is the fraction of neurons pinned on the floor or the ceiling of the
% nonlinearity. RIGHT: one dot per calibration state and seed: the
% participation ratio of the squared gains, <d^2>^2 / <d^4> (the density at
% the rim of the bulk relative to a uniform disk, Ahmadian et al. 2015,
% Eq. 2.18), against the fraction of neurons with nonzero gain; dotted:
% the two are equal when every nonzero gain is the same, so the distance
% below it is what the spread of the nonzero gains costs.
%
% No text on the panels; the numbers (medians and ranges by condition and
% by cell type) are in Fig_Eigdens_Gains<suffix>_table.md. Options as
% fig_eigdens_radial (network, conditions, run_dir, data_file, out_dir,
% save, visible, verbose).
%
% See also: fig_eigdens_radial, eigdens_run_measures, eigdens_figure_data

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
    cfg.fig_size    (1,2) double  = [NaN NaN]   % inches [w h]; NaN = 26 x 11 cm (fig_supp6_grouped lays it out at the canvas width)
    cfg.conditions                = {}
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = eigdens_figure_data(cfg, 'fig_eigdens_gains');
st = F.st;
n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];
ax_style = {'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out'};
pooled = cellfun(@(rs) double(reshape(cat(2, rs.gains), [], 1)), F.by_cond, 'UniformOutput', false);
g_max = max(1, max(cellfun(@max, pooled)));

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 26 11]);
if all(isfinite(cfg.fig_size)); fig.Units = 'inches'; fig.Position(3:4) = cfg.fig_size; end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
for c = 1:n_c
    g = sort(pooled{c});
    stairs(ax, [0; g; g_max], [0; (1:numel(g))' / numel(g); 1], '-', 'Color', F.colors{c}, 'LineWidth', 2.5);
end
set(ax, ax_style{:}, 'XLim', [0 g_max], 'YLim', [0 1], 'YTick', 0:0.2:1);
xlabel(ax, 'Gain d', 'FontSize', st.label_fs);
ylabel(ax, {'Fraction of neurons', 'with gain \leq d'}, 'FontSize', st.label_fs);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, [0 1], [0 1], ':', 'Color', gray, 'LineWidth', 1.8);
for c = 1:n_c
    rs = F.by_cond{c};
    plot(ax, 1 - [rs.p0], [rs.pr], '.', 'Color', F.colors{c}, 'MarkerSize', 12);
end
set(ax, ax_style{:}, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.5:1, 'YTick', 0:0.2:1, 'XTickLabelRotation', 0);   % three x ticks: a panel of S6 is about 1.5 in wide
xlabel(ax, {'Fraction of neurons', 'with nonzero gain'}, 'FontSize', st.label_fs);
ylabel(ax, '\langle d^2 \rangle^2 / \langle d^4 \rangle', 'FontSize', st.label_fs);

fig_tag = ['Fig_Eigdens_Gains' F.suffix];
rows = cell(1, n_c);
for c = 1:n_c
    rs = F.by_cond{c};
    bt = cat(3, rs.gain_by_type);                % C x 3 x runs: p0, mean and cv of the nonzero gains
    by_type = strjoin(arrayfun(@(b) sprintf('%.3f / %.3f / %.2f', median(bt(b, 1, :), 'omitnan'), ...
        median(bt(b, 2, :), 'omitnan'), median(bt(b, 3, :), 'omitnan')), 1:size(bt, 1), 'UniformOutput', false), ' | ');
    rows{c} = sprintf('| %s | %d | %s | %s | %s | %s | %s |', F.titles{c}, numel(rs), ...
        mmm([rs.p0]), mmm([rs.cv_nz]), mmm([rs.pr]), mmm([rs.rms_gain]), by_type);
end
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        n_types = size(F.net.s2, 1);
        fprintf(fid, '# %s\n\n', fig_tag);
        fprintf(fid, ['The gains behind the eigenvalue densities; network `%s`. The gain d of a neuron is the slope of ' ...
            'the nonlinearity times the available vesicle fraction, read at the calibration states of the stage ' ...
            'outlier_spectra. LEFT: cumulative distribution of d over every neuron, state and seed, per condition in ' ...
            'its color; the step at d = 0 is the fraction pinned on the floor or the ceiling of the nonlinearity. ' ...
            'RIGHT: one dot per state and seed: the participation ratio of the squared gains <d^2>^2/<d^4> (the density ' ...
            'at the rim of the bulk relative to a uniform disk, Ahmadian et al. 2015, Eq. 2.18) against the fraction of ' ...
            'neurons with nonzero gain; dotted: equality, reached when every nonzero gain is the same. Data: `%s`.\n\n'], ...
            F.network, F.file);
        fprintf(fid, ['Median [min, max] over states and seeds. `p0`: fraction with zero gain; `cv`: coefficient of ' ...
            'variation of the nonzero gains; `PR`: <d^2>^2/<d^4>; `rms`: root mean square gain. By type (%d types, in ' ...
            'the model''s order): p0 / mean of the nonzero gains / cv of the nonzero gains, all states of a run pooled, ' ...
            'median over runs.\n\n'], n_types);
        fprintf(fid, '| Condition | runs | p0 | cv | PR | rms | %s |\n', strjoin(arrayfun(@(b) sprintf('type %d', b), ...
            1:n_types, 'UniformOutput', false), ' | '));
        fprintf(fid, '|%s\n', repmat('---|', 1, 6 + n_types));
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
end

%% ------------------------------------------------------------------------
function s = mmm(v)
v = v(isfinite(v));
if isempty(v); s = 'NaN'; return; end
s = sprintf('%.3f [%.3f, %.3f]', median(v), min(v), max(v));
end

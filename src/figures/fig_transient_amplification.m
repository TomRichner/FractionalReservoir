function out = fig_transient_amplification(cfg)
% FIG_TRANSIENT_AMPLIFICATION Numerical vs spectral abscissa of the Jacobian.
%
%   out = FIG_TRANSIENT_AMPLIFICATION('run_dir', d)
%
% Reads the eig-heatmap stage's data (run_eig_heatmap samples the Jacobian at
% ~150 states per condition and records, per state, the NUMERICAL ABSCISSA
% omega = max eig((J_xx + J_xx')/2) and the SPECTRAL ABSCISSA
% alpha = max real eig(J_xx) of the DENDRITIC BLOCK J_xx -- the rate-network
% Jacobian with adaptation and depression held fixed; on the full J the
% numerical abscissa reflects the state coordinates' units, not dynamics,
% see sample_eigenvalues in run_eig_heatmap). Two panels:
%
%   (a) omega(t) (solid) and alpha(t) (dashed) per condition, condition
%       colours from manuscript_style. omega bounds the instantaneous growth
%       rate of any perturbation, alpha the growth of the eigen-directions;
%       the gap is non-normal transient amplification.
%   (b) the per-condition distribution of omega - alpha (box plots), with the
%       finite-time lambda_1 of each condition marked for reference.
%
% This is the Jacobian-side "transient divergence" measure, the companion of
% the trajectory-side one the sweeps store (fraction of time the local
% Lyapunov rate is positive; see sweep_metrics 'fpos', 'p95', 'exc'). A
% stable network (lambda_1 < 0, alpha < 0 most of the time) with omega > 0
% supports transient growth; a normal Jacobian would have omega = alpha.
%
% See also: run_eig_heatmap, fig_eig_heatmap, sweep_metrics

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'eig_heatmap')}, ...
    'eig_heatmap_data.mat', ...
    'Run run_eig_heatmap first');
D = load(data_file);
if ~isfield(D, 'num_abscissa_by_cond')
    error('fig_transient_amplification:OldData', ...
        ['%s predates the numerical-abscissa fields; rerun run_eig_heatmap ' ...
         '(2026-09-12 or later).'], data_file);
end
n_cond = numel(D.cond_names);
colors = cellfun(@(n) st.condition_color(n), D.cond_names, 'UniformOutput', false);
titles = D.condition_titles;

fig = figure('Color', 'w', 'Position', [100 100 1100 420]);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) time series
ax1 = nexttile(tl); hold(ax1, 'on');
h_leg = gobjects(1, n_cond);
for i = 1:n_cond
    t = D.J_times_by_cond{i};
    h_leg(i) = plot(ax1, t, D.num_abscissa_by_cond{i}, '-', 'Color', colors{i}, 'LineWidth', 1.8, ...
        'DisplayName', titles{i});
    plot(ax1, t, D.spec_abscissa_by_cond{i}, '--', 'Color', colors{i}, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
yline(ax1, 0, ':', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
hold(ax1, 'off');
xlabel(ax1, 'time (s)', 'FontSize', st.label_fs);
ylabel(ax1, 'abscissa of J_{xx} (1/s)', 'FontSize', st.label_fs);
title(ax1, '\omega(J_{xx}) solid, \alpha(J_{xx}) dashed', 'FontWeight', 'normal', 'FontSize', st.title_fs);
legend(ax1, h_leg, 'Location', 'best', 'FontSize', 11);
set(ax1, 'FontSize', st.tick_fs); box(ax1, 'off');

% (b) distributions of the non-normal margin
ax2 = nexttile(tl); hold(ax2, 'on');
for i = 1:n_cond
    m = D.num_abscissa_by_cond{i} - D.spec_abscissa_by_cond{i};
    q = prctile(m, [5 25 50 75 95]);
    plot(ax2, [i i], q([1 5]), '-', 'Color', colors{i}, 'LineWidth', 1.2);
    patch(ax2, i + 0.25 * [-1 1 1 -1], q([2 2 4 4]), colors{i}, 'FaceAlpha', 0.35, 'EdgeColor', colors{i});
    plot(ax2, i + 0.25 * [-1 1], q([3 3]), '-', 'Color', colors{i}, 'LineWidth', 2.5);
end
hold(ax2, 'off');
short = cellfun(@(n) short_of(st, n), D.cond_names, 'UniformOutput', false);
labels = arrayfun(@(i) sprintf('%s\\newline\\lambda_1 = %+.2f', short{i}, D.lle_by_cond(i)), 1:n_cond, 'UniformOutput', false);
set(ax2, 'XTick', 1:n_cond, 'XTickLabel', labels, 'FontSize', st.tick_fs - 2, 'XTickLabelRotation', 0);
xlim(ax2, [0.3, n_cond + 0.7]);
ylabel(ax2, '\omega(J_{xx}) - \alpha(J_{xx})  (1/s)', 'FontSize', st.label_fs);
title(ax2, 'non-normal margin per state (5-95%, IQR, median)', 'FontWeight', 'normal', 'FontSize', st.title_fs);
box(ax2, 'off');

if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_Transient_Amplification';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

function s = short_of(st, name)
if st.condition_short.isKey(name); s = st.condition_short(name); else; s = strrep(name, '_', ' '); end
end

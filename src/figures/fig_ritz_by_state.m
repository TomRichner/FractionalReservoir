function out = fig_ritz_by_state(cfg)
% FIG_RITZ_BY_STATE Ritz values of each sampled state against its time. A diagnostic, not for the paper.
%
%   out = FIG_RITZ_BY_STATE('run_dir', d)
%   out = FIG_RITZ_BY_STATE('run_dir', d, 'ritz_scaling', 'type_sd', 'ritz_k', 25)
%
% The Ritz overlay of fig_eig_heatmap pools the Ritz values of every sampled
% state. Here they are kept apart (FD-033; Tom, 2026-11-19: are they per
% snapshot, and do the points right of zero come from a few states or from
% all of them?). One row of panels per condition:
%   left    the real part of every Ritz value of a state against the
%           state's time (k points per state, in the condition's color;
%           right of zero in red), with, per state, the largest real part of
%           the FULL Jacobian's spectrum (black line), the largest
%           eigenvalue of (H + H')/2 (the largest instantaneous growth rate
%           inside the tracked subspace, which bounds the Ritz real parts;
%           gray line) and the pass's largest Lyapunov exponent (dashed)
%   middle  per state, how many of the k Ritz values lie right of zero
%           (bars), how many DIAGONAL entries of H, the instantaneous local
%           Lyapunov rates, are positive (black), and how many eigenvalues of
%           the full Jacobian have positive real part (blue)
%   right   the Ritz residual ||J Q - Q H||_F / ||J Q||_F per state: 0 if
%           the tracked subspace were invariant under the Jacobian
% Needs a run of run_eig_heatmap with the Ritz pass (ritz_K > 0).
%
% See also: fig_eig_heatmap, run_eig_heatmap, ritz_diagnostics

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.conditions                = {}
    cfg.ritz_k       (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    cfg.ritz_scaling (1,:) char {mustBeMember(cfg.ritz_scaling, {'raw', 'type_sd'})} = 'raw'
end
setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
data_file = resolve_data_file(cfg.data_file, cfg.run_dir, {fullfile(cfg.run_dir, 'eig_heatmap')}, ...
    'eig_heatmap_data.mat', 'Run run_eig_heatmap first');
D = load(data_file);
if ~isfield(D, 'ritz_by_cond')
    error('fig_ritz_by_state:NoRitz', '%s holds no Ritz pass (run_eig_heatmap with ''ritz_K'' > 0).', data_file);
end
[cond_names, idx] = select_conditions(D.cond_names, cfg.conditions);
titles = titles_by_name(cond_names, D.condition_titles(idx));
st = manuscript_style();
n_cond = numel(cond_names);
fig = figure('Color', 'w', 'Position', [100, 100, 1250, 250 * n_cond]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, n_cond, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:n_cond
    R = D.ritz_by_cond{idx(i)};
    sc = find(strcmp(R.scalings, cfg.ritz_scaling), 1);
    if isempty(sc)
        error('fig_ritz_by_state:NoRitzScaling', 'The run holds no ''%s'' Ritz pass.', cfg.ritz_scaling);
    end
    k = cfg.ritz_k; if k == 0; k = R.K; end
    k = min(k, R.K);
    states = find(R.valid(:))'; t = R.t(states);
    col = [0 0 0];
    if st.condition_color.isKey(cond_names{i}); col = st.condition_color(cond_names{i}); end
    re = nan(k, numel(states)); n_diag = nan(1, numel(states)); res = n_diag; om = n_diag;
    for s = 1:numel(states)
        Dg = ritz_diagnostics(R.H(:, :, states(s), sc), R.JQ_gram(:, :, states(s), sc), k);
        re(:, s) = real(Dg.values); n_diag(s) = nnz(Dg.diag > 0); res(s) = Dg.residual; om(s) = Dg.num_abscissa;
    end
    T = repmat(t(:)', k, 1);

    ax = nexttile(tl); hold(ax, 'on');
    scatter(ax, T(re <= 0), re(re <= 0), 5, col, 'filled', 'MarkerFaceAlpha', 0.35, 'Tag', 'ritz_by_state');
    scatter(ax, T(re > 0), re(re > 0), 7, [0.85 0.1 0.1], 'filled', 'MarkerFaceAlpha', 0.7, 'Tag', 'ritz_by_state_positive');
    h1 = plot(ax, t, R.full_max_re(states), '-', 'Color', [0 0 0], 'LineWidth', 1);
    h2 = plot(ax, t, om, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    h3 = yline(ax, R.LE_spectrum(1, sc), '--', 'Color', col, 'LineWidth', 1);
    yline(ax, 0, ':', 'Color', st.zeroline_color, 'LineWidth', 1.2);
    ylabel(ax, {titles{i}, 'Re(\theta) (1/s)'}); box(ax, 'off');
    if i == 1
        title(ax, sprintf('Ritz values of each state (k = %d, %s norm)', k, strrep(cfg.ritz_scaling, '_', '\_')), 'FontWeight', 'normal');
        legend(ax, [h1 h2 h3], {'full J: max Re', 'max eig (H+H'')/2', '\lambda_1 of the pass'}, ...
            'Location', 'best', 'Box', 'off', 'FontSize', 8);
    end
    if i == n_cond; xlabel(ax, 'time of the sampled state (s)'); end

    ax = nexttile(tl); hold(ax, 'on');
    bar(ax, t, sum(re > 0, 1), 1, 'FaceColor', [0.85 0.1 0.1], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(ax, t, n_diag, '-o', 'Color', [0 0 0], 'MarkerSize', 3, 'LineWidth', 0.8);
    plot(ax, t, R.full_n_positive(states), '-', 'Color', [0.2 0.4 0.9], 'LineWidth', 0.8);
    ylabel(ax, 'count per state'); box(ax, 'off');
    if i == 1
        title(ax, 'right of zero, per state', 'FontWeight', 'normal');
        legend(ax, {sprintf('Ritz values (of %d)', k), 'local rates (diag H)', 'eigenvalues of the full J'}, ...
            'Location', 'best', 'Box', 'off', 'FontSize', 8);
    end
    if i == n_cond; xlabel(ax, 'time of the sampled state (s)'); end

    ax = nexttile(tl);
    plot(ax, t, res, '-o', 'Color', col, 'MarkerSize', 3, 'LineWidth', 0.8);
    ylim(ax, [0 1]); ylabel(ax, 'residual'); box(ax, 'off');
    if i == 1; title(ax, '||JQ - QH||_F / ||JQ||_F (0 = invariant subspace)', 'FontWeight', 'normal'); end
    if i == n_cond; xlabel(ax, 'time of the sampled state (s)'); end
end
fig_tag = 'fig_ritz_by_state';
if strcmp(cfg.ritz_scaling, 'type_sd'); fig_tag = [fig_tag '_typesd']; end
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

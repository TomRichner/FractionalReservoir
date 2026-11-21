function out = fig_ostruct_h1(cfg)
% FIG_OSTRUCT_H1 Hypothesis H1 of the outlier_structures stage: the observed shift of the intended mode against the first-order prediction from the gains of its carriers.
%
%   out = FIG_OSTRUCT_H1('run_dir', d)
%
% Draws what the opt-in stage outlier_structures saved (FD-043). Reads one
% run directory or the file given; simulates nothing. Supplementary.
%
% H1. Gains scale the columns of W, so to first order an eigenvalue moves
% by the participation-weighted gain of the neurons that carry it
% (ospec_mode_reference). The relative shift of the mode from the
% uniform-gain picture, lam / (d_rms lam0) - 1 (complex; 0: the mode moved
% exactly as the bulk did; negative real part: it lost more than the bulk),
% is drawn OBSERVED (vertical: the followed eigenvalue of the frozen-gain
% matrix) against PREDICTED (horizontal: l0' M r0 / (l0' r0), from the
% structural left and right eigenvectors and the state's gains), real
% parts left and imaginary parts right. One point per RUN (every network,
% condition and seed; the median over the followed states of the run):
% color the condition, marker the group (circles added structure,
% squares local outliers and the reference, triangles block means). Gray:
% equality; dotted: zero. Equal axes, one range per panel covering every
% point. No text on the panels but axis labels.
% Fig_Ostruct_H1_table.md: the caption and the numbers (median and range
% over seeds per network and condition).
%
% Options: run_dir, data_file, out_dir, save, visible, verbose,
% preset_name (unused), conditions.
%
% See also: ostruct_figure_data, ostruct_mode_measures, ospec_snapshot

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

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = ostruct_figure_data(cfg, 'fig_ostruct_h1');
st = F.st;
T = F.T;
n_net = numel(F.keys); n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];
marks = {'o', 's', '^'};

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 9.5 4.6], 'Visible', vis);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
parts = {'re', 'im'};
plab = {'real part', 'imaginary part'};
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    xo = T.(['shift_pred_' parts{p}]); yo = T.(['shift_obs_' parts{p}]);
    sel = ismember(T.condition, F.cond_names);
    v = [xo(sel); yo(sel); 0];
    v = v(isfinite(v));
    lim = [min(v), max(v)] + [-1 1] * max(0.05, 0.05 * (max(v) - min(v)));
    plot(ax, lim, lim, '-', 'Color', gray, 'LineWidth', 1.5);
    xline(ax, 0, ':', 'Color', gray); yline(ax, 0, ':', 'Color', gray);
    for w = 1:n_net
        for c = 1:n_c
            k = F.row_of(w, c);
            if isempty(k); continue; end
            plot(ax, xo(k), yo(k), marks{F.group(w)}, 'MarkerSize', 8, 'MarkerFaceColor', F.colors{c}, ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.75, 'LineStyle', 'none');
        end
    end
    axis(ax, 'equal'); xlim(ax, lim); ylim(ax, lim);
    set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    xlabel(ax, sprintf('predicted shift, %s', plab{p}), 'FontSize', st.label_fs - 2);
    ylabel(ax, sprintf('observed shift, %s', plab{p}), 'FontSize', st.label_fs - 2);
end

fig_tag = 'Fig_Ostruct_H1';
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    sel = ismember(T.condition, F.cond_names);
    r_re = corr_(T.shift_pred_re(sel), T.shift_obs_re(sel));
    r_im = corr_(T.shift_pred_im(sel), T.shift_obs_im(sel));
    caption = sprintf(['H1 across every network, condition and seed of the outlier_structures stage. The shift of the ' ...
        'intended mode''s eigenvalue from the uniform-gain picture, lambda / (d_rms lambda_0) - 1 (complex; 0: the ' ...
        'mode moved exactly as the bulk did; negative real part: it lost more than the bulk), observed (vertical; the ' ...
        'eigenvalue of the frozen-gain matrix followed by its eigenvector) against the first-order prediction from the ' ...
        'gains of the neurons that carry it, l0'' M r0 / (l0'' r0) / (d_rms lambda_0) - 1 with l0 and r0 the ' ...
        'structural left and right eigenvectors (horizontal). LEFT: real parts. RIGHT: imaginary parts. One point per ' ...
        'run (%d runs, %s per network; median over the followed states of the run); color: condition; circles: added ' ...
        'low-rank structure; squares: local outliers and the reference (no structure put in); triangles: block-mean ' ...
        'outliers. Gray: equality; dotted: zero. Correlation of observed with predicted over the runs: %.3f (real), ' ...
        '%.3f (imaginary). mean_over_rms is the mean gain over the rms gain (what a mode spread evenly over all ' ...
        'neurons loses against the bulk from the spread of the gains alone); gain_ratio the participation-weighted ' ...
        'gain over the mean gain (what the carriers lose beyond the average neuron); mode_ov the overlap of the ' ...
        'followed eigenvector with the structural one (where it is low the observed value is only the best match).'], ...
        nnz(sel), F.seeds_text, r_re, r_im);
    cols = {'lam0_re_R', 'lam0_im_R', 'cond0', 'p_E', 'mode_ov', 'shift_obs_re', 'shift_pred_re', 'shift_obs_im', ...
        'shift_pred_im', 'mean_over_rms', 'gain_ratio', 'gain_top_rest', 'on_top', 'on_rest', 'frac_on_E', 'frac_on_I', ...
        'gain_ratio_EI'};
    name = ostruct_write_table(out_dir, fig_tag, caption, F, cols, 1:n_net);
    if ~isempty(name); out.files{end + 1} = name; end
end
vprintf(cfg.verbose, 'verbose', '%s: %d runs\n', fig_tag, height(T));
end

%% ------------------------------------------------------------------------
function r = corr_(a, b)
ok = isfinite(a) & isfinite(b);
a = a(ok) - mean(a(ok)); b = b(ok) - mean(b(ok));
den = sqrt(sum(a.^2) * sum(b.^2));
if den == 0; r = NaN; else; r = sum(a .* b) / den; end
end

function out = fig_ostruct_rotation(cfg)
% FIG_OSTRUCT_ROTATION The rotation sweep of the outlier_structures stage: is a rotating mode reined in as a real one is?
%
%   out = FIG_OSTRUCT_ROTATION('run_dir', d)
%
% Draws what the opt-in stage outlier_structures saved (FD-043), for the
% networks with ADDED low-rank structure (group 'rotation'). Reads one run
% directory or the file given; simulates nothing. Supplementary.
%
% THE HYPOTHESIS (H2 of the exploration): adaptation and depression follow
% the rate, which is slow and one-signed, so a real or slowly rotating
% mode is reined in and a fast-rotating one is not.
%
% WHAT IS DRAWN, against the rotation rate of the mode in the uniform-gain
% picture, Im(lam0) d_rms / tau_d (rad/s; it differs between conditions
% because the rms gain does), one line per condition, circles for the
% networks at about 1.4 bulk radii and squares for those at about 2.5.
% Every marker is the median over seeds, the vertical bar the range over
% seeds, the horizontal position the median over seeds.
%   a, b  the fraction of the mode's structural detachment from the bulk
%         radius that the actual gains removed (1, gray line: back at the
%         radius; more: inside the bulk). WITHOUT ADAPTATION ALONE IN a,
%         WITH ADAPTATION IN b, each panel with its own vertical range:
%         without adaptation the network is pinned, the mode is pushed far
%         inside the bulk and the fraction reaches 3 and more, which on a
%         shared axis flattened the adapted lines into the line at 1 (the
%         exploration's figure). Nothing is clipped.
%   c     the real part of the mode as a rate, frozen gains (open,
%         dashed) and in the full Jacobian (filled, solid); dashed black
%         zero.
%   d     the rotation rate of the mode's continuation in the full
%         Jacobian; gray: equal to the uniform-gain rate.
% No text on the panels but axis labels. Fig_Ostruct_Rotation_table.md:
% the caption and the numbers (median and range over seeds).
%
% Options: run_dir, data_file, out_dir, save, visible, verbose,
% preset_name (unused), conditions.
%
% See also: ostruct_figure_data, fig_ostruct_summary, fig_ostruct_plane

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
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
F = ostruct_figure_data(cfg, 'fig_ostruct_rotation');
st = F.st;
w_idx = find(F.group == find(strcmp(F.group_names, 'rotation')));
if isempty(w_idx)
    error('fig_ostruct_rotation:NoNetwork', 'The stage in %s has no finished run of an added-structure network.', F.folder);
end
n_c = numel(F.cond_names);
gray = [0.55 0.55 0.55];
is_none = strcmp(F.cond_names, 'no_adaptation');
strong = F.med.lam0_re_R(w_idx, 1) > 2;                 % the 2.5 R networks
for c = 2:n_c
    strong = strong | F.med.lam0_re_R(w_idx, c) > 2;
end

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 16 4.4], 'Visible', vis);
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = {'reined_none', 'reined_adapted', 'rate', 'rot'};
ylab = {'Fraction of detachment removed', 'Fraction of detachment removed', 'Real part of the mode (s^{-1})', ...
    'Rotation rate, full Jacobian (rad s^{-1})'};
axs = gobjects(1, 4);
for p = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    switch panels{p}
        case 'reined_none';    conds = find(is_none);
        case 'reined_adapted'; conds = find(~is_none);
        otherwise;             conds = 1:n_c;
    end
    for c = conds
        col = F.colors{c};
        for s = [false true]
            w = w_idx(strong == s);
            w = w(F.n_seed(w, c) > 0);
            if isempty(w); continue; end
            x = F.med.rot_struct(w, c);
            [x, o] = sort(x); w = w(o);
            mk = 'o'; if s; mk = 's'; end
            switch panels{p}
                case {'reined_none', 'reined_adapted'}
                    series(ax, x, w, c, F, 'mode_reined', ['-' mk], col, true);
                case 'rate'
                    series(ax, x, w, c, F, 'rate_frozen', ['--' mk], col, false);
                    series(ax, x, w, c, F, 'rate_full', ['-' mk], col, true);
                otherwise
                    series(ax, x, w, c, F, 'rot_full', ['-' mk], col, true);
            end
        end
    end
    switch panels{p}
        case {'reined_none', 'reined_adapted'}
            yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
        case 'rate'
            yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
    end
    set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    ylabel(ax, ylab{p}, 'FontSize', st.label_fs - 2);
    axs(p) = ax;
end
x_hi = max(cellfun(@(l) l(2), get(axs, 'XLim')));
set(axs, 'XLim', [0, x_hi]);                            % one horizontal range in every panel
y4 = ylim(axs(4)); hi4 = max(y4(2), x_hi);
plot(axs(4), [0 hi4], [0 hi4], '-', 'Color', gray, 'LineWidth', 1.5);   % equality
set(axs(4), 'YLim', [min(0, y4(1)), hi4]);
xlabel(tl, 'Rotation rate of the mode, uniform gain (rad s^{-1})', 'FontSize', st.label_fs - 2);

fig_tag = 'Fig_Ostruct_Rotation';
out = struct('figs', fig, 'files', {{}}, 'source', F.file);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    caption = sprintf(['The rotation sweep: low-rank structure added to the weights with a real part of about 1.4 ' ...
        'structural bulk radii and an imaginary part from 0 to about one bulk radius (circles; the rank-one, real case ' ...
        'at 0), and two networks at about 2.5 bulk radii, real and with an imaginary part of one radius (squares). ' ...
        'Horizontal axis: the rotation rate of the intended mode in the uniform-gain picture, Im(lambda_0) d_rms / ' ...
        'tau_d, which differs between conditions because the rms gain does. One line per condition; markers the median ' ...
        'over %s, vertical bars the range over seeds. (a, b) The fraction of the mode''s structural detachment from ' ...
        'the theoretical bulk radius that the actual gains removed (1, gray line: back at the radius; more: inside the ' ...
        'bulk): without adaptation in (a), with adaptation in (b), each on its own vertical range so that the ' ...
        'no-adaptation values (a pinned network pushes the mode far inside the bulk) do not flatten the adapted ones; ' ...
        'nothing is clipped. (c) The real part of the mode as a rate with the slow variables frozen (open, dashed) and ' ...
        'free, its continuation in the full Jacobian (filled, solid); dashed black: zero. (d) The rotation rate of the ' ...
        'continuation in the full Jacobian; gray: equal to the uniform-gain rate. Full-Jacobian values are medians over ' ...
        'the calibration states of a run, frozen ones over the followed states.'], F.seeds_text);
    cols = {'lam0_re_R', 'lam0_im_R', 'mode_rho', 'mode_rho_uniform', 'mode_reined', 'mode_reined_abs', 'mode_ov', ...
        'rate_struct', 'rate_frozen', 'rate_full', 'slow_shift', 'rot_struct', 'rot_frozen', 'rot_full', 'cont_ov', ...
        'gain_ratio', 'mean_over_rms', 'shift_obs_re', 'shift_pred_re', 'frac_on_E', 'frac_on_I'};
    name = ostruct_write_table(out_dir, fig_tag, caption, F, cols, w_idx);
    if ~isempty(name); out.files{end + 1} = name; end
end
vprintf(cfg.verbose, 'verbose', '%s: %d networks, %d conditions\n', fig_tag, numel(w_idx), n_c);
end

%% ------------------------------------------------------------------------
function series(ax, x, w, c, F, col_name, spec, col, filled)
% One line over networks w of condition c: medians, with the range over seeds.
y = F.med.(col_name)(w, c);
lo = F.lo.(col_name)(w, c); hi = F.hi.(col_name)(w, c);
for k = 1:numel(w)
    if isfinite(lo(k)) && hi(k) > lo(k)
        plot(ax, [x(k) x(k)], [lo(k) hi(k)], '-', 'Color', col, 'LineWidth', 1.0);
    end
end
face = 'w'; if filled; face = col; end
edge = col; if filled; edge = 'w'; end
plot(ax, x, y, spec, 'Color', col, 'LineWidth', 1.5 + 0.5 * filled, 'MarkerSize', 7 + filled, ...
    'MarkerFaceColor', face, 'MarkerEdgeColor', edge);
end

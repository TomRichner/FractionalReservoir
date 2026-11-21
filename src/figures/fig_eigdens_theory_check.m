function out = fig_eigdens_theory_check(cfg)
% FIG_EIGDENS_THEORY_CHECK The radial-density theory against sampled random matrices, five cases.
%
%   out = FIG_EIGDENS_THEORY_CHECK()
%   out = FIG_EIGDENS_THEORY_CHECK('out_dir', o, 'n', 500, 'n_draws', 20)
%
% The check behind the eigenvalue-density figures (FD-042):
% eigdens_check_theory samples n_draws Gaussian matrices of size n per case
% (uniform gains, broad log-normal gains, half the gains zero, Harris et
% al.'s two populations, a non-separable block profile with broad gains)
% and compares their cumulative radial distribution with the theory
% (eigdens_radial_theory_blocks). One panel per case: x the radius over the
% predicted bulk radius R, y the fraction of eigenvalues within it; thick
% blue: sampled; dashed black: theory; dotted gray: a uniform disk.
% Supplementary.
%
% HOW IT FITS THE FIGURE SIGNATURE. It takes the standard cfg (run_dir,
% out_dir, save, visible, verbose, preset_name, conditions) so the figure
% registry can call it like every other entry, but its data are not a
% run's: they are drawn here, from a private seeded stream
% (cfg.seed), so they are the same on every call and need no stage.
% run_dir is accepted and IGNORED (as fig_std_1TS_vs_MTS, which draws from
% the preset's route models). At the default size (500, 20 draws) the
% sampling takes about ten seconds.
%
% No text on the panels beyond one title line (the case's name); the
% numbers (largest differences against the sampling error, the agreement
% of the two solvers and of Harris's closed form) are in
% Fig_Eigdens_Theory_Check_table.md.
%
% See also: eigdens_check_theory, eigdens_radial_theory_blocks,
%           fig_eigdens_radial

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.images      (1,1) logical = true    % false with save: write the _table.md only (the composite figures)
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''    % accepted, unused: the data are sampled here
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}    % unused
    cfg.n           (1,1) double  = 500
    cfg.n_draws     (1,1) double  = 20
    cfg.seed        (1,1) double  = 1
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st = manuscript_style();
C = eigdens_check_theory('n', cfg.n, 'n_draws', cfg.n_draws, 'seed', cfg.seed);

vis = 'on'; if ~cfg.visible; vis = 'off'; end
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 34 8.5]);
tl = tiledlayout(fig, 1, numel(C), 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:numel(C)
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, C(c).u, C(c).u.^2, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    plot(ax, C(c).u, C(c).sampled, '-', 'Color', [0.000 0.447 0.698], 'LineWidth', 3);
    plot(ax, C(c).u, C(c).theory, '--', 'Color', 'k', 'LineWidth', 1.2);
    set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', ...
        'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.5:1, 'YTick', 0:0.5:1);
    title(ax, C(c).name, 'FontSize', st.title_fs, 'FontWeight', 'normal');
    if c == 1; ylabel(ax, 'Fraction within r', 'FontSize', st.label_fs); end
    xlabel(ax, 'r / R', 'FontSize', st.label_fs);
end

fig_tag = 'Fig_Eigdens_Theory_Check';
rows = arrayfun(@(k) sprintf('| %s | %.4g | %.4f | %.4f | %.1e | %.1e | %.3f |', C(k).name, C(k).R, C(k).max_diff, ...
    C(k).sampling_error, C(k).max_diff_solvers, C(k).max_diff_harris, C(k).edge_mean_max), 1:numel(C), ...
    'UniformOutput', false);
out = struct('figs', fig, 'files', {{}}, 'source', 'sampled (eigdens_check_theory)', 'check', C);
if cfg.save
    if cfg.images; save_figure_stable(out_dir, fig_tag, fig); end
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        fprintf(fid, '# %s\n\n', fig_tag);
        fprintf(fid, ['The radial-density theory against sampled random matrices (`eigdens_check_theory`; %d matrices ' ...
            'of size %d per case, Gaussian entries, no mean, seed %g). Each panel: the fraction of eigenvalues within ' ...
            'radius r of the origin, r in units of the predicted bulk radius R. Thick blue: sampled. Dashed black: ' ...
            'theory (Ahmadian et al. 2015, Eqs. 2.13 to 2.16, solved for block variances by ' ...
            '`eigdens_radial_theory_blocks`). Dotted gray: a uniform disk, (r/R)^2. A curve above the dotted one is a ' ...
            'spectrum concentrated toward the center.\n\n'], cfg.n_draws, cfg.n, cfg.seed);
        fprintf(fid, ['Cases: `uniform` every gain 1; `broad` log-normal gains with coefficient of variation 1; ' ...
            '`zeros` half the gains 0 and half 1 (the curve starts at 0.5: the zero-gain columns'' eigenvalues are ' ...
            'exactly at the origin); `two types` the E and I column variances of Harris et al. (2023) Fig. 2(c) at ' ...
            'alpha = 0.3; `blocks` variances by (post, pre) type [0.10 0.046; 0.046 0.03] with broad gains, a fifth of ' ...
            'them zero (the non-separable case).\n\n']);
        fprintf(fid, ['| case | R | largest difference, sampled vs theory | sampling error 1/sqrt(n draws) | column solver ' ...
            'vs block solver | vs Harris closed form (density / uniform) | mean largest modulus / R |\n' ...
            '|---|---|---|---|---|---|---|\n']);
        fprintf(fid, '%s\n', rows{:});
        fprintf(fid, '\nNaN: not applicable.\n');
        fclose(fid);
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
end

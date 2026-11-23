function out = fig_supp7_grouped(cfg)
% FIG_SUPP7_GROUPED Supplementary Figure S7: the reined-in mode once the adaptive variables respond, and the same test for other structures.
%
%   out = FIG_SUPP7_GROUPED('run_dir', d)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md). Reads the
% opt-in stages outlier_spectra and outlier_structures of ONE run
% directory; simulates nothing; errors (resolve_data_file, in the source
% figures) when the run lacks one of them.
%
%   (A)  fig_outlier_slow_variables: the full Jacobian near zero, the
%        frozen-gain mode and its continuation, networks by conditions.
%        One seed (the stage's first; cfg.seed picks another).
%   (B)  fig_ostruct_summary: for every structure of the outlier_structures
%        stage, the detachment of the intended mode (linear axis, the bulk
%        shaded), the responsive shares of E and I, and the mode's rate
%        with the slow variables frozen and free.
%
% The panels are the source figures' own drawings (composite_figure).
% Fig_Supp7_Grouped_table.md: a caption paragraph and their sidecars.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions, networks (the rows of A), seed, fig_tag.
%
% See also: fig_outlier_slow_variables, fig_ostruct_summary, composite_figure

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
    cfg.networks                  = {'reference', 'EEx1p5'}
    cfg.seed double               = []
    cfg.fig_tag     (1,:) char    = 'Fig_Supp7_Grouped'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
% A source figure that is already drawn when a later one fails is closed,
% not left behind.
before = findall(groot, 'Type', 'figure');
try
    out = build(cfg, out_dir);
catch err
    delete(setdiff(findall(groot, 'Type', 'figure'), before));
    rethrow(err);
end
end

function out = build(cfg, out_dir)

[f_slow, t_slow, o_slow] = composite_part(@fig_outlier_slow_variables, cfg, 'networks', cfg.networks, 'seed', cfg.seed);
[f_sum, t_sum, o_sum] = composite_part(@fig_ostruct_summary, cfg);
[fig, B] = composite_figure({struct('fig', f_slow, 'width', []), struct('fig', f_sum, 'width', [])}, 'visible', cfg.visible);
composite_letter(fig, B{1}(1), 'A');
composite_letter(fig, B{2}(1), 'B');

caption = sprintf(['Supplementary Figure S7. (A) The reined-in mode once the adaptive variables respond: the full ' ...
    'Jacobian''s eigenvalues near zero as rates, rows the networks %s, columns the conditions; the frozen-gain mode ' ...
    '(open markers) and the eigenvalue of the full Jacobian it continues into (filled), one seed. (B) The same test ' ...
    'for other structures: every network of the stage outlier_structures, in three groups separated by dotted lines ' ...
    '(low-rank structure added to the weights; local outliers of the sparse matrix and the reference network; ' ...
    'outliers of the block means). Top: detachment of the intended mode, its real part over the bulk radius, on a ' ...
    'linear axis with the bulk (-1 to 1) shaded; open diamonds: the same with every gain equal. Middle: responsive ' ...
    'share of the E (filled) and I (open) neurons. Bottom: the mode''s rate with the slow variables frozen (open) ' ...
    'and free (filled). Conditions in their colors.'], strjoin(cfg.networks, ', '));
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S7: the reined-in mode with the adaptive variables free; other structures', ...
    caption, {'Panel A (fig_outlier_slow_variables)', t_slow; 'Panel B (fig_ostruct_summary)', t_sum}, {o_slow.source, o_sum.source});
end

function out = fig_main5_mechanism(cfg)
% FIG_MAIN5_MECHANISM Main Figure 5: what adaptation does to the effective connectivity, and how that holds as the E-to-E imbalance grows.
%
%   out = FIG_MAIN5_MECHANISM('run_dir', d)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md), registry name
% fig_main5_grouped. Reads the opt-in stage outlier_spectra of ONE run
% directory; simulates nothing; errors (through resolve_data_file, in the
% source figures) when the run does not hold the stage.
%
%   (A)  fig_outlier_plane: the eigenvalues of the frozen-gain effective
%        connectivity W D as rates, reference network (top row) and the
%        E-to-E mean raised by half (bottom row), one column per
%        condition; the theoretical bulk radius; the equal-gain position of
%        the mode as a diamond, beyond a break in the real axis when it is
%        far out; in each panel the neurons on their response curve, E
%        above the curve and I below it, marked E and I. ONE SEED is drawn:
%        the first seed of the stage's index (cfg.seed picks another); the
%        sidecar names it.
%   (B)  fig_outlier_dose, left: detachment of the rightmost frozen-gain
%        eigenvalue over the bulk radius against the structural outlier,
%        with the equal-gain diamonds.
%   (C)  the same figure, middle: the responsive share of the E neurons.
%   (D)  right: mean E gain over mean I gain.
%
% Both source figures are drawn by their own functions and copied here
% (composite_figure), so the stand-alone appendix versions and this figure
% cannot drift apart. Fig_Main5_Grouped_table.md: a caption paragraph with
% the seed drawn, then the two source figures' own sidecars.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions (select_conditions), networks (the rows of A), seed, fig_tag.
%
% See also: fig_outlier_plane, fig_outlier_dose, composite_figure,
%           fig_main6_amplification

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the stage's files
    cfg.conditions                = {}    % names to draw, in order; empty = all the stage ran
    cfg.networks                  = {'reference', 'EEx1p5'}   % the rows of (A)
    cfg.seed double               = []    % the seed drawn in (A); [] = the first of the stage's index
    cfg.fig_tag     (1,:) char    = 'Fig_Main5_Grouped'
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

[f_plane, t_plane, o_plane] = composite_part(@fig_outlier_plane, cfg, 'networks', cfg.networks, 'seed', cfg.seed);
[f_dose, t_dose] = composite_part(@fig_outlier_dose, cfg);
width = f_plane.Position(3);                         % inches: fig_outlier_plane lays itself out in inches
[fig, B] = composite_figure({struct('fig', f_plane, 'width', []), struct('fig', f_dose, 'width', width)}, ...
    'visible', cfg.visible, 'margin', [0.15 0.15 -0.2 0.15]);
% (A) at the height of the condition titles, beside the first panel: the
% plane figure keeps half an inch for its titles above the panels, and a
% letter in the strip above that sat far from the panel (polish list,
% 2026-11-22). The top margin is taken back by as much.
composite_letter(fig, B{1}(1), 'A', 'dy', -0.36);
dose_ax = B{2}(1).ax;
letters = 'BCD';
for k = 1:min(3, numel(dose_ax))
    composite_letter(fig, B{2}(1), letters(k), dose_ax(k), 'dy', -0.12);
end

caption = sprintf(['Main Figure 5. (A) What adaptation does to the effective connectivity: eigenvalues of W D (each ' ...
    'neuron''s outgoing weights times its gain, slow variables held) drawn as rates in 1/s, so the imaginary axis is ' ...
    'the stability boundary; rows: the networks %s; columns: the conditions in their colors. ONE SEED IS DRAWN in (A): ' ...
    'network seed [%d %d], the first of the stage''s index unless another was asked for. Solid gray circle: the ' ...
    'theoretical bulk; open diamond: where the mode would be with every gain equal, beyond a break in the real axis ' ...
    'when far out. Inset of each panel: the neurons on their response curve, excitatory above the curve (E) and ' ...
    'inhibitory below it (I). (B) to (D), every seed of the stage: against the structural outlier over the bulk ' ...
    'radius, (B) the detachment of the rightmost frozen-gain eigenvalue with the equal-gain diamonds, (C) the ' ...
    'responsive share of the excitatory neurons, (D) the mean excitatory gain over the mean inhibitory gain. The ' ...
    'sections below are the sidecars of the two source figures, with every number.'], ...
    strjoin(cfg.networks, ', '), o_plane.seed);
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Main Figure 5: mechanism', caption, ...
    {'Panel A (fig_outlier_plane)', t_plane; 'Panels B to D (fig_outlier_dose)', t_dose}, {o_plane.source});
out.seed = o_plane.seed;
end

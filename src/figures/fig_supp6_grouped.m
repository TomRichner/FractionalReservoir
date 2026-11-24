function out = fig_supp6_grouped(cfg)
% FIG_SUPP6_GROUPED Supplementary Figure S6: gains, the density of the effective connectivity, and the frequency response of one neuron.
%
%   out = FIG_SUPP6_GROUPED('run_dir', d)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md), with S8's two
% panels added as D and E (Plan 3, question 11, 2026-11-22). Reads the
% opt-in stages eigenvalue_density and frequency_response of ONE run
% directory; simulates nothing; errors (resolve_data_file, in the source
% figures) when the run lacks them.
%
%   (A)  the per-neuron gains by condition (fig_eigdens_gains, reference
%        network): cumulative distribution, and the participation ratio of
%        the squared gains against the share of neurons with nonzero gain.
%   (B)  the fraction of eigenvalues within r, measured (colored) and from
%        the random-matrix theory with the saved gains (dashed), one panel
%        per network of cfg.networks, titled by the network, one y axis for
%        them (the later panels' axis name and tick labels left out): the
%        reference network and the E-to-E mean raised by half.
%   (C)  the measures by condition (fig_eigdens_measures, reference
%        network): a box and whisker over the networks of each condition,
%        the reference values as open diamonds beside it (Tom's notes of
%        2026-11-22, S6 note 3).
%   (D)  the gain of one neuron, the excitatory neuron of median rate, with
%        its slow variables free, over the frozen gain, against frequency
%        (fig_freqresp_single_neuron 'columns' 1; S8 (A) until now).
%   (E)  its phase lead (S8 (B) until now). One key under D and E.
%
% LAYOUT ON THE 9 IN CANVAS: three rows, A and B side by side (four
% panels), C (four panels), D and E with their key. Before, A, B and C had
% a row each and the figure was 13.2 in tall.
%
% The panels are the source figures' own drawings (composite_figure).
% Fig_Supp6_Grouped_table.md: a caption paragraph and their sidecars.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions, networks (the panels of B; the first is the network of A and
% C), fig_tag.
%
% See also: fig_eigdens_gains, fig_eigdens_radial, fig_eigdens_measures,
%           fig_freqresp_single_neuron, fig_supp8_grouped, composite_figure

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
    cfg.networks    (1,:) cell    = {'reference', 'EEx1p5'}
    cfg.fig_tag     (1,:) char    = 'Fig_Supp6_Grouped'
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
st = manuscript_style();
net = cfg.networks{1};

% The common canvas (Plan 2, canvas standard), 9 in since Plan 3: every
% row is canvas_w less the margins; the blocks are drawn at their size in
% the composite (not scaled). Type by role at the end (canvas_type).
MARGIN = [0.15 0.15 0.1 0.15]; COL_GAP = 0.2;
row_w = st.canvas_w - MARGIN(1) - MARGIN(2);
n_b = numel(cfg.networks);
H_AB = 2.75; H_C = 2.85; H_DE = 3.45;                % inches: the rows' heights
% (B)'s first panel carries the y axis (about 0.85 in of it); the others
% are as wide as their axes need.
w_b1 = 2.3; w_b = 1.5;
w_a = row_w - w_b1 - (n_b - 1) * w_b - n_b * COL_GAP;
[f_gain, t_gain, o_gain] = composite_part(@fig_eigdens_gains, cfg, 'network', net, 'fig_size', [w_a, H_AB]);
radial = struct('fig', cell(1, n_b), 'width', []);
parts = {sprintf('Panel A (fig_eigdens_gains, %s)', net), t_gain};
for k = 1:n_b
    w = w_b; if k == 1; w = w_b1; end
    [radial(k).fig, t_rad] = composite_part(@fig_eigdens_radial, cfg, 'network', cfg.networks{k}, 'panels', {'cumulative'}, ...
        'fig_size', [w, H_AB]);
    % The one thing that tells the panels of (B) apart: the network.
    ax = findall(radial(k).fig, 'Type', 'axes');
    title(ax(1), network_title(cfg.networks{k}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
    if k > 1                                         % one y axis for (B): the first panel's
        ylabel(ax(1), ''); set(ax(1), 'YTickLabel', []);
        drawnow;
    end
    parts(end + 1, :) = {sprintf('Panel B, panel %d (fig_eigdens_radial, %s; the cumulative panel is drawn)', k, cfg.networks{k}), t_rad}; %#ok<AGROW>
end
[f_meas, t_meas] = composite_part(@fig_eigdens_measures, cfg, 'network', net, 'fig_size', [row_w, H_C]);
parts(end + 1, :) = {sprintf('Panel C (fig_eigdens_measures, %s)', net), t_meas};
[f_one, t_one] = composite_part(@fig_freqresp_single_neuron, cfg, 'columns', 1, 'key', 'south', 'fig_size', [row_w, H_DE]);
parts(end + 1, :) = {'Panels D and E (fig_freqresp_single_neuron, neuron 1; the table has all three)', t_one};

row1 = [struct('fig', f_gain, 'width', []), radial];
[fig, B] = composite_figure({row1, struct('fig', f_meas, 'width', []), struct('fig', f_one, 'width', [])}, ...
    'visible', cfg.visible, 'margin', MARGIN, 'col_gap', COL_GAP, 'row_gap', 0.2, 'width', st.canvas_w);
composite_letter(fig, B{1}(1), 'A');
composite_letter(fig, B{1}(2), 'B', B{1}(2).ax(1));
composite_letter(fig, B{2}(1), 'C');
de = B{3}(1).ax;                                     % gain, then phase
composite_letter(fig, B{3}(1), 'D', de(1));
composite_letter(fig, B{3}(1), 'E', de(2));
canvas_type(fig);

caption = sprintf(['Supplementary Figure S6. (A) Gains of the neurons by condition, network `%s`: left, cumulative ' ...
    'distribution of the gain d (slope of the response curve times available vesicle fraction) over every neuron, ' ...
    'calibration state and seed, the step at zero being the share pinned at the floor or ceiling; right, per state ' ...
    'and seed, the participation ratio of the squared gains against the share of neurons with nonzero gain (dotted: ' ...
    'equality). (B) Fraction of the eigenvalues of the frozen-gain effective connectivity within radius r, in units ' ...
    'of the theoretical bulk radius, one panel per network (%s): thick colored lines the median over seeds and ' ...
    'states with the interquartile range shaded, thin dashed black the random-matrix theory with the saved gains, ' ...
    'dotted gray a uniform disk. (C) Measures by condition, network `%s`, over its networks (seeds; each the mean ' ...
    'over its calibration states): box the interquartile range, bar the median, whiskers to the furthest network ' ...
    'within 1.5 interquartile ranges, a dot a network beyond; open diamonds beside the box, one per network, the ' ...
    'reference value: rim occupancy relative to a uniform disk (diamonds: theory), eigenvalues of the frozen-gain ' ...
    'dendritic block in the band just below the boundary (diamonds: uniform disk) and above it (diamonds: the full ' ...
    'Jacobian), and the fast modes of the full Jacobian in that band (no reference). (D) One neuron''s frequency ' ...
    'response, the excitatory neuron of median rate in the first network: the gain from its dendritic input to its ' ...
    'synaptic output with its slow variables free, relative to the frozen gain, against frequency. (E) Its phase ' ...
    'lead. In D and E, lines as in the key: solid, spike-frequency adaptation and depression both free (median ' ...
    'over seeds, quartiles as a band); dashed, spike-frequency adaptation alone; dotted, depression alone; the ' ...
    'inhibitory neuron of median rate and a hard-driven excitatory neuron are in the table below. Conditions in ' ...
    'their colors.'], net, strjoin(cfg.networks, ', '), net);
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S6: gains and the density of the effective connectivity', ...
    caption, parts, {o_gain.source});
end

%% ------------------------------------------------------------------------
function s = network_title(name)
% 'reference' -> 'Reference'; 'EEx1p5' -> 'E-to-E × 1.5'; anything else as it is.
if strcmp(name, 'reference')
    s = 'Reference';
    return
end
tok = regexp(name, '^([EI])([EI])x(\d+)p(\d+)$', 'tokens', 'once');
if isempty(tok)
    s = name;
else
    s = sprintf('%s-to-%s %s %s.%s', tok{1}, tok{2}, char(215), tok{3}, tok{4});   % the sign itself, as Figure 5, S3 and S7 A write it (TeX's times set it with wide spaces)
end
end

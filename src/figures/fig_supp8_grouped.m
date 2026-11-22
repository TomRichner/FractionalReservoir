function out = fig_supp8_grouped(cfg)
% FIG_SUPP8_GROUPED Supplementary Figure S8: the frequency response of one neuron, gain (A) and phase lead (B) at three operating points.
%
%   out = FIG_SUPP8_GROUPED('run_dir', d)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md: "A: gain against
% frequency, adaptation alone, depression alone, both. B: phase lead").
% Reads the opt-in stage frequency_response of ONE run directory; simulates
% nothing; errors (resolve_data_file, in the source figure) when the run
% lacks it.
%
% Both panels are fig_freqresp_single_neuron's own drawing
% (composite_figure): three neurons in columns, titled, with one key at the
% right for the conditions' colors and the three line styles (adaptation
% and depression both free, solid; adaptation alone, dashed; depression
% alone, dotted).
%   (A)  its top row: gain over the frozen gain against frequency
%   (B)  its bottom row: phase lead against frequency
%
% UNTIL 2026-11-22 (A) was that whole figure, gain and phase, and (B)
% fig_freqresp_phase, the phase lead once more at a matched operating
% point, at half the width: the phase was shown twice and the gain shared
% a letter with it. Now each is drawn once. fig_freqresp_phase (the
% matched operating point, with the dendrite's phase lag in gray) stays a
% figure of its own in the appendix set and is no part of S8.
%
% Fig_Supp8_Grouped_table.md: a caption paragraph and the source's sidecar.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions, fig_tag.
%
% See also: fig_freqresp_single_neuron, fig_freqresp_phase, composite_figure

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
    cfg.fig_tag     (1,:) char    = 'Fig_Supp8_Grouped'
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

[f_one, t_one, o_one] = composite_part(@fig_freqresp_single_neuron, cfg);
[fig, B] = composite_figure({struct('fig', f_one, 'width', [])}, 'visible', cfg.visible, 'margin', [0.45 0.15 0.1 0.15]);
blk = B{1}(1);
assert(numel(blk.ax) == 6, 'fig_supp8_grouped:Panels', ...
    'Expected the 2 x 3 panels of fig_freqresp_single_neuron, found %d axes.', numel(blk.ax));
% (A) in the strip above the block; (B) at the same x, level with the top
% of the phase row (the first axes of the second row), left of its y label.
hA = composite_letter(fig, blk, 'A', [], 'dx', -0.4);
u = blk.ax(4).Units; blk.ax(4).Units = 'inches';
p = blk.ax(4).Position;
blk.ax(4).Units = u;
hB = composite_letter(fig, blk, 'B', [], 'dx', -0.4, 'dy', (p(2) + p(4)) - blk.top - 0.12);
out_letters = [hA hB]; %#ok<NASGU>

caption = ['Supplementary Figure S8, one neuron''s frequency response: the gain from its dendritic input to its ' ...
    'synaptic output with its slow variables free, relative to the frozen gain. (A) Its magnitude against ' ...
    'frequency. (B) Its phase lead. Columns, left to right, as titled: the excitatory neuron of median rate in ' ...
    'the first network of the stage, the inhibitory neuron of median rate there, and a hard-driven excitatory ' ...
    'neuron (90th percentile of rate) in the last network. Colors: the conditions with slow variables. Lines, as ' ...
    'in the key: solid, adaptation and depression both free (quartiles over seeds as a band); dashed, adaptation ' ...
    'alone; dotted, depression alone. The phase lead at a matched operating point, with the dendrite''s phase lag, ' ...
    'is fig_freqresp_phase, a figure of its own.'];
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S8: frequency response of one neuron', ...
    caption, {'Panels A and B (fig_freqresp_single_neuron)', t_one}, {o_one.source});
end

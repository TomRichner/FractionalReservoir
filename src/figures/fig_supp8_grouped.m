function out = fig_supp8_grouped(cfg)
% FIG_SUPP8_GROUPED Supplementary Figure S8: the frequency response of one neuron, and the phase lead at a matched operating point.
%
%   out = FIG_SUPP8_GROUPED('run_dir', d)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md). Reads the
% opt-in stage frequency_response of ONE run directory; simulates nothing;
% errors (resolve_data_file, in the source figures) when the run lacks it.
%
%   (A)  fig_freqresp_single_neuron: gain over the frozen gain (top) and
%        phase lead (bottom) against frequency, for three neurons
%        (columns), with adaptation and depression both free (solid),
%        adaptation alone (dashed) and depression alone (dotted).
%   (B)  fig_freqresp_phase: the phase lead at a matched operating point,
%        the same three line styles; gray: the dendrite's phase lag.
%
% The panels are the source figures' own drawings (composite_figure), one
% above the other. Fig_Supp8_Grouped_table.md: a caption paragraph and their
% sidecars.
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
[f_ph, t_ph] = composite_part(@fig_freqresp_phase, cfg);
[fig, B] = composite_figure({struct('fig', f_one, 'width', []), struct('fig', f_ph, 'width', [])}, 'visible', cfg.visible);
composite_letter(fig, B{1}(1), 'A');
composite_letter(fig, B{2}(1), 'B');

caption = ['Supplementary Figure S8. (A) One neuron''s frequency response: the gain from its dendritic input to its ' ...
    'synaptic output with its slow variables free, over the frozen gain (top), and its phase lead (bottom). Columns, ' ...
    'left to right: the excitatory neuron of median rate in the first network of the stage, the inhibitory neuron of ' ...
    'median rate there, and a hard-driven excitatory neuron in the last network. Solid: adaptation and depression ' ...
    'both free (quartiles over seeds as a band); dashed: adaptation alone; dotted: depression alone. (B) The phase ' ...
    'lead at a matched operating point (every excitatory neuron put at one rate on the linear part of the response ' ...
    'curve, slow variables at their steady state), median over neurons and seeds, the same line styles; gray: the ' ...
    'phase lag of the dendrite, which a mode''s gain must cancel for the mode to be marginal at that frequency. ' ...
    'Conditions in their colors; a condition without slow variables has no response to draw.'];
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S8: frequency response of one neuron', ...
    caption, {'Panel A (fig_freqresp_single_neuron)', t_one; 'Panel B (fig_freqresp_phase)', t_ph}, {o_one.source});
end

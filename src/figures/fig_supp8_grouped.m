function out = fig_supp8_grouped(cfg)
% FIG_SUPP8_GROUPED Supplementary Figure S8: the frequency response of one neuron, gain (A) and phase lead (B).
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
% ONE NEURON, the excitatory neuron of median rate in the stage's first
% network (Tom's notes of 2026-11-22: separate columns for the E, the I
% and a hard-driven E neuron are not needed; Plan 3, question 9), median
% and quartiles over seeds. Both panels are fig_freqresp_single_neuron's
% own drawing with 'columns' 1 (composite_figure), side by side, the key
% under them: the three line styles (spike-frequency adaptation and
% depression together, solid; spike-frequency adaptation alone, dashed;
% depression alone, dotted) and the conditions' colors.
%   (A)  gain over the frozen gain against frequency
%   (B)  phase lead against frequency
% The other two neurons are in the table (the source's sidecar keeps all
% three).
%
% THE SAME TWO PANELS ARE PANELS D AND E OF S6 (Plan 3, question 11:
% S8's panels go into S6). This figure draws them on their own until the
% registry drops it.
%
% UNTIL 2026-11-22 (A) and (B) were three columns each (the three neurons)
% and the key said "adaptation" for spike-frequency adaptation.
%
% Fig_Supp8_Grouped_table.md: a caption paragraph and the source's sidecar.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions, fig_tag.
%
% See also: fig_freqresp_single_neuron, fig_supp6_grouped, composite_figure

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
% The common canvas (Plan 2, canvas standard): the source is laid out at
% canvas_w less the margins. Type by role at the end (canvas_type).
% The tiled layout keeps about half an inch free right of the last tick
% label; the right margin is negative by part of that, so that the drawing
% reaches the right edge of the canvas (the cropped PNG within canvas_tol).
st = manuscript_style(); MARGIN = [0.15 -0.3 0.1 0.15];
row_w = st.canvas_w - MARGIN(1) - MARGIN(2);
[f_one, t_one, o_one] = composite_part(@fig_freqresp_single_neuron, cfg, 'columns', 1, 'key', 'south', ...
    'fig_size', [row_w, 4.3]);
[fig, B] = composite_figure({struct('fig', f_one, 'width', [])}, 'visible', cfg.visible, 'margin', MARGIN, 'width', st.canvas_w);
blk = B{1}(1);
assert(numel(blk.ax) == 2, 'fig_supp8_grouped:Panels', ...
    'Expected the gain and the phase panel of fig_freqresp_single_neuron, found %d axes.', numel(blk.ax));
composite_letter(fig, blk, 'A', blk.ax(1));
composite_letter(fig, blk, 'B', blk.ax(2));
canvas_type(fig);

caption = ['Supplementary Figure S8, one neuron''s frequency response: the gain from its dendritic input to its ' ...
    'synaptic output with its slow variables free, relative to the frozen gain, for the excitatory neuron of ' ...
    'median rate in the first network of the stage. (A) Its magnitude against frequency. (B) Its phase lead. ' ...
    'Colors: the conditions with slow variables. Lines, as in the key: solid, spike-frequency adaptation and ' ...
    'depression both free (median over seeds, quartiles as a band); dashed, spike-frequency adaptation alone; ' ...
    'dotted, depression alone. The inhibitory neuron of median rate and a hard-driven excitatory neuron (90th ' ...
    'percentile of rate, last network) are in the table below. The phase lead at a matched operating point, with ' ...
    'the dendrite''s phase lag, is fig_freqresp_phase, a figure of its own.'];
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S8: frequency response of one neuron', ...
    caption, {'Panels A and B (fig_freqresp_single_neuron, neuron 1; the table has all three)', t_one}, {o_one.source});
end

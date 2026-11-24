function args = canvas_export_args()
%CANVAS_EXPORT_ARGS exportgraphics options for the 18 manuscript figures: the whole canvas, untrimmed.
%
%   exportgraphics(fig, file, canvas_export_args{:})
%   save_figure_stable(out_dir, tag, fig, canvas_export_args())
%
% WHY (2026-11-23). exportgraphics trims a PNG or SVG to its ink by default
% (Padding "tight"), so the 18 figures, every one drawn on the 9 in canvas
% (manuscript_style().canvas_w), came out 8.6 to 9.0 in wide, and the trim
% differed between two sessions on the same figure (S6 was 8.64 in in one
% and failed test_canvas_standard's lower limit). Padding "figure" (R2026a's
% help exportgraphics: Padding "tight" (default) | "figure" | number) saves
% the figure's own size, so every image is canvas_w wide and the same in
% every session. Used for PNG and SVG by composite_save, fig_grouped_main
% (Figures 2 and 3), fig_main1_concept_networks,
% fig_supp3_connectivity_gallery, fig_supp10_grouped and stim_fig8_compose.
%
% See also: save_figure_stable, save_some_figs_to_folder_2, manuscript_style
args = {'Padding', 'figure'};
end

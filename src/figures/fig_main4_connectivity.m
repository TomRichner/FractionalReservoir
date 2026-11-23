function out = fig_main4_connectivity(cfg)
% FIG_MAIN4_CONNECTIVITY Main Figure 4: sensitivity (A), rate (B) and what decays slowly (C).
%
%   out = FIG_MAIN4_CONNECTIVITY('run_dir', d)
%
% The paper's Figure 4 of the figure plan of 2026-11-21. Panels (A) and (B)
% are the grouped Figure 4 as it was (fig_grouped_main, group 4: the
% sensitivity medians over the 1-D sweeps and lambda_1 against mean rate
% over the joint sample), unchanged in size and style. Panel (C), new, is
% fig_leading_direction in its 'column' layout, placed to the right of
% them: for each condition, the dendritic share of the leading perturbation
% against lambda_1 over the joint sample. No line at the uncoupled rate
% (Tom, 2026-11-21); its numbers are in the table when the run holds the
% uncoupled_reference stage.
%
% WHY (C) IS A COLUMN ON THE RIGHT, NOT A ROW BELOW. A third row would make
% the canvas about 1040 px tall. A window taller than the workstation's
% screen (1073 px, 2026-11-21) is shrunk by MATLAB at a moment of its
% choosing, so the canvas keeps the height of the present Figure 4
% (807.5 px) and grows in width by the column. It fits the screen, and is
% exported by the ordinary path.
%
% HOW IT IS BUILT. With package D's composite helpers (composite_figure,
% composite_save): one row of two blocks, each a source function's own
% figure at its own size, no margins and no letter strip, because the
% blocks carry their letters inside. fig_grouped_main(4, ...) is called with
% 'save', false and returns its styled figure; its two block letters are
% annotations, which a copy of the axes does not take along, so they are
% redrawn at the same places. fig_leading_direction's column is made as tall
% as the canvas below the level of (A)'s letter, moved to the bottom of the
% row, and lettered (C) at that level. Nothing is simulated or recomputed here, and
% fig_grouped_main is not edited: whatever changes in panels (A) and (B)
% (condition colors, sentence-case titles) arrives here with it.
% The first version copied the axes with a helper of its own and saved with
% print; that was a second mechanism beside package D's and is gone.
%
% OPTIONS
%   conditions       passed to both inner figures (FD-006)
%   reference_file, require_reference   passed to fig_leading_direction
%                    (the uncoupled reference; see there)
%   reference_line   passed to fig_leading_direction: the dotted line at
%                    the uncoupled rate (default false)
%   fig_tag          the file name stem, default 'Fig_Main4_Connectivity'.
%                    The registry passes 'Fig_Main4_Grouped' so that the
%                    manuscript's figure_sources.json needs no repointing.
%
% Saves <tag>.png, .svg, .fig and <tag>_table.md (panel (C)'s numbers, and
% where the numbers of (A) and (B) are). out.figs, out.files, out.source,
% out.stats (fig_leading_direction's, per condition).
%
% See also: fig_grouped_main, style_main4_grouped, fig_leading_direction,
%           composite_figure, composite_save

arguments
    cfg.verbose        (1,:) char    = 'minimal'
    cfg.data_file      (1,:) char    = ''     % unused; accepted for a uniform call
    cfg.out_dir        (1,:) char    = ''
    cfg.save           (1,1) logical = true
    cfg.visible        (1,1) logical = false
    cfg.run_dir        (1,:) char    = ''
    cfg.preset_name    (1,:) char    = ''
    cfg.conditions                   = {}
    cfg.reference_file (1,:) char    = ''
    cfg.require_reference (1,1) logical = false
    cfg.reference_line (1,1) logical = false
    cfg.fig_tag        (1,:) char    = 'Fig_Main4_Connectivity'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
C_WIDTH = 340;     % px, the column of panel (C)

%% Panels (A) and (B): the grouped Figure 4, as it is
ab = fig_grouped_main(4, 'run_dir', cfg.run_dir, 'preset_name', cfg.preset_name, 'save', false, ...
    'visible', false, 'verbose', cfg.verbose, 'conditions', cfg.conditions);
drawnow;
ab_size = ab.figs.Position(3:4);
letter_a = findall(ab.figs, 'Tag', 'main4_block_A');
letter_b = findall(ab.figs, 'Tag', 'main4_block_B');
assert(isscalar(letter_a) && isscalar(letter_b), 'fig_main4_connectivity:Letters', ...
    'The grouped Figure 4 no longer carries its two block letters (main4_block_A, main4_block_B).');
a_box = letter_a.Position; b_box = letter_b.Position;      % normalized in the grouped figure
top = a_box(2) + a_box(4);                                   % the top of (A)'s letter

%% Panel (C): the column, as tall as the canvas below the level of (A)
c = fig_leading_direction('run_dir', cfg.run_dir, 'save', false, 'visible', false, 'verbose', cfg.verbose, ...
    'conditions', cfg.conditions, 'layout', 'column', 'reference_file', cfg.reference_file, ...
    'require_reference', cfg.require_reference, 'reference_line', cfg.reference_line);
c.figs.Position(3:4) = [C_WIDTH, round(top * ab_size(2))];
drawnow;
sources = [reshape(cellstr(ab.source), 1, []), reshape(c.source, 1, [])];
sources = sources(~cellfun(@isempty, sources));

%% One row, two blocks, at their own size
blocks = struct('fig', {ab.figs, c.figs});
[fig, B] = composite_figure({blocks}, 'visible', cfg.visible, 'margin', [0 0 0 0], 'strip', 0, 'col_gap', 0);
fp = fig.Position;                                           % inches
% composite_figure aligns the blocks of a row at their top edges; the column
% is the shorter block and belongs at the bottom, level with (B)'s axis names.
% (Giving the column's tiled layout an OuterPosition on a full-height figure
% was tried first: the tiles did not move.)
drop = B{1}(1).rect(4) - B{1}(2).rect(4);
for k = 1:numel(B{1}(2).ax)
    B{1}(2).ax(k).Position(2) = B{1}(2).ax(k).Position(2) - drop;
end
r = B{1}(1).rect;
letters = {'(A)', a_box, 'main4_block_A'; '(B)', b_box, 'main4_block_B'};
for k = 1:size(letters, 1)
    p = letters{k, 2};
    annotation(fig, 'textbox', [(r(1) + p(1) * r(3)) / fp(3), (r(2) + p(2) * r(4)) / fp(4), p(3) * r(3) / fp(3), p(4) * r(4) / fp(4)], ...
        'String', letters{k, 1}, 'LineStyle', 'none', 'FontSize', 14, 'Margin', 0, 'Tag', letters{k, 3});
end
rc = B{1}(2).rect;
annotation(fig, 'textbox', [(rc(1) + 0.02) / fp(3), (r(2) + a_box(2) * r(4)) / fp(4), a_box(3) * r(3) / fp(3), a_box(4) * r(4) / fp(4)], ...
    'String', '(C)', 'LineStyle', 'none', 'FontSize', 14, 'Margin', 0, 'Tag', 'main4_block_C');

caption = ['(A) Leading growth rate against each connectivity parameter, median and interquartile range per condition ' ...
    '(1-D sweeps). (B) Leading growth rate against mean firing rate over the joint sample, per condition, colored by ' ...
    'the E:I weight balance. (C) For each condition, the share of the leading perturbation that lies in the dendritic ' ...
    'coordinates against the leading growth rate, one point per configuration of the joint sample; dashed, zero ' ...
    'growth'];
if cfg.reference_line
    caption = [caption '; dotted, the growth rate of the same network with its recurrent weights scaled to near zero'];
end
caption = [caption '.'];
parts = {'Panels (A) and (B)', ['the grouped Figure 4 (fig_grouped_main, group 4). Their numbers are in the tables of ' ...
    'fig_sensitivity_medians and fig_lle_vs_rate (Fig_LLE_vs_Rate_table.md).']; ...
    'Panel (C)', c.table_text};
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Main Figure 4: connectivity', caption, parts, sources);
out.stats = c.stats;
end

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
% against lambda_1 over the joint sample, with the uncoupled reference.
%
% WHY (C) IS A COLUMN ON THE RIGHT, NOT A ROW BELOW. A third row would make
% the canvas about 1040 px tall. A figure taller than the workstation's
% screen is shrunk by MATLAB at creation or at export, not every time (seen
% on the four-row Figure 3, 2026-11-21), so the canvas keeps the height of
% the present Figure 4 (807.5 px) and grows in width by the column.
%
% HOW IT IS BUILT. fig_grouped_main(4, ...) is called with 'save', false and
% returns its styled figure; every axes, legend and colorbar is copied to
% the left part of a wider canvas at the same pixel size, and the two block
% letters, which are annotations, are redrawn at the same places. The
% column of fig_leading_direction is copied to the right part and lettered
% (C) at the height of (A). Nothing is simulated or recomputed here.
% fig_grouped_main is not edited, so whatever package A changes in panels
% (A) and (B) (condition colors, sentence-case titles) arrives here with it.
%
% OPTIONS
%   conditions       passed to both inner figures (FD-006)
%   reference_file, require_reference   passed to fig_leading_direction
%                    (the uncoupled reference; see there)
%   fig_tag          the file name stem, default 'Fig_Main4_Connectivity'.
%                    The registry passes 'Fig_Main4_Grouped' so that the
%                    manuscript's figure_sources.json needs no repointing.
%
% Saves <tag>.png, <tag>.fig and <tag>_table.md (panel (C)'s numbers, and
% where the numbers of (A) and (B) are). out.figs, out.files, out.source,
% out.table_text.
%
% See also: fig_grouped_main, style_main4_grouped, fig_leading_direction,
%           fig_sensitivity_medians, fig_lle_vs_rate

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

%% Panel (C)
c = fig_leading_direction('run_dir', cfg.run_dir, 'save', false, 'visible', false, 'verbose', cfg.verbose, ...
    'conditions', cfg.conditions, 'layout', 'column', 'reference_file', cfg.reference_file, ...
    'require_reference', cfg.require_reference);
% The column takes the height of the canvas, less the room (A)'s legend
% takes above its panels, so that (C) starts level with (A).
letter_a = findall(ab.figs, 'Tag', 'main4_block_A');
letter_b = findall(ab.figs, 'Tag', 'main4_block_B');
assert(isscalar(letter_a) && isscalar(letter_b), 'fig_main4_connectivity:Letters', ...
    'The grouped Figure 4 no longer carries its two block letters (main4_block_A, main4_block_B).');
a_box = letter_a.Position; b_box = letter_b.Position;
top = a_box(2) + a_box(4);                       % normalized: the top of (A)'s letter
c.figs.Position(3:4) = [C_WIDTH, round(top * ab_size(2))];
drawnow;

%% The wider canvas
W = ab_size(1) + C_WIDTH; H = ab_size(2);
fig = figure('Color', 'w', 'Position', [40 40 W H]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
copy_into(ab.figs, fig, [0, 0, ab_size(1) / W, 1]);
copy_into(c.figs, fig, [ab_size(1) / W, 0, C_WIDTH / W, top]);
letters = {'(A)', a_box, 'main4_block_A'; '(B)', b_box, 'main4_block_B'};
for k = 1:size(letters, 1)
    p = letters{k, 2}; p([1 3]) = p([1 3]) * ab_size(1) / W;
    annotation(fig, 'textbox', p, 'String', letters{k, 1}, 'LineStyle', 'none', 'FontSize', 14, 'Margin', 0, ...
        'Tag', letters{k, 3});
end
annotation(fig, 'textbox', [ab_size(1) / W + 0.004, a_box(2), a_box(3) * ab_size(1) / W, a_box(4)], 'String', '(C)', ...
    'LineStyle', 'none', 'FontSize', 14, 'Margin', 0, 'Tag', 'main4_block_C');
sources = [reshape(cellstr(ab.source), 1, []), reshape(c.source, 1, [])];
close(ab.figs); close(c.figs);

table_text = sprintf(['# Main Figure 4: connectivity\n\nRun: `%s`\n\nPanels (A) and (B) are the grouped Figure 4 ' ...
    '(fig_grouped_main, group 4): their numbers are in the tables of fig_sensitivity_medians and fig_lle_vs_rate ' ...
    '(Fig_LLE_vs_Rate_table.md). Panel (C) is fig_leading_direction; its table follows.\n\n%s\n'], ...
    cfg.run_dir, regexprep(c.table_text, '^# ', '## Panel (C): ', 'once'));

out = struct('figs', fig, 'files', {{}}, 'source', {sources}, 'table_text', table_text, 'stats', c.stats);
if cfg.save
    for ext = {'.png', '.fig', '.svg'}
        old = fullfile(out_dir, [cfg.fig_tag ext{1}]);
        if isfile(old); delete(old); end
    end
    % print at the canvas size, not exportgraphics: exportgraphics was seen
    % to resize the window during the export (fig_grouped_main's note of
    % 2026-11-12, and again on 2026-11-21), which print does not depend on.
    drawnow;
    px = fig.Position(3:4);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 px / 96], 'PaperSize', px / 96);
    print(fig, fullfile(out_dir, [cfg.fig_tag '.png']), '-dpng', '-r200');
    fig.Position(3:4) = px;
    savefig(fig, fullfile(out_dir, [cfg.fig_tag '.fig']));
    out.files = existing_outputs(out_dir, cfg.fig_tag);
    fid = fopen(fullfile(out_dir, [cfg.fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n', table_text);
        fclose(fid);
    end
end
end

function copy_into(src, dst, region)
% Copy every axes, legend and colorbar of src into the normalized REGION of
% dst, each at its place within src. One copyobj call, so a legend or
% colorbar keeps its axes.
drawnow;
h = findall(src, '-isa', 'matlab.graphics.axis.Axes', '-or', 'Type', 'legend', '-or', 'Type', 'colorbar');
if isempty(h); return; end
rect = zeros(numel(h), 4);
for k = 1:numel(h); rect(k, :) = getpixelposition(h(k), true); end
sz = getpixelposition(src); rect = rect ./ [sz(3:4) sz(3:4)];
new = copyobj(h, dst);
for k = 1:numel(new)
    new(k).Units = 'normalized';
    new(k).Position = [region(1:2) + rect(k, 1:2) .* region(3:4), rect(k, 3:4) .* region(3:4)];
    if isa(new(k), 'matlab.graphics.axis.Axes'); new(k).PositionConstraint = 'innerposition'; end
end
end

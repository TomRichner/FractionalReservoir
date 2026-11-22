function out = fig_supp7_grouped(cfg)
% FIG_SUPP7_GROUPED Supplementary Figure S7: the reined-in mode once the adaptive variables respond, and the same test for the structures of the gallery S3.
%
%   out = FIG_SUPP7_GROUPED('run_dir', d)
%
% The figure plan of 2026-11-21 (manuscript repository,
% docs/Advancing_the_narrative/figures_main_supp_appen.md). Reads the
% opt-in stages outlier_spectra and outlier_structures of ONE run
% directory, and the gallery of its concept_examples stage when it has
% one; simulates nothing; errors (resolve_data_file, in the source
% figures) when the run lacks one of the first two.
%
%   (A)  fig_outlier_slow_variables: the full Jacobian near zero, the
%        frozen-gain mode and its continuation, networks by conditions.
%        One seed (the stage's first; cfg.seed picks another). The
%        eigenvalues its dendritic-share filter leaves out are drawn
%        faintly in gray (see that function).
%   (B)  fig_ostruct_summary for THE STRUCTURES OF THE GALLERY S3, in the
%        gallery's order and under its tile titles (cx_gallery_names): the
%        detachment of the intended mode (linear axis, the bulk shaded),
%        the responsive shares of E and I, and the mode's rate with the
%        slow variables frozen and free. So panel B summarizes exactly what
%        S3 shows. A run without a gallery (no concept_examples stage)
%        draws every network of the outlier_structures stage under the
%        stage's own names, as this panel did until 2026-11-22; the
%        networks of that stage which the gallery does not show stay in
%        the stand-alone fig_ostruct_summary (appendix). A gallery
%        structure the outlier_structures stage did not run is left out of
%        B and named in the table.
%
% LAYOUT. A and B side by side, each at a width of its own, so that the
% figure stays under 0.9 of the workstation screen in both dimensions and
% is exported cropped as usual (until 2026-11-22 they were stacked at their
% own sizes: 19 by 20 inches, a third of A's row empty, and the text small
% at any printed width). cfg.width_a and cfg.width_b are those widths.
%
% The panels are the source figures' own drawings (composite_figure).
% Fig_Supp7_Grouped_table.md: a caption paragraph and their sidecars.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions, networks (the rows of A), seed, fig_tag, width_a, width_b,
% gallery (true: B draws the gallery's structures when the run has them).
%
% See also: fig_outlier_slow_variables, fig_ostruct_summary, composite_figure,
%           cx_gallery_names, fig_supp3_connectivity_gallery

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
    cfg.width_a     (1,1) double  = 8.45  % inches
    cfg.width_b     (1,1) double  = 6.85
    cfg.gallery     (1,1) logical = true
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
T = composite_tools();

% Smaller marks than the stand-alone figure's 7 pt: at this width they
% would reach the tick labels beside the origin.
[f_slow, t_slow, o_slow] = composite_part(@fig_outlier_slow_variables, cfg, 'networks', cfg.networks, 'seed', cfg.seed, ...
    'mark_size', 5.5);
drawnow;
px = getpixelposition(f_slow);
h_a = px(4) / px(3) * cfg.width_a;                  % the height A has at its width

% (B): the gallery's structures, when the run has a gallery.
[keys, labels, missing, gallery_file] = gallery_networks(cfg);
b_args = {'fig_size', [cfg.width_b, min(max(h_a, 7.0), 8.6)], 'compact', true};
if ~isempty(keys); b_args = [b_args, {'networks', keys, 'labels', labels}]; end
[f_sum, t_sum, o_sum] = composite_part(@fig_ostruct_summary, cfg, b_args{:});

[fig, B] = composite_figure({[struct('fig', f_slow, 'width', cfg.width_a), struct('fig', f_sum, 'width', cfg.width_b)]}, ...
    'visible', cfg.visible, 'strip', 0.6);
% (A): the conditions over its columns on two lines (at this width the
% full names on one line run into each other).
top = B{1}(1).ax;
top = top(arrayfun(@(a) ~isempty(strtrim(strjoin(cellstr(a.Title.String), ' '))), top));
for k = 1:numel(top)
    try
        name = grouped_condition_label({strjoin(cellstr(top(k).Title.String), ' ')}, 'name_of_title');
        top(k).Title.String = T.condition_title(name{1}, true);
    catch err
        if ~any(strcmp(err.identifier, {'grouped_condition_label:UnknownTitle', 'grouped_condition_label:AmbiguousTitle'}))
            rethrow(err);
        end
    end
end
composite_letter(fig, B{1}(1), 'A', [], 'dy', 0.3);
composite_letter(fig, B{1}(2), 'B', [], 'dy', 0.3);

if isempty(keys)
    what_b = ['every network of the stage outlier_structures, in three groups separated by dotted lines (low-rank ' ...
        'structure added to the weights; local outliers of the sparse matrix and the reference network; outliers ' ...
        'of the block means)'];
else
    what_b = sprintf(['the structures of the gallery S3, in its order and under its titles (%s); dotted lines where the ' ...
        'kind of structure changes (the balanced reference network; outliers of the block means; a local outlier of ' ...
        'the sparse matrix; low-rank structure added to the weights)'], strjoin(labels, ', '));
end
caption = sprintf(['Supplementary Figure S7. (A) The reined-in mode once the adaptive variables respond: the full ' ...
    'Jacobian''s eigenvalues near zero as rates, rows the networks %s, columns the conditions; the frozen-gain mode ' ...
    '(open markers) and the eigenvalue of the full Jacobian it continues into (filled), one seed; faint gray: the ' ...
    'eigenvalues with little weight in the dendritic variables (slow modes of single neurons). (B) The same test ' ...
    'for other structures: %s. Top: detachment of the intended mode, its real part over the bulk radius, on a ' ...
    'linear axis with the bulk (-1 to 1) shaded; open diamonds: the same with every gain equal. Middle: responsive ' ...
    'share of the E (filled) and I (open) neurons. Bottom: the mode''s rate with the slow variables frozen (open) ' ...
    'and free (filled). Conditions in their colors. Figure %.2f in wide (A %.2f, B %.2f).'], ...
    strjoin(cfg.networks, ', '), what_b, fig.Position(3), cfg.width_a, cfg.width_b);
note = '';
if ~isempty(gallery_file)
    note = sprintf('The structures of panel B are those of the gallery in `%s`, keys in the outlier_structures stage: %s.', ...
        gallery_file, strjoin(keys, ', '));
    if ~isempty(missing)
        note = sprintf('%s NOT IN PANEL B, because the outlier_structures stage did not run them: %s.', note, strjoin(missing, ', '));
    end
    note = [note newline newline];
end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S7: the reined-in mode with the adaptive variables free; other structures', ...
    caption, {'Panel A (fig_outlier_slow_variables)', t_slow; 'Panel B (fig_ostruct_summary)', [note t_sum]}, {o_slow.source, o_sum.source});
out.gallery_keys = keys;
out.gallery_labels = labels;
out.gallery_missing = missing;
end

%% ------------------------------------------------------------------------
function [keys, labels, missing, file] = gallery_networks(cfg)
% The gallery's structures as keys of the outlier_structures stage, with
% S3's titles; the titles of those that stage did not run in missing.
% Everything empty when B is to draw the whole stage (cfg.gallery false, or
% a run without a gallery).
keys = {}; labels = {}; missing = {}; file = '';
f = fullfile(cfg.run_dir, 'concept_examples', 'concept_examples.mat');
if ~cfg.gallery || ~isfile(f); return; end
C = cx_load(cfg.run_dir, '', 'fig_supp7_grouped');
if isempty(C.gallery); return; end
[titles, okeys] = cx_gallery_names(C.gallery, C.settings.seed);
have = ostruct_stage_keys(cfg);
in_stage = ismember(okeys, have);
if ~any(in_stage)
    error('fig_supp7_grouped:GalleryNotTested', ['None of the gallery''s structures (%s) is a network of the ' ...
        'outlier_structures stage of %s (%s). Pass ''gallery'', false to draw that stage''s own list.'], ...
        strjoin(okeys, ', '), cfg.run_dir, strjoin(have, ', '));
end
keys = okeys(in_stage); labels = titles(in_stage); missing = titles(~in_stage); file = C.file;
end

function keys = ostruct_stage_keys(cfg)
F = ostruct_figure_data(struct('run_dir', cfg.run_dir, 'data_file', '', 'conditions', {cfg.conditions}, ...
    'verbose', cfg.verbose), 'fig_supp7_grouped');
keys = F.keys(:)';
end

function out = fig_supp7_grouped(cfg)
% FIG_SUPP7_GROUPED Supplementary Figure S7: the reined-in mode once the adaptive variables respond, the same test for the structures of the gallery S3, and the dwell time beyond zero of the reference network.
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
%   (C)  fig_outlier_dwell for THE REFERENCE NETWORK (Tom, 2026-11-23,
%        FD-057: "I want the dwell-time figure promoted to supplement as
%        part of S7"): how long the rightmost eigenvalue stays beyond zero,
%        censored survival curves, median across seeds, dashed thin with the
%        gains frozen, solid thick for the full Jacobian; a key of the two
%        line styles in gray to the right of the axes (the condition colors
%        are A's titles). It reads cfg.dwell_run_dir when given (a variant
%        drawing C from a longer run while A and B read run_dir), else
%        run_dir. The stand-alone fig_outlier_dwell stays, and the spectral
%        supplement still draws it as its own figure.
%
% LAYOUT WITH PANEL C (2026-11-23). To fit 11.8 in (canvas_max_h) with C,
% A's planes are cfg.panel_w_a wide (1.95 in, centered in their columns,
% which keep the condition titles' room; they filled the columns, 2.6 in,
% before) with cfg.row_gap_a between the rows; B is cfg.height_b tall (4.7
% in, from 5.7) with its names turned 22 degrees (from 28), which gives the
% three rows back most of the height their two-line names need; C is
% cfg.height_c (2.4 in) tall with axes cfg.axes_w_c (4.3 in, about half
% the canvas: one panel of survival curves needs no more), its key in the
% other half. C's source keeps its title strip, where its letter goes, so
% its row has no letter strip of its own. The figure is 11.7 in tall on the
% fast run (A's height follows its data: the planes are 1:1). Type is
% manuscript_style's throughout; nothing was shrunk.
%
% LAYOUT (Plan 3, the 9 in canvas, 2026-11-22). The figure is canvas_w wide
% and both blocks are DRAWN at the width of their place (A through its
% fig_w, B through its fig_size), so nothing is scaled when they are
% copied: A on top across the width, so that its three columns hold the
% condition titles side by side; B under it across the width,
% cfg.height_b tall, its structure names turned cfg.label_rotation degrees
% (less than the stand-alone figure's 35, so that the three rows keep the
% height their two-line axis names need). Under Plan 2 the layout made for
% 13 in was shrunk onto 9 in (A scaled to 0.88 of the row): the condition
% titles ran into each other and the frozen modes of the reference network
% sat on the tick label "5". Before that A and B stood side by side
% (15.4 in wide).
%
% The panels are the source figures' own drawings (composite_figure).
% Fig_Supp7_Grouped_table.md: a caption paragraph and their sidecars.
%
% Options: run_dir, out_dir, save, visible, verbose, preset_name (unused),
% conditions, networks (the rows of A), seed, fig_tag, height_b,
% label_rotation, gallery (true: B draws the gallery's structures when the
% run has them), panel_w_a, row_gap_a, height_c, axes_w_c, dwell_run_dir.
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
    cfg.height_b    (1,1) double  = 4.7   % inches, B with its structure names (5.7 before panel C)
    cfg.label_rotation (1,1) double = 22  % degrees, the structure names under B (28 before panel C)
    cfg.gallery     (1,1) logical = true
    cfg.panel_w_a   (1,1) double  = 1.95  % inches, each plane of A (they filled their columns, 2.6 in, before panel C)
    cfg.row_gap_a   (1,1) double  = 0.18  % inches between A's two rows (0.25 before)
    cfg.height_c    (1,1) double  = 2.4   % inches, panel C
    cfg.axes_w_c    (1,1) double  = 4.3   % inches, C's axes (about half the canvas; its key to the right)
    cfg.dwell_run_dir (1,:) char  = ''    % the run directory C reads; '' = run_dir (a variant: C from a longer run)
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

% The common canvas (manuscript_style canvas_w). Both source figures are
% DRAWN at the width of their block, so nothing is scaled when they are
% copied (as fig_main5_mechanism). A's source keeps room for its titles
% above the panels, where the letter goes: the negative top margin takes
% the letters' strip back.
st = manuscript_style(); MARGIN = [0.15 0.1 -0.31 0.05];
width = st.canvas_w - MARGIN(1) - MARGIN(2);
width_a = width; width_b = width;

% Smaller marks than the stand-alone figure's 7 pt: at this width they
% would reach the tick labels beside the origin.
[f_slow, t_slow, o_slow] = composite_part(@fig_outlier_slow_variables, cfg, 'networks', cfg.networks, 'seed', cfg.seed, ...
    'mark_size', 5.5, 'fig_w', width_a, 'panel_w', cfg.panel_w_a, 'row_gap', cfg.row_gap_a);

% (B): the gallery's structures, when the run has a gallery; drawn at its
% place's size (its three rows share the structures' axis).
[keys, labels, missing, gallery_file] = gallery_networks(cfg);
b_args = {'fig_size', [width_b, cfg.height_b], 'label_rotation', cfg.label_rotation};
if ~isempty(keys); b_args = [b_args, {'networks', keys, 'labels', labels}]; end
[f_sum, t_sum, o_sum] = composite_part(@fig_ostruct_summary, cfg, b_args{:});

% (C): the dwell time of the reference network (Tom, 2026-11-23, FD-057),
% one panel about half the canvas wide, its key of the line styles to the
% right; read from cfg.dwell_run_dir when given (a longer run), else run_dir.
cfg_c = cfg;
if ~isempty(cfg.dwell_run_dir); cfg_c.run_dir = cfg.dwell_run_dir; end
[f_dw, t_dw, o_dw] = composite_part(@fig_outlier_dwell, cfg_c, 'networks', {'reference'}, ...
    'fig_size', [width, cfg.height_c], 'axes_w', cfg.axes_w_c, 'key', true);

[fig, B] = composite_figure({struct('fig', f_slow, 'width', width_a), struct('fig', f_sum, 'width', width_b), ...
    struct('fig', f_dw, 'width', width, 'strip', 0.05)}, 'visible', cfg.visible, 'margin', MARGIN, 'width', st.canvas_w, 'row_gap', -0.1);
composite_letter(fig, B{1}(1), 'A', [], 'dy', -0.27);   % at the height of the condition titles
composite_letter(fig, B{2}(1), 'B', [], 'dy', -0.1);
composite_letter(fig, B{3}(1), 'C', [], 'dy', -0.27);   % at the height of its title, as A
canvas_type(fig, supp7_type_exceptions());

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
    'and free (filled). Conditions in their colors. %s Figure %.2f by %.2f in.'], ...
    strjoin(cfg.networks, ', '), what_b, dwell_caption(o_dw), fig.Position(3), fig.Position(4));
note = '';
if ~isempty(gallery_file)
    note = sprintf('The structures of panel B are those of the gallery in `%s`, keys in the outlier_structures stage: %s.', ...
        gallery_file, strjoin(keys, ', '));
    if ~isempty(missing)
        note = sprintf('%s NOT IN PANEL B, because the outlier_structures stage did not run them: %s.', note, strjoin(missing, ', '));
    end
    note = [note newline newline];
end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S7: the reined-in mode with the adaptive variables free; other structures; the dwell time beyond zero', ...
    caption, {'Panel A (fig_outlier_slow_variables)', t_slow; 'Panel B (fig_ostruct_summary)', [note t_sum]; ...
    'Panel C (fig_outlier_dwell)', t_dw}, {o_slow.source, o_sum.source, o_dw.source});
out.gallery_keys = keys;
out.gallery_labels = labels;
out.gallery_missing = missing;
end

%% ------------------------------------------------------------------------
function s = dwell_caption(o)
% Panel C's caption paragraph, with the window, seeds and steps of the data drawn.
q = o.curves(~cellfun(@isempty, o.curves));
q = [q{:}];
win = unique([q.window_s]); n_seeds = max([q.n_seeds]);
fr = o.curves(1, :, 1); fr = fr(~cellfun(@isempty, fr));   % frozen gains
fq = o.curves(1, :, 2); fq = fq(~cellfun(@isempty, fq));   % full Jacobian
dt_fr = unique(cellfun(@(c) c.dt, fr));
dt_full = unique(cellfun(@(c) c.dt, fq));
s = sprintf(['(C) How long the rightmost eigenvalue of the reference network stays beyond zero. An excursion is a run ' ...
    'of consecutive sampled time points at which the rightmost eigenvalue has a positive real part; the curves give the ' ...
    'fraction of excursions lasting longer than the time on the horizontal axis (log), the median over %d seeds of each ' ...
    'seed''s Kaplan-Meier estimate, the band the 25th to 75th percentile across seeds. Dashed, thin: the rightmost ' ...
    'eigenvalue with the gains frozen (slow variables held); solid, thick: the full Jacobian''s rightmost network mode ' ...
    '(slow variables free); colors the conditions as in (A). An excursion cut off by the start or the end of the %s s ' ...
    'window is counted as still going at its observed length, not as ended there; a curve that ends above zero ends ' ...
    'in a triangle: the excursions still going there were cut off by the window, and nothing is known beyond. Sampling ' ...
    'step %s s (gains frozen) and %s s (full Jacobian): no excursion shorter than one step is seen.'], n_seeds, ...
    strjoin(arrayfun(@(v) sprintf('%g', v), win, 'UniformOutput', false), ' and '), ...
    strjoin(arrayfun(@(v) sprintf('%g', v), dt_fr, 'UniformOutput', false), ' and '), ...
    strjoin(arrayfun(@(v) sprintf('%g', v), dt_full, 'UniformOutput', false), ' and '));
end

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

function ex = supp7_type_exceptions()
% The type that may differ from its role's size (canvas_audit), with why. None so far.
ex = cell(0, 2);
end

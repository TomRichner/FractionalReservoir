function out = fig_main_mc_grouped(cfg)
% FIG_MAIN_MC_GROUPED Main figure on memory capacity, in the grouped house style; saved data only.
%
%   out = FIG_MAIN_MC_GROUPED('run_dir', d, 'out_dir', o)
%
% Tom (2026-11-17): the memory-capacity figure of the paper, between group 6
% (fig_main6_grouped) and the stimulation figure. Assembled from two block C
% figures of a mu7chain_gainbias run (FD-029), each drawn from the run's
% saved tables by its own function and its axes copied onto one canvas; no
% simulation runs.
%   (A) the gain x self-coupling maps as the paper draws them
%       (fig_mc_gain_bias_maps with the registry entry
%       'fig_mc_maps_selfcoupling' of mu7chain_gainbias_config): memory
%       capacity (top row) and lambda_1 with the task input on (bottom row,
%       symmetric log, guppy), one column per condition, star at the preset.
%   (B) memory capacity against the input (Wiener) noise sigma_u at the
%       preset (fig_mc_noise_sweeps on the grid noise1d_wiener only, without
%       its lambda_1 row; the release-noise sweep is left out).
% The source figures' header text (seeds, star, symlog and clipping, the
% broken noise axis) is dropped: it belongs in the caption, and the
% provenance file beside the PNG repeats what it said. The colour-bar
% labels are shortened to the quantity, for the same reason.
%
% STYLE as fig_grouped_main and style_main6_grouped: 14-point type
% throughout, axes line width 1, '(A)' panel letters in normal weight,
% condition titles from grouped_condition_label in the condition colours,
% PNG at 200 dpi with 20 px padding, plus the .fig and an .svg.
%
% Options
%   run_dir      the run directory (data/mu7chain_gainbias_<mode>)
%   run_mode     the mode whose registry entry styles the maps ('fastmedium')
%   figure_set   the figure set whose conditions are drawn ('3cond')
%   conditions   names to draw, in order; empty = the figure set's
%   out_dir      where the files go (default: beside this function)
%   save, visible
%
% See also: fig_grouped_main, style_main6_grouped, fig_mc_gain_bias_maps,
%           fig_mc_noise_sweeps, mu7chain_gainbias_config

arguments
    cfg.run_dir    (1,:) char    = ''
    cfg.run_mode   (1,:) char    = 'fastmedium'
    cfg.figure_set (1,:) char    = '3cond'
    cfg.conditions               = {}
    cfg.out_dir    (1,:) char    = ''
    cfg.save       (1,1) logical = true
    cfg.visible    (1,1) logical = false
    cfg.verbose    (1,:) char    = 'minimal' %#ok<INUSA> accepted for a uniform call
end
setup_paths();
root = fileparts(which('setup_paths'));
bundle = mu7chain_gainbias_config(cfg.run_mode, cfg.figure_set);
if isempty(cfg.run_dir); cfg.run_dir = fullfile(root, bundle.run_dir); end
if isempty(cfg.conditions); cfg.conditions = bundle.conditions; end
names = cellstr(cfg.conditions);
nc = numel(names);
entry = registry_entry(bundle.figures, 'fig_mc_maps_selfcoupling');

%% The two source figures, drawn invisibly from saved data
A = fig_mc_gain_bias_maps('run_dir', cfg.run_dir, 'conditions', names, ...
    'save', false, 'visible', false, entry.args{:});
B = fig_mc_noise_sweeps('run_dir', cfg.run_dir, 'grids', {'noise1d_wiener'}, ...
    'show_lle', false, 'conditions', names, 'save', false, 'visible', false);
src_a = A.figs; src_b = B.figs;
cleanup = onCleanup(@() close_valid([src_a src_b]));

%% The canvas (pixels; normalized at the end)
W = 1300; H = 780;
fig = figure('Visible', on_off(cfg.visible), 'Color', 'w', 'Units', 'pixels', ...
    'Position', [40 40 W H]);
drawnow;
% (A): the six maps and their two colour bars, copied together so each
% colour bar keeps its axes.
ha = findall(src_a, '-isa', 'matlab.graphics.axis.Axes', '-or', 'Type', 'colorbar');
maps_src = findall(src_a, 'Type', 'axes');
assert(numel(maps_src) == 2 * nc, 'fig_main_mc_grouped:Maps', ...
    'Expected %d map axes (two rows, %d conditions), found %d.', 2 * nc, nc, numel(maps_src));
% Each row is its own inner tiledlayout, so a position is only comparable
% within a row: the row is read from the colour limits (memory capacity
% from 0, lambda_1 symmetric about 0), the column from x within the row.
% The source axes are tagged, and the tags survive the copy.
is_lle_row = arrayfun(@(h) h.CLim(1) < 0, maps_src);
assert(nnz(is_lle_row) == nc, 'fig_main_mc_grouped:Rows', 'Expected %d lambda_1 maps.', nc);
for r = 1:2
    row = maps_src(is_lle_row == (r == 2));
    [~, ord] = sort(arrayfun(@(h) h.Position(1), row));
    for c = 1:nc; row(ord(c)).Tag = sprintf('main_mc_r%d_c%d', r, c); end
end
copied = copyobj(ha, fig);
cbs = copied(arrayfun(@(h) isa(h, 'matlab.graphics.illustration.ColorBar'), copied));
assert(numel(cbs) == 2, 'fig_main_mc_grouped:Colorbars', 'Expected two colour bars.');
maps = gobjects(2, nc);                    % maps(r, c); row 1 memory capacity
for r = 1:2
    for c = 1:nc; maps(r, c) = findobj(copied, 'flat', 'Tag', sprintf('main_mc_r%d_c%d', r, c)); end
end
is_lle_cb = arrayfun(@(cb) contains(cb.Label.String, 'lambda'), cbs);
cb_mc = cbs(~is_lle_cb); cb_l = cbs(is_lle_cb);
% (B): the memory-capacity panel and its legend.
hb = findall(src_b, '-isa', 'matlab.graphics.axis.Axes', '-or', 'Type', 'legend');
copied_b = copyobj(hb, fig);
axb = copied_b(arrayfun(@(h) isa(h, 'matlab.graphics.axis.Axes'), copied_b));
lgb = copied_b(arrayfun(@(h) isa(h, 'matlab.graphics.illustration.Legend'), copied_b));
assert(isscalar(axb) && isscalar(lgb), 'fig_main_mc_grouped:Noise', ...
    'Expected one noise-sweep axes and its legend.');

%% Geometry (pixels)
mw = 200; mh = 225; gap = 16;            % one map; gap between columns
x0 = 120;                                % left edge of the first map column
y_top = 470; y_bot = 205;                % bottom edges of the two map rows
for r = 1:2
    for c = 1:nc
        set(maps(r, c), 'Units', 'pixels', 'PositionConstraint', 'innerposition', ...
            'Position', [x0 + (c - 1) * (mw + gap), (r == 1) * y_top + (r == 2) * y_bot, mw, mh]);
    end
end
x_cb = x0 + nc * mw + (nc - 1) * gap + 14;
set(cb_mc, 'Units', 'pixels', 'Position', [x_cb, y_top, 14, mh]);
set(cb_l, 'Units', 'pixels', 'Position', [x_cb, y_bot, 14, mh]);
xb = x_cb + 210; wb = W - xb - 30;
set(axb, 'Units', 'pixels', 'PositionConstraint', 'innerposition', ...
    'Position', [xb, y_top - 75, wb, mh + 75]);

%% Style
all_ax = [maps(:); axb];
set(all_ax, 'FontSize', 14, 'LineWidth', 1, 'TitleFontSizeMultiplier', 1, ...
    'LabelFontSizeMultiplier', 1);
set([cb_mc cb_l], 'FontSize', 14, 'LineWidth', 1);
titles = grouped_condition_label(names, 'two_line');
gain_lbl = maps(2, 1).XLabel.String;
s_lbl = maps(2, 1).YLabel.String;
if iscell(s_lbl); s_lbl = s_lbl{end}; end
for c = 1:nc
    col = maps(1, c).Title.Color;           % the condition colour the source used
    title(maps(1, c), titles{c}, 'FontWeight', 'normal', 'FontSize', 14, 'Color', col, ...
        'Interpreter', 'tex');
    for r = 1:2
        ax = maps(r, c);
        % Every grid level keeps its tick; every other one is labelled, so
        % eleven levels fit a map at 14 points.
        thin_labels(ax, 'X', label_levels(ax.XTickLabel));
        if c == 1
            thin_labels(ax, 'Y', label_levels(ax.YTickLabel));
            ylabel(ax, title_case(s_lbl), 'FontSize', 14);
        else
            ax.YLabel.String = '';
        end
        if r == 2
            xlabel(ax, title_case(gain_lbl), 'FontSize', 14);
        end
    end
end
cb_mc.Label.String = 'Memory Capacity';
cb_l.Label.String = 'Lyapunov Exp., \lambda_1 (s^{-1})';
set([cb_mc.Label cb_l.Label], 'FontSize', 14, 'Interpreter', 'tex');
% (B)
noise_lbl = axb.XLabel.String;
title(axb, '');
xlabel(axb, title_case(noise_lbl), 'FontSize', 14);
ylabel(axb, 'Memory Capacity', 'FontSize', 14);
thin_labels(axb, 'X', label_levels(axb.XTickLabel));
lgb.String = grouped_condition_label(names, 'title');
set(lgb, 'FontSize', 14, 'Box', 'off', 'Orientation', 'vertical', 'NumColumns', 1, ...
    'Units', 'pixels');
drawnow;
lp = lgb.Position;
lgb.Position = [xb + wb / 2 - lp(3) / 2, y_bot - 15, lp(3), lp(4)];
% Panel letters, as style_main6_grouped places them.
delete(findall(fig, 'Tag', 'main_mc_panel_letter'));
letters = {maps(1, 1), axb};
for k = 1:2
    ax = letters{k}; p = ax.Position;
    text(ax, -95 / p(3), 1 + 55 / p(4), sprintf('(%c)', 'A' + k - 1), ...
        'Units', 'normalized', 'Clipping', 'off', 'FontSize', 14, 'FontWeight', 'normal', ...
        'Interpreter', 'none', 'Tag', 'main_mc_panel_letter');
end
set([maps(:); axb; cb_mc; cb_l; lgb], 'Units', 'normalized');

%% Save
tag = 'Fig_Main_MC_Grouped';
out = struct('figs', fig, 'files', {{}}, 'source', {[{A.source}, B.source]}, 'layout', ...
    struct('A', 'gain x self-coupling maps: memory capacity (top), lambda_1 (bottom)', ...
           'B', 'memory capacity against input noise sigma_u'));
if cfg.save
    outdir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    if ~isfolder(outdir); mkdir(outdir); end
    for ext = {'.png', '.fig', '.svg'}
        old = fullfile(outdir, [tag ext{1}]);
        if isfile(old); delete(old); end
    end
    drawnow;
    % The export-resize guard of fig_grouped_main.
    pos = fig.Position;
    for attempt = 1:3
        exportgraphics(fig, fullfile(outdir, [tag '.png']), 'Resolution', 200, 'Padding', 20);
        if isequal(fig.Position(3:4), pos(3:4)); break; end
        fig.Position = pos; drawnow;
        if attempt == 3
            warning('fig_main_mc_grouped:ExportResized', ...
                'exportgraphics resized %s on every attempt; the PNG may be on the wrong canvas.', tag);
        end
    end
    print(fig, fullfile(outdir, [tag '.svg']), '-dsvg', '-vector');
    savefig(fig, fullfile(outdir, [tag '.fig']));
    out.files = existing_outputs(outdir, tag);
    write_provenance(fullfile(outdir, [tag '_provenance.md']), cfg, root, out.source, names, entry);
    out.files{end+1} = [tag '_provenance.md'];
end
end

%% ------------------------------------------------------------------------
function e = registry_entry(figures, name)
e = [];
for k = 1:numel(figures)
    if strcmp(figures{k}.name, name); e = figures{k}; return; end
end
error('fig_main_mc_grouped:NoEntry', 'The bundle registers no figure %s.', name);
end

function keep = label_levels(labels)
% Every other level labelled, ending with the last, so the top of the
% range always reads (with an odd count the first reads too).
n = numel(labels);
keep = false(1, n);
keep(n:-2:1) = true;
end

function thin_labels(ax, which, keep)
lbl = cellstr(ax.([which 'TickLabel']));
lbl(~keep) = {''};
ax.([which 'TickLabel']) = lbl;
end

function s = title_case(s)
% Capitalise each word of an axis label (house style: 'Longest SFA Time
% Constant (s)'), leaving symbols such as \sigma_u alone.
if iscell(s); s = strjoin(s, ' '); end
w = strsplit(s, ' ');
for k = 1:numel(w)
    % A one-letter word is a symbol (the self-coupling s): left as it is.
    if numel(w{k}) > 1 && isletter(w{k}(1)); w{k}(1) = upper(w{k}(1)); end
end
s = strjoin(w, ' ');
s = regexprep(s, '(\w)-([a-z])', '$1-${upper($2)}');
end

function close_valid(h)
h = h(isgraphics(h));
if ~isempty(h); close(h); end
end

function v = on_off(b)
if b; v = 'on'; else; v = 'off'; end
end

function write_provenance(f, cfg, root, sources, names, entry)
[~, head] = system(sprintf('git -C "%s" rev-parse HEAD', root));
[~, dirty] = system(sprintf('git -C "%s" status --porcelain -- src/figures/fig_main_mc_grouped.m', root));
fid = fopen(f, 'w'); guard = onCleanup(@() fclose(fid));
fprintf(fid, '# Fig_Main_MC_Grouped\n\n');
fprintf(fid, 'Main memory-capacity figure, drawn by `src/figures/fig_main_mc_grouped.m` at commit `%s`%s.\n\n', ...
    strtrim(head), ternary(isempty(strtrim(dirty)), '', ' (with uncommitted changes to the function)'));
fprintf(fid, 'Run: `%s`\n\n', cfg.run_dir);
fprintf(fid, 'Conditions: %s\n\n', strjoin(names, ', '));
fprintf(fid, '## Sources (saved data only; no simulation)\n\n');
for k = 1:numel(sources); fprintf(fid, '- `%s`\n', string(sources{k})); end
fprintf(fid, ['\n## Panels\n\n' ...
    '- (A) `fig_mc_gain_bias_maps` with the registry entry `fig_mc_maps_selfcoupling` of ' ...
    '`mu7chain_gainbias_config(''%s'', ''%s'')`: memory capacity (median over seeds, linear, white to ' ...
    'the guppy center colour; top row) and lambda_1 with the task input on (median over seeds, ' ...
    'symmetric log, linear within +-%g 1/s, clipped at +-%g 1/s, guppy, darkest = 0; bottom row), ' ...
    'gain by self-coupling s, one column per condition. White star: the preset. Grey: no successful job.\n' ...
    '- (B) `fig_mc_noise_sweeps` on the grid `noise1d_wiener` only, `show_lle` false: memory capacity ' ...
    'of the preset network against the input (Wiener) noise sigma_u; line median over networks, band ' ...
    'smallest to largest network, points each network. Noise axis: 0, then log-spaced levels, broken ' ...
    'between 0 and the lowest positive level.\n\n'], cfg.run_mode, cfg.figure_set, ...
    arg_of(entry.args, 'lle_lin'), arg_of(entry.args, 'lle_clip'));
fprintf(fid, ['## Dropped from the source figures\n\n' ...
    '- The header text of both source figures (seeds, star, darkest colour, the broken noise axis): ' ...
    'it goes in the caption (stated above).\n' ...
    '- The colour-bar labels'' scale wording ("symmetric log, clipped at +-5"): the bars now name ' ...
    'the quantity only; the scale is stated above.\n' ...
    '- The noise-sweep figure''s lambda_1 row and its release-noise (sigma_gamma) column.\n' ...
    '- Every other tick label on the map and noise axes (every level keeps its tick mark).\n']);
end

function v = arg_of(args, name)
k = find(strcmp(args(1:2:end), name), 1);
v = args{2 * k};
end

function s = ternary(c, a, b)
if c; s = a; else; s = b; end
end

function out = fig_supp2_grouped(cfg)
% FIG_SUPP2_GROUPED Supplementary Figure S2: numerical checks; saved data only.
%
%   out = FIG_SUPP2_GROUPED('run_dir', d, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S2), for the
% Methods' "Numerical integration", as revised after Tom's review of
% 2026-11-22 (Plan2_figure_feedback.md, S2). One row per check, one column
% per condition, each row taken from the numerics figure that owns it
% (fig_numerics_verification, the numerics_verification stage of run_dir),
% drawn invisibly and unsaved and copied here (composite_tools):
%   (A) step refinement of the solver: the bottom row of the 'solver'
%       variant (rms reshooting error against the step, noise-free and
%       with noise, with the reference slopes; the slope-2 reference is a
%       pale band behind the noise-free data, which lie on it).
%   (B) the two Lyapunov estimators: the top row of the 'lya_method'
%       variant, with FOUR traces: Benettin's local rate and exponent
%       (thick, the condition's colour, behind) and the QR method's local
%       lambda_1 and lambda_1 (thin, on top). The lower QR exponents are
%       not drawn (n_qr_show = 1). The y axis is zoomed (local_ylim_*), so
%       the transient at the start of the run is clipped.
%   (C) the analytic Jacobian against central finite differences: the
%       'jacobian' variant (relative Frobenius error and the mismatch of
%       the leading eigenvalues at the sampled time points, the snapshots),
%       every column on jacobian_ylim. The acceptance threshold the source
%       draws is taken out: the panel reports what the error was, not what
%       was allowed.
% The other rows of the source figures (example trajectories, the
% reshooting error in time, the full QR spectrum) are not drawn; the
% numerics figures themselves stay in the appendix set.
%
% LAYOUT. Made on the common canvas itself (canvas_w of manuscript_style, in
% pixels), not on a wider sheet squeezed onto it: three columns across the
% whole width and each row's key in a line under the row. (Until 2026-11-22
% the keys stood at the right of the rows, laid out on 1100 px; squeezed to
% 9 in they reached past the canvas.)
%
% TEXT. The source panels carry their numbers in multi-line titles and, for
% (A), the fitted slopes in every legend. Here the columns are titled by
% condition once, each row has one key, and those numbers are in
% Fig_Supp2_Grouped_table.md, read from the source figures. The table also
% carries two sentences for the caption (the measured slopes; the network
% without adaptation at its fixed point).
%
% Options: run_dir, conditions (empty = all the run holds), n_qr_show,
% local_ylim_no_adaptation, local_ylim_adaptation, jacobian_ylim, out_dir,
% save, visible, fig_tag; preset_name is accepted and unused (the stage
% records its preset).
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, composite_save, fig_numerics_verification

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
    cfg.n_qr_show   (1,1) double  = 1     % QR local exponents drawn in (B): lambda_1 alone
    cfg.local_ylim_no_adaptation (1,2) double = [-15 5]        % (B), the column without adaptation (1/s)
    cfg.local_ylim_adaptation    (1,2) double = [-4 4]         % (B), every other column
    cfg.jacobian_ylim            (1,2) double = [1e-12 1e-8]   % (C), every column
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.fig_tag     (1,:) char    = 'Fig_Supp2_Grouped'
end
setup_paths();
T = composite_tools();
st = manuscript_style();
variants = {'solver', 'lya_method', 'jacobian'};
source_row = {'last', 'first', 'first'};       % which row of each variant is drawn

%% The source figures
S = cell(1, 3);
cleanup = onCleanup(@() T.close(findall(groot, 'Type', 'figure', 'Tag', 'fig_supp2_source')));
for v = 1:3
    S{v} = fig_numerics_verification('variant', variants{v}, 'run_dir', cfg.run_dir, ...
        'conditions', cfg.conditions, 'n_qr_show', cfg.n_qr_show, 'verbose', cfg.verbose, 'save', false, 'visible', false);
    set(S{v}.figs, 'Tag', 'fig_supp2_source');
end

%% The canvas: the common width itself, so that T.fit has nothing to stretch
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
X0 = 98; RM = 16; GAP = 62;     % left of the first plot box; right margin; between columns (the next column's tick labels)
RH = 160;                       % a row's plot boxes
TOP = 64;                       % column titles on two lines, and the letter (A)
BELOW = 58;                     % tick labels and the axis name under a row
KEY_H = [50 28 28];             % room for each row's key (two lines, one, one)
BETWEEN = 34;                   % from a key to the next row (its letter)
H = TOP + 3 * RH + 3 * BELOW + sum(KEY_H) + 2 * BETWEEN + 10;
fig = T.canvas(W, H, cfg.visible);
ROW_Y = zeros(1, 3);
y = H - TOP;
for v = 1:3
    ROW_Y(v) = y - RH;
    y = ROW_Y(v) - BELOW - KEY_H(v) - BETWEEN;
end
ax = cell(1, 3); lg = cell(1, 3); notes = cell(1, 3);
nc = 0; CW = 0;
for v = 1:3
    G = T.grid(findall(S{v}.figs, 'Type', 'axes'));
    r = 1;
    if strcmp(source_row{v}, 'last'); r = size(G, 1); end
    if v == 1
        nc = size(G, 2);
        CW = (W - X0 - RM - (nc - 1) * GAP) / nc;
        % The conditions drawn, from the first row's titles (their first line).
        cond_titles = arrayfun(@(a) first_line(a.Title.String), G(1, :), 'UniformOutput', false);
    end
    assert(size(G, 2) == nc, 'fig_supp2_grouped:Columns', ...
        'The ''%s'' variant draws %d columns, the ''solver'' variant %d.', variants{v}, size(G, 2), nc);
    notes{v} = panel_text(G, T);
    Cp = T.copy(reshape(G(r, :), [], 1), fig);
    ax{v} = Cp.axes;
    lg{v} = Cp.legends;
    for c = 1:nc
        ax{v}(c).Position = [X0 + (c - 1) * (CW + GAP), ROW_Y(v), CW, RH];
    end
end

%% Style
T.style(fig);
names = condition_names(cond_titles);
for v = 1:3
    for c = 1:nc
        title(ax{v}(c), '');
        if c > 1; ylabel(ax{v}(c), ''); end
    end
    % One key per row: the first column's, the others go.
    keep = lg{v}(arrayfun(@(l) isequal(l.Axes, ax{v}(1)), lg{v}));
    delete(lg{v}(~ismember(lg{v}, keep)));
    lg{v} = keep;
end
if all(~cellfun(@isempty, names))
    T.retitle(ax{1}, names, true);
else
    for c = 1:nc; title(ax{1}(c), cond_titles{c}, 'FontWeight', 'normal', 'FontSize', st.title_fs); end
end
% (A): the fitted slopes differ by condition; the key names the lines.
if ~isempty(lg{1})
    lg{1}.String = regexprep(lg{1}.String, ',\s*slope\s*[-+0-9.]+', '');
end
ylabel(ax{1}(1), {'RMS |Error|', 'over 0.02 s'});
for c = 1:nc
    xlabel(ax{1}(c), 'Step (s)');
    % The end ticks a little inside the panel, so that their labels stay level
    % and clear of the next column.
    xl = sort(ax{1}(c).XLim);
    if strcmp(ax{1}(c).XScale, 'log') && xl(1) > 0; ax{1}(c).XLim = [xl(1) / 1.12, xl(2) * 1.12]; end
    ax{1}(c).XTickLabelRotation = 0;
end
% (B): zoomed, even where that clips the first transient.
ylabel(ax{2}(1), {'Local Lyapunov', 'Exponent (s^{-1})'}, 'Interpreter', 'tex');
for c = 1:nc
    xlabel(ax{2}(c), 'Time (s)');
    yl = cfg.local_ylim_adaptation;
    if strcmp(names{c}, 'no_adaptation'); yl = cfg.local_ylim_no_adaptation; end
    ax{2}(c).YLim = yl;
    ax{2}(c).YTick = zoom_ticks(yl);
end
% (C): one y axis for every column, labelled on the first; no threshold; the
% counts in the source's x label are in the table.
delete(findall([ax{3}], 'Tag', 'numerics_jacobian_threshold'));
ylabel(ax{3}(1), 'Error');
n_clipped = 0;
for c = 1:nc
    xlabel(ax{3}(c), 'Snapshot');
    set(ax{3}(c), 'YLim', cfg.jacobian_ylim, 'YTick', decade_ticks(cfg.jacobian_ylim), 'YMinorTick', 'off');
    if c > 1; ax{3}(c).YTickLabel = {}; end
    ln = findall(ax{3}(c), 'Type', 'line');
    for k = 1:numel(ln)
        n_clipped = n_clipped + nnz(ln(k).YData < cfg.jacobian_ylim(1) | ln(k).YData > cfg.jacobian_ylim(2));
    end
end
% The keys: in sentence case, in one line where that fits the canvas, centred under their rows.
for v = 1:3
    if isempty(lg{v}); continue; end
    lg{v}.String = T.sentence(lg{v}.String);
    set(lg{v}, 'Units', 'pixels', 'Location', 'none', 'Orientation', 'horizontal', 'NumColumns', numel(lg{v}.String));
    lg{v}.ItemTokenSize = [22 18];
    drawnow;
    if lg{v}.Position(3) > W - 2 * RM; lg{v}.NumColumns = ceil(numel(lg{v}.String) / 2); end
    for pass = 1:2
        drawnow;
        p = lg{v}.Position;
        lg{v}.Position = [(W - p(3)) / 2, ROW_Y(v) - BELOW - p(4) - 2, p(3), p(4)];
    end
end

%% The common canvas (the layout is already on it; fit leaves it as it is)
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);

%% Panel letters
letters = 'ABC';
dy = [40 10 10];
for v = 1:3; T.letter(ax{v}(1), letters(v), -(X0 - 6), dy(v)); end

%% Save
panels = struct( ...
    'A', 'fig_numerics_verification, variant ''solver'' (its bottom row): convergence of the reshooting error with the step', ...
    'B', 'fig_numerics_verification, variant ''lya_method'' (its top row): local rates of Benettin''s and the QR method', ...
    'C', 'fig_numerics_verification, variant ''jacobian'': analytic Jacobian against central finite differences');
sources = cellfun(@(s) s.source, S, 'UniformOutput', false);
parts = cell(0, 2);
if cfg.save
    f = fieldnames(panels);
    for v = 1:3
        body = [{'What the source figure writes on its panels (titles and legends, every row of the source):', ''}, notes{v}, {''}];
        if v == 2
            body = [body, {sprintf(['Drawn on y limits %s (no adaptation) and %s (the other columns): the first transient of ' ...
                'the local rates is clipped. The table rows below give the whole range of each trace.'], ...
                mat2str(cfg.local_ylim_no_adaptation), mat2str(cfg.local_ylim_adaptation)), ''}]; %#ok<AGROW>
        elseif v == 3
            body = [body, {sprintf(['Every column on y limits %g to %g; %d point(s) fall outside them. The acceptance ' ...
                'threshold the source draws (see the legends above) is not drawn here.'], cfg.jacobian_ylim, n_clipped), ''}]; %#ok<AGROW>
        end
        for c = 1:nc
            body{end + 1} = T.axes_table(ax{v}(c), sprintf('Column %d, %s', c, cond_titles{c})); %#ok<AGROW>
        end
        parts(end + 1, :) = {sprintf('Panel %s, %s', f{v}, panels.(f{v})), strjoin(body, newline)}; %#ok<AGROW>
    end
    parts(end + 1, :) = {'For the caption (Tom''s review of 2026-11-22, S2 notes 1 and 9)', strjoin({ ...
        ['- (A) The order of the solver with noise is at least 1.5 (what SRA1 guarantees for additive noise, the lower bound ' ...
         'drawn) and close to 2 over the steps tested, where the error of the drift (order 2) dominates. The measured slopes ' ...
         'are in the legends listed under Panel A.'], ...
        ['- (B, C) Without adaptation this small network (the reduced network of the checks) settles to a fixed point, so its ' ...
         'exponents equal -1/tau_d and the sampled time points coincide: the errors of (C) are the same at every snapshot ' ...
         'of that column. This is expected.'], ''}, newline)};
end
caption = sprintf(['Supplementary Figure S2, numerical checks (drawn by fig_supp2_grouped; run `%s`), one column per ' ...
    'condition. (A) Step refinement of the solver: rms error of the fixed-step SRA1 solution reshot from a tight ' ...
    'ode45 reference over 0.02 s, against the step, without noise (filled circles) and with the input noise (open ' ...
    'squares), with the reference slopes: 2 (the pale band behind the noise-free data, which lie on it) and 1.5 ' ...
    '(the lower bound of SRA1 with noise; the measured slopes are steeper, close to 2). (B) The leading exponent by ' ...
    'two estimators on the reduced network: the local rate of Benettin''s method (thick, the condition''s colour) with ' ...
    'the local lambda_1 of the QR method on top of it (thin, dark gray), and both accumulated values (horizontal ' ...
    'lines: Benettin dashed, QR dotted). The y axis is zoomed and clips the first transient. (C) The analytic Jacobian ' ...
    'against central finite differences at %d sampled time points (snapshots) along the run: relative Frobenius error ' ...
    'and mismatch of the leading eigenvalues. Without adaptation the reduced network settles to a fixed point, so its ' ...
    'exponents equal -1/tau_d and its snapshots coincide. The numbers the source figures print in their titles and ' ...
    'legends are in this table. Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, n_snapshots(ax{3}(1)), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S2: numerical checks', caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function t = zoom_ticks(yl)
% Ticks of a zoomed axis of (B): the two limits and 0 when the axis is
% symmetric, else every 5 between the limits.
if yl(1) == -yl(2)
    t = [yl(1), 0, yl(2)];
else
    t = ceil(yl(1) / 5) * 5 : 5 : floor(yl(2) / 5) * 5;
end
end

function t = decade_ticks(yl)
% Every second decade between the limits of a log axis, from the lower limit.
e = ceil(log10(yl(1)) - 1e-9) : 2 : floor(log10(yl(2)) + 1e-9);
t = 10 .^ e;
end

function n = n_snapshots(ax)
% The number of points of the longest line of a panel of (C).
ln = findall(ax, 'Type', 'line');
n = max([0, arrayfun(@(h) numel(h.XData), ln(:)')]);
end

function s = first_line(str)
c = cellstr(str);
s = strtrim(c{1});
end

function names = condition_names(titles)
% Condition names from the titles the source painted; '' where a title is
% no condition's (then the source's own text is kept).
names = cell(size(titles));
for k = 1:numel(titles)
    try
        n = grouped_condition_label(titles(k), 'name_of_title');
        names{k} = n{1};
    catch
        names{k} = '';
    end
end
end

function lines = panel_text(G, T) %#ok<INUSD>
% Titles and legend entries of every panel of a source figure, row by row.
lines = {};
src = ancestor(G(1), 'figure');
lgs = findall(src, 'Type', 'legend');
for r = 1:size(G, 1)
    for c = 1:size(G, 2)
        if ~isgraphics(G(r, c)); continue; end
        t = strjoin(strtrim(cellstr(G(r, c).Title.String)), '; ');
        lines{end + 1} = sprintf('- row %d, column %d: %s', r, c, strrep(t, '|', '/')); %#ok<AGROW>
        own = lgs(arrayfun(@(l) isequal(l.Axes, G(r, c)), lgs));
        if ~isempty(own)
            lines{end + 1} = sprintf('  - legend: %s', strrep(strjoin(own(1).String, '; '), '|', '/')); %#ok<AGROW>
        end
    end
end
end

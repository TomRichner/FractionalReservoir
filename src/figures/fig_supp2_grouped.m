function out = fig_supp2_grouped(cfg)
% FIG_SUPP2_GROUPED Supplementary Figure S2: numerical checks; saved data only.
%
%   out = FIG_SUPP2_GROUPED('run_dir', d, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S2), for the
% Methods' "Numerical integration". One row per check, one column per
% condition, each row taken from the numerics figure that owns it
% (fig_numerics_verification, the numerics_verification stage of run_dir),
% drawn invisibly and unsaved and copied here (composite_tools):
%   (A) step refinement of the solver: the bottom row of the 'solver'
%       variant (rms reshooting error against the step, noise-free and
%       with noise, with the reference slopes).
%   (B) the two Lyapunov estimators: the top row of the 'lya_method'
%       variant (the local rates of Benettin's method and of the QR method
%       on the reduced network, with both accumulated values).
%   (C) the analytic Jacobian against central finite differences: the
%       'jacobian' variant (relative Frobenius error and the mismatch of
%       the leading eigenvalues at the sampled states, with the threshold).
% The other rows of the source figures (example trajectories, the
% reshooting error in time, the full QR spectrum) are not drawn; the
% numerics figures themselves stay in the appendix set.
%
% TEXT. The source panels carry their numbers in multi-line titles and, for
% (A), the fitted slopes in every legend. Here the columns are titled by
% condition once, each row has one legend at the right, and those numbers
% are in Fig_Supp2_Grouped_table.md, read from the source figures.
%
% Options: run_dir, conditions (empty = all the run holds), out_dir, save,
% visible, fig_tag; preset_name is accepted and unused (the stage records
% its preset).
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, composite_save, fig_numerics_verification

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused
    cfg.conditions                = {}
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
        'conditions', cfg.conditions, 'verbose', cfg.verbose, 'save', false, 'visible', false);
    set(S{v}.figs, 'Tag', 'fig_supp2_source');
end

%% The canvas
W = 1100; H = 915;
fig = T.canvas(W, H, cfg.visible);
X0 = 100; CW = 195; GAP = 62; RH = 190;     % CW: the keys at the right end near the canvas edge
ROW_Y = [650 365 80];
ax = cell(1, 3); lg = cell(1, 3); notes = cell(1, 3);
nc = 0;
for v = 1:3
    G = T.grid(findall(S{v}.figs, 'Type', 'axes'));
    r = 1;
    if strcmp(source_row{v}, 'last'); r = size(G, 1); end
    if v == 1
        nc = size(G, 2);
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
    % One legend per row, at the right: the first column's, the others go.
    keep = lg{v}(arrayfun(@(l) isequal(l.Axes, ax{v}(1)), lg{v}));
    delete(lg{v}(~ismember(lg{v}, keep)));
    lg{v} = keep;
    if ~isempty(keep); keep(1).String = T.sentence(keep(1).String); end   % legend entries in sentence case
end
if all(~cellfun(@isempty, names))
    T.retitle(ax{1}, names, true);
else
    for c = 1:nc; title(ax{1}(c), cond_titles{c}, 'FontWeight', 'normal', 'FontSize', st.title_fs); end
end
% (A): the fitted slopes differ by condition; the legend names the lines.
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
% (B)
ylabel(ax{2}(1), {'Local Lyapunov', 'Exponent (s^{-1})'}, 'Interpreter', 'tex');
for c = 1:nc; xlabel(ax{2}(c), 'Time (s)'); end
% (C): the counts in the source's x label are in the table.
ylabel(ax{3}(1), 'Error');
for c = 1:nc; xlabel(ax{3}(c), 'Sampled State'); end
x_lg = X0 + nc * CW + (nc - 1) * GAP + 12;
for v = 1:3
    if isempty(lg{v}); continue; end
    set(lg{v}, 'Units', 'pixels', 'Location', 'none');
    lg{v}.ItemTokenSize = [22 18];
    for pass = 1:2
        drawnow;
        p = lg{v}.Position;
        lg{v}.Position = [x_lg, ROW_Y(v) + RH - p(4), p(3), p(4)];
    end
end

%% The common canvas: the layout above, made on W px, stretched across canvas_w
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);

%% Panel letters
letters = 'ABC';
dy = [42 22 22];
for v = 1:3; T.letter(ax{v}(1), letters(v), -90, dy(v)); end

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
        for c = 1:nc
            body{end + 1} = T.axes_table(ax{v}(c), sprintf('Column %d, %s', c, cond_titles{c})); %#ok<AGROW>
        end
        parts(end + 1, :) = {sprintf('Panel %s, %s', f{v}, panels.(f{v})), strjoin(body, newline)}; %#ok<AGROW>
    end
end
caption = sprintf(['Supplementary Figure S2, numerical checks (drawn by fig_supp2_grouped; run `%s`), one column per ' ...
    'condition. (A) Step refinement of the solver: rms error of the fixed-step SRA1 solution reshot from a tight ' ...
    'ode45 reference over 0.02 s, against the step, without noise (filled circles) and with the input noise (open ' ...
    'squares), with the reference slopes 2 and 1.5. (B) The leading exponent by two estimators on the reduced ' ...
    'network: the local rates of the QR method (thin, gray, the first exponent dark to the last light) and of Benettin''s method (thick), and both accumulated ' ...
    'values (horizontal lines). (C) The analytic Jacobian against central finite differences at sampled states: ' ...
    'relative Frobenius error and mismatch of the leading eigenvalues, with the acceptance threshold (dotted). ' ...
    'The numbers the source figures print in their titles and legends are in this table. Canvas %d x %d px (canvas_w of manuscript_style), ' ...
    'type in points by role.'], cfg.run_dir, getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S2: numerical checks', caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
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

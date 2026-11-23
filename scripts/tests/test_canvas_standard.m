function T = test_canvas_standard(names, opts)
%TEST_CANVAS_STANDARD The 18 main and supplementary figures on the common canvas, with their type by role.
%
%   test_canvas_standard                              % all 18
%   test_canvas_standard({'fig_main5_grouped', 'fig_supp6_grouped'})
%   test_canvas_standard('fig_main2_grouped', 'out_dir', 'figs/p3')
%   T = test_canvas_standard(...)                     % the table, no error
%
% THE STANDARD (the manuscript repository's
% docs/Advancing_the_narrative/Plan2_canvas_standard.md; Tom, 2026-11-22:
% "The figures should have similar canvas widths so that the letters appear
% proportionally the same", and later the same day "we should probably get
% all figures down to a canvas width of 9 inches"). The numbers are
% manuscript_style's (canvas_w, canvas_max_h, canvas_tol, tick_fs,
% label_fs, title_fs, letter_fs, axis_lw); none is typed here. Every figure
% is built as its bundle's registry builds it, from the `fast` run, and
% checked:
%
%   width    the figure is canvas_w wide (to 0.01 in), measured on the
%            figure as drawn (MATLAB may shrink a window larger than the
%            screen, so the size asked for is not trusted);
%   height   at most canvas_max_h;
%   type     every axes (its tick labels), axis label, title, legend,
%            colour bar, free text and panel letter at its role's size
%            (canvas_audit: ticks, legends and free text tick_fs; axis and
%            colour-bar labels label_fs; titles title_fs; panel letters
%            letter_fs), except the items named in EXCEPTIONS below, each
%            with its reason. An exception the figure no longer uses fails
%            too, so the list stays true;
%   lines    the axis lines of every visible axes are axis_lw;
%   inside   no axes (with its tick labels, axis names and title: its tight
%            inset) and no legend reaches past the canvas;
%   png      the exported PNG is within canvas_tol of canvas_w wide (the
%            cropped image fills the width) and no taller than canvas_max_h.
%
% For each figure a row is printed: its size, the PNG's size, the sizes
% found per role, the exceptions used, PASS or FAIL, and for a FAIL every
% reason. A figure whose data are not on this machine is SKIP, with why.
% At the end, the failing figures by name. Called without an output, the
% test throws (test_canvas_standard:Failed) when any figure fails or is
% skipped; with an output it returns the table and throws nothing.
%
% names    registry names (FIGURES below), one or several; empty (the
%          default): all 18. A package runs it on its own figures.
% out_dir  a folder: each figure is saved into <out_dir>/<name>/ and kept
%          (to look at). Default: a temporary folder, deleted after.
%          Figure 8's images hold clinical recordings: they stay under
%          figs/ and are never added to git.
%
% THE DATA. data/mu7chain_AC_fast and data/IED_model/ied_stim_fast of this
% checkout, else of the main checkout. The mu7chain run is COPIED to a
% temporary folder first, because fig_param_space_allStd writes a prep
% folder into the run directory it is given. Figure 8 also needs its model
% figure, figs/IED_model/ied_stim_fast/fig_ied_stim_blocks/
% fig_ied_stim_blocks.fig of this checkout, and the clinical panels named
% by ied_stim_config.
%
% Drawing only, nothing is simulated; about 15 minutes for all 18, figures
% invisible. Assumes setup_paths has run (pool_set(0) will do).
%
% Replaces test_canvas_standard_a and _b (Plan 2's groups A and B, which
% had seven and eleven figures between them and their own checks).
%
% See also: canvas_audit, canvas_type, manuscript_style, composite_tools

arguments
    names = {}
    opts.out_dir (1,:) char = ''
end

FIGURES = { ...
    'fig_main1_grouped', 'fig_main2_grouped', 'fig_main3_grouped', 'fig_main4_grouped', ...
    'fig_main5_grouped', 'fig_main6_grouped', 'fig_main7_grouped', 'Fig_Main8_Stim', ...
    'fig_supp1_grouped', 'fig_supp2_grouped', 'fig_supp3_grouped', 'fig_supp4_grouped', ...
    'fig_supp5_grouped', 'fig_supp6_grouped', 'fig_supp7_grouped', 'fig_supp8_grouped', ...
    'fig_supp9_grouped', 'fig_supp10_grouped'};

% The type that may depart from its role's size, by figure: {pattern on
% canvas_audit's name, reason}. Every figure not listed has none.
EXCEPTIONS = struct();
EXCEPTIONS.fig_main5_grouped = {'^[EI] \| axes fig_outlier_plane_inset', ...
    'E and I marks inside the response-curve insets (inset_mark_fs): at 14 pt they cover the points they label'};

if isempty(names); names = FIGURES; end
names = cellstr(names);
unknown = setdiff(names, FIGURES);
assert(isempty(unknown), 'test_canvas_standard:Unknown', 'Not one of the 18 figures: %s', strjoin(unknown, ', '));
names = FIGURES(ismember(FIGURES, names));          % in the paper's order

st = manuscript_style();
root = fileparts(which('setup_paths'));
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));
out_root = opts.out_dir;
if isempty(out_root); out_root = fullfile(tmp, 'figs'); end
if ~is_absolute(out_root); out_root = fullfile(root, out_root); end
figs_before = findall(groot, 'Type', 'figure');

ied_names = {'Fig_Main8_Stim', 'fig_supp10_grouped'};
run_dir = ''; ied_dir = find_run(root, fullfile('data', 'IED_model', 'ied_stim_fast'));
if any(~ismember(names, ied_names))
    src = find_run(root, fullfile('data', 'mu7chain_AC_fast'));
    if ~isempty(src)
        run_dir = fullfile(tmp, 'mu7chain_AC_fast');
        copyfile(src, run_dir);
        B = mu7chain_AC_config('fast', '3cond');
        reg_names = cellfun(@(e) e.name, B.figures, 'UniformOutput', false);
    end
end

n = numel(names);
T = table(names(:), strings(n, 1), nan(n, 1), nan(n, 1), nan(n, 1), nan(n, 1), strings(n, 1), strings(n, 1), ...
    zeros(n, 1), strings(n, 1), 'VariableNames', ...
    {'name', 'status', 'width', 'height', 'png_width', 'png_height', 'sizes', 'exceptions', 'seconds', 'why'});

for j = 1:n
    name = names{j};
    out_dir = fullfile(out_root, name);
    t0 = tic;
    try
        if strcmp(name, 'Fig_Main8_Stim')
            o = draw_fig8(root, ied_dir, out_dir);
        elseif strcmp(name, 'fig_supp10_grouped')
            o = draw_supp10(ied_dir, out_dir);
        else
            if isempty(run_dir); error('test_canvas_standard:Skip', 'no data/mu7chain_AC_fast'); end
            e = B.figures{strcmp(reg_names, name)};
            args = [{'run_dir', run_dir, 'out_dir', out_dir, 'save', true, 'visible', false}, e.args];
            if ~ismember(name, {'fig_main1_grouped', 'fig_supp3_grouped'})   % the two without network conditions
                args = [args, {'conditions', B.conditions}]; %#ok<AGROW>
            end
            if ismember(name, {'fig_main2_grouped', 'fig_main3_grouped'})
                % The archived source figures are read only without the saved
                % trajectories, which the fast run holds: a folder of nothing
                % keeps the main checkout out of it.
                at = find(strcmp(args(1:2:end), 'source_fig_root'));
                args{2 * at} = fullfile(tmp, 'no_source_figs');
            end
            o = e.fn(args{:});
        end
    catch err
        if strcmp(err.identifier, 'test_canvas_standard:Skip')
            T.status(j) = "SKIP"; T.why(j) = err.message;
        else
            T.status(j) = "FAIL"; T.why(j) = "does not build: " + err.message;
        end
        T.seconds(j) = toc(t0);
        close(setdiff(findall(groot, 'Type', 'figure'), figs_before));
        continue
    end
    fig = o.figs(1);
    ex = cell(0, 2);
    if isfield(EXCEPTIONS, name); ex = EXCEPTIONS.(name); end
    why = {};

    % The canvas, measured on the figure as drawn.
    drawnow;
    [w, h] = figure_inches(fig);
    T.width(j) = w; T.height(j) = h;
    if abs(w - st.canvas_w) >= 0.01
        why{end + 1} = sprintf('%.2f in wide, not canvas_w %g', w, st.canvas_w); %#ok<AGROW>
    end
    if h > st.canvas_max_h + 0.01
        why{end + 1} = sprintf('%.2f in tall, more than canvas_max_h %g', h, st.canvas_max_h); %#ok<AGROW>
    end

    % The type by role.
    A = canvas_audit(fig, ex);
    T.sizes(j) = sizes_by_role(A);
    bad = A(~A.ok, :);
    if ~isempty(bad)
        why{end + 1} = sprintf('%d piece(s) of type off their role''s size: %s', height(bad), list_some( ...
            strcat(bad.kind, {' "'}, cellfun(@(s) shorten(s, 50), bad.name, 'UniformOutput', false), {'" '}, ...
            string(bad.size), {' pt ('}, bad.role, {' '}, string(bad.want), {')'}))); %#ok<AGROW>
    end
    used = unique(A.reason(~cellfun(@isempty, A.reason)));
    T.exceptions(j) = strjoin(string(used), '; ');
    for x = 1:size(ex, 1)
        if ~any(strcmp(used, ex{x, 2}))
            why{end + 1} = sprintf('the exception "%s" is not used any more: take it out of EXCEPTIONS', ex{x, 1}); %#ok<AGROW>
        end
    end

    % Axis lines.
    lw = axis_widths(fig);
    if any(abs(lw - st.axis_lw) > 1e-9)
        why{end + 1} = sprintf('axis lines of %s pt (axis_lw %g)', num2str(unique(lw(:))'), st.axis_lw); %#ok<AGROW>
    end

    % Nothing past the canvas.
    out_of = outside(fig);
    if ~isempty(out_of)
        why{end + 1} = sprintf('reaches past the canvas: %s', list_some(out_of)); %#ok<AGROW>
    end

    % The PNG.
    png = the_png(o, out_dir);
    if isempty(png)
        why{end + 1} = 'no PNG written'; %#ok<AGROW>
    else
        info = imfinfo(png);
        dpi = png_dpi(info, fig, name);
        T.png_width(j) = info.Width / dpi; T.png_height(j) = info.Height / dpi;
        if T.png_width(j) < st.canvas_w - st.canvas_tol - 0.01 || T.png_width(j) > st.canvas_w + 0.01
            why{end + 1} = sprintf('the PNG is %.2f in wide (canvas_w %g, tolerance %g)', T.png_width(j), st.canvas_w, st.canvas_tol); %#ok<AGROW>
        end
        if T.png_height(j) > st.canvas_max_h + 0.01
            why{end + 1} = sprintf('the PNG is %.2f in tall (canvas_max_h %g)', T.png_height(j), st.canvas_max_h); %#ok<AGROW>
        end
    end

    if isempty(why); T.status(j) = "PASS"; else; T.status(j) = "FAIL"; T.why(j) = strjoin(why, ' | '); end
    T.seconds(j) = toc(t0);
    fprintf('%-4s %-19s %5.2f x %5.2f in  PNG %5.2f x %5.2f in  (%.0f s)  %s\n', T.status(j), name, w, h, ...
        T.png_width(j), T.png_height(j), T.seconds(j), T.sizes(j));
    close(setdiff(findall(groot, 'Type', 'figure'), figs_before));
end
close(setdiff(findall(groot, 'Type', 'figure'), figs_before));

%% The report
fprintf('\ntest_canvas_standard: canvas_w %g in (tolerance %g), canvas_max_h %g in; type tick %g, label %g, title %g, letter %g pt; axis_lw %g\n', ...
    st.canvas_w, st.canvas_tol, st.canvas_max_h, st.tick_fs, st.label_fs, st.title_fs, st.letter_fs, st.axis_lw);
fprintf('%-19s %-5s %-13s %-13s %s\n', 'figure', '', 'W x H (in)', 'PNG (in)', 'type by role (pt); exceptions');
for j = 1:n
    fprintf('%-19s %-5s %5.2f x %5.2f %5.2f x %5.2f %s', T.name{j}, T.status(j), T.width(j), T.height(j), ...
        T.png_width(j), T.png_height(j), T.sizes(j));
    if strlength(T.exceptions(j)) > 0; fprintf('; exception: %s', T.exceptions(j)); end
    fprintf('\n');
end
not_pass = find(T.status ~= "PASS");
if isempty(not_pass)
    fprintf('PASS test_canvas_standard: all %d figures on the common canvas\n', n);
else
    fprintf('\n%d of %d figures do not pass:\n', numel(not_pass), n);
    for j = not_pass'
        fprintf('  %s %s: %s\n', T.status(j), T.name{j}, T.why(j));
    end
end
if nargout == 0
    clear T
    if ~isempty(not_pass)
        error('test_canvas_standard:Failed', '%d of %d figures not on the canvas standard: %s', numel(not_pass), n, ...
            strjoin(cellstr(names(not_pass)), ', '));
    end
end
end

%% ------------------------------------------------------------------------
function o = draw_fig8(root, ied_dir, out_dir)
% Figure 8 through ied_stim_fig_combined, as the ied_stim bundle draws it
% (invisible; it closes its figure), then its saved .fig reopened invisibly.
if isempty(ied_dir); error('test_canvas_standard:Skip', 'no data/IED_model/ied_stim_fast'); end
c = ied_stim_config('fast', 'run_dir', ied_dir, 'fig_root', fullfile(root, 'figs', 'IED_model', 'ied_stim_fast'));
model_fig = fullfile(c.fig_root, 'fig_ied_stim_blocks', 'fig_ied_stim_blocks.fig');
if ~isfolder(fullfile(c.clinical_dir, 'panels')) || ~isfile(model_fig)
    error('test_canvas_standard:Skip', 'no clinical panels in %s or no model figure %s', c.clinical_dir, model_fig);
end
r = ied_stim_fig_combined('out_dir', out_dir, 'model_fig', model_fig, 'clinical_dir', c.clinical_dir, ...
    'clinical_files', c.clinical_files, 'run_name', c.run_name, 'sim', c.sim, 'burst_mode_dir', ied_dir, ...
    'tag', 'Fig_Main8_Stim', 'visible', false);
fig_file = r.files(endsWith(r.files, '.fig'));
o = struct('figs', openfig(fig_file{1}, 'invisible'), 'files', {r.files});
end

function o = draw_supp10(ied_dir, out_dir)
if isempty(ied_dir); error('test_canvas_standard:Skip', 'no data/IED_model/ied_stim_fast'); end
o = fig_supp10_grouped('run_dir', ied_dir, 'ied_run_dir', ied_dir, ...
    'out_dir', out_dir, 'save', true, 'visible', false);
end

function [w, h] = figure_inches(fig)
% The figure's size in inches as drawn. composite_save and composite_figure
% record the size they asked for ('composite_size', in inches or pixels);
% when the window was shrunk by MATLAB it is put back first.
ppi = get(groot, 'ScreenPixelsPerInch');
sz = getappdata(fig, 'composite_size');
if numel(sz) == 2
    u = fig.Units;
    if all(sz < 100); fig.Units = 'inches'; else; fig.Units = 'pixels'; end
    if any(abs(fig.Position(3:4) - sz) > 1e-6); fig.Position(3:4) = sz; drawnow; end
    fig.Units = u;
end
px = getpixelposition(fig);
w = px(3) / ppi; h = px(4) / ppi;
end

function s = sizes_by_role(A)
roles = {'tick', 'label', 'title', 'letter'};
parts = {};
for r = roles
    v = unique(A.size(strcmp(A.role, r{1})))';
    if isempty(v); continue; end
    parts{end + 1} = sprintf('%s %s', r{1}, strjoin(arrayfun(@(x) sprintf('%g', x), v, 'UniformOutput', false), '/')); %#ok<AGROW>
end
s = string(strjoin(parts, ', '));
end

function lw = axis_widths(fig)
ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes', 'Visible', 'on');
ax = ax(~strcmp(get(ax, {'Tag'}), 'print_whole_check'));
lw = arrayfun(@(a) a.LineWidth, ax);
end

function bad = outside(fig)
% Axes (with their tight inset) and legends that reach past the figure's
% edge by more than a pixel: they would be cut in the export.
bad = {};
fp = getpixelposition(fig);
for a = findall(fig, '-isa', 'matlab.graphics.axis.Axes', 'Visible', 'on')'
    if strcmp(a.Tag, 'print_whole_check'); continue; end
    p = getpixelposition(a, true);
    u = a.Units; a.Units = 'pixels'; ti = a.TightInset; a.Units = u;
    box = [p(1) - ti(1), p(2) - ti(2), p(1) + p(3) + ti(3), p(2) + p(4) + ti(4)];
    if box(1) < -1 || box(2) < -1 || box(3) > fp(3) + 1 || box(4) > fp(4) + 1
        bad{end + 1} = sprintf('axes %s "%s" spans x %.0f to %.0f, y %.0f to %.0f px of %.0f x %.0f', a.Tag, ...
            shorten(a.Title.String, 30), box([1 3 2 4]), fp(3:4)); %#ok<AGROW>
    end
end
for l = findall(fig, 'Type', 'legend', 'Visible', 'on')'
    p = getpixelposition(l, true);
    if p(1) < -1 || p(2) < -1 || p(1) + p(3) > fp(3) + 1 || p(2) + p(4) > fp(4) + 1
        bad{end + 1} = sprintf('legend "%s" spans x %.0f to %.0f px of %.0f', shorten(l.String, 30), ...
            p(1), p(1) + p(3), fp(3)); %#ok<AGROW>
    end
end
end

function png = the_png(o, out_dir)
% The figure's PNG: the one it reports writing, else the one PNG in its
% folder that is not a panel of it.
png = '';
if isfield(o, 'files') && ~isempty(o.files)
    f = cellstr(o.files);
    f = f(endsWith(f, '.png', 'IgnoreCase', true) & ~contains(f, '_panel'));
    for k = 1:numel(f)
        if ~is_absolute(f{k}); f{k} = fullfile(out_dir, f{k}); end
    end
    f = f(cellfun(@isfile, f));
    if ~isempty(f); png = f{1}; return; end
end
d = dir(fullfile(out_dir, '*.png'));
d = d(~contains({d.name}, '_panel'));
if isscalar(d); png = fullfile(d.folder, d.name); end
end

function dpi = png_dpi(info, fig, name)
% The resolution the PNG was written at: its own record where it has one,
% else what the writer chose (fig_grouped_main writes at 200;
% save_some_figs_to_folder_2's png_resolution otherwise).
dpi = NaN;
if isfield(info, 'XResolution') && ~isempty(info.XResolution)
    dpi = info.XResolution;
    if strcmpi(info.ResolutionUnit, 'meter'); dpi = dpi * 0.0254; end
end
if isnan(dpi) || dpi < 10
    if ismember(name, {'fig_main2_grouped', 'fig_main3_grouped'})
        dpi = 200;
    else
        px = getpixelposition(fig);
        dpi = min(600, floor(4000 * get(groot, 'ScreenPixelsPerInch') / max(px(3:4))));
    end
end
end

function d = find_run(root, rel)
% <project root>/<rel>, else the main checkout's; '' when neither has it.
d = fullfile(root, rel);
if isfolder(d); return; end
d = '';
[s, common] = system(sprintf('git -C "%s" rev-parse --path-format=absolute --git-common-dir', root));
if s == 0
    main = fileparts(strtrim(common));
    if isfolder(fullfile(main, rel)); d = fullfile(main, rel); end
end
end

function tf = is_absolute(p)
tf = startsWith(p, {'/', '\'}) || ~isempty(regexp(p, '^[A-Za-z]:', 'once'));
end

function s = shorten(x, n)
if iscell(x) || isstring(x); x = strjoin(cellstr(x), ' '); end
if size(x, 1) > 1; x = strjoin(cellstr(x), ' '); end
s = strtrim(char(x));
if numel(s) > n; s = [s(1:n - 3) '...']; end
end

function s = list_some(c)
c = cellstr(c);
k = numel(c);
s = strjoin(c(1:min(k, 8)), '; ');
if k > 8; s = sprintf('%s; and %d more', s, k - 8); end
end

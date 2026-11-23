%% TEST_CANVAS_STANDARD_B Figures 1, 2, 3, 4, 7 and S1 to S5, S9 are on the common canvas, with their type by role.
%
% Group B of Plan 2's canvas standard (docs Plan2_canvas_standard.md, Tom,
% 2026-11-22: "The figures should have similar canvas widths so that the
% letters appear proportionally the same"). Every composite here is built
% as the fast bundle's registry builds it (mu7chain_AC_config('fast',
% '3cond'), its entry's function and arguments) and checked against the
% numbers of manuscript_style:
%
%   1. the figure is canvas_w wide (to a hundredth of an inch) and at most
%      canvas_max_h tall, measured on the figure after it is drawn (MATLAB
%      shrinks a window larger than the screen, so the size asked for is
%      not trusted);
%   2. type by role, in points: axes (tick labels) tick_fs; axis labels
%      and colour-bar labels label_fs; titles title_fs; legends and free
%      text tick_fs; panel letters letter_fs. Only what is drawn is
%      checked (visible, non-empty). The exceptions allowed are listed by
%      name in EXCEPTIONS below, each with its reason;
%   3. axis lines of every visible axes are axis_lw wide;
%   4. the PNG that composite_save (or fig_grouped_main) writes is within
%      canvas_tol of canvas_w wide, at the resolution it was written at.
%
% THE DATA. The fast run, data/mu7chain_AC_fast of this checkout, else of
% the main checkout; it is COPIED to a temporary folder first, because
% fig_param_space_allStd writes a prep folder into the run directory it is
% given. Without the run the test is SKIPPED (it prints SKIP and passes).
%
% Two environment variables, for drawing the figures at real size by hand:
%   CANVAS_B_OUT    a folder: the PNGs are written to <folder>/<entry>/
%                   instead of a temporary folder, and kept
%   CANVAS_B_ONLY   comma-separated entry names (fig_supp4_grouped,...):
%                   only those are built
% A table of what was measured (width x height, and every size found per
% role) is printed for every figure.
%
% About ten minutes for all eleven. Figures invisible. No simulation.
%
% See also: manuscript_style, composite_tools, composite_save

st = manuscript_style();
n_fail = 0;
figs_before = findall(groot, 'Type', 'figure');

src_run = find_fast_run();
if isempty(src_run)
    fprintf('SKIP: no fast run (data/mu7chain_AC_fast); test_canvas_standard_b built nothing\n');
    return
end
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));
run_dir = fullfile(tmp, 'mu7chain_AC_fast');
copyfile(src_run, run_dir);

out_root = getenv('CANVAS_B_OUT');
if isempty(out_root); out_root = fullfile(tmp, 'figs'); end
entries = {'fig_main1_grouped', 'fig_main2_grouped', 'fig_main3_grouped', 'fig_main4_grouped', ...
    'fig_main7_grouped', 'fig_supp1_grouped', 'fig_supp2_grouped', 'fig_supp3_grouped', ...
    'fig_supp4_grouped', 'fig_supp5_grouped', 'fig_supp9_grouped'};
only = strtrim(strsplit(getenv('CANVAS_B_ONLY'), ','));
only = only(~cellfun(@isempty, only));
if ~isempty(only)
    assert(all(ismember(only, entries)), 'test_canvas_standard_b:Only', 'CANVAS_B_ONLY names an entry not of group B.');
    entries = entries(ismember(entries, only));
end

B = mu7chain_AC_config('fast', '3cond');
names = cellfun(@(e) e.name, B.figures, 'UniformOutput', false);
X = exceptions();
ppi = get(groot, 'ScreenPixelsPerInch');
for j = 1:numel(entries)
    e = B.figures{strcmp(names, entries{j})};
    out_dir = fullfile(out_root, e.name);
    args = [{'run_dir', run_dir, 'out_dir', out_dir, 'save', true, 'visible', false}, e.args];
    if ~strcmp(e.name, 'fig_main1_grouped') && ~strcmp(e.name, 'fig_supp3_grouped')
        args = [args, {'conditions', B.conditions}]; %#ok<AGROW>
    end
    if ismember(e.name, {'fig_main2_grouped', 'fig_main3_grouped'})
        % The source figure root is read only without the saved trajectories,
        % which the fast run holds; a folder of this checkout keeps the
        % main checkout out of it.
        at = find(strcmp(args(1:2:end), 'source_fig_root'));
        args{2 * at} = fullfile(tmp, 'no_source_figs');
    end
    t0 = tic;
    try
        o = e.fn(args{:});
    catch err
        n_fail = n_fail + check(false, sprintf('%s builds: %s', e.name, err.message));
        continue
    end
    fig = o.figs(1);
    drawnow;
    px = getpixelposition(fig);
    w = px(3) / ppi; h = px(4) / ppi;
    n_fail = n_fail + check(abs(w - st.canvas_w) < 0.01, sprintf('%s is %.2f in wide (canvas_w %g)', e.name, w, st.canvas_w));
    n_fail = n_fail + check(h <= st.canvas_max_h + 0.01, sprintf('%s is %.2f in tall (at most %g)', e.name, h, st.canvas_max_h));
    want = getappdata(fig, 'composite_size');
    [bad, seen] = type_by_role(fig, st, X(strcmp({X.entry}, e.name)));
    n_fail = n_fail + check(isempty(bad), sprintf('%s: type by role%s', e.name, list_bad(bad)));
    out_of = outside(fig);
    n_fail = n_fail + check(isempty(out_of), sprintf('%s: nothing reaches past the canvas%s', e.name, list_bad(out_of)));
    lw = axis_widths(fig);
    n_fail = n_fail + check(all(abs(lw - st.axis_lw) < 1e-9), sprintf('%s: axis lines %s pt (axis_lw %g)', ...
        e.name, num2str(unique(lw(:))'), st.axis_lw));
    png = dir(fullfile(out_dir, '*.png'));
    png = png(~contains({png.name}, '_panel'));
    ok = isscalar(png);
    pw = NaN;
    if ok
        info = imfinfo(fullfile(out_dir, png(1).name));
        pw = info.Width / png_dpi(info, px, e.name);
        ok = pw >= st.canvas_w - st.canvas_tol - 0.01 && pw <= st.canvas_w + 0.01;
    end
    n_fail = n_fail + check(ok, sprintf('%s: the PNG is %.2f in wide (within %g of %g)', e.name, pw, st.canvas_tol, st.canvas_w));
    fprintf('  %-18s %5.2f x %5.2f in (asked %s), PNG %.2f in wide, %.0f s; sizes by role: %s\n', e.name, w, h, ...
        mat2str(round(want, 2)), pw, toc(t0), seen);
    close(o.figs(isgraphics(o.figs)));
end
close(setdiff(findall(groot, 'Type', 'figure'), figs_before));
if n_fail > 0
    error('test_canvas_standard_b:Failed', '%d check(s) failed.', n_fail);
end
fprintf('test_canvas_standard_b: all checks PASS\n');

%% ------------------------------------------------------------------------
function X = exceptions()
% Type that may depart from its role's size: entry, a test on the object,
% the size allowed, and why. Every one is named here; nothing else may.
X = struct('entry', {}, 'what', {}, 'match', {}, 'size', {}, 'why', {});
end

function [bad, seen] = type_by_role(fig, st, X)
% Every drawn piece of type against the size of its role. bad: one line per
% departure; seen: the sizes found, per role, as one line.
roles = struct('name', {'tick', 'label', 'title', 'legend', 'text', 'letter'}, ...
    'size', {st.tick_fs, st.label_fs, st.title_fs, st.tick_fs, st.tick_fs, st.letter_fs}, 'found', {[]});
bad = {};
owned = gobjects(0);
ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes');
ax = ax(~strcmp(get(ax, 'Tag'), 'print_whole_check'));
for k = 1:numel(ax)
    a = ax(k);
    parts = {a.Title, a.XLabel, a.YLabel, a.ZLabel};
    if isprop(a, 'Subtitle'); parts{end+1} = a.Subtitle; end %#ok<AGROW>
    for p = 1:numel(parts); owned(end+1) = parts{p}; end %#ok<AGROW>
    if strcmp(a.Visible, 'on') && (~isempty(a.XTick) || ~isempty(a.YTick))
        judge(a, 'tick', a.FontSize, describe(a));
    end
    judge_text(a.Title, 'title', a);
    if isprop(a, 'Subtitle'); judge_text(a.Subtitle, 'title', a); end
    if strcmp(a.Visible, 'on')
        judge_text(a.XLabel, 'label', a); judge_text(a.YLabel, 'label', a);
        if numel(a.ZLabel.String) > 0 && ~isempty(a.ZTick); judge_text(a.ZLabel, 'label', a); end
    end
end
for c = findall(fig, 'Type', 'colorbar')'
    owned(end+1) = c.Label; %#ok<AGROW>
    if strcmp(c.Visible, 'off'); continue; end
    judge(c, 'tick', c.FontSize, 'a colour bar');
    judge_text(c.Label, 'label', c);
end
for l = findall(fig, 'Type', 'legend')'
    if strcmp(l.Visible, 'off'); continue; end
    judge(l, 'legend', l.FontSize, sprintf('the legend "%s"', strjoin(cellstr(l.String), '/')));
    owned(end+1) = l.Title; %#ok<AGROW>
    judge_text(l.Title, 'legend', l);
end
for t = findall(fig, 'Type', 'text')'
    if any(owned == t) || ~isempty(ancestor(t, 'legend')); continue; end
    if strcmp(t.Tag, 'composite_panel_letter') || is_letter(t.String)
        judge_text(t, 'letter', []);
    else
        judge_text(t, 'text', []);
    end
end
for t = findall(fig, 'Type', 'textboxshape', '-or', 'Type', 'textarrowshape')'
    if startsWith(t.Tag, 'panel_letter') || is_letter(t.String)
        judge_text(t, 'letter', []);
    else
        judge_text(t, 'text', []);
    end
end
seen = strjoin(arrayfun(@(r) sprintf('%s %s', r.name, mat2str(unique(r.found(:))')), roles, 'UniformOutput', false), '; ');

    function judge_text(t, role, owner)
        s = t.String;
        if isempty(s) || strcmp(t.Visible, 'off') || all(cellfun(@(x) isempty(strtrim(x)), cellstr(s))); return; end
        if ~isempty(owner); where = sprintf('%s of %s', role, describe(owner)); else; where = role; end
        judge(t, role, t.FontSize, sprintf('%s "%s"', where, short(s)));
    end

    function judge(h, role, fs, what)
        r = strcmp({roles.name}, role);
        roles(r).found(end+1) = fs;
        if abs(fs - roles(r).size) < 1e-6; return; end
        for x = 1:numel(X)
            if X(x).match(h) && abs(fs - X(x).size) < 1e-6; return; end
        end
        bad{end+1} = sprintf('%s at %g pt (%s %g)', what, fs, role, roles(r).size);
    end
end

function tf = is_letter(s)
tf = (ischar(s) || isstring(s)) && ~isempty(regexp(char(s), '^\([A-Z]\)$', 'once'));
end

function s = describe(a)
% An axes by its title, else its labels, else its tag.
s = short(a.Title.String);
if isempty(s) && isprop(a, 'YLabel'); s = short(a.YLabel.String); end
if isempty(s) && isprop(a, 'XLabel'); s = short(a.XLabel.String); end
if isempty(s); s = a.Tag; end
s = sprintf('axes "%s"', s);
end

function s = short(x)
if iscell(x) || isstring(x); x = strjoin(cellstr(x), ' '); end
if size(x, 1) > 1; x = strjoin(cellstr(x), ' '); end
s = strtrim(char(x));
if numel(s) > 40; s = [s(1:37) '...']; end
end

function bad = outside(fig)
% Axes (with their tick labels, axis names and title: the tight inset) and
% legends that reach past the figure's edge by more than a pixel: they
% would be cut in the export.
bad = {};
fp = getpixelposition(fig);
for a = findall(fig, '-isa', 'matlab.graphics.axis.Axes', 'Visible', 'on')'
    if strcmp(a.Tag, 'print_whole_check'); continue; end
    p = getpixelposition(a, true);
    u = a.Units; a.Units = 'pixels'; ti = a.TightInset; a.Units = u;
    box = [p(1) - ti(1), p(2) - ti(2), p(1) + p(3) + ti(3), p(2) + p(4) + ti(4)];
    if box(1) < -1 || box(2) < -1 || box(3) > fp(3) + 1 || box(4) > fp(4) + 1
        bad{end+1} = sprintf('%s spans x %.0f to %.0f, y %.0f to %.0f px of %.0f x %.0f', describe(a), box([1 3 2 4]), fp(3:4)); %#ok<AGROW>
    end
end
for l = findall(fig, 'Type', 'legend', 'Visible', 'on')'
    p = getpixelposition(l, true);
    if p(1) < -1 || p(2) < -1 || p(1) + p(3) > fp(3) + 1 || p(2) + p(4) > fp(4) + 1
        bad{end+1} = sprintf('the legend "%s" spans x %.0f to %.0f px of %.0f', short(l.String), p(1), p(1) + p(3), fp(3)); %#ok<AGROW>
    end
end
end

function lw = axis_widths(fig)
ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes', 'Visible', 'on');
ax = ax(~strcmp(get(ax, 'Tag'), 'print_whole_check'));
lw = arrayfun(@(a) a.LineWidth, ax);
end

function dpi = png_dpi(info, px, name)
% The resolution the PNG was written at: its own record where it has one,
% else what the writer chose (save_some_figs_to_folder_2's png_resolution;
% fig_grouped_main writes at 200).
dpi = NaN;
if isfield(info, 'XResolution') && ~isempty(info.XResolution) && strcmpi(info.ResolutionUnit, 'meter')
    dpi = info.XResolution * 0.0254;
end
if isnan(dpi) || dpi < 10
    if ismember(name, {'fig_main2_grouped', 'fig_main3_grouped'})
        dpi = 200;
    else
        dpi = min(600, floor(4000 * get(groot, 'ScreenPixelsPerInch') / max(px(3:4))));
    end
end
end

function s = list_bad(bad)
if isempty(bad); s = ''; return; end
n = numel(bad);
s = sprintf(' (%d departures): %s', n, strjoin(bad(1:min(n, 25)), '; '));
if n > 25; s = [s sprintf('; and %d more', n - 25)]; end
end

function d = find_fast_run()
% <project root>/data/mu7chain_AC_fast, else the main checkout's.
rel = fullfile('data', 'mu7chain_AC_fast');
root = fileparts(which('setup_paths'));
d = fullfile(root, rel);
if isfolder(d); return; end
[s, common] = system(sprintf('git -C "%s" rev-parse --path-format=absolute --git-common-dir', root));
d = '';
if s == 0
    main = fileparts(strtrim(common));
    if isfolder(fullfile(main, rel)); d = fullfile(main, rel); end
end
end

function f = check(ok, what)
if ok
    fprintf('PASS: %s\n', what); f = 0;
else
    fprintf('FAIL: %s\n', what); f = 1;
end
end

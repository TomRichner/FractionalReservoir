function path = write_figure_report(fig_root)
% WRITE_FIGURE_REPORT One scrollable sheet of every PNG a figure run produced.
%
%   path = WRITE_FIGURE_REPORT(fig_root)   % writes <fig_root>/report.md
%
%   write_figure_report(fullfile(project_root, 'figs', 'single_multi_TS'))
%
% Reviewing a run used to mean opening twenty folders one at a time. This walks
% the directory and emits a single Markdown file that renders in VS Code's
% previewer (Ctrl+Shift+V), or converts to one shareable self-contained page:
%
%   pandoc report.md -o report.html --standalone --embed-resources
%
% IMAGES ARE REFERENCED, NOT COPIED. The report is a few KB whatever the run
% weighs, regenerates instantly, and cannot hold a stale copy of a figure --
% which a report that embedded its images necessarily could.
%
% IT SCANS THE DIRECTORY rather than reading manifest.md, deliberately. The
% manifest records what a run INTENDED and reports; this reports what is on
% disk. Scanning also makes it work on any figs directory, including ones
% predating the manifest and ones whose run died partway -- exactly the cases
% worth looking at. manifest.md is linked, and a few of its header values are
% echoed for convenience, but nothing here depends on parsing it.
%
% SECTION ORDER MATCHES VS CODE'S FILE EXPLORER, which is the whole point of the
% ordering code below: scrolling the report and scrolling the tree agree, so a
% figure seen here can be found there without searching. See natural_order.
%
% Not called by make_all_paper_figures -- run it yourself against a fig_root.
%
% See also: write_figure_manifest, make_all_paper_figures, existing_outputs

arguments
    fig_root (1,:) char
end

if ~isfolder(fig_root)
    error('write_figure_report:NoSuchRoot', ...
        'Figure root does not exist: %s', fig_root);
end

path = fullfile(fig_root, 'report.md');

% Sections are the subdirectories, in explorer order. Root-level PNGs become a
% trailing section because VS Code lists folders before files.
sections = collect_sections(fig_root);
root_pngs = natural_order(list_files(fig_root, '*.png', false));

n_images = sum(cellfun(@(s) numel(s.pngs), sections)) + numel(root_pngs);

fid = fopen(path, 'w');
if fid < 0
    error('write_figure_report:CannotOpen', 'Could not open %s', path);
end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

write_header(fid, fig_root, numel(sections) + ~isempty(root_pngs), n_images);
write_provenance(fid, fig_root);
write_toc(fid, sections, root_pngs);

for k = 1:numel(sections)
    write_section(fid, sections{k}.name, sections{k}.pngs, sections{k}.mds, fig_root);
end
if ~isempty(root_pngs)
    write_section(fid, '(root)', root_pngs, {}, fig_root);
end

% fig_root goes in verbatim: fprintf processes escapes in the FORMAT string
% only, never in a %s argument, so doubling the backslashes here would print
% them doubled.
fprintf(fid, ['\n> Regenerate with `write_figure_report(''%s'')`. `figs/` is ' ...
    'gitignored; this report references the images in place rather than ' ...
    'copying them, so it is only ever as current as its generation stamp ' ...
    'above.\n'], fig_root);
end

%% ------------------------------------------------------------------------
function sections = collect_sections(fig_root)
% One section per subdirectory, in explorer order, each carrying the PNGs found
% anywhere beneath it and the .md files sitting directly in it.
d = dir(fig_root);
names = {d([d.isdir]).name};
names = names(~ismember(names, {'.', '..'}));
names = natural_order(names);

sections = cell(1, numel(names));
for k = 1:numel(names)
    sub = fullfile(fig_root, names{k});
    sections{k} = struct( ...
        'name', names{k}, ...
        'pngs', {natural_order(list_files(sub, '*.png', true))}, ...
        'mds',  {natural_order(list_files(sub, '*.md', false))});
end
end

function files = list_files(base, pattern, recurse)
% Paths RELATIVE TO base, with forward slashes. Markdown links need '/' on every
% platform -- a Windows filesep here renders as a broken image in the previewer.
if recurse
    d = dir(fullfile(base, '**', pattern));
else
    d = dir(fullfile(base, pattern));
end
d = d(~[d.isdir]);
files = cell(1, numel(d));
prefix = [base filesep];
for k = 1:numel(d)
    full = fullfile(d(k).folder, d(k).name);
    % strncmp, not erase: erase strips EVERY occurrence, so a path that repeats
    % the root name deeper down would be mangled.
    if strncmp(full, prefix, numel(prefix))
        rel = full(numel(prefix)+1:end);
    else
        rel = d(k).name;
    end
    files{k} = strrep(rel, '\', '/');
end
end

%% ------------------------------------------------------------------------
function out = natural_order(names)
% Sort as VS CODE'S EXPLORER does, which MATLAB's sort does not.
%
% Three differences, all of which bite here:
%
%   - sort() compares by CODE POINT, so 'fig_EI_param_space' and
%     'fig_SFA_steady_state' land before 'fig_adaptation_methods' ('E' is 69,
%     'a' is 97). The explorer is case-insensitive.
%   - digit runs compare NUMERICALLY there, so fig_panel2 precedes fig_panel10.
%     No folder in the paper's runs contains a digit today, so lower() alone
%     would look correct and silently diverge the first time a figure writes
%     numbered panels. That is why this is here up front.
%   - names differing only in case fall back to a case-sensitive comparison,
%     which keeps the order stable rather than arbitrary.
%
% Folders-before-files is handled by the caller: subdirectories become sections
% and root-level PNGs are appended after all of them.
if isempty(names); out = names; return; end

keys = cell(1, numel(names));
for k = 1:numel(names)
    keys{k} = split_runs(names{k});
end

idx = 1:numel(names);
for a = 2:numel(idx)                       % insertion sort; n is ~20
    j = a;
    while j > 1 && precedes(keys{idx(j)}, names{idx(j)}, keys{idx(j-1)}, names{idx(j-1)})
        idx([j-1 j]) = idx([j j-1]);
        j = j - 1;
    end
end
out = names(idx);
end

function parts = split_runs(name)
% Alternating text/number runs, e.g. 'fig_panel10' -> {'fig_panel', 10}.
tok = regexp(lower(name), '\d+|\D+', 'match');
parts = cell(1, numel(tok));
for k = 1:numel(tok)
    if all(isstrprop(tok{k}, 'digit'))
        parts{k} = str2double(tok{k});
    else
        parts{k} = tok{k};
    end
end
end

function tf = precedes(pa, na, pb, nb)
for k = 1:min(numel(pa), numel(pb))
    a = pa{k}; b = pb{k};
    if isnumeric(a) && isnumeric(b)
        if a ~= b; tf = a < b; return; end
    elseif isnumeric(a) ~= isnumeric(b)
        tf = isnumeric(a);          % a number sorts before text at the same slot
        return;
    else
        if ~strcmp(a, b)
            c = sort({a, b});
            tf = strcmp(c{1}, a);
            return;
        end
    end
end
if numel(pa) ~= numel(pb)
    tf = numel(pa) < numel(pb);
    return;
end
c = sort({na, nb});                 % identical case-insensitively; keep it stable
tf = strcmp(c{1}, na) && ~strcmp(na, nb);
end

%% ------------------------------------------------------------------------
function write_header(fid, fig_root, n_sections, n_images)
fprintf(fid, '# Figure report\n\n');
fprintf(fid, '_Generated %s on %s, MATLAB %s._\n\n', ...
    char(datetime('now')), hostname(), version());
fprintf(fid, '| | |\n|---|---|\n');
fprintf(fid, '| Figure root | `%s` |\n', fig_root);
fprintf(fid, '| Sections | %d |\n', n_sections);
fprintf(fid, '| Images | %d |\n', n_images);
end

function write_provenance(fid, fig_root)
% Tolerant on purpose: every row is optional and a parse failure omits it rather
% than erroring. This has to work on a figs directory that predates the manifest
% or whose run died before writing one.
rows = {};

git_file = fullfile(fig_root, 'git_provenance.txt');
if isfile(git_file)
    txt = fileread(git_file);
    for key = {'commit_short', 'branch', 'dirty'}
        val = regexp(txt, ['^' key{1} ':\s*(\S.*?)\s*$'], 'tokens', 'once', 'lineanchors');
        if ~isempty(val)
            rows(end+1, :) = {key{1}, ['`' strtrim(val{1}) '`']}; %#ok<AGROW>
        end
    end
end

man_file = fullfile(fig_root, 'manifest.md');
if isfile(man_file)
    txt = fileread(man_file);
    for key = {'Run directory', 'Preset', 'Run mode', 'Result', 'Elapsed'}
        val = regexp(txt, ['^\|\s*' key{1} '\s*\|\s*(.*?)\s*\|\s*$'], ...
            'tokens', 'once', 'lineanchors');
        if ~isempty(val)
            rows(end+1, :) = {key{1}, strtrim(val{1})}; %#ok<AGROW>
        end
    end
end

if isempty(rows); fprintf(fid, '\n'); return; end

for k = 1:size(rows, 1)
    fprintf(fid, '| %s | %s |\n', rows{k, 1}, rows{k, 2});
end
if isfile(man_file)
    fprintf(fid, '\n_Values above are echoed from [manifest.md](manifest.md), which is the record of the run._\n');
end
fprintf(fid, '\n');
end

function write_toc(fid, sections, root_pngs)
fprintf(fid, '## Contents\n\n');
for k = 1:numel(sections)
    fprintf(fid, '- [%s](#%s) — %s\n', ...
        sections{k}.name, anchor(sections{k}.name), ...
        count_label(numel(sections{k}.pngs), numel(sections{k}.mds)));
end
if ~isempty(root_pngs)
    fprintf(fid, '- [(root)](#root) — %s\n', count_label(numel(root_pngs), 0));
end
fprintf(fid, '\n');
end

function s = count_label(n_png, n_md)
parts = {};
if n_png > 0; parts{end+1} = sprintf('%d image%s', n_png, plural(n_png)); end
if n_md  > 0; parts{end+1} = sprintf('%d table%s', n_md,  plural(n_md));  end
if isempty(parts); s = '**empty**'; else; s = strjoin(parts, ', '); end
end

function s = plural(n)
if n == 1; s = ''; else; s = 's'; end
end

function write_section(fid, name, pngs, mds, fig_root)
% NOT escaped. The underscores in these names are intraword, which CommonMark
% explicitly does not read as emphasis, so escaping is unnecessary -- and it
% would put a backslash in the rendered heading text that the TOC anchors,
% computed from the raw name, would then fail to match.
fprintf(fid, '## %s\n\n', name);

for k = 1:numel(pngs)
    if strcmp(name, '(root)'); rel = pngs{k}; else; rel = [name '/' pngs{k}]; end
    fprintf(fid, '`%s`\n\n', rel);
    fprintf(fid, '![%s](%s)\n\n', rel, url_encode(rel));
end

for k = 1:numel(mds)
    rel = [name '/' mds{k}];
    fprintf(fid, '_Inlined from [`%s`](%s)._\n\n', rel, url_encode(rel));
    fprintf(fid, '%s\n\n', demote_headings(fileread(fullfile(fig_root, rel))));
end

if isempty(pngs) && isempty(mds)
    % A figure that ran and wrote NOTHING is exactly what this report exists to
    % surface -- this pipeline has more than once reported success for a figure
    % that produced no files.
    fprintf(fid, '> **No PNG or Markdown files in this folder.**\n\n');
end
end

%% ------------------------------------------------------------------------
function txt = demote_headings(txt)
% Push an inlined file's headings down two levels so it nests under this
% report's '##' section instead of restarting the hierarchy at '#'.
%
% Textual, not a Markdown parser: it tracks ``` fences so code samples are left
% alone, but a '#' opening a line inside a $$...$$ math block would be demoted.
% No generated table does that today.
lines = regexp(txt, '\r\n|\r|\n', 'split');
in_fence = false;
for k = 1:numel(lines)
    if ~isempty(regexp(strtrim(lines{k}), '^(```|~~~)', 'once'))
        in_fence = ~in_fence;
        continue
    end
    if in_fence; continue; end
    hashes = regexp(lines{k}, '^(#{1,6})\s', 'tokens', 'once');
    if ~isempty(hashes)
        lines{k} = [repmat('#', 1, min(6, numel(hashes{1}) + 2)) ...
            lines{k}(numel(hashes{1})+1:end)];
    end
end
txt = strjoin(lines, newline);
end

function a = anchor(name)
% GitHub/VS Code heading anchors: lowercase, spaces to hyphens, punctuation
% other than hyphen and underscore dropped.
a = lower(name);
a = regexprep(a, '\s+', '-');
a = regexprep(a, '[^a-z0-9_\-]', '');
end

function s = url_encode(s)
% Only spaces need handling for these paths, and %20 is what a previewer wants.
s = strrep(s, ' ', '%20');
end

function h = hostname()
[st, h] = system('hostname');
if st ~= 0; h = 'unknown'; end
h = strtrim(h);
end

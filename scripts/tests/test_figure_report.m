% TEST_FIGURE_REPORT Checks write_figure_report against a synthetic figs tree.
%
% Aimed at the ORDERING, which is the part that can be wrong without looking
% wrong: the report only earns its keep if scrolling it and scrolling VS Code's
% file tree agree, and MATLAB's sort does not order names the way the explorer
% does. The rest of the checks cover the ways a link can render broken.
%
% Assumes setup_paths has already run. Prints PASS/FAIL per check.

fprintf('\n=== test_figure_report ===\n\n');
n_pass = 0; n_fail = 0;

root = fullfile(tempdir, sprintf('figrep_%s', datestr(now, 'yyyymmdd_HHMMSSFFF'))); %#ok<TNOW1,DATST>
% Cleaned up explicitly at the end rather than with onCleanup: in a SCRIPT that
% has local functions, the destructor cannot resolve an anonymous handle back to
% the script's scope and warns instead of deleting. Leaving the tree behind on
% failure is the useful behaviour anyway -- the report can then be inspected.

% --- a tree that separates the three ways sort() differs from the explorer ---
% 'Bravo' vs 'alpha' is the CASE difference (code point puts 'B' first);
% panel2/panel10 is the NUMERIC one; 'fig_Panel3' checks the two interact.
mkdir(fullfile(root, 'Bravo'));
mkdir(fullfile(root, 'alpha'));
mkdir(fullfile(root, 'fig_panel2'));
mkdir(fullfile(root, 'fig_panel10'));
mkdir(fullfile(root, 'fig_Panel3'));
mkdir(fullfile(root, 'empty_folder'));
mkdir(fullfile(root, 'nested', 'inner'));
mkdir(fullfile(root, 'tables'));

write_png(fullfile(root, 'Bravo', 'b.png'));
write_png(fullfile(root, 'alpha', 'a.png'));
write_png(fullfile(root, 'fig_panel2', 'p.png'));
write_png(fullfile(root, 'fig_panel10', 'p.png'));
write_png(fullfile(root, 'fig_Panel3', 'p.png'));
write_png(fullfile(root, 'nested', 'inner', 'deep.png'));
fid = fopen(fullfile(root, 'tables', 'tbl.md'), 'w');
fprintf(fid, '# Top\n\n## Sub\n\n```\n# not a heading\n```\n\ntext\n');
fclose(fid);

path = write_figure_report(root);
txt = fileread(path);
lines = regexp(txt, '\r\n|\r|\n', 'split');

% --- 1. section order matches VS Code's explorer -------------------------
% [^\r\n]+ rather than .+ : MATLAB's '.' matches newlines by DEFAULT, unlike
% most regex flavours, so '^## (.+)$' with lineanchors swallows the whole file
% from the first heading onwards.
heads = section_headings(txt);
expected = {'alpha', 'Bravo', 'empty_folder', 'fig_panel2', 'fig_Panel3', ...
    'fig_panel10', 'nested', 'tables'};
[n_pass, n_fail] = check(isequal(heads, expected), ...
    'section order matches the explorer', n_pass, n_fail);
if ~isequal(heads, expected)
    fprintf('       got: %s\n', strjoin(heads, ', '));
    fprintf('  expected: %s\n', strjoin(expected, ', '));
end

% --- 2. a folder that produced nothing is called out --------------------
[n_pass, n_fail] = check(contains(txt, '**No PNG or Markdown files in this folder.**'), ...
    'empty folder is flagged, not skipped', n_pass, n_fail);

% --- 3. nested images link with forward slashes -------------------------
% A Windows filesep here renders as a broken image in the previewer.
[n_pass, n_fail] = check(contains(txt, '(nested/inner/deep.png)'), ...
    'nested image links use forward slashes', n_pass, n_fail);
[n_pass, n_fail] = check(~contains(txt, '\inner\'), ...
    'no backslashes in image links', n_pass, n_fail);

% --- 4. inlined markdown is demoted, code fences untouched --------------
[n_pass, n_fail] = check(contains(txt, '### Top') && contains(txt, '#### Sub'), ...
    'inlined headings demoted two levels', n_pass, n_fail);
[n_pass, n_fail] = check(contains(txt, sprintf('\n# not a heading')), ...
    'heading inside a code fence left alone', n_pass, n_fail);

% --- 5. every linked path exists ----------------------------------------
links = regexp(txt, '^!\[[^\]]*\]\(([^)]+)\)$', 'tokens', 'lineanchors');
links = cellfun(@(c) strrep(c{1}, '%20', ' '), links, 'UniformOutput', false);
all_exist = ~isempty(links) && all(cellfun(@(p) isfile(fullfile(root, p)), links));
[n_pass, n_fail] = check(all_exist, ...
    sprintf('all %d image links resolve to files', numel(links)), n_pass, n_fail);
[n_pass, n_fail] = check(numel(links) == 6, ...
    'one link per PNG on disk (6)', n_pass, n_fail);

% --- 6. TOC anchors match the headings they point at --------------------
toc = regexp(txt, '^- \[([^\]]+)\]\(#([^)]+)\)', 'tokens', 'lineanchors');
ok_anchor = true;
for k = 1:numel(toc)
    want = lower(regexprep(toc{k}{1}, '\s+', '-'));
    want = regexprep(want, '[^a-z0-9_\-]', '');
    ok_anchor = ok_anchor && strcmp(want, toc{k}{2});
end
[n_pass, n_fail] = check(ok_anchor && numel(toc) == numel(expected), ...
    'TOC has one correct anchor per section', n_pass, n_fail);

% --- 7. works with no manifest.md and no git_provenance.txt -------------
% The interesting figs directories are often the broken ones; the report must
% still be produced for a run that died before writing a manifest.
[n_pass, n_fail] = check(isfile(path) && ~isempty(txt), ...
    'writes with no manifest or provenance present', n_pass, n_fail);

% --- 8. root-level PNGs come last, matching folders-before-files --------
write_png(fullfile(root, 'loose.png'));
h2 = section_headings(fileread(write_figure_report(root)));
[n_pass, n_fail] = check(strcmp(h2{end}, '(root)'), ...
    'root-level PNGs form the last section', n_pass, n_fail);

% --- 9. a missing root errors -------------------------------------------
try
    write_figure_report(fullfile(root, 'does_not_exist'));
    [n_pass, n_fail] = check(false, 'missing root errors', n_pass, n_fail);
catch err
    [n_pass, n_fail] = check(strcmp(err.identifier, 'write_figure_report:NoSuchRoot'), ...
        'missing root errors with NoSuchRoot', n_pass, n_fail);
end

fprintf('\n--- %d passed, %d failed ---\n', n_pass, n_fail);
if n_fail > 0
    fprintf('Temp tree left for inspection: %s\n', root);
    error('test_figure_report:Failed', '%d check(s) failed.', n_fail);
end
rmdir(root, 's');
fprintf('ALL PASS\n\n');

%% ------------------------------------------------------------------------
function [n_pass, n_fail] = check(tf, label, n_pass, n_fail)
if tf
    fprintf('  PASS  %s\n', label);
    n_pass = n_pass + 1;
else
    fprintf('  FAIL  %s\n', label);
    n_fail = n_fail + 1;
end
end

function write_png(p)
imwrite(uint8(zeros(4, 4, 3)), p);
end

function heads = section_headings(txt)
h = regexp(txt, '^## ([^\r\n]+)$', 'tokens', 'lineanchors');
heads = cellfun(@(c) c{1}, h, 'UniformOutput', false);
heads = heads(~strcmp(heads, 'Contents'));
end

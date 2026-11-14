function M = read_figure_manifest(fig_root)
% READ_FIGURE_MANIFEST The per-entry results a figure pass recorded in <fig_root>/manifest.md.
%
%   M = READ_FIGURE_MANIFEST(fig_root)
%
% The inverse of write_figure_manifest, for the one caller that needs it: a
% partial pass of make_all_paper_figures (cfg.only), which redraws some
% entries of a figure set and must write a manifest for the WHOLE set. The
% rows of the entries it did not touch come from here, exactly as the
% earlier pass recorded them.
%
%   M.results    struct array with the fields make_all_paper_figures builds:
%                name, ok, in_paper, n_files, seconds, files, err,
%                conditions_applied. err is the manifest's one-line text
%                (cut at 120 characters there); seconds has the manifest's
%                one decimal.
%   M.generated  the text of the "Generated ..." line, without the markup
%   M.notes      the lines of the "Partial passes" section (cellstr, empty
%                when the manifest records a full pass only)
%
% Errors when there is no manifest or it holds no entry row: a partial pass
% has nothing to add to then, and the caller should make a full pass.
%
% See also: write_figure_manifest, make_all_paper_figures

arguments
    fig_root (1,:) char
end

path = fullfile(fig_root, 'manifest.md');
if ~isfile(path)
    error('read_figure_manifest:NoManifest', ...
        'No manifest.md in %s. A partial pass adds to a full one; make the full pass first.', fig_root);
end
lines = regexp(fileread(path), '\r?\n', 'split');

M = struct('results', [], 'generated', '', 'notes', {{}});
g = regexp(lines, '^_Generated (.*?)\.?_$', 'tokens', 'once');
g = g(~cellfun(@isempty, g));
if ~isempty(g); M.generated = g{1}{1}; end

results = struct('name', {}, 'ok', {}, 'in_paper', {}, 'n_files', {}, ...
                 'seconds', {}, 'files', {}, 'err', {}, 'conditions_applied', {});
row = '^\| `([^`]+)` \| (yes|--) \| (\d+) \| ([\d.]+) \| (.*?) \| (.*) \|$';
section = '';
current = 0;                                   % the entry whose files are being listed
for k = 1:numel(lines)
    ln = lines{k};
    if startsWith(ln, '## ')
        section = strtrim(ln(4:end)); current = 0; continue;
    end
    switch section
        case 'Entries'
            t = regexp(ln, row, 'tokens', 'once');
            if isempty(t); continue; end
            ok = strcmp(t{6}, 'ok');
            err = '';
            if ~ok; err = regexprep(t{6}, '^\*\*FAILED\*\* -- ', ''); end
            if strcmp(t{5}, 'yes'); applied = 'yes';
            elseif startsWith(t{5}, '**no**'); applied = 'no';
            else; applied = '';
            end
            results(end+1) = struct('name', t{1}, 'ok', ok, 'in_paper', strcmp(t{2}, 'yes'), ...
                'n_files', str2double(t{3}), 'seconds', str2double(t{4}), 'files', {{}}, ...
                'err', err, 'conditions_applied', applied); %#ok<AGROW>
        case 'Files'
            t = regexp(ln, '^- `(.*)/`$', 'tokens', 'once');
            if ~isempty(t)
                current = find(strcmp({results.name}, t{1}), 1);
                if isempty(current); current = 0; end
                continue;
            end
            t = regexp(ln, '^  - `(.*)`$', 'tokens', 'once');
            if ~isempty(t) && current > 0
                results(current).files{end+1} = t{1};
            end
        case 'Partial passes'
            t = regexp(ln, '^- (.*)$', 'tokens', 'once');
            if ~isempty(t); M.notes{end+1} = t{1}; end
    end
end
if isempty(results)
    error('read_figure_manifest:NoEntries', 'The manifest %s lists no entry.', path);
end
M.results = results;
end

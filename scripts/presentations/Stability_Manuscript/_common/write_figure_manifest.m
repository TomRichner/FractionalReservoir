function path = write_figure_manifest(fig_root, cfg, run_dir, root_mode, results, minutes)
% WRITE_FIGURE_MANIFEST One record of what a figure run produced, at the root.
%
%   path = WRITE_FIGURE_MANIFEST(fig_root, cfg, run_dir, root_mode, results, minutes)
%
% Writes <fig_root>/manifest.md. Replaces the per-figure README_*.txt and the
% fifteen near-identical git_provenance.txt files that used to sit one per
% figure folder.
%
% THE LINE THAT MATTERS IS run_dir. Two figures used to read a .mat frozen beside
% their own source while the other sixteen read the run they were handed; on
% 2026-08-26 that put Aug 22 data in a sheet whose siblings were all Aug 25, and
% nothing in fifteen provenance files said so, because each recorded the COMMIT
% and none recorded the DATA. A manifest naming one run directory is what makes
% that visible.
%
% Written after the loop regardless of failures, so a partial run still describes
% itself -- a manifest listing what did NOT build is more useful than no manifest.
%
% Git state (commit, branch, dirty flag, working_changes.patch) is captured
% separately by capture_git_provenance, called once at this same root.
%
% See also: make_all_paper_figures, capture_git_provenance, paper_config

arguments
    fig_root  (1,:) char
    cfg       struct
    run_dir   (1,:) char
    root_mode (1,:) char
    results   struct
    minutes   (1,1) double
end

path = fullfile(fig_root, 'manifest.md');
fid = fopen(path, 'w');
if fid < 0
    error('write_figure_manifest:CannotOpen', 'Could not open %s', path);
end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

n_ok    = sum([results.ok]);
n_paper = sum([results.in_paper]);
n_pok   = sum([results.ok] & [results.in_paper]);

fprintf(fid, '# Paper figures\n\n');
fprintf(fid, '_Generated %s on %s, MATLAB %s._\n\n', ...
    char(datetime('now')), hostname(), version());

fprintf(fid, '| | |\n|---|---|\n');
fprintf(fid, '| Run directory | `%s` |\n', or_none(run_dir));
fprintf(fid, '| Preset | `%s` |\n', cfg.preset_name);
fprintf(fid, '| Run mode | `%s` |\n', cfg.run_mode);
fprintf(fid, '| Figure root | `%s` (%s) |\n', fig_root, root_mode);
fprintf(fid, '| Result | %d/%d succeeded, %d/%d in-paper |\n', n_ok, numel(results), n_pok, n_paper);
fprintf(fid, '| Elapsed | %.1f min |\n\n', minutes);

if isempty(run_dir)
    fprintf(fid, ['> **No run directory resolved.** Replot figures will have ' ...
        'failed; analytic and self-simulating ones will not.\n\n']);
end

fprintf(fid, '## Entries\n\n');
fprintf(fid, '| Entry | In paper | Files | Seconds | Result |\n|---|---|---:|---:|---|\n');
for k = 1:numel(results)
    r = results(k);
    if r.ok; verdict = 'ok'; else; verdict = ['**FAILED** -- ' one_line(r.err)]; end
    fprintf(fid, '| `%s` | %s | %d | %.1f | %s |\n', ...
        r.name, tick(r.in_paper), r.n_files, r.seconds, verdict);
end
fprintf(fid, '\n');

% Name the files, so the manifest says what is actually on disk rather than what
% was intended. existing_outputs feeds these, and it globs the directory.
fprintf(fid, '## Files\n\n');
for k = 1:numel(results)
    if isempty(results(k).files); continue; end
    fprintf(fid, '- `%s/`\n', results(k).name);
    fprintf(fid, '  - `%s`\n', results(k).files{:});
end
fprintf(fid, '\n');

fprintf(fid, ['> Figures are regenerable: this directory is gitignored, and ' ...
    '`make_all_paper_figures` rebuilds it from the run directory and commit ' ...
    'named above. Git state is in `git_provenance.txt` beside this file.\n']);
end

%% ------------------------------------------------------------------------
function s = tick(tf)
if tf; s = 'yes'; else; s = '--'; end
end

function s = or_none(s)
if isempty(s); s = '(none resolved)'; end
end

function s = one_line(s)
% Error messages can be multi-line; a table cell cannot.
if isempty(s); s = ''; return; end
s = regexprep(char(s), '\s+', ' ');
if numel(s) > 120; s = [s(1:117) '...']; end
end

function h = hostname()
[st, h] = system('hostname');
if st ~= 0; h = 'unknown'; end
h = strtrim(h);
end

function path = write_figure_manifest(fig_root, cfg, run_dir, root_mode, results, minutes, notes)
% WRITE_FIGURE_MANIFEST One record of what a figure run produced, at the root.
%
%   path = WRITE_FIGURE_MANIFEST(fig_root, cfg, run_dir, root_mode, results, minutes)
%   path = WRITE_FIGURE_MANIFEST(..., notes)     % a partial pass; see below
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
% THE CONDITIONS (FD-006). A "Conditions" row gives cfg.conditions, the
% selection requested for the pass (empty = all the run holds), so a _3cond and
% a _4cond folder built from one run say which they are; the entries table says
% per entry whether the selection was applied (results.conditions_applied), and
% a note names every entry that drew all conditions because its function does
% not declare the option.
%
% A PARTIAL PASS (make_all_paper_figures with cfg.only) writes the manifest of
% the whole set: `results` then holds every entry, the redrawn ones from this
% pass and the others as read_figure_manifest read them from the earlier
% manifest, and `notes` holds one line per pass, written as the section
% "Partial passes" (which entries, when, where that pass's git state is).
% The "Elapsed" row is then the partial pass's own time.
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
    notes     cell = {}     % a partial pass's record, one line each; see the header
end

path = fullfile(fig_root, 'manifest.md');
fid = fopen(path, 'w');
if fid < 0
    error('write_figure_manifest:CannotOpen', 'Could not open %s', path);
end
closer = onCleanup(@() fclose(fid));

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
fprintf(fid, '| Conditions | %s |\n', conditions_text(cfg));
fprintf(fid, '| Figure root | `%s` (%s) |\n', fig_root, root_mode);
fprintf(fid, '| Result | %d/%d succeeded, %d/%d in-paper |\n', n_ok, numel(results), n_pok, n_paper);
fprintf(fid, '| Elapsed | %.1f min |\n\n', minutes);

if isempty(run_dir)
    fprintf(fid, ['> **No run directory resolved.** Replot figures will have ' ...
        'failed; analytic and self-simulating ones will not.\n\n']);
end

% Per entry, whether the condition selection was applied: "yes", "**no**"
% (the function does not declare 'conditions' and drew every condition), or
% "--" (nothing selected, or the entry failed first). Results from callers
% predating the field read as "--".
has_sel = isfield(results, 'conditions_applied');
not_applied = {};
fprintf(fid, '## Entries\n\n');
fprintf(fid, '| Entry | In paper | Files | Seconds | Conditions applied | Result |\n|---|---|---:|---:|---|---|\n');
for k = 1:numel(results)
    r = results(k);
    if r.ok; verdict = 'ok'; else; verdict = ['**FAILED** -- ' one_line(r.err)]; end
    ca = '--';
    if has_sel
        switch r.conditions_applied
            case 'yes', ca = 'yes';
            case 'no',  ca = '**no** (all drawn)'; not_applied{end+1} = r.name; %#ok<AGROW>
        end
    end
    fprintf(fid, '| `%s` | %s | %d | %.1f | %s | %s |\n', ...
        r.name, tick(r.in_paper), r.n_files, r.seconds, ca, verdict);
end
fprintf(fid, '\n');
if ~isempty(not_applied)
    fprintf(fid, ['> **The condition selection was NOT applied to %d entr%s** (the function ' ...
        'does not declare `conditions`, so it drew every condition the run holds): %s.\n\n'], ...
        numel(not_applied), plural_y(numel(not_applied)), ...
        strjoin(strcat('`', not_applied, '`'), ', '));
end

% Name the files, so the manifest says what is actually on disk rather than what
% was intended. existing_outputs feeds these, and it globs the directory.
fprintf(fid, '## Files\n\n');
for k = 1:numel(results)
    if isempty(results(k).files); continue; end
    fprintf(fid, '- `%s/`\n', results(k).name);
    fprintf(fid, '  - `%s`\n', results(k).files{:});
end
fprintf(fid, '\n');

% read_figure_manifest reads this section back by its heading and its "- "
% lines, so a later partial pass carries the earlier ones forward.
if ~isempty(notes)
    fprintf(fid, '## Partial passes\n\n');
    fprintf(fid, ['Some entries were redrawn after the full pass (`make_all_paper_figures` with `cfg.only`). ' ...
        'The rows of the other entries are as the earlier manifest recorded them; "Elapsed" above is the last pass''s.\n\n']);
    fprintf(fid, '- %s\n', notes{:});
    fprintf(fid, '\n');
end

fprintf(fid, ['> Figures are regenerable: this directory is gitignored, and ' ...
    '`make_all_paper_figures` rebuilds it from the run directory and commit ' ...
    'named above. Git state is in `git_provenance.txt` beside this file.\n']);
end

%% ------------------------------------------------------------------------
function s = tick(tf)
if tf; s = 'yes'; else; s = '--'; end
end

function s = conditions_text(cfg)
% The selection requested for the pass (cfg.conditions; empty = all).
c = {};
if isfield(cfg, 'conditions') && ~isempty(cfg.conditions)
    c = cellstr(cfg.conditions);
end
if isempty(c) || isequal(c, {''})
    s = 'all the run holds (none selected)';
else
    s = strjoin(strcat('`', c, '`'), ', ');
    s = sprintf('%s (in this order)', s);
end
if isfield(cfg, 'figure_set') && ~isempty(cfg.figure_set)
    s = sprintf('%s; figure set `%s`', s, cfg.figure_set);
end
end

function s = plural_y(k)
if k == 1; s = 'y'; else; s = 'ies'; end
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

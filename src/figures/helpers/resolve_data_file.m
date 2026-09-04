function mat_file = resolve_data_file(explicit, run_dir, searched, pattern, hint)
% RESOLVE_DATA_FILE The newest match for a data file, from INSIDE one run directory.
%
%   f = RESOLVE_DATA_FILE(explicit, run_dir, searched, pattern, hint)
%
%   explicit  a caller-supplied path to the file itself, or '' to search
%   run_dir   the run directory. REQUIRED when explicit is '': the search is
%             confined to it, and every entry of `searched` must lie within it
%   searched  cellstr of directories to try, MOST SPECIFIC FIRST, all under run_dir
%   pattern   'eig_heatmap_data.mat', or a glob like '*_results.mat'
%   hint      what to run to produce the file, quoted in the error
%
% ONE RUN, OR AN EXPLICIT FILE. NOTHING ELSE.
%
% This used to search the run directory and THEN fall back to a standalone
% location under data/. That fallback caused the failure it was meant to prevent.
% On 2026-09-03 a single_multi_TS run lost its memory_capacity stage (an
% unsupported run_mode), so fig_memory_capacity fell through to
% data/memory_capacity/ and plotted a .mat from 2026-08-22 -- a different
% network, twelve days old -- and reported SUCCESS with three files written. A
% stale figure was indistinguishable from a good one, and the run's manifest said
% nothing.
%
% The fallback was never a separate capability, which is what makes removing it
% cheap: `<project_root>/data/memory_capacity` is just `<run_dir>/memory_capacity`
% with run_dir = <project_root>/data. So a standalone analysis is still plottable,
% by SAYING SO:
%
%   fig_memory_capacity('run_dir', fullfile(project_root, 'data'))   % standalone
%   fig_memory_capacity('mat_file', '/path/to/that_one.mat')         % one file
%
% Same answer as before; the difference is that the caller chose it. A figure that
% silently reaches outside the run it was handed cannot notice it is doing so.
%
% EVERY SEARCHED DIRECTORY MUST BE UNDER RUN_DIR, and that is enforced rather
% than documented: a future edit re-adding a data/ tier fails loudly here instead
% of quietly plotting another experiment.
%
% RUN_DIR MUST BE NON-EMPTY. The previous version tried to guard this by dropping
% empty entries from `searched`, which never fired: fullfile('', 'memory_capacity')
% is 'memory_capacity', not '', so an unresolved run silently searched a RELATIVE
% path against whatever the current directory happened to be.
%
% See also: fig_memory_capacity, fig_memory_capacity_example, fig_eig_heatmap,
%           fig_dc_lle, resolve_run_dir

arguments
    explicit (1,:) char
    run_dir  (1,:) char
    searched (1,:) cell
    pattern  (1,:) char
    hint     (1,:) char
end

%% An explicit file wins, and is the escape hatch for anything unusual
if ~isempty(explicit)
    if ~isfile(explicit)
        error('resolve_data_file:NoSuchFile', ...
            'The file given explicitly does not exist:\n  %s', explicit);
    end
    mat_file = explicit;
    return
end

%% Otherwise a run directory is required
if isempty(strtrim(run_dir))
    error('resolve_data_file:NoRunDir', ...
        ['No run directory, so there is nothing to search for ''%s''.\n\n' ...
         'This figure reads data produced by an analysis stage, and will only ' ...
         'look inside ONE run directory -- it will not fall back to a ' ...
         'standalone location, because doing so once made a figure plot ' ...
         'twelve-day-old data from a different network and report success.\n\n' ...
         'Give it one of:\n' ...
         '  ''run_dir'', <a run directory>\n' ...
         '  ''run_dir'', fullfile(project_root, ''data'')  (a standalone analysis)\n' ...
         '  the file itself, via the mat_file / data_file argument'], pattern);
end
if ~isfolder(run_dir)
    error('resolve_data_file:NoSuchRunDir', ...
        'The run directory does not exist:\n  %s', run_dir);
end

%% The search may not leave the run directory
searched = searched(~cellfun(@isempty, searched));
outside  = searched(~cellfun(@(d) is_within(d, run_dir), searched));
if ~isempty(outside)
    error('resolve_data_file:OutsideRunDir', ...
        ['These search directories are not inside the run directory:\n%s\n' ...
         'run_dir: %s\n\n' ...
         'A figure must read one run. If you meant the standalone location, ' ...
         'pass it AS run_dir rather than adding it as a second tier.'], ...
        sprintf('    %s\n', outside{:}), run_dir);
end

for k = 1:numel(searched)
    hits = dir(fullfile(searched{k}, pattern));
    hits = hits(~[hits.isdir]);
    if ~isempty(hits)
        [~, newest] = max([hits.datenum]);
        mat_file = fullfile(hits(newest).folder, hits(newest).name);
        return
    end
end

error('resolve_data_file:NotFound', ...
    ['No file matching ''%s'' in this run. Looked in:\n%s' ...
     '%s\n' ...
     'The search does NOT extend outside the run directory, so a stale copy ' ...
     'elsewhere cannot be picked up by accident. Point run_dir at the run that ' ...
     'has this stage, or pass the file explicitly.\n' ...
     '(the .mat files are gitignored, so a fresh clone has none).'], ...
    pattern, sprintf('    %s\n', searched{:}), hint);
end

%% ------------------------------------------------------------------------
function tf = is_within(child, parent)
% Is `child` at or below `parent`? Compared as normalised text, not by resolving
% the filesystem: the directory need not exist yet, and a missing stage folder is
% a normal outcome here rather than an error.
c = normalise(child);
p = normalise(parent);
tf = strcmp(c, p) || startsWith(c, [p filesep]);
end

function p = normalise(p)
p = strrep(p, '/', filesep);
p = strrep(p, '\', filesep);
while ~isempty(p) && endsWith(p, filesep)
    p = p(1:end-1);
end
if ispc; p = lower(p); end     % NTFS is case-insensitive; comparisons must be too
end

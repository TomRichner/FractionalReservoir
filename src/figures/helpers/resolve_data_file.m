function mat_file = resolve_data_file(explicit, searched, pattern, hint)
% RESOLVE_DATA_FILE The newest match for a data file, from the most specific place that has one.
%
%   f = RESOLVE_DATA_FILE(explicit, searched, pattern, hint)
%
%   explicit  a caller-supplied path, or '' to search
%   searched  cellstr of directories, MOST SPECIFIC FIRST (typically the run
%             directory's stage folder, then the standalone data/ location)
%   pattern   'eig_heatmap_data.mat', or a glob like '*_results.mat'
%   hint      what to run to produce the file, quoted in the error
%
% Errors naming EVERY directory it looked in, so a miss says where to put the
% file rather than only that it is absent.
%
% WHY THIS EXISTS. fig_eig_heatmap and fig_memory_capacity_example each loaded a
% single hardcoded path -- a .mat sitting beside their own .m -- and declared
% cfg.run_dir "unused". In the 2026-08-26 figure run they read data from Aug 22
% while the other sixteen figures used the Aug 25 sweep, and every folder's
% git_provenance.txt claimed the same commit, so nothing flagged it. A figure
% that reads one fixed location cannot notice that the run it was handed
% contains something newer.
%
% Extracted from fig_memory_capacity's local resolve_mc_results, which had this
% right, so the three now share one implementation rather than three copies
% drifting apart.
%
% See also: fig_eig_heatmap, fig_memory_capacity, fig_memory_capacity_example

arguments
    explicit (1,:) char
    searched (1,:) cell
    pattern  (1,:) char
    hint     (1,:) char
end

if ~isempty(explicit)
    if ~isfile(explicit)
        error('resolve_data_file:NoSuchFile', ...
            'The file given explicitly does not exist:\n  %s', explicit);
    end
    mat_file = explicit;
    return
end

% Drop empties: run_dir is '' when no run could be resolved, and an empty
% fullfile would otherwise search the current directory by accident.
searched = searched(~cellfun(@isempty, searched));

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
    ['No file matching ''%s'' found. Looked in:\n%s' ...
     '%s\n' ...
     '(the .mat files are gitignored, so a fresh clone has none).'], ...
    pattern, sprintf('    %s\n', searched{:}), hint);
end

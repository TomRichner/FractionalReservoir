function B = bm_load_runs(run_dir, index_file, which_runs)
%BM_LOAD_RUNS What the burst_mode stage wrote into one run directory: its index, and the networks asked for.
%
%   B = bm_load_runs(run_dir)                 index, and the first network
%   B = bm_load_runs(run_dir, '', 'all')      index, and every network
%   B = bm_load_runs('', index_file, ...)     an explicit index file
%
% B.runs      the index table (run_burst_mode_stage: one row per network)
% B.settings  the stage's plan without its runs
% B.folder    the folder holding the files
% B.R         1 x k cell of bm_run_network structs, in index order: the
%             first network only (default, which_runs 'first') or all of
%             them ('all')
% B.index_file
%
% Errors through resolve_data_file when the stage has not run in run_dir
% (no fallback to another run), and when a network of the index failed.
%
% See also: run_burst_mode_stage, fig_burst_mode, fig_supp10_grouped

arguments
    run_dir (1,:) char = ''
    index_file (1,:) char = ''
    which_runs (1,:) char {mustBeMember(which_runs, {'first', 'all'})} = 'first'
end
searched = {};
if ~isempty(run_dir); searched = {fullfile(run_dir, 'burst_mode')}; end
f = resolve_data_file(index_file, run_dir, searched, 'burst_mode_index.mat', ...
    'the opt-in stage burst_mode (run_all_paper_analyses with ''burst_mode'' in cfg.stages)');
L = load(f, 'runs', 'settings');
if ~all(L.runs.ok)
    error('bm_load_runs:FailedRuns', 'The burst_mode index %s lists failed network(s): %s', f, ...
        strjoin(L.runs.file(~L.runs.ok)', ', '));
end
B = struct('runs', L.runs, 'settings', L.settings, 'folder', fileparts(f), 'index_file', f);
n = height(L.runs);
if strcmp(which_runs, 'first'); n = 1; end
B.R = cell(1, n);
for k = 1:n
    S = load(fullfile(B.folder, L.runs.file{k}), 'R');
    B.R{k} = S.R;
end
end

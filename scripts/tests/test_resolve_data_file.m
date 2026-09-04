% TEST_RESOLVE_DATA_FILE A figure reads ONE run directory, or errors.
%
% THE BUG THIS GUARDS. resolve_data_file used to search the run directory and
% then fall back to a standalone location under data/. On 2026-09-03 a
% single_multi_TS run lost its memory_capacity stage to an unsupported run_mode,
% so fig_memory_capacity fell through to data/memory_capacity/ and plotted a
% .mat from 2026-08-22 -- a different network, twelve days old -- and reported
% SUCCESS with three files written. Nothing distinguished that figure from a good
% one, and the manifest recorded no problem.
%
% The fallback was never a separate capability: <root>/data/memory_capacity is
% just <run_dir>/memory_capacity with run_dir = <root>/data. So it is removed
% rather than warned about, and a standalone analysis is plotted by naming its
% location AS run_dir, or by passing the file itself.
%
% Three properties are asserted here, all of which failed in some form before:
%   1. a missing file ERRORS instead of reaching outside the run
%   2. an empty run_dir errors, rather than searching a RELATIVE path against the
%      current directory -- which the old "drop empty entries" guard did not
%      prevent, since fullfile('', 'x') is 'x', not ''
%   3. a search directory outside run_dir is refused, so the fallback cannot be
%      reintroduced by a later edit
%
% Uses a temporary directory tree; touches nothing in the repo.
%
% See also: resolve_data_file, fig_memory_capacity, fig_eig_heatmap

all_passed = true;

root = fullfile(tempdir, 'rdf_test');
if isfolder(root); rmdir(root, 's'); end
run_a   = fullfile(root, 'run_A');
run_b   = fullfile(root, 'run_B');
stage_a = fullfile(run_a, 'memory_capacity');
stage_b = fullfile(run_b, 'memory_capacity');
mkdir(stage_a); mkdir(stage_b);

% run_A has the file; run_B's stage exists but is EMPTY -- a stage that failed.
file_a = fullfile(stage_a, 'MC_run_A_results.mat');
x = 1; save(file_a, 'x');

fprintf('\n-- finds the file inside the run it is given --\n');
got = resolve_data_file('', run_a, {stage_a}, '*_results.mat', 'run it');
all_passed = check('run_A resolves to run_A''s file', strcmp(got, file_a)) && all_passed;

fprintf('\n-- a run whose stage produced nothing ERRORS --\n');
% THE CENTRAL CASE. run_B has no results, and another run next door does. The old
% code would have found run_A's (via the data/ tier); this must error instead.
all_passed = check('run_B (empty stage) errors rather than borrowing run_A''s', ...
    throws_id(@() resolve_data_file('', run_b, {stage_b}, '*_results.mat', 'run it'), ...
              'resolve_data_file:NotFound')) && all_passed;

fprintf('\n-- an empty run_dir errors --\n');
% Was silently searching a RELATIVE path against cwd: fullfile('', 'x') is 'x'.
all_passed = check('empty run_dir -> NoRunDir', ...
    throws_id(@() resolve_data_file('', '', {fullfile('', 'memory_capacity')}, ...
                                    '*_results.mat', 'run it'), ...
              'resolve_data_file:NoRunDir')) && all_passed;
all_passed = check('whitespace run_dir -> NoRunDir', ...
    throws_id(@() resolve_data_file('', '   ', {stage_a}, '*_results.mat', 'run it'), ...
              'resolve_data_file:NoRunDir')) && all_passed;
all_passed = check('a nonexistent run_dir -> NoSuchRunDir', ...
    throws_id(@() resolve_data_file('', fullfile(root, 'nope'), {stage_a}, ...
                                    '*_results.mat', 'run it'), ...
              'resolve_data_file:NoSuchRunDir')) && all_passed;

fprintf('\n-- the fallback cannot be reintroduced --\n');
% Naming run_A's stage while claiming to search run_B is exactly the shape of the
% old data/ tier, and is refused.
all_passed = check('a search dir outside run_dir -> OutsideRunDir', ...
    throws_id(@() resolve_data_file('', run_b, {stage_b, stage_a}, ...
                                    '*_results.mat', 'run it'), ...
              'resolve_data_file:OutsideRunDir')) && all_passed;
% run_dir itself, and nested subfolders, are inside and must be accepted --
% fig_dc_lle legitimately walks <run_dir>/dc_lle/dc_lle_nSeeds_*.
all_passed = check('run_dir itself counts as inside', ...
    strcmp(resolve_data_file('', run_a, {run_a, stage_a}, '*_results.mat', 'run it'), ...
           file_a)) && all_passed;

fprintf('\n-- an explicit file still wins, and is checked --\n');
all_passed = check('explicit path is used as-is', ...
    strcmp(resolve_data_file(file_a, '', {}, '*_results.mat', 'run it'), file_a)) ...
    && all_passed;
all_passed = check('a nonexistent explicit path -> NoSuchFile', ...
    throws_id(@() resolve_data_file(fullfile(root, 'ghost.mat'), run_a, {stage_a}, ...
                                    '*_results.mat', 'run it'), ...
              'resolve_data_file:NoSuchFile')) && all_passed;

fprintf('\n-- the standalone location works when NAMED as run_dir --\n');
% This is what replaces the fallback: same directory, chosen explicitly.
data_root  = fullfile(root, 'data');
data_stage = fullfile(data_root, 'memory_capacity');
mkdir(data_stage);
file_d = fullfile(data_stage, 'MC_standalone_results.mat');
save(file_d, 'x');
all_passed = check('run_dir = <data> finds the standalone file', ...
    strcmp(resolve_data_file('', data_root, {data_stage}, '*_results.mat', 'run it'), ...
           file_d)) && all_passed;

rmdir(root, 's');

%% ------------------------------------------------------------------------
if all_passed
    fprintf('\ntest_resolve_data_file: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_resolve_data_file: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function tf = throws_id(fn, want)
try
    fn();
    tf = false;
catch ME
    tf = strcmp(ME.identifier, want);
    if ~tf; fprintf(2, '       got %s\n', ME.identifier); end
end
end

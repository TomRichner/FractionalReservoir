% TEST_POOL_ENTRY_POINTS The two entry points refuse to start with no pool allocation.
%
% run_all_paper_analyses and make_all_paper_figures check pool_config().is_set
% before creating anything (FD-017): a session that was never told how many
% pool workers it may hold (pool_set(N), or pool_set(0)) gets an error with
% the entry point's own identifier, and no run or figure directory.
%
% The pool module itself (pool_set, pool_ensure, pool_restart, pool_release,
% pool_status) lives in the shared MatlabAgenticOrg repository since
% 2026-11-19, with its own test: run('<MatlabAgenticOrg>/tests/test_pool_module.m').
% These checks were part of that test while the module lived here.
%
% It opens and closes no pool. The session's allocation is unset for the
% two calls and put back afterwards.
%
% See also: run_all_paper_analyses, make_all_paper_figures, pool_config, use_matlab_agentic_org

c = pool_config();
saved_workers = getenv(c.env_workers);
restore = onCleanup(@() setenv(c.env_workers, saved_workers));
all_passed = true;

copies = which('pool_set', '-all');
all_passed = check('exactly one pool_set is on the path, and it is not in this repository', ...
    isscalar(copies) && ~startsWith(copies{1}, fileparts(which('setup_paths')))) && all_passed;

setenv(c.env_workers, '');
root = fileparts(which('setup_paths'));
probe_run = fullfile(root, 'data', '__pool_probe__');
probe_fig = fullfile(root, 'figs', '__pool_probe__');
all_passed = check('run_all_paper_analyses refuses (NoPoolAllocation)', ...
    throws_id(@() run_all_paper_analyses(paper_config('run_mode', 'fast', 'run_dir', 'data/__pool_probe__')), ...
    'run_all_paper_analyses:NoPoolAllocation')) && all_passed;
all_passed = check('...and creates no run directory', ~isfolder(probe_run)) && all_passed;
all_passed = check('make_all_paper_figures refuses (NoPoolAllocation)', ...
    throws_id(@() make_all_paper_figures(paper_config('fig_root', 'figs/__pool_probe__')), ...
    'make_all_paper_figures:NoPoolAllocation')) && all_passed;
all_passed = check('...and creates no figure directory', ~isfolder(probe_fig)) && all_passed;
if isfolder(probe_run); rmdir(probe_run, 's'); end
if isfolder(probe_fig); rmdir(probe_fig, 's'); end

setenv(c.env_workers, saved_workers);
all_passed = check('the session''s allocation is back as it was found', ...
    strcmp(getenv(c.env_workers), saved_workers)) && all_passed;

if all_passed
    fprintf('\ntest_pool_entry_points: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_pool_entry_points: FAILURES ABOVE\n');
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
    if ~tf; fprintf(2, '     got %s: %s\n', ME.identifier, ME.message); end
end
end

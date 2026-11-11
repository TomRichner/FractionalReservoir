function [ok, diffs, info] = regression_check(opts)
% REGRESSION_CHECK Rerun the FULL regression bundle and compare it with the saved reference (about 50 minutes).
%
% This is the full tier: every analysis stage and every figure of the
% manuscript's bundle, on the 100-neuron network. For the check to run after
% every edit (about 20 s, the three classes only) see regression_rapid_check.
%
%   ok = regression_check()                         % full: analyses and figures
%   ok = regression_check('figures', false)         % analyses only (faster)
%   ok = regression_check('tag', 'existing_run')    % compare a run already on disk
%   [ok, diffs, info] = regression_check('keep', true, 'rtol', 1e-9)
%
% The reference is scripts/tests/fixtures/regression_fast_digest.mat, a
% digest (regression_digest) of the bundle run on the code as it stood when
% the reference was made. This runs the bundle again under a new tag, digests
% it, and compares (regression_compare). Use it before and after any change
% to the model, the analysis stages or the figures:
%
%   - OK true: nothing the pipeline saves has moved.
%   - OK false: DIFFS lists what moved. Either the change was meant to move
%     it (then say so, and make a new reference deliberately with
%     make_regression_reference), or it was not (then it is a regression).
%
% Options
%   tag      name of the run under data/regression_fast/. Default: a new
%            name from the clock ('check_<yyyymmdd_HHMMSS>'), which is run.
%            If the named run directory already exists it is NOT rerun; it
%            is digested as it is.
%   figures  also run make_all_paper_figures and compare the figure file
%            list (default true). Ignored when comparing an existing run.
%   keep     keep the run directory and its figures afterwards. Default
%            false: a passing check deletes them (they are several tens of
%            megabytes). A FAILING check always keeps them for inspection.
%   rtol, atol, ignore   passed to regression_compare.
%
% INFO has the run directory, the wall time in minutes [analyses, figures],
% the digest of this run, and the figure results.
%
% See also: regression_fast_run, regression_digest, regression_compare,
%           make_regression_reference, test_regression_bundle
arguments
    opts.tag (1,:) char = ''
    opts.figures (1,1) logical = true
    opts.keep (1,1) logical = false
    opts.rtol (1,1) double = 0
    opts.atol (1,1) double = 0
    opts.ignore cell = {}
end
project_root = fileparts(which('setup_paths'));
ref_file = fullfile(project_root, 'scripts', 'tests', 'fixtures', 'regression_fast_digest.mat');
if ~isfile(ref_file)
    error('regression_check:NoReference', ...
        ['No reference at %s.\nMake one from a clean checkout with ' ...
         'make_regression_reference(''full'').'], ref_file);
end
R = load(ref_file);
ref = R.digest;

tag = opts.tag;
if isempty(tag); tag = ['check_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]; end
cfg = regression_fast_config(tag);
run_dir = fullfile(project_root, cfg.run_dir);
fig_root = fullfile(project_root, cfg.fig_root);

info = struct('run_dir', run_dir, 'fig_root', fig_root, 'minutes', [0 0], ...
    'ran', false, 'fig_results', [], 'digest', []);
if isfolder(run_dir) && numel(dir(run_dir)) > 2
    fprintf('regression_check: using the existing run %s (not rerun)\n', run_dir);
else
    fprintf('regression_check: running the regression bundle as tag ''%s''\n', tag);
    [~, info.fig_results, info.minutes] = regression_fast_run(tag, 'figures', opts.figures);
    info.ran = true;
end

if isfolder(fig_root)
    new = regression_digest(run_dir, 'fig_root', fig_root);
else
    new = regression_digest(run_dir);
end
info.digest = new;

fprintf('\nReference: %s\n  made %s from commit %s\n', ref_file, R.made_on, R.commit);
[ok, diffs] = regression_compare(ref, new, 'rtol', opts.rtol, 'atol', opts.atol, 'ignore', opts.ignore);

% The outcome also goes to a small text file beside the run directories, so
% a run launched through the MATLAB MCP server can be followed by reading a
% file: a call that stays silent for 30 minutes is dropped by the client
% although the run carries on (2026-11-10).
info.result_file = fullfile(fileparts(run_dir), [tag '.result.txt']);
live = diffs(~diffs.ignored, :);
fid = fopen(info.result_file, 'w');
if fid > 0
    closer = onCleanup(@() fclose(fid));
    if ok; verdict = 'PASS'; else; verdict = 'DIFFERENCES'; end
    fprintf(fid, 'regression_check %s\n', verdict);
    fprintf(fid, 'tag %s\nrun_dir %s\n', tag, run_dir);
    fprintf(fid, 'reference commit %s made %s\n', R.commit, R.made_on);
    fprintf(fid, 'minutes analyses %.1f figures %.1f\n', info.minutes(1), info.minutes(2));
    if ~isempty(info.fig_results)
        fprintf(fid, 'figures ok %d of %d\n', sum([info.fig_results.ok]), numel(info.fig_results));
    end
    fprintf(fid, 'differences %d (ignored %d)\n', height(live), nnz(diffs.ignored));
    for k = 1:min(height(live), 200)
        fprintf(fid, '  %-7s %s | %s\n', live.kind{k}, live.key{k}, live.detail{k});
    end
end

if ok && info.ran && ~opts.keep
    rmdir(run_dir, 's');
    if isfolder(fig_root); rmdir(fig_root, 's'); end
    fprintf('regression_check: PASS; removed %s\n', run_dir);
elseif ok
    fprintf('regression_check: PASS; kept %s\n', run_dir);
else
    fprintf('regression_check: DIFFERENCES FOUND; kept %s for inspection\n', run_dir);
end
fprintf('regression_check: outcome written to %s\n', info.result_file);
end

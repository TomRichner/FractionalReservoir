function [ok, diffs, info] = regression_rapid_chain_check(opts)
% REGRESSION_RAPID_CHAIN_CHECK Run the chain's ultrarapid regression and compare it with its saved reference (about half a minute).
%
%   ok = regression_rapid_chain_check()
%   [ok, diffs, info] = regression_rapid_chain_check('rtol', 1e-9, 'keep', true)
%
% regression_rapid_check's twin for networks with chained vesicle-pool
% depression (regression_rapid_chain_run): chain networks with noise off
% and on, a four-cell-type network with a chain beside every other kind of
% route, the three Lyapunov estimators, a clearly chaotic chain network,
% transient gain, a small sweep and a short memory-capacity run. Its
% reference is scripts/tests/fixtures/regression_rapid_chain_digest.mat.
% Run it, with regression_rapid_check, after every edit to
% SRNNCellTypePairs, the route models, SRNN_ESN_reservoir or
% ParamSpaceAnalysis2 that should move no number.
%
%   - OK true: nothing those classes produce on a chain network here has
%     moved, to the last bit (exact mode, the default).
%   - OK false: DIFFS lists what moved. If the change was meant to move it,
%     make a new reference deliberately
%     (make_regression_reference('rapid_chain')) and say why in the commit
%     message. A rewrite of the chain's arithmetic is gated by the
%     tolerance tier's chain suite (regression_tol_check('chain')), not by
%     this.
%
% Options: rtol, atol, ignore (passed to regression_compare); keep (keep the
% output folder; by default it is removed when the check passes).
%
% See also: regression_rapid_chain_run, regression_rapid_check,
%           regression_tol_check, make_regression_reference
arguments
    opts.rtol (1,1) double = 0
    opts.atol (1,1) double = 0
    opts.ignore cell = {}
    opts.keep (1,1) logical = false
end
project_root = fileparts(which('setup_paths'));
ref_file = fullfile(project_root, 'scripts', 'tests', 'fixtures', 'regression_rapid_chain_digest.mat');
if ~isfile(ref_file)
    error('regression_rapid_chain_check:NoReference', ...
        'No reference at %s.\nMake one from a clean checkout with make_regression_reference(''rapid_chain'').', ref_file);
end
R = load(ref_file);

[out_dir, seconds] = regression_rapid_chain_run();
new = regression_digest(out_dir, 'verbose', false);
fprintf('Reference: %s\n  made %s from commit %s\n', ref_file, R.made_on, R.commit);
[ok, diffs] = regression_compare(R.digest, new, 'rtol', opts.rtol, 'atol', opts.atol, 'ignore', opts.ignore);
info = struct('out_dir', out_dir, 'seconds', seconds, 'digest', new);

if ok && ~opts.keep
    rmdir(out_dir, 's');
    fprintf('regression_rapid_chain_check: PASS (%.0f s)\n', seconds);
elseif ok
    fprintf('regression_rapid_chain_check: PASS (%.0f s); kept %s\n', seconds, out_dir);
else
    fprintf('regression_rapid_chain_check: DIFFERENCES FOUND; kept %s for inspection\n', out_dir);
end
end

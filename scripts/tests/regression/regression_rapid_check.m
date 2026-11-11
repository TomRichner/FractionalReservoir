function [ok, diffs, info] = regression_rapid_check(opts)
% REGRESSION_RAPID_CHECK Run the ultrarapid regression and compare it with its saved reference (about 20 s).
%
%   ok = regression_rapid_check()
%   [ok, diffs, info] = regression_rapid_check('rtol', 1e-9, 'keep', true)
%
% The ultrarapid tier (regression_rapid_run): the model in the three
% conditions with noise off and on, a small sweep, and a short
% memory-capacity run. Its reference is
% scripts/tests/fixtures/regression_rapid_digest.mat. Run this after every
% edit to SRNNCellTypePairs, SRNN_ESN_reservoir or ParamSpaceAnalysis2.
%
%   - OK true: nothing those three classes produce here has moved, to the
%     last bit (exact mode, the default).
%   - OK false: DIFFS lists what moved. If the change was meant to move it,
%     make a new reference deliberately (make_regression_reference('rapid'))
%     and say why in the commit message.
%
% It does not exercise the analysis stages or the figures; the full tier
% (regression_check) does, and takes about 50 minutes.
%
% Options: rtol, atol, ignore (passed to regression_compare); keep (keep the
% output folder; by default it is removed when the check passes).
%
% See also: regression_rapid_run, regression_check, make_regression_reference
arguments
    opts.rtol (1,1) double = 0
    opts.atol (1,1) double = 0
    opts.ignore cell = {}
    opts.keep (1,1) logical = false
end
project_root = fileparts(which('setup_paths'));
ref_file = fullfile(project_root, 'scripts', 'tests', 'fixtures', 'regression_rapid_digest.mat');
if ~isfile(ref_file)
    error('regression_rapid_check:NoReference', ...
        'No reference at %s.\nMake one from a clean checkout with make_regression_reference(''rapid'').', ref_file);
end
R = load(ref_file);

[out_dir, seconds] = regression_rapid_run();
new = regression_digest(out_dir, 'verbose', false);
fprintf('Reference: %s\n  made %s from commit %s\n', ref_file, R.made_on, R.commit);
[ok, diffs] = regression_compare(R.digest, new, 'rtol', opts.rtol, 'atol', opts.atol, 'ignore', opts.ignore);
info = struct('out_dir', out_dir, 'seconds', seconds, 'digest', new);

if ok && ~opts.keep
    rmdir(out_dir, 's');
    fprintf('regression_rapid_check: PASS (%.0f s)\n', seconds);
elseif ok
    fprintf('regression_rapid_check: PASS (%.0f s); kept %s\n', seconds, out_dir);
else
    fprintf('regression_rapid_check: DIFFERENCES FOUND; kept %s for inspection\n', out_dir);
end
end

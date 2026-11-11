function [ok, diffs, info] = regression_stages_check(opts)
% REGRESSION_STAGES_CHECK Run the middle regression tier and compare it with its saved reference (a few minutes).
%
%   ok = regression_stages_check()
%   [ok, diffs, info] = regression_stages_check('rtol', 1e-9, 'keep', true)
%
% The middle tier (regression_stages_run): three analysis stages called
% directly on the 100-neuron regression network, and the figures that read
% them. Its reference is scripts/tests/fixtures/regression_stages_digest.mat.
% Run it before committing a change that touches an analysis stage, a figure,
% or anything the stages call in the model class.
%
%   - OK true: nothing those stages save has moved, and the same figure files
%     were written.
%   - OK false: DIFFS lists what moved. If the change was meant to move it,
%     make a new reference deliberately (make_regression_reference('stages'))
%     and say why in the commit message.
%
% Options: rtol, atol, ignore (passed to regression_compare); keep (keep the
% output folder; by default it is removed when the check passes).
%
% See also: regression_stages_run, regression_rapid_check, regression_check,
%           make_regression_reference
arguments
    opts.rtol (1,1) double = 0
    opts.atol (1,1) double = 0
    opts.ignore cell = {}
    opts.keep (1,1) logical = false
end
project_root = fileparts(which('setup_paths'));
ref_file = fullfile(project_root, 'scripts', 'tests', 'fixtures', 'regression_stages_digest.mat');
if ~isfile(ref_file)
    error('regression_stages_check:NoReference', ...
        'No reference at %s.\nMake one from a clean checkout with make_regression_reference(''stages'').', ref_file);
end
R = load(ref_file);

[out_dir, seconds] = regression_stages_run();
new = regression_digest(out_dir, 'verbose', false);
fprintf('Reference: %s\n  made %s from commit %s\n', ref_file, R.made_on, R.commit);
[ok, diffs] = regression_compare(R.digest, new, 'rtol', opts.rtol, 'atol', opts.atol, 'ignore', opts.ignore);
info = struct('out_dir', out_dir, 'seconds', seconds, 'digest', new);

if ok && ~opts.keep
    rmdir(out_dir, 's');
    fprintf('regression_stages_check: PASS (%.0f s)\n', seconds);
elseif ok
    fprintf('regression_stages_check: PASS (%.0f s); kept %s\n', seconds, out_dir);
else
    fprintf('regression_stages_check: DIFFERENCES FOUND; kept %s for inspection\n', out_dir);
end
end

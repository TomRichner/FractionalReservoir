function [ok, leaves, info] = regression_tol_check(name, opts)
% REGRESSION_TOL_CHECK Run the tolerance regression tier and compare it with a named reference.
%
%   ok = regression_tol_check()                     % against the anchor
%   ok = regression_tol_check('before_jacobian_edit')
%   [ok, leaves, info] = regression_tol_check('anchor', 'rtol', 0)
%
% The gate for a change to the model's ARITHMETIC (FD-020; Plan A, sections
% 2.2 and 6): the right-hand side, both Jacobians and the packed physics at
% stored states (layer 1, tight), and small stable networks end to end
% (layer 2, looser). About 20 s, in the client; no pool is needed and none
% is used. Run it after every edit to SRNNCellTypePairs,
% SRNN_ESN_reservoir or ParamSpaceAnalysis2 that is NOT expected to be
% bit-identical; for edits that should move no number the byte-exact
% regression_rapid_check is the stricter gate.
%
%   NAME   'anchor' (default): scripts/tests/fixtures/regression_tol_anchor.mat,
%          the model class before Plan A touched it. Any other name: a
%          pre-edit reference made by make_regression_tol_reference(NAME),
%          in data/regression_tol/.
%
%   - OK true: every leaf of the reference is within tolerance. The
%     printout says how close the worst leaf of each layer was.
%   - OK false: LEAVES lists what moved, by how much and where. THAT IS A
%     DEFECT TO FIND. It is never a reason to loosen a tolerance, and the
%     anchor is never re-captured to make this pass.
%
% The pointwise layer is evaluated at the states STORED IN THE REFERENCE,
% so nothing chaotic is compared along a trajectory.
%
% Options: rtol (a scalar replacing every tolerance; 0 demands exact
% equality, which holds for the same code on the same thread count),
% settings (default regression_tol_settings()), verbose.
%
% A case that throws does not stop the check: its error is printed, its
% leaves count as missing, and OK is false (see regression_tol_run).
%
% INFO has the reference's meta, the fresh run (info.data), its duration
% and the per-layer summary of regression_tol_compare.
%
% See also: regression_tol_run, regression_tol_compare,
%           make_regression_tol_reference, regression_tol_settings,
%           regression_tol_band_check
arguments
    name (1,:) char = 'anchor'
    opts.rtol double {mustBeScalarOrEmpty, mustBeNonnegative} = []
    opts.settings (1,1) struct = regression_tol_settings()
    opts.verbose (1,1) logical = true
end
ref_file = regression_tol_reference_file(name);
if ~isfile(ref_file)
    error('regression_tol_check:NoReference', ...
        'No reference at %s.\nMake one with make_regression_tol_reference(''%s'').', ref_file, name);
end
ref = load(ref_file);

[data, seconds] = regression_tol_run(ref.data);

if opts.verbose
    fprintf('Reference: %s\n  made %s from commit %s%s, %d thread(s), MATLAB %s\n', ref_file, ...
        ref.meta.made_on, ref.meta.commit(1:min(7, end)), ...
        repmat(' (dirty tree)', 1, ~isempty(ref.meta.dirty)), ref.meta.threads, ref.meta.matlab_release);
    if data.info.threads ~= ref.meta.threads
        fprintf('  this run uses %d thread(s): expect round-off differences, not exact equality\n', ...
            data.info.threads);
    end
    if ~strcmp(version, ref.meta.matlab_release)
        fprintf('  this is MATLAB %s: library round-off may differ\n', version);
    end
end
[ok, leaves, summary] = regression_tol_compare(ref.data, data, 'settings', opts.settings, ...
    'rtol', opts.rtol, 'verbose', opts.verbose);
info = struct('ref_file', ref_file, 'meta', ref.meta, 'data', data, 'seconds', seconds, ...
    'summary', summary);
% Cases that threw. Their reference leaves already count as missing; this
% says why, and keeps OK false even for a case the reference does not hold.
if ~isempty(data.info.errors)
    ok = false;
    if opts.verbose
        fprintf('  %d case(s) did not run:\n', numel(data.info.errors));
        fprintf('    ERROR   %s\n', data.info.errors{:});
    end
end
if opts.verbose
    if ok; verdict = 'PASS'; else; verdict = 'DIFFERENCES FOUND'; end
    fprintf('regression_tol_check(''%s''): %s (%.0f s)\n', name, verdict, seconds);
end
end

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
%          the model class before Plan A touched it. 'chain':
%          .../fixtures/regression_tol_chain.mat, networks with chained
%          vesicle-pool depression as first merged (the gate for work on
%          route_pool: with facilitation off, the chain is unchanged).
%          'stf_chain': .../fixtures/regression_tol_stf_chain.mat, networks
%          with release-coupled facilitation as first built (the gate for
%          the per-neuron spread and any later edit to it). 'spread':
%          .../fixtures/regression_tol_spread.mat, networks whose
%          depression and facilitation parameters differ per neuron
%          (std_spread, stf_spread) as first built. Any
%          other name: a pre-edit reference made by
%          make_regression_tol_reference(NAME), in data/regression_tol/.
%          The reference says which suite of cases it holds
%          (regression_tol_run: 'base', 'chain', 'stf_chain' or 'spread')
%          and that suite is run; 'anchor', 'chain', 'stf_chain' and
%          'spread' are run one after the other by
%          test_regression_tolerance, test_regression_tolerance_chain,
%          test_regression_tolerance_stf_chain and
%          test_regression_tolerance_spread.
%
%   - OK true: every leaf of the reference is within tolerance. The
%     printout says how close the worst leaf of each layer was.
%   - OK false: LEAVES lists what moved, by how much and where. THAT IS A
%     DEFECT TO LOOK FOR FIRST. A tolerance is not adjusted, and the
%     anchor not re-captured, while something may be wrong. A tolerance
%     set once is not absolute, though (TR, 2026-11-11): when a thorough
%     investigation explains why it cannot be met and shows the code is
%     right, the lead adjusts it in regression_tol_settings and writes
%     the cause, the evidence and the remaining margin beside the value.
%
% SUPERSEDED LEAVES OF THE ANCHOR. The anchor is a historical reference: it
% is never re-captured, so it keeps tying the code to the model class as it
% was before Plan A for every leaf that has not moved. When a deliberate
% change of physics moves some of its leaves, those leaves are listed in
% regression_tol_settings().anchor_superseded ({pattern, reason} rows),
% and ONLY when NAME is 'anchor' they are taken out of pass / fail: they
% get status 'superseded' in LEAVES (their differences still computed and
% kept in the table), the printout says "N leaves superseded, not
% compared: <reason>", and info.superseded lists each pattern with its
% reason, count and keys. No tolerance changes and no other leaf is
% treated differently. The suites 'chain', 'stf_chain' and 'spread' and
% every pre-edit capture are compared in full. The first entry is F3
% (commit cce78c9, 2026-11-12): with one SFA timescale a neuron's time
% constant reads the fast-end draw of its per-neuron ladder, which moved
% the 176 leaves of the condition sfa1_std1 that depend on it; those cases
% are guarded since by the exact ultrarapid tier, re-captured at bdf211d.
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
% Which suite of cases the reference holds (regression_tol_run). The anchor
% and references made before the chain suite existed do not say: 'base'.
suite = 'base';
if isfield(ref.meta, 'suite'); suite = ref.meta.suite; end

[data, seconds] = regression_tol_run(ref.data, suite);

if opts.verbose
    fprintf('Reference: %s\n  suite ''%s'', made %s from commit %s%s, %d thread(s), MATLAB %s\n', ref_file, ...
        suite, ref.meta.made_on, ref.meta.commit(1:min(7, end)), ...
        repmat(' (dirty tree)', 1, ~isempty(ref.meta.dirty)), ref.meta.threads, ref.meta.matlab_release);
    if data.info.threads ~= ref.meta.threads
        fprintf('  this run uses %d thread(s): expect round-off differences, not exact equality\n', ...
            data.info.threads);
    end
    if ~strcmp(version, ref.meta.matlab_release)
        fprintf('  this is MATLAB %s: library round-off may differ\n', version);
    end
end
% Superseded leaves apply to the anchor alone (see the header).
superseded = cell(0, 2);
if strcmp(name, 'anchor') && isfield(opts.settings, 'anchor_superseded')
    superseded = opts.settings.anchor_superseded;
end
[ok, leaves, summary] = regression_tol_compare(ref.data, data, 'settings', opts.settings, ...
    'rtol', opts.rtol, 'verbose', opts.verbose, 'superseded', superseded);
sup = struct('pattern', superseded(:, 1), 'reason', superseded(:, 2), 'n', [], 'keys', []);
left = strcmp(leaves.status, 'superseded');     % first matching pattern wins, as in the compare
for k = 1:numel(sup)
    hit = left & ~cellfun(@isempty, regexp(leaves.key, sup(k).pattern, 'once'));
    sup(k).keys = leaves.key(hit);
    sup(k).n = nnz(hit);
    left = left & ~hit;
end
info = struct('ref_file', ref_file, 'meta', ref.meta, 'data', data, 'seconds', seconds, ...
    'summary', summary);
info.superseded = sup;
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

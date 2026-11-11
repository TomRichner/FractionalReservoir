% test_regression_bundle.m
%
% Reruns the regression bundle (the manuscript's pipeline on a 100-neuron
% network; regression_fast_config) and compares everything it saves with the
% reference in scripts/tests/fixtures/regression_fast_digest.mat.
%
% SLOW: it runs every analysis stage and every figure. It is on
% run_all_tests' slow list. Run it before and after a change to the model
% class, an analysis stage or a figure:
%
%   run_all_tests('test_regression_bundle')
%
% or, with more control (analyses only, keep the run, a tolerance):
%
%   [ok, diffs, info] = regression_check('figures', false, 'keep', true);
%
% If it fails because a change was MEANT to move results, the reference is
% replaced deliberately with make_regression_reference, and the commit
% message says why. See docs/ISSUES.md, FD-003.

[ok, diffs, info] = regression_check();

if info.ran
    fprintf('wall time: analyses %.1f min, figures %.1f min\n', info.minutes(1), info.minutes(2));
    if ~isempty(info.fig_results)
        n_ok = sum([info.fig_results.ok]);
        fprintf('figures: %d of %d succeeded\n', n_ok, numel(info.fig_results));
        assert(n_ok == numel(info.fig_results), 'test_regression_bundle:FigureFailed', ...
            '%d figure(s) failed to produce their files.', numel(info.fig_results) - n_ok);
    end
end
assert(ok, 'test_regression_bundle:Differences', ...
    '%d difference(s) from the regression reference; see the table above and regression_check.', ...
    nnz(~diffs.ignored));
fprintf('PASS  regression bundle matches the reference\n');

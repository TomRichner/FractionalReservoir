% test_regression_tol_superseded.m
%
% Checks that the anchor's SUPERSEDED leaves (regression_tol_settings().
% anchor_superseded, applied by regression_tol_check to the anchor only)
% hide nothing they were not written for. A few seconds; runs nothing.
%
%   (a) On the anchor's own key list, the patterns match only leaves of the
%       condition sfa1_std1, and exactly the 176 that F3 (commit cce78c9,
%       2026-11-12) moved. Comparing the anchor with itself flattens its
%       keys the way regression_tol_check does.
%   (b) On two small structs, a leaf that differs and lies OUTSIDE the
%       patterns still fails, including a leaf of sfa1_std1 that the
%       patterns leave compared (the stored state S), while the matching
%       leaf is reported as superseded and does not fail.
%   (c) Without the option (every reference other than the anchor), the
%       same differences all fail.

T = regression_tol_settings();
sup = T.anchor_superseded;
assert(iscell(sup) && size(sup, 2) == 2 && ~isempty(sup), ...
    'anchor_superseded must be a non-empty n x 2 cell {pattern, reason}');

%% (a) the patterns on the anchor's keys
ref = load(fullfile(fileparts(which('setup_paths')), 'scripts', 'tests', 'fixtures', ...
    'regression_tol_anchor.mat'));
[ok_a, L] = regression_tol_compare(ref.data, ref.data, 'verbose', false, 'superseded', sup);
assert(ok_a, 'the anchor compared with itself must pass');
hit = strcmp(L.status, 'superseded');
assert(all(contains(L.key(hit), '.sfa1_std1')), ...
    'a superseded pattern matches a leaf outside sfa1_std1: %s', ...
    strjoin(L.key(hit & ~contains(L.key, '.sfa1_std1')), ', '));
assert(nnz(hit) == 176, 'the superseded patterns match %d anchor leaves, expected 176 (F3)', nnz(hit));
fprintf('PASS  (a) %d of %d anchor leaves superseded, all in sfa1_std1 (%d sfa1_std1 leaves stay compared)\n', ...
    nnz(hit), height(L), nnz(contains(L.key, '.sfa1_std1') & ~hit));

%% (b) a differing leaf outside the patterns still fails
case_on = struct('dS', [1; 2; 3], 'S', [4; 5; 6]);
r = struct();
r.pointwise.pairs_steps.sfa1_std1.rk4_noise_off = case_on;
r.pointwise.pairs_steps.sfa3_std2.rk4_noise_off = case_on;
r.stable.mc.sfa1_std1 = struct('MC', 2.5, 'R2', [0.9 0.5]);
r.stable.mc.sfa3_std2 = struct('MC', 2.5, 'R2', [0.9 0.5]);
n = r;
n.pointwise.pairs_steps.sfa1_std1.rk4_noise_off.dS(2) = 2.1;   % matched: superseded
n.pointwise.pairs_steps.sfa1_std1.rk4_noise_off.S(2) = 5.1;    % sfa1_std1, not matched: fails
n.pointwise.pairs_steps.sfa3_std2.rk4_noise_off.dS(2) = 2.1;   % other condition: fails
n.stable.mc.sfa1_std1.MC = 3;                                  % matched: superseded
n.stable.mc.sfa3_std2.MC = 3;                                  % other condition: fails

[ok_b, Lb] = regression_tol_compare(r, n, 'verbose', false, 'superseded', sup);
status = containers.Map(Lb.key, Lb.status);
assert(~ok_b, 'differences outside the superseded patterns must fail the comparison');
assert(strcmp(status('pointwise.pairs_steps.sfa1_std1.rk4_noise_off.dS'), 'superseded'));
assert(strcmp(status('stable.mc.sfa1_std1.MC'), 'superseded'));
assert(strcmp(status('pointwise.pairs_steps.sfa1_std1.rk4_noise_off.S'), 'FAIL'));
assert(strcmp(status('pointwise.pairs_steps.sfa3_std2.rk4_noise_off.dS'), 'FAIL'));
assert(strcmp(status('stable.mc.sfa3_std2.MC'), 'FAIL'));
% Superseded: dS and MC, which moved, and R2, which matches but did not.
assert(nnz(strcmp(Lb.status, 'FAIL')) == 3 && nnz(strcmp(Lb.status, 'superseded')) == 3);
% A superseded leaf keeps what the comparison found, so it stays visible.
row = Lb(strcmp(Lb.key, 'stable.mc.sfa1_std1.MC'), :);
assert(row.max_abs == 0.5 && contains(row.detail{1}, 'FAIL'), ...
    'a superseded leaf must keep its measured difference and say it would have failed');
fprintf('PASS  (b) differing leaves outside the patterns fail; the matched ones are superseded\n');

%% (c) without the option, every difference fails
[ok_c, Lc] = regression_tol_compare(r, n, 'verbose', false);
assert(~ok_c && nnz(strcmp(Lc.status, 'FAIL')) == 5 && ~any(strcmp(Lc.status, 'superseded')), ...
    'without the superseded option all five differences must fail');
fprintf('PASS  (c) without the option all five differences fail\n');

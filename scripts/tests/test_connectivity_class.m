% test_connectivity_class.m - verification of the SRNNModel2 connectivity
% diagnostics: the static connectivity_class / strong_connectivity_probability
% and the instance method checkConnectivityClass.
%
% Assumes setup_paths() has already been run in this session (see CLAUDE.md).
%
% Covers:
%   1. Four-way ground truth on hand-built adjacency matrices
%   2. Transpose (edge-reversal) invariance of the class label
%   3. Insensitivity to self-loops
%   4. Calibration of the Erdos-Renyi P(strong) estimate against Monte Carlo
%   5. A live SRNNModel2 build at the fig_stim_engages_adaptation settings

clc;
fprintf('=== test_connectivity_class ===\n\n');

n_pass = 0; n_fail = 0;

%% ------------------------------------------------------------------
%  1. Four-way ground truth
%  ------------------------------------------------------------------
% NOTE: these are written in the MODEL's convention, A(i,j) ~= 0 means j -> i,
% so each matrix below is the transpose of the usual i->j adjacency. The class
% labels are reversal-invariant, so the orientation does not change the answer
% (test 2 makes that explicit), but writing them this way exercises the same
% code path the model uses.

cases = struct('name', {}, 'A', {}, 'expected', {});

% Directed 3-cycle 1->2->3->1 : strongly connected
A_cycle = zeros(3); A_cycle(2,1) = 1; A_cycle(3,2) = 1; A_cycle(1,3) = 1;
cases(end+1) = struct('name', 'directed 3-cycle', 'A', A_cycle, 'expected', 'strong');

% Path 1->2->3 : every pair comparable one way, but not strongly connected
A_path = zeros(3); A_path(2,1) = 1; A_path(3,2) = 1;
cases(end+1) = struct('name', 'path 1->2->3', 'A', A_path, 'expected', 'unilateral');

% Out-star 1->2, 1->3 : 2 and 3 mutually unreachable, underlying graph connected
A_star = zeros(3); A_star(2,1) = 1; A_star(3,1) = 1;
cases(end+1) = struct('name', 'out-star 1->{2,3}', 'A', A_star, 'expected', 'weak');

% Two disjoint edges : underlying graph is disconnected
A_split = zeros(4); A_split(2,1) = 1; A_split(4,3) = 1;
cases(end+1) = struct('name', 'two disjoint edges', 'A', A_split, 'expected', 'disconnected');

fprintf('--- 1. Four-way ground truth ---\n');
for k = 1:numel(cases)
    cls = SRNNModel2.connectivity_class(cases(k).A);
    ok  = strcmp(cls, cases(k).expected);
    n_pass = n_pass + ok; n_fail = n_fail + ~ok;
    fprintf('  [%s] %-20s -> %-13s (expected %s)\n', ...
        pass_str(ok), cases(k).name, cls, cases(k).expected);
end

%% ------------------------------------------------------------------
%  2. Transpose (edge-reversal) invariance
%  ------------------------------------------------------------------
fprintf('\n--- 2. Transpose invariance ---\n');
for k = 1:numel(cases)
    cls_a = SRNNModel2.connectivity_class(cases(k).A);
    cls_b = SRNNModel2.connectivity_class(cases(k).A.');
    ok = strcmp(cls_a, cls_b);
    n_pass = n_pass + ok; n_fail = n_fail + ~ok;
    fprintf('  [%s] %-20s A: %-13s A'': %s\n', ...
        pass_str(ok), cases(k).name, cls_a, cls_b);
end

%% ------------------------------------------------------------------
%  3. Self-loops must not change the class
%  ------------------------------------------------------------------
fprintf('\n--- 3. Self-loop insensitivity ---\n');
for k = 1:numel(cases)
    A_loops = cases(k).A;
    A_loops(1:size(A_loops,1)+1:end) = 1;      % add a full diagonal
    [cls_loops, info_loops] = SRNNModel2.connectivity_class(A_loops);
    ok = strcmp(cls_loops, cases(k).expected) && ...
         info_loops.n_selfloops == size(A_loops, 1);
    n_pass = n_pass + ok; n_fail = n_fail + ~ok;
    fprintf('  [%s] %-20s -> %-13s (%d self-loops counted, excluded from degrees)\n', ...
        pass_str(ok), cases(k).name, cls_loops, info_loops.n_selfloops);
end

%% ------------------------------------------------------------------
%  4. Calibration of the Erdos-Renyi P(strong) estimate
%  ------------------------------------------------------------------
fprintf('\n--- 4. P(strong) estimate vs Monte Carlo (n = 50, indegree = 4) ---\n');
n_mc      = 200;
n_nodes   = 50;
indeg     = 4;
alpha_mc  = indeg / n_nodes;

rng(12345);
is_strong_mc = false(n_mc, 1);
for k = 1:n_mc
    A_mc = double(rand(n_nodes) < alpha_mc);   % same mask construction as RMTMatrix
    is_strong_mc(k) = strcmp(SRNNModel2.connectivity_class(A_mc), 'strong');
end
p_emp = mean(is_strong_mc);
[p_est, lam_est] = SRNNModel2.strong_connectivity_probability(n_nodes, indeg);

% The estimate only rules out zero-degree nodes and assumes those events are
% independent, so it is an approximation rather than a bound (it can land on
% either side of the truth). Require agreement within Monte Carlo noise plus a
% modest modelling slack.
se = sqrt(p_emp * (1 - p_emp) / n_mc);
ok = abs(p_emp - p_est) < (3*se + 0.05);
n_pass = n_pass + ok; n_fail = n_fail + ~ok;
fprintf('  empirical  P(strong) = %.3f  (%d/%d draws, se = %.3f)\n', ...
    p_emp, sum(is_strong_mc), n_mc, se);
fprintf('  estimated  P(strong) = %.3f  (expected %.2f source/sink nodes)\n', ...
    p_est, lam_est);
fprintf('  [%s] estimate agrees with empirical within Monte Carlo noise\n', pass_str(ok));

fprintf('\n  Headroom above indegree = 4 (n = 50):\n');
fprintf('    %-10s %-12s %s\n', 'indegree', 'P(strong)', 'E[sources+sinks]');
for d = 4:2:16
    [pd, lamd] = SRNNModel2.strong_connectivity_probability(n_nodes, d);
    fprintf('    %-10g %-12.3f %.3f\n', d, pd, lamd);
end

%% ------------------------------------------------------------------
%  5. Live model at the bursting-figure settings
%  ------------------------------------------------------------------
fprintf('\n--- 5. Live SRNNModel2 (n = 50, indegree = 4, rng_seeds = [19 20]) ---\n');

% Not-built guard first.
model_unbuilt = SRNNModel2('n', 50, 'indegree', 4);
try
    model_unbuilt.checkConnectivityClass();
    ok = false;
    fprintf('  [FAIL] expected SRNNModel2:NotBuilt before build()\n');
catch ME
    ok = strcmp(ME.identifier, 'SRNNModel2:NotBuilt');
    fprintf('  [%s] errors with %s before build()\n', pass_str(ok), ME.identifier);
end
n_pass = n_pass + ok; n_fail = n_fail + ~ok;

% Build with the check flag on: the connectivity line should follow the
% spectral-radius line.
fprintf('\n  build() with check_connectivity = true:\n');
model = SRNNModel2('n', 50, 'f', 0.5, 'indegree', 4, ...
    'n_a_E', 3, 'n_b_E', 1, 'tau_d', 0.025, ...
    'T_range', [0 1], 'rng_seeds', [19 20], ...
    'lya_method', 'none', 'check_connectivity', true);
model.build();

fprintf('\n  model.checkConnectivityClass() summary:\n');
model.checkConnectivityClass();

[cls_live, info_live] = model.checkConnectivityClass();
ok = ismember(cls_live, {'strong', 'unilateral', 'weak', 'disconnected'}) && ...
     numel(info_live.scc_bins) == model.n && ...
     info_live.largest_scc_frac > 0 && info_live.largest_scc_frac <= 1 && ...
     info_live.n_sources == numel(info_live.zero_indegree);
n_pass = n_pass + ok; n_fail = n_fail + ~ok;
fprintf('\n  [%s] info struct is self-consistent\n', pass_str(ok));
fprintf('      class = %s, SCC sizes = %s\n', cls_live, ...
    mat2str(info_live.scc_sizes(1:min(5, end))));
fprintf('      realized mean indegree = %.2f (nominal %g), alpha_hat = %.4f (nominal %.4f)\n', ...
    info_live.mean_indegree, info_live.indegree_nominal, ...
    info_live.alpha_hat, info_live.alpha_nominal);
fprintf('      zero-indegree nodes: %s\n', mat2str(info_live.zero_indegree));
fprintf('      zero-outdegree nodes: %s\n', mat2str(info_live.zero_outdegree));

% A tolerance high enough to delete every edge must give 'disconnected'.
cls_tol = model.checkConnectivityClass('Tolerance', max(abs(model.W(:))));
ok = strcmp(cls_tol, 'disconnected');
n_pass = n_pass + ok; n_fail = n_fail + ~ok;
fprintf('  [%s] Tolerance above max|W| yields %s\n', pass_str(ok), cls_tol);

%% ------------------------------------------------------------------
fprintf('\n=== %d passed, %d failed ===\n', n_pass, n_fail);

function s = pass_str(ok)
    if ok
        s = 'PASS';
    else
        s = 'FAIL';
    end
end

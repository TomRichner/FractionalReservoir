% test_multi_timescale_std.m - Verify multi-timescale STD in SRNNModel2.
%
% Checks:
%   1. A model with n_b_E=3, n_b_I=2 (vector tau_b_*_rec) builds and runs.
%   2. The analytic Jacobian (SRNNModel2.compute_Jacobian_fast) matches a
%      central finite-difference Jacobian of the RHS (dynamics_fast), with
%      std_zero_floor both off and on.
%   3. Regression: the n_b_E=1 path still matches finite differences.
%   4. Validation: n_b_E=2 with a scalar tau_b_E_rec errors clearly.
%
% Uses SRNN2TestAccess (a test-only subclass) to reach the protected RHS.

close all; clear; clc;
setup_paths();

rng(0);
tol = 1e-6;         % max abs Jacobian error tolerance
h   = 1e-6;         % central-difference step
all_pass = true;

%% Test 1 & 2: n_b_E = 3, floor off and on
fprintf('=== Multi-timescale STD (n_b_E=3, n_b_I=2) ===\n');
params3 = make_model(3, [0.3, 1.0, 3.0]);
fprintf('  N_sys_eqs = %d\n', params3.N_sys_eqs);

err_off = compare_jacobian(params3, false, h);
fprintf('  Jacobian max abs err (floor OFF): %.3e ... %s\n', err_off, pass_str(err_off < tol));
all_pass = all_pass && (err_off < tol);

err_on = compare_jacobian(params3, true, h);
fprintf('  Jacobian max abs err (floor ON):  %.3e ... %s\n', err_on, pass_str(err_on < tol));
all_pass = all_pass && (err_on < tol);

%% Test 3: regression n_b_E = 1
fprintf('=== Regression (n_b_E=1) ===\n');
params1 = make_model(1, 1.0);
err1_off = compare_jacobian(params1, false, h);
fprintf('  Jacobian max abs err (floor OFF): %.3e ... %s\n', err1_off, pass_str(err1_off < tol));
all_pass = all_pass && (err1_off < tol);
err1_on = compare_jacobian(params1, true, h);
fprintf('  Jacobian max abs err (floor ON):  %.3e ... %s\n', err1_on, pass_str(err1_on < tol));
all_pass = all_pass && (err1_on < tol);

%% Test 4: validation error for scalar tau_b_E_rec with n_b_E=2
fprintf('=== Validation (scalar tau_b_E_rec with n_b_E=2 must error) ===\n');
threw = false;
try
    bad = SRNN2TestAccess('n', 60, 'indegree', 20, 'n_b_E', 2, 'tau_b_E_rec', 1.0);
    bad.build();
catch ME
    threw = strcmp(ME.identifier, 'SRNNModel:InvalidParams');
    fprintf('  Caught: %s\n', ME.message);
end
fprintf('  Errored as expected ... %s\n', pass_str(threw));
all_pass = all_pass && threw;

%% Summary
fprintf('\n=== %s ===\n', ternary(all_pass, 'ALL TESTS PASSED', 'TESTS FAILED'));
assert(all_pass, 'test_multi_timescale_std: one or more checks failed.');

%% ---- Local functions ----

function params = make_model(n_b_E, tau_b_E_rec)
    % Build a small model exercising all adaptation/STD blocks.
    model = SRNN2TestAccess( ...
        'n', 60, 'indegree', 20, ...
        'n_a_E', 2, 'n_a_I', 1, ...
        'n_b_E', n_b_E, 'n_b_I', 2, ...
        'tau_b_E_rec', tau_b_E_rec, 'tau_b_E_rel', 0.5, ...
        'tau_b_I_rec', [0.5, 2.0], 'tau_b_I_rel', 0.4, ...
        'T_range', [0, 1]);
    model.build();
    params = model.get_params();
    params.u_interpolant = model.u_interpolant;   % needed by dynamics_fast
end

function S = random_state(params)
    % Random valid state (b in (0,1], a >= 0, x moderate).
    nE = params.n_E; nI = params.n_I; n = params.n;
    a_E = 0.5 * rand(nE * params.n_a_E, 1);
    a_I = 0.5 * rand(nI * params.n_a_I, 1);
    b_E = 0.4 + 0.6 * rand(nE * params.n_b_E, 1);
    b_I = 0.4 + 0.6 * rand(nI * params.n_b_I, 1);
    x   = 0.8 * (2 * rand(n, 1) - 1);
    S = [a_E; a_I; b_E; b_I; x];
end

function J = fd_jacobian(S, params, h)
    % Central finite-difference Jacobian of the RHS.
    N = numel(S);
    J = zeros(N, N);
    for k = 1:N
        Sp = S; Sp(k) = S(k) + h;
        Sm = S; Sm(k) = S(k) - h;
        J(:, k) = (SRNN2TestAccess.rhs(0, Sp, params) - SRNN2TestAccess.rhs(0, Sm, params)) / (2 * h);
    end
end

function err = compare_jacobian(params, floor_on, h)
    % One FD-vs-analytic comparison; returns max abs error.
    params.std_zero_floor = floor_on;
    S = random_state(params);
    J_an = full(SRNNModel2.compute_Jacobian_fast(S, params));
    J_fd = fd_jacobian(S, params, h);
    err = max(abs(J_an(:) - J_fd(:)));
end

function s = pass_str(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function s = ternary(tf, a, b)
    if tf, s = a; else, s = b; end
end

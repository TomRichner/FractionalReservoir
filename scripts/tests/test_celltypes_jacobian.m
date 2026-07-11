% test_celltypes_jacobian.m
% Correctness gate for SRNNModelCellTypes: compare the analytic Jacobian
% (compute_Jacobian_fast_ct) against a central finite-difference of the RHS
% (dynamics_fast_ct), at a random state, for two configurations:
%   (A) SFA + STD on, STF off
%   (B) SFA + STD + STF all on
% External input is overridden to a constant (0) so it cancels in the difference
% and is state-independent. Prints CELLTYPES_JAC_PASS / CELLTYPES_JAC_FAIL.

setup_paths();

tol = 1e-6;
h   = 1e-6;
all_ok = true;

configs = struct( ...
    'name', {'SFA+STD (STF off)', 'SFA+STD+STF'}, ...
    'n_a',  {1, 1}, ...
    'n_b',  {1, 1}, ...
    'n_u',  {0, 1});

for ci = 1:numel(configs)
    cfg = configs(ci);
    model = SRNNModelCellTypes('n', 12, 'n_a', cfg.n_a, 'n_b', cfg.n_b, 'n_u', cfg.n_u, ...
        'T_range', [0, 2], 'fs', 100, 'rng_seeds', [11 12], 'level_of_chaos', 1.2, ...
        'rescale_by_abscissa', true, 'lya_method', 'none');
    model.build();

    p = model.cached_params;
    p.u_interpolant = @(t) zeros(1, p.n);      % constant, state-independent input

    N = p.N_sys_eqs;
    rng(7);
    S = 0.4 * randn(N, 1);                      % arbitrary interior state

    Jan = full(SRNNModelCellTypes.compute_Jacobian_fast_ct(S, p));
    Jnum = zeros(N);
    t0 = 1.0;
    for i = 1:N
        Sp = S; Sp(i) = Sp(i) + h;
        Sm = S; Sm(i) = Sm(i) - h;
        fp = SRNNModelCellTypes.dynamics_fast_ct(t0, Sp, p);
        fm = SRNNModelCellTypes.dynamics_fast_ct(t0, Sm, p);
        Jnum(:, i) = (fp - fm) / (2 * h);
    end

    err = max(abs(Jan - Jnum), [], 'all');
    ok = err <= tol;
    all_ok = all_ok && ok;
    fprintf('  [%s] %-20s N=%d  max|Jan-Jnum| = %.3e\n', ...
        ternary(ok, 'ok', 'FAIL'), cfg.name, N, err);

    if ~ok
        [~, lin] = max(abs(Jan(:) - Jnum(:)));
        [ri, cj] = ind2sub([N N], lin);
        fprintf('      worst entry (%d,%d): analytic=%.6g numeric=%.6g\n', ri, cj, Jan(ri, cj), Jnum(ri, cj));
    end
end

if all_ok
    disp('CELLTYPES_JAC_PASS');
else
    disp('CELLTYPES_JAC_FAIL');
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

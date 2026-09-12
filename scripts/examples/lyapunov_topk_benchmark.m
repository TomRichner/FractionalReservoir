% LYAPUNOV_TOPK_BENCHMARK The top-K Lyapunov spectrum on the paper's network.
%
% One seed of the paper preset at n = 500 (N ~ 4000 states), each of the
% three adaptation regimes, noise off: Benettin for lambda_1, then 'topk'
% at K = 10, 50 and 200 (and 400 if K = 200 took under ten minutes),
% timing each and printing lambda_1..lambda_5, the number of positive
% exponents, h_KS in bit/s, D_KY and whether it resolved within K, and the
% conditioning diagnostic. This is the cost/benefit input for choosing K in
% the ensemble stage, and the first look at the full network's spectrum.
%
% Everything is built on the client, so 'twister' seeds; ~15-40 min.
% Assumes setup_paths has run.
%
% See also: lyapunov_topk, test_lyapunov_topk, run_numerics_verification

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
[~, ~, conditions] = srnn_param_preset(P);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
T = 20;
base = {'sigma_u_noise', 0, 'ode_solver', 'rk4', 'rng_seeds', [1 2], 'fs', 400, ...
    'T_range', [0 T], 'lya_T_interval', [T/2 T], 'lya_warmup', 5, 'lya_dt', 0.1, ...
    'store_full_state', true};

fprintf('== top-K benchmark: %s, n = 500, T = %d s, window [%g %g], noise off ==\n', P, T, T/2, T);
for i = 1:numel(names)
    fprintf('\n-- %s --\n', names{i});
    evalc('m = build_from_preset(P, names{i}, base{:}, ''lya_method'', ''benettin'');');
    t0 = tic; evalc('m.run();'); tb = toc(t0);
    fprintf('  benettin      lambda_1 %+.4f                       %6.0f s (N = %d)\n', m.lya_results.LLE, tb, m.N_sys_eqs);
    Ks = [10 50 200];
    k = 1;
    while k <= numel(Ks)
        K = Ks(k);
        evalc('m = build_from_preset(P, names{i}, base{:}, ''lya_method'', ''topk'', ''lya_K'', K);');
        t0 = tic; evalc('m.run();'); tt = toc(t0);
        r = m.lya_results;
        if r.D_KY_resolved; dky = sprintf('%.1f', r.D_KY); else; dky = sprintf('> %d (unresolved)', K); end
        fprintf('  topk K = %-4d lambda_1..5 %s  n_pos %d  h_KS %.2f bit/s  D_KY %s  cond_max %.1e  %6.0f s (lya %.0f s)\n', ...
            K, mat2str(r.LE_spectrum(1:min(5, K))', 3), r.n_positive, r.h_KS_bits, dky, r.cond_max, tt, r.seconds);
        if K == 200 && r.seconds < 600 && ~r.D_KY_resolved
            Ks(end + 1) = 400; %#ok<SAGROW>
        end
        k = k + 1;
    end
end

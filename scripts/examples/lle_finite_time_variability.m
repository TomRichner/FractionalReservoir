% LLE_FINITE_TIME_VARIABILITY Why one Benettin run in the single-timescale
% regime is not a number to quote on its own.
%
% Companion to figs/numerics_verification_test_med/interpretation_report_medium.md
% section 7. The numerics-verification stage found the full-network LLE of
% the sfa1_std1 regime differing between ode45 (1e-10) and SRA1 (400 Hz) on
% one seed, 0.32 vs 0.58 over a 10 s window, while the two integrators'
% trajectory error had been shown to be 1e-6 per segment. This script asks
% three questions of that gap, all Benettin, noise off, seed [1 2]:
%
%   A. Does it survive changing ONLY the initial condition (x0 rescaled by
%      1-4%, W and stimulus fixed)? -> both integrators scatter over the same
%      0.3-0.7 range; five-IC means agree to 0.01.
%   B. Does the renormalisation interval lya_dt matter? -> no, to 4 digits.
%   C. Does it shrink with a longer window? -> SRA1 sits at 0.567 at 30 and
%      50 s; ode45 goes 0.32 -> 0.70 -> 0.66. The gap reverses and shrinks.
%
% Runs about 15 minutes. Assumes setup_paths has run.
%
% See also: run_numerics_verification, build_from_preset

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
solvers = {'sra1', 'ode45'};

fprintf('\n== A. same network, initial condition scaled by (1 + k/100), T = [0 20], window [10 20], lya_dt 0.02 ==\n');
fprintf('%6s %8s %8s\n', 'k', 'sra1', 'ode45');
for k = 0:4
    v = zeros(1, 2);
    for s = 1:2
        v(s) = lle_of(P, solvers{s}, 'T_range', [0 20], 'lya_T_interval', [10 20], 'x0_std', 0.1 * (1 + k / 100));
    end
    fprintf('%6d %+8.4f %+8.4f\n', k, v);
end

fprintf('\n== B. lya_dt, seed [1 2], T = [0 20], window [10 20] ==\n');
fprintf('%6s %8s %8s\n', 'lya_dt', 'sra1', 'ode45');
for d = [0.01 0.02 0.05 0.1]
    v = zeros(1, 2);
    for s = 1:2
        v(s) = lle_of(P, solvers{s}, 'T_range', [0 20], 'lya_T_interval', [10 20], 'lya_dt', d);
    end
    fprintf('%6.2f %+8.4f %+8.4f\n', d, v);
end

fprintf('\n== C. longer accumulation, seed [1 2], lya_dt 0.02 ==\n');
fprintf('%10s %8s %8s\n', 'window', 'sra1', 'ode45');
for T = [40 60]
    v = zeros(1, 2);
    for s = 1:2
        v(s) = lle_of(P, solvers{s}, 'T_range', [0 T], 'lya_T_interval', [10 T]);
    end
    fprintf('%10s %+8.4f %+8.4f\n', sprintf('[10 %d]', T), v);
end

function L = lle_of(P, solver, varargin)
evalc('m = build_from_preset(P, ''sfa1_std1'', ''sigma_u_noise'', 0, ''ode_solver'', solver, ''fs'', 400, ''lya_method'', ''benettin'', varargin{:});');
if strcmp(solver, 'ode45')
    m.ode_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'MaxStep', 1 / 400);
end
evalc('m.run();');
L = m.lya_results.LLE;
end

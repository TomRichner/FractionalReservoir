% TAU_SPREAD_ALIGNMENT Does a per-neuron SFA ladder spread change how fast
% Benettin (K = 1) and the top-K method converge on lambda_1 in the stable
% regime?
%
% Theory (docs/notes/Lyapunov_estimation_methods.md sec. 2): the finite-time
% bias of a single vector in the stable regime is the decay of its components
% in the faster-contracting directions (the STD band at -0.25 to -0.5 /s),
% which the spread does not touch; the spread only widens the slow SFA band
% itself. So the prediction is: Benettin's bias vs window length is
% UNCHANGED by the spread; top-K's likewise; the band width lambda_1 -
% lambda_10 grows in proportion to the spread; and lambda_1 walks toward
% -1/tau of the slowest drawn neuron. Alignment WITHIN the band gets slower,
% not faster (the top-two spacing is ~0.3 sigma / 10 s).
%
% Network: the tauSpread preset's physics (sfaEI, S_c spread 0.1, mu 8.25)
% in the sfa3_std2 regime at n = 100, F_tracks_network so the spectral
% radius matches the full network; noise off, rk4, seed [1 2], 70 s, with
% tau_a_spread overridden to 0 / 0.02 / 0.05 / 0.1. Per spread, three runs:
%   ref   top-10, warmup 40 s, accumulate [40 70]  (the best-aligned estimate)
%   topk  top-10, warmup 5 s,  accumulate [5 70]   (finite-time values read
%                                                  off at 15, 25, 45, 70 s)
%   ben   Benettin, same window                   (same read-outs)
% Printed, not asserted. ~5 min. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs.tau_a_spread, test_SRNNCellTypePairs_tau_a_spread,
%           lyapunov_topk

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p05_noise0p025_dualStd_3cond_mu8p25';
cond = 'sfa3_std2';
spreads = [0 0.02 0.05 0.1];
T = 70; reads = [15 25 45 70];
base = {'n', 100, 'indegree', 20, 'F_tracks_network', true, 'sigma_u_noise', 0, ...
    'ode_solver', 'rk4', 'rng_seeds', [1 2], 'T_range', [0 T], 'lya_dt', 0.1, ...
    'store_full_state', true};

fprintf('== tau_a_spread vs alignment: %s / %s, n = 100, %d s, noise off ==\n', P, cond, T);
fprintf('%-7s %-9s %-9s %-9s | %-31s | %-31s\n', 'spread', '-1/tau_sl', 'ref l1', 'l1-l10', ...
    'Benettin  l1 - ref  at 15/25/45/70 s', 'top-10    l1 - ref  at 15/25/45/70 s');
t_all = tic;
for s = spreads
    args = [base, {'tau_a_spread', s}];
    ref  = run_it(P, cond, [args, {'lya_method', 'topk', 'lya_K', 10, 'lya_T_interval', [40 T], 'lya_warmup', 40}]);
    topk = run_it(P, cond, [args, {'lya_method', 'topk', 'lya_K', 10, 'lya_T_interval', [5 T],  'lya_warmup', 5}]);
    ben  = run_it(P, cond, [args, {'lya_method', 'benettin',          'lya_T_interval', [5 T],  'lya_warmup', 5}]);
    if isempty(ref.tau_a_matrix)
        tau_sl = max(cellfun(@(x) max(x), ref.tau_a));
    else
        tau_sl = max(cellfun(@(x) max([x(:); 0]), ref.tau_a_matrix));
    end
    l1_ref = ref.lya_results.LLE;
    band = ref.lya_results.LE_spectrum(1) - ref.lya_results.LE_spectrum(10);
    fb = arrayfun(@(t) at_time(ben.lya_results.t_lya,  ben.lya_results.finite_lya, t), reads) - l1_ref;
    ft = arrayfun(@(t) at_time(topk.lya_results.t_lya, topk.lya_results.finite_LE_spectrum_t(:, 1), t), reads) - l1_ref;
    fprintf('%-7.3f %-9.4f %-9.4f %-9.4f | %s | %s\n', s, -1 / tau_sl, l1_ref, band, ...
        sprintf('%+8.4f', fb), sprintf('%+8.4f', ft));
end
fprintf('(%.0f s wall)\n', toc(t_all));
fprintf(['\nReading: "l1 - ref" is the finite-time lambda_1 accumulated over [5, t] minus the\n' ...
    '40 s-warmup top-10 value; negative = still biased toward more negative. Compare the\n' ...
    'ROWS: if the spread changed alignment, the Benettin column would shrink faster with\n' ...
    'spread. -1/tau_sl is what lambda_1 would be if the slowest neuron''s adaptation\n' ...
    'direction were the leading one and decoupled.\n']);

%% ------------------------------------------------------------------------
function m = run_it(P, cond, args) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, args{:});');
evalc('m.run();');
end

function v = at_time(t_lya, finite, t)
% Finite-time value accumulated up to (the segment ending at) t.
k = find(t_lya <= t - 0.1 + 1e-9, 1, 'last');
v = finite(k);
end

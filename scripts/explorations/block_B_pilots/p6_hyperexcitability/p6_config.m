function P = p6_config()
% P6_CONFIG Paths, thresholds and the networks chosen for block B pilot 6.
%
% One place for the numbers every p6 script shares, so the table, the figures
% and the rebuilt examples cannot disagree about a threshold.

main_root = 'C:\Users\m218089\Desktop\github_repos\FractionalReservoir';
P.run_dir   = fullfile(main_root, 'data', 'mu7chain_essential');
P.out_dir   = fullfile(main_root, 'figs', 'explore_B', 'p6_hyperexcitability');
P.conds     = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
P.cond_lbl  = {'No adaptation', 'Single-timescale', 'Chained MTS', 'Product MTS'};
P.cond_col  = [0 0 0; 0.72 0.49 0; 0 0.45 0.70; 0 0.7 0];

% Thresholds. The mean rate is the mean of r in [0, 1] over all neurons (E and
% I) and the whole 20-s run, as ParamSpaceAnalysis2 stores it.
%   rate_hi = 0.40 sits in the widest gap of the rate distribution in three of
%             the four conditions (0.35-0.60, 0.22-0.53, 0.24-0.52); the
%             product-MTS condition has no gap there.
%   rate_lo = 0.15 sits in the next gap down (0.12-0.24, 0.14-0.22, 0.12-0.21).
%   lle_thr = 0; lle_margin = 0.1 1/s is reported as the sensitivity of (a):
%             a stable network with adaptation sits at about -0.1 1/s (its
%             slowest adaptation rate), so +-0.1 is the scale of "marginal".
P.rate_hi    = 0.40;
P.rate_lo    = 0.15;
P.lle_thr    = 0;
P.lle_margin = 0.1;

% The six networks rebuilt for question 3 (config_idx in the joint sample),
% classed by their stored values under chained MTS.
P.rebuild_idx   = [240 969 1400 1879 284 761];
P.rebuild_class = {'chaotic + high rate', 'chaotic + high rate', ...
    'weakly chaotic, mid rate', ...
    'high rate, not chaotic', ...
    'chaotic, low rate', 'stable, low rate'};
P.rebuild_conds = {'no_adaptation', 'sfa3_stdc2'};
end

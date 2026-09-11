function cfg = numerics_verification_trials_config()
% NUMERICS_VERIFICATION_TRIALS_CONFIG The numerics stage at 'medium' with five network seeds.
%
%   cfg = numerics_verification_trials_config();
%   wait_for_parpool(13);                  % first: poll the licence, hold a seat
%   numerics_verification_trials_run       % then: stage (~45 min on 13 workers), three figures
%
% The ensemble version of numerics_verification_test_med_config. Every
% sub-experiment (noise-free reshoot, noisy reshoot on a shared Brownian
% path, Benettin LLE with ode45 vs SRA1, Benettin vs QR on the reduced
% network) is repeated on several networks, rng_seeds = [k, k+1]: five for the
% reshoot experiments, whose error is a local quantity with hundreds of
% restarts per seed, and twenty-five for the two Lyapunov comparisons. The
% point is the finite-time scatter of the LLE in the intermittent
% single-timescale regime (about +-0.2 over 10 s with either integrator), which
% a single seed cannot separate from integrator bias and a paired ensemble can.
%
% EVERY SETTING IS STATED HERE, like every *_config.m since 2026-09-09. Its
% own roots, so neither the 'fast' nor the single-seed 'medium' run is touched.
% The stage does not require an empty run_dir, so this can be rerun in place.
%
% See also: numerics_verification_trials_run, run_numerics_verification,
%           fig_numerics_verification, numerics_verification_test_med_config

cfg = struct();

cfg.preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
cfg.run_mode    = 'medium';
cfg.n_trials_reshoot = 5;      % A and B, passed to the stage by the _run script
cfg.n_trials_lle     = 25;     % L and C
cfg.n_workers        = 13;     % of 14 cores: 25 LLE seeds go in two batches, not three; ~1.5 GB per worker

cfg.run_dir  = 'data/numerics_verification_trials';  % the stage writes <run_dir>/numerics_verification
cfg.fig_root = 'figs/numerics_verification_trials';  % overwritten in place

cfg.visible_figures = false;

% The three variants of one figure function; names are the output folders and
% the three fig_tags are not prefixes of one another.
F = {};
F = add(F, 'fig_numerics_solver',     @fig_numerics_verification, false, {'variant', 'solver'});
F = add(F, 'fig_numerics_lya_method', @fig_numerics_verification, false, {'variant', 'lya_method'});
F = add(F, 'fig_numerics_ensemble',   @fig_numerics_verification, false, {'variant', 'ensemble'});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end

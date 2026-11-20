function detail = run_lyapunov_spectrum_n80_stage(cfg, run_dir, plan)
% RUN_LYAPUNOV_SPECTRUM_N80_STAGE A bundle's lyapunov_spectrum_n80 stage: the top-240 spectrum of the 80-neuron network into <run_dir>/lyapunov_spectrum_n80.
%
%   detail = RUN_LYAPUNOV_SPECTRUM_N80_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (its verbosity; the plan already holds
%            what the config asked for)
%   run_dir  the bundle's run directory
%   plan     from lyapunov_spectrum_n80_stage_plan(cfg)
%   detail   the file written, and what was run, on one line
%
% The run half of the opt-in stage 'lyapunov_spectrum_n80'
% (paper_opt_in_stages; FD-032). It calls run_lyapunov_spectrum once, with
% the plan's arguments, into <run_dir>/<plan.folder>, where
% fig_lyapunov_ranked_spectrum reads it. The seeds of each condition run in
% the session's pool (pool_ensure); the conditions run one after another.
%
% MEMORY. run() holds the whole trajectory (store_full_state): at 140 s,
% 500 Hz and 640 state variables (both multiple-timescale conditions at
% n = 80; one timescale has 320, no adaptation 80) about 360 MB per run;
% one run per pool worker at a time: a few GB on 12 workers, fine on the
% workstation.
%
% TIME, measured 2026-11-19 at 'essential' (one seed) on a MacBook Air M3:
% no adaptation 16 s, one timescale 49 s, chained MTS 81 s, product MTS
% 74 s; the stage 3.9 min. The conditions run one after another and a
% condition's seeds in parallel, so from 'fastmedium' up (three seeds) it
% takes about as long given three free workers.
%
% See also: lyapunov_spectrum_n80_stage_plan, run_lyapunov_spectrum,
%           fig_lyapunov_ranked_spectrum, run_all_paper_analyses

arguments
    cfg     (1,1) struct
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
vprintf(verbose, 'minimal', ...
    '  lyapunov_spectrum_n80: %d condition(s) x %d seed(s), K %d, T %g s, window [%g %g] s, warm-up %g s, %g Hz\n', ...
    numel(plan.conditions), plan.n_seeds, plan.K, plan.T, plan.lya_T_interval, plan.lya_warmup, plan.fs);
mat_file = run_lyapunov_spectrum(plan.args{:}, 'out_dir', fullfile(run_dir, plan.folder));
detail = sprintf('%s (%d conditions x %d seeds, K %d)', mat_file, ...
    numel(plan.conditions), plan.n_seeds, plan.K);
end

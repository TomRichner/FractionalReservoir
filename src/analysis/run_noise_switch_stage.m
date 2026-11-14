function detail = run_noise_switch_stage(cfg, run_dir, plan)
% RUN_NOISE_SWITCH_STAGE A bundle's noise_switch stage: the noise-switch rollouts into <run_dir>/noise_switch.
%
%   detail = RUN_NOISE_SWITCH_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (its verbosity; the plan already holds
%            what the config asked for)
%   run_dir  the bundle's run directory
%   plan     from noise_switch_stage_plan(cfg)
%   detail   the folder written, and what was run, on one line
%
% The run half of the opt-in stage 'noise_switch' (paper_opt_in_stages).
% It calls run_noise_switch_rollouts once, with the plan's arguments, into
% <run_dir>/<plan.folder>, where fig_noise_switch_rollouts reads it. The
% copies of each condition run in the session's pool (pool_ensure).
%
% See also: noise_switch_stage_plan, run_noise_switch_rollouts,
%           fig_noise_switch_rollouts, run_all_paper_analyses

arguments
    cfg     (1,1) struct
    run_dir (1,:) char
    plan    (1,1) struct
end

vprintf(cfg.verbose, 'minimal', ...
    '  noise_switch: %d condition(s) x %d copies, shared %g s then %g s on own seeds\n', ...
    numel(plan.conditions), plan.n_copies, plan.t_shared, plan.t_rollout);
out_dir = run_noise_switch_rollouts(plan.args{:}, 'out_dir', fullfile(run_dir, plan.folder));
detail = sprintf('%s (%d conditions x %d copies, rollout %g s)', out_dir, ...
    numel(plan.conditions), plan.n_copies, plan.t_rollout);
end

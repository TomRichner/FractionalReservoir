function registry = paper_opt_in_stages()
% PAPER_OPT_IN_STAGES The stages run_all_paper_analyses runs only for a bundle that names them.
%
%   registry = PAPER_OPT_IN_STAGES()
%
% A struct array, one element per opt-in stage, in run order:
%   name   the stage's name, as a bundle writes it in cfg.stages
%   plan   @(cfg) -> plan. Validates the bundle's settings for this stage
%          and returns whatever its run needs. No compute and no file
%          access: the entry point calls it in its pre-flight, before the
%          run directory exists, so a bad setting fails in milliseconds.
%   run    @(cfg, run_dir, plan) -> char. Runs the stage into run_dir and
%          returns a one-line description of what it wrote. It throws on
%          failure; the entry point wraps it like every stage.
%
% WHY OPT-IN. A stage that only some bundles want (block C's sweep is hours
% of compute) must not start running in every bundle that names no stages.
% So "every stage" means paper_stage_names' first list, and the stages here
% run only when cfg.stages names them, after the others.
%
% ADDING ONE is a row here and its two functions. paper_stage_names,
% paper_stage_plan and run_all_paper_analyses read this list and name no
% opt-in stage themselves.
%
% See also: paper_stage_names, paper_stage_plan, run_all_paper_analyses

% mc_gain_sweep  block C: memory capacity over gain and a second axis
% noise_switch   copies of one network that share the noise, then stop
%                sharing it (run_noise_switch_rollouts; Tom, 2026-11-14)
registry = struct( ...
    'name', {'mc_gain_sweep',            'noise_switch'}, ...
    'plan', {@mc_gain_sweep_stage_plan,  @noise_switch_stage_plan}, ...
    'run',  {@run_mc_gain_sweep_stage,   @run_noise_switch_stage});
end

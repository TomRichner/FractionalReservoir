function [solver, problem] = ospec_bundle_solver(cfg, who)
% OSPEC_BUNDLE_SOLVER The integrator a bundle hands the outlier_spectra runs: cfg.ode_solver, checked for ospec_run.
%
%   [solver, problem] = OSPEC_BUNDLE_SOLVER(cfg, who)
%
% For the plans of the stages whose runs are ospec_run (outlier_spectra,
% outlier_structures). solver is the bundle's cfg.ode_solver (stage_solver;
% Tom, 2026-11-17: every paper simulation at SRA1, stated once in the
% bundle config), or '' when the bundle states none, which leaves
% ospec_build_model its own rule (sra1 for a noisy preset, rk4 otherwise).
% who names the caller in stage_solver's error for an unknown solver name.
%
% problem is '' or a sentence for the plan's own bad-config error: ospec_run
% integrates in pieces (OspecNet.run_sampled) and so needs a FIXED-STEP
% integrator, 'rk4' or one of stochastic_solver_names; 'ode45' and the like
% are refused here, in the pre-flight, instead of in the first run.
%
% See also: outlier_spectra_stage_plan, outlier_structures_stage_plan,
%           ospec_build_model, stage_solver, stage_fs

arguments
    cfg (1,1) struct
    who (1,:) char
end
solver = '';
problem = '';
if isfield(cfg, 'ode_solver'); solver = stage_solver(cfg.ode_solver, '', who); end
if ~isempty(solver) && ~ismember(solver, [{'rk4'}, stochastic_solver_names()])
    problem = sprintf(['The runs of this stage integrate in pieces and need a fixed-step integrator ' ...
        '(rk4, %s); the bundle''s cfg.ode_solver is ''%s''.'], strjoin(stochastic_solver_names(), ', '), solver);
end
end

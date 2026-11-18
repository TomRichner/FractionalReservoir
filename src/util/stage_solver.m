function solver = stage_solver(solver_opt, default_solver, who)
% STAGE_SOLVER The integrator a stage runs with: the caller's, or the stage's own.
%
%   solver = STAGE_SOLVER(solver_opt, default_solver, who)
%
% solver_opt is the 'ode_solver' option a stage was given. Empty means "the
% stage's own choice" (default_solver: its run-mode table's, or the one
% build_from_preset picks from the preset's noise), so a stage called
% without the option behaves exactly as before the option existed. Anything
% else must be a valid ode_solver name (check_ode_solver) and is returned as
% is. who names the caller for the error identifier (<who>:BadOdeSolver).
%
% Why (Tom, 2026-11-17): every simulation of every paper analysis uses SRA1
% at fs = 500 Hz, so the Methods state one scheme. The bundle config states
% the integrator once (cfg.ode_solver, 'sra1' in mu7chain_config), and
% run_all_paper_analyses hands it to every stage beside cfg.fs (stage_fs).
% run_numerics_verification is the one exception: it compares solvers.
% Without noise SRA1 is Ralston's second-order Runge-Kutta scheme, so the
% noise-free stages (memory capacity, the noise-off Lyapunov spectrum, the
% single-neuron illustration) run it legally at sigma_u_noise = 0.
%
% See also: stage_fs, run_all_paper_analyses, check_ode_solver

if isempty(solver_opt)
    solver = default_solver;
    return
end
if ~(ischar(solver_opt) || (isstring(solver_opt) && isscalar(solver_opt)))
    error([who ':BadOdeSolver'], 'ode_solver must be empty (the stage''s own) or a solver name.');
end
solver = char(solver_opt);
check_ode_solver(solver, who);
end

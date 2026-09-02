function names = stochastic_solver_names()
%STOCHASTIC_SOLVER_NAMES Fixed-step SDE schemes (see sde_fixed_step).
%
% Required when sigma_u_noise > 0; usable at sigma = 0 too, where they
% degenerate to their deterministic parents -- sde_fixed_step reads
% has_noise = ~isempty(noise) && noise.sigma ~= 0, so an absent noise tensor
% simply drops the diffusion term. That is what the convergence tests rely on,
% and it is what lets run_memory_capacity offer 'sra1' at sigma = 0 for parity
% with the sweeps.
%
% See also: solver_names, deterministic_solver_names, sde_fixed_step

names = {'euler', 'heun', 'sra1'};
end

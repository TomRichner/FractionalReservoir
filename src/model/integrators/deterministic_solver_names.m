function names = deterministic_solver_names()
%DETERMINISTIC_SOLVER_NAMES Solvers that ignore sigma_u_noise.
%
% The only ones usable at sigma = 0 -- and the only ones usable AT ALL at
% sigma = 0 if you want bit-identical results with pre-SDE runs.
%
% See also: solver_names, stochastic_solver_names, check_noise_settings

names = {'ode45', 'ode15s', 'rk4'};
end

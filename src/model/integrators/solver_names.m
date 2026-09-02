function names = solver_names()
%SOLVER_NAMES The valid values of the `ode_solver` property, both classes.
%
%   names = SOLVER_NAMES()
%
% One list, shared. This and its two halves used to be statics duplicated
% VERBATIM on SRNNModel2 and SRNNCellTypePairs, which is why
% test_ode_solver_name carried an assertion whose only job was checking the two
% copies still agreed. A shared function makes that assertion vacuous.
%
% See also: deterministic_solver_names, stochastic_solver_names,
%           resolve_solver, check_ode_solver, check_noise_settings

names = [deterministic_solver_names(), stochastic_solver_names()];
end

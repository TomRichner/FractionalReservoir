function check_noise_settings(sigma_u_noise, ode_solver, err_id_prefix)
%CHECK_NOISE_SETTINGS Validate sigma_u_noise and its pairing with the integrator.
%
%   CHECK_NOISE_SETTINGS(sigma_u_noise, ode_solver, err_id_prefix)
%
% Shared by both model classes, by ParamSpaceAnalysis2's pre-flight, and by the
% memory-capacity and DC-LLE analyses -- so a swept or preset sigma with a
% deterministic solver fails BEFORE a run starts rather than at the first
% nonzero grid point.
%
% This was a static on SRNNModel2, which meant every SRNNCellTypePairs
% construction in the repo depended on SRNNModel2 existing: validate() calls it,
% and validate() runs from both the constructor and build(). Moving it here is
% what lets the paper pipeline run without SRNNModel2 on the path.
%
% err_id_prefix names the caller in the error identifier, so the message reads
% as coming from the thing the user invoked.
%
% See also: solver_names, stochastic_solver_names, check_ode_solver

if nargin < 3 || isempty(err_id_prefix), err_id_prefix = 'check_noise_settings'; end

if ~isscalar(sigma_u_noise) || ~isnumeric(sigma_u_noise) || ...
        ~isreal(sigma_u_noise) || ~isfinite(sigma_u_noise) || sigma_u_noise < 0
    error([err_id_prefix ':InvalidParams'], ...
        'sigma_u_noise must be a finite non-negative real scalar.');
end
if sigma_u_noise > 0 && ...
        ~ismember(lower(char(string(ode_solver))), stochastic_solver_names())
    error([err_id_prefix ':InvalidParams'], ...
        ['sigma_u_noise = %g requires a stochastic integrator; ' ...
         'ode_solver is ''%s''. Set ode_solver to one of %s. ' ...
         '(The adaptive solvers cannot step an SDE at all, and ' ...
         '''rk4'' is kept deterministic so sigma = 0 work stays ' ...
         'bit-identical to earlier runs.)'], ...
        sigma_u_noise, char(string(ode_solver)), ...
        strjoin(stochastic_solver_names(), ', '));
end
end

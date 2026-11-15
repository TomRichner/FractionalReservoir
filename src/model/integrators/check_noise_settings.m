function check_noise_settings(sigma_u_noise, ode_solver, err_id_prefix, release_sigma)
%CHECK_NOISE_SETTINGS Validate the noise settings and their pairing with the integrator.
%
%   CHECK_NOISE_SETTINGS(sigma_u_noise, ode_solver, err_id_prefix)
%   CHECK_NOISE_SETTINGS(sigma_u_noise, ode_solver, err_id_prefix, release_sigma)
%
% Shared by both model classes, by ParamSpaceAnalysis2's pre-flight, and by the
% memory-capacity and DC-LLE analyses -- so a swept or preset sigma with a
% deterministic solver fails BEFORE a run starts rather than at the first
% nonzero grid point.
%
% NOISE IS ON WHEN EITHER SOURCE IS ON: the Wiener term on x
% (sigma_u_noise > 0) or the release noise (any release_sigma > 0, the
% gamma increments of SRNNCellTypePairs.release_sigma; omitted or empty
% means none). Either needs a stochastic integrator. A spontaneous release
% RATE with release_sigma = 0 is a steady leak in the drift, not noise,
% and runs on any integrator. The values of release_sigma themselves are
% validated by the class that owns them.
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
if nargin < 4, release_sigma = []; end

if ~isscalar(sigma_u_noise) || ~isnumeric(sigma_u_noise) || ...
        ~isreal(sigma_u_noise) || ~isfinite(sigma_u_noise) || sigma_u_noise < 0
    error([err_id_prefix ':InvalidParams'], ...
        'sigma_u_noise must be a finite non-negative real scalar.');
end
is_stochastic = ismember(lower(char(string(ode_solver))), stochastic_solver_names());
if sigma_u_noise > 0 && ~is_stochastic
    error([err_id_prefix ':InvalidParams'], ...
        ['sigma_u_noise = %g requires a stochastic integrator; ' ...
         'ode_solver is ''%s''. Set ode_solver to one of %s. ' ...
         '(The adaptive solvers cannot step an SDE at all, and ' ...
         '''rk4'' is kept deterministic so sigma = 0 work stays ' ...
         'bit-identical to earlier runs.)'], ...
        sigma_u_noise, char(string(ode_solver)), ...
        strjoin(stochastic_solver_names(), ', '));
end
if isnumeric(release_sigma) && any(release_sigma(:) > 0) && ~is_stochastic
    error([err_id_prefix ':InvalidParams'], ...
        ['release_sigma = %s (release noise) requires a stochastic integrator; ' ...
         'ode_solver is ''%s''. Set ode_solver to one of %s: the release kicks ' ...
         'are applied by the fixed-step stochastic schemes, and a deterministic ' ...
         'solver would silently run the leak alone.'], ...
        mat2str(release_sigma), char(string(ode_solver)), ...
        strjoin(stochastic_solver_names(), ', '));
end
end

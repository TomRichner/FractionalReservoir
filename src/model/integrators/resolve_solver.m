function fn = resolve_solver(name, noise, err_id_prefix)
%RESOLVE_SOLVER Map an ode_solver name to a callable.
%
%   fn = RESOLVE_SOLVER(name)
%   fn = RESOLVE_SOLVER(name, noise)
%   fn = RESOLVE_SOLVER(name, noise, err_id_prefix)
%
% Returns a callable with the solver(odefun, tspan, y0, opts) signature.
%
% The stochastic schemes close over the pre-generated noise. A closure is safe
% here where it would not be on the property itself: this handle is built per
% call and never stored or compared -- what is persisted and compared is the
% NAME. That is the whole reason ode_solver is a string rather than a handle.
%
% err_id_prefix names the caller in the error identifier, so a model class
% raises its own ID rather than one naming this file. It defaults to
% 'resolve_solver' for a direct call. NOTE SRNNModel2 passes 'SRNNModel', not
% its own class name -- a pre-existing inconsistency that test_ode_solver_name
% asserts, so it is preserved deliberately rather than tidied.
%
% See also: solver_names, check_ode_solver, sde_fixed_step, ode_rk4

if nargin < 2, noise = []; end
if nargin < 3 || isempty(err_id_prefix), err_id_prefix = 'resolve_solver'; end

switch lower(name)
    case 'ode45'
        fn = @ode45;
    case 'ode15s'
        fn = @ode15s;
    case 'rk4'
        fn = @ode_rk4;
    case {'euler', 'heun', 'sra1'}
        scheme = lower(name);
        fn = @(f, tsp, y0, o) sde_fixed_step(f, tsp, y0, o, noise, scheme);
    otherwise
        error([err_id_prefix ':InvalidParams'], ...
            'Unknown ode_solver ''%s''. Valid: %s.', ...
            char(string(name)), strjoin(solver_names(), ', '));
end
end

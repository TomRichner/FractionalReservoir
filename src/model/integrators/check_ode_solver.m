function check_ode_solver(value, err_id_prefix)
%CHECK_ODE_SOLVER Reject handles by name and validate the ode_solver string.
%
%   CHECK_ODE_SOLVER(value)
%   CHECK_ODE_SOLVER(value, err_id_prefix)
%
% Called from a model class's constructor (so a bad name fails at the point of
% the mistake) and from its validate() (so a later assignment is caught too).
% Deliberately NOT a set method: ParamSpaceAnalysis2 assigns properties
% weakest-first and may set one before the ones that give it meaning, so both
% model classes have no setters at all.
%
% err_id_prefix names the caller in the error identifier. SRNNModel2 passes
% 'SRNNModel' -- not its own class name -- which test_ode_solver_name asserts
% at lines 67, 72 and 84; SRNNCellTypePairs passes its own name, asserted at
% line 89. The inconsistency is pre-existing and preserved on purpose: fixing
% it is a separate decision with its own test update.
%
% See also: solver_names, resolve_solver, check_noise_settings

if nargin < 2 || isempty(err_id_prefix), err_id_prefix = 'check_ode_solver'; end

if isa(value, 'function_handle')
    legacy = struct('ode45', 'ode45', 'ode15s', 'ode15s', ...
        'ode_rk4', 'rk4', 'ode23', 'ode45', 'ode113', 'ode45');
    as_str = func2str(value);
    if isfield(legacy, as_str)
        hint = sprintf('use ''%s'' instead of @%s', legacy.(as_str), as_str);
    else
        hint = sprintf('valid names are %s', strjoin(solver_names(), ', '));
    end
    error([err_id_prefix ':RenamedProperty'], ...
        ['''ode_solver'' now takes a NAME rather than a function ' ...
         'handle, so the choice survives into resolved_defaults and ' ...
         'compares cleanly across runs: %s.'], hint);
end
if ~(ischar(value) || isstring(value)) || ...
        ~ismember(lower(char(string(value))), solver_names())
    error([err_id_prefix ':InvalidParams'], ...
        'Unknown ode_solver ''%s''. Valid: %s.', ...
        char(string(value)), strjoin(solver_names(), ', '));
end
end

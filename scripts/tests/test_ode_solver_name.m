% test_ode_solver_name.m
% Verify that the integrator is chosen BY NAME, in the mould of `activation`.
%
% ode_solver used to hold a function handle. A handle does not survive into
% resolved_defaults as comparable data (same_config had to special-case it via
% func2str), it cannot be carried by a preset without smuggling code through a
% struct, and it let a caller pass any solver at all -- including ones the
% Lyapunov paths cannot use. It is now a string, and passing a handle raises
% RenamedProperty naming the replacement.
%
% Also checked here: lya_method = 'qr' now works with a fixed-step trajectory
% integrator. It never could before, because the QR variational solve was
% handed the model's own solver and integrates on a 2-point span, which
% ode_rk4 rejects outright. The variational equation is deterministic and
% independent of how the fiducial trajectory was produced, so it pins ode45.
%
% See also: SRNNModel2, SRNNCellTypePairs, ode_rk4, test_SRNN2_activation

fprintf('=== Testing the named ode_solver property ===\n\n');
all_passed = true;

%% The names, and the default
all_passed = check('the default is ''ode45''', ...
    strcmp(SRNNModel2().ode_solver, 'ode45')) && all_passed;
all_passed = check('the deterministic names are ode45/ode15s/rk4', ...
    isequal(SRNNModel2.deterministic_solver_names(), {'ode45', 'ode15s', 'rk4'})) && all_passed;
all_passed = check('the stochastic names are euler/heun/sra1', ...
    isequal(SRNNModel2.stochastic_solver_names(), {'euler', 'heun', 'sra1'})) && all_passed;
all_passed = check('solver_names is the union of the two', ...
    isequal(SRNNModel2.solver_names(), ...
        [SRNNModel2.deterministic_solver_names(), ...
         SRNNModel2.stochastic_solver_names()])) && all_passed;
all_passed = check('both classes offer the same names', ...
    isequal(SRNNModel2.solver_names(), SRNNCellTypePairs.solver_names())) && all_passed;

%% resolve_solver maps names to the expected callables
all_passed = check('''ode45'' resolves to @ode45', ...
    strcmp(func2str(SRNNModel2.resolve_solver('ode45')), 'ode45')) && all_passed;
all_passed = check('''rk4'' resolves to @ode_rk4', ...
    strcmp(func2str(SRNNModel2.resolve_solver('rk4')), 'ode_rk4')) && all_passed;
all_passed = check('the name is case-insensitive', ...
    strcmp(func2str(SRNNModel2.resolve_solver('RK4')), 'ode_rk4')) && all_passed;

%% Every name actually integrates, on the requested time grid
% The stochastic schemes are included at sigma = 0, where they degenerate to
% their deterministic parents.
for name = {'ode45', 'ode15s', 'rk4', 'euler', 'heun', 'sra1'}
    m = tiny(name{1});
    all_passed = check(sprintf('''%s'' runs and returns the requested grid', name{1}), ...
        m.has_run && numel(m.t_out) == numel(m.t_ex) && ...
        max(abs(m.t_out(:) - m.t_ex(:))) < 1e-9) && all_passed;
end

%% ode45 and rk4 agree on the trajectory without being identical
m45 = tiny('ode45');
mrk = tiny('rk4');
traj_diff = max(abs(m45.S_out(:) - mrk.S_out(:)));
fprintf('  (max |ode45 - rk4| over the trajectory: %.2e)\n', traj_diff);
all_passed = check('ode45 and rk4 agree to solver tolerance', ...
    traj_diff < 1e-4) && all_passed;
all_passed = check('...but are not the same integrator', ...
    traj_diff > 0) && all_passed;

%% Handles are rejected, by name, wherever they are supplied
[threw, err] = capture(@() SRNNModel2('ode_solver', @ode_rk4));
all_passed = check('a @ode_rk4 handle errors and names ''rk4''', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty') && ...
    contains(err.message, '''rk4''')) && all_passed;

[threw, err] = capture(@() SRNNModel2('ode_solver', @ode45));
all_passed = check('a @ode45 handle errors and names ''ode45''', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty') && ...
    contains(err.message, '''ode45''')) && all_passed;

[threw, err] = capture(@() SRNNModel2('ode_solver', 'rk45'));
all_passed = check('an unknown name errors and lists the valid ones', ...
    threw && contains(err.message, 'ode15s')) && all_passed;

% Assigned after construction rather than passed in: caught at build().
m_late = SRNNModel2('n', 12, 'indegree', 6);
m_late.ode_solver = @ode_rk4;
[threw, err] = capture(@() m_late.build());
all_passed = check('a handle assigned after construction is caught at build', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty')) && all_passed;

%% The Pairs class rejects handles the same way
[threw, err] = capture(@() pairs_model('ode_solver', @ode_rk4));
all_passed = check('SRNNCellTypePairs rejects a handle and names ''rk4''', ...
    threw && strcmp(err.identifier, 'SRNNCellTypePairs:RenamedProperty') && ...
    contains(err.message, '''rk4''')) && all_passed;
p = pairs_model('ode_solver', 'rk4');
all_passed = check('SRNNCellTypePairs runs from a name', ...
    p.has_run && numel(p.t_out) == numel(p.t_ex)) && all_passed;

%% QR now works with a fixed-step trajectory integrator
% ode_rk4 hard-errors on the 2-point span the variational solve uses, so this
% combination used to be impossible; the variational solve pins ode45 itself.
q = SRNNModel2('n', 8, 'indegree', 4, 'n_a_E', 0, 'n_a_I', 0, ...
    'n_b_E', 0, 'n_b_I', 0, 'T_range', [0 3], 'fs', 200, ...
    'ode_solver', 'rk4', 'lya_method', 'qr', 'lya_T_interval', [1 3], ...
    'lya_warmup', 1, 'store_full_state', true);
q.build();
evalc('q.run();');
all_passed = check('lya_method ''qr'' runs with ode_solver ''rk4''', ...
    numel(q.lya_results.LE_spectrum) == q.N_sys_eqs && ...
    all(isfinite(q.lya_results.LE_spectrum))) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function m = tiny(solver_name)
m = SRNNModel2('n', 16, 'indegree', 8, 'T_range', [0 2], 'fs', 200, ...
    'ode_solver', solver_name, 'lya_method', 'none', 'store_full_state', true);
m.build();
evalc('m.run();');
end

function p = pairs_model(varargin)
p = SRNNCellTypePairs('n', 16, 'indegree', 8, 'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], ...
    'n_a', [2 1], 'T_range', [0 2], 'fs', 200, 'lya_method', 'none', ...
    'store_full_state', true, varargin{:});
p.build();
evalc('p.run();');
end

function [threw, err] = capture(fcn)
threw = false; err = [];
try
    fcn();
catch err
    threw = true;
end
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end

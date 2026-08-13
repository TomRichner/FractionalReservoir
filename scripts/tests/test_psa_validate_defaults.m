% test_psa_validate_defaults.m
% Verify ParamSpaceAnalysis2.validate_model_defaults rejects model_defaults
% fields that can never take effect, and warns on fields that are shadowed by
% the grid or by a condition.
%
% No simulations are run -- validation is pure introspection over
% SRNNModel2's property list.
%
% See also: ParamSpaceAnalysis2, SRNNModel2

fprintf('=== Testing ParamSpaceAnalysis2.validate_model_defaults ===\n\n');

all_passed = true;

%% Helper: build a minimal PSA with one grid parameter
make_psa = @() build_test_psa();

%% 1. Unknown property (typo) -> error naming the intended property
psa = make_psa();
psa.model_defaults.tau_D = 1;
[threw, err] = capture_error(@() psa.validate_model_defaults());
all_passed = check('typo tau_D errors', threw && ...
    strcmp(err.identifier, 'ParamSpaceAnalysis2:InvalidModelDefaults')) && all_passed;
all_passed = check('  ... and suggests tau_d', threw && contains(err.message, '''tau_d''')) && all_passed;

%% 2. Dependent property -> error
psa = make_psa();
psa.model_defaults.alpha = 0.3;
[threw, err] = capture_error(@() psa.validate_model_defaults());
all_passed = check('Dependent alpha errors', threw && contains(err.message, 'Dependent')) && all_passed;

%% 3. Non-public property -> error
psa = make_psa();
psa.model_defaults.W = eye(3);
[threw, err] = capture_error(@() psa.validate_model_defaults());
all_passed = check('protected W errors', threw && contains(err.message, 'not publicly settable')) && all_passed;

%% 4. Several bad fields at once -> ONE error listing all of them
psa = make_psa();
psa.model_defaults.tau_D = 1;
psa.model_defaults.alpha = 0.3;
psa.model_defaults.W = eye(3);
[threw, err] = capture_error(@() psa.validate_model_defaults());
all_passed = check('accumulates into one error', threw && contains(err.message, '3 field(s)')) && all_passed;

%% 5. Grid-shadowed field -> silent. A parameter preset carries a value for
%    every axis and each sweep varies a different subset, so this is expected.
psa = make_psa();  % grid parameter is 'f'
psa.model_defaults.f = 0.6;
[threw, ~, warn_id] = capture_error(@() psa.validate_model_defaults());
all_passed = check('grid-shadowed f is silent', ~threw && isempty(warn_id)) && all_passed;

%% 6. Condition-SET field -> warning (it can never take effect)
psa = make_psa();
psa.model_defaults.n_a_E = 3;
[threw, ~, warn_id] = capture_error(@() psa.validate_model_defaults());
all_passed = check('condition-set n_a_E warns, does not error', ...
    ~threw && strcmp(warn_id, 'ParamSpaceAnalysis2:ConditionParamShadowed')) && all_passed;

%% 6b. A condition field NO condition sets is not shadowed -- it takes effect.
%     The default conditions set n_a_E / n_b_E only.
psa = make_psa();
psa.model_defaults.n_a_I = 2;
[threw, ~, warn_id] = capture_error(@() psa.validate_model_defaults());
all_passed = check('n_a_I does not warn (no condition sets it)', ...
    ~threw && isempty(warn_id)) && all_passed;

%% 6c. Condition fields may never be swept
psa = make_psa();
[threw, err] = capture_error(@() psa.add_grid_parameter('n_a_E', [0, 3]));
all_passed = check('add_grid_parameter(''n_a_E'') errors', threw && ...
    strcmp(err.identifier, 'ParamSpaceAnalysis2:ConditionFieldAsGridParam')) && all_passed;

psa = make_psa();
[threw, err] = capture_error(@() psa.add_vector_parameter('n_b_E', ...
    'fixed_value', 1, 'vary_range', [1, 2], 'n_elements', 2));
all_passed = check('add_vector_parameter(''n_b_E'') errors', threw && ...
    strcmp(err.identifier, 'ParamSpaceAnalysis2:ConditionFieldAsGridParam')) && all_passed;

psa = make_psa();
psa.grid_params{end+1} = 'n_b_I';   % bypass the setter
[threw, err] = capture_error(@() psa.run());
all_passed = check('run() catches a condition field assigned into grid_params', ...
    threw && strcmp(err.identifier, 'ParamSpaceAnalysis2:ConditionFieldAsGridParam')) && all_passed;

%% 7. Clean defaults -> passes silently (no error, no warning)
psa = make_psa();
psa.model_defaults.tau_d = 0.1;
psa.model_defaults.c_E = 0.05;
psa.model_defaults.ode_solver = 'rk4';
psa.model_defaults.input_config = struct('n_steps', 3, 'no_stim_pattern', [true false true]);
[threw, ~, warn_id] = capture_error(@() psa.validate_model_defaults());
all_passed = check('clean defaults pass silently', ~threw && isempty(warn_id)) && all_passed;

%% 8. Empty model_defaults -> no-op
psa = make_psa();
[threw, ~, warn_id] = capture_error(@() psa.validate_model_defaults());
all_passed = check('empty model_defaults is a no-op', ~threw && isempty(warn_id)) && all_passed;

%% 9. run() rejects before creating the output directory
psa = make_psa();
psa.model_defaults.tau_D = 1;
dir_before = psa.output_dir;
[threw, ~] = capture_error(@() psa.run());
all_passed = check('run() rejects bad defaults', threw) && all_passed;
all_passed = check('  ... and leaves no dated output folder', ...
    strcmp(psa.output_dir, dir_before) && ~exist(psa.output_dir, 'dir')) && all_passed;

%% 10. A swept sigma_u_noise with a deterministic integrator is caught up front
% The model would catch this too, but only at the first NONZERO grid point --
% which for a sweep starting at sigma = 0 can be an hour in, after a partial
% run directory exists. run() therefore pre-flights it beside the defaults
% check, before the output directory is created.
psa = make_psa();
psa.add_grid_parameter('sigma_u_noise', [0, 0.5]);
psa.model_defaults.ode_solver = 'rk4';
dir_before = psa.output_dir;
[threw, err] = capture_error(@() psa.run());
all_passed = check('a swept sigma with a deterministic solver is rejected', ...
    threw && contains(err.message, 'sra1')) && all_passed;
all_passed = check('  ... before any output folder is created', ...
    strcmp(psa.output_dir, dir_before) && ~exist(psa.output_dir, 'dir')) && all_passed;

% The same sweep with a stochastic integrator passes the pre-flight.
psa = make_psa();
psa.add_grid_parameter('sigma_u_noise', [0, 0.5]);
psa.model_defaults.ode_solver = 'sra1';
[threw, ~] = capture_error(@() psa.validate_noise_settings());
all_passed = check('the same sweep with ''sra1'' passes', ~threw) && all_passed;

% A sweep that never leaves zero needs no stochastic solver at all.
psa = make_psa();
psa.add_grid_parameter('sigma_u_noise', [0, 0]);
psa.model_defaults.ode_solver = 'rk4';
[threw, ~] = capture_error(@() psa.validate_noise_settings());
all_passed = check('an all-zero sigma sweep is left alone', ~threw) && all_passed;

% A scalar default (not a grid axis) is checked the same way.
psa = make_psa();
psa.model_defaults.sigma_u_noise = 0.2;
psa.model_defaults.ode_solver = 'ode45';
[threw, err] = capture_error(@() psa.validate_noise_settings());
all_passed = check('a scalar sigma default with ode45 is rejected too', ...
    threw && contains(err.message, 'sra1')) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helper functions
function psa = build_test_psa()
% Minimal PSA: one grid parameter ('f') and the default conditions, which
% carry n_a_E / n_b_E. Nothing is run.
psa = ParamSpaceAnalysis2('n_levels', 2, 'verbose', false);
psa.add_grid_parameter('f', [0.4, 0.6]);
psa.output_dir = fullfile(tempdir, 'psa_validate_defaults_test');
end

function [threw, err, warn_id] = capture_error(fcn)
% Run FCN, reporting whether it errored and the identifier of the last
% warning it raised (empty when it warned about nothing). Warnings are left
% enabled on purpose: a disabled warning does not update lastwarn, and seeing
% the expected text in the console is useful when a case fails.
threw = false;
err = [];
lastwarn('', '');
try
    fcn();
catch err
    threw = true;
end
[~, warn_id] = lastwarn();
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

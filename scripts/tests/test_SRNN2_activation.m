% test_SRNN2_activation.m
% Verify the named-activation scheme: SRNNModel2 chooses its nonlinearity by
% NAME ('logistic' | 'piecewise' | 'tanh') parameterised by S_a / S_c, and the
% activation_function / _derivative handles are Dependent, rebuilt from those.
%
% The point of the scheme is that S_a and S_c always describe the function
% actually in use, so the tests below check the PARAMETERS BITE (or correctly
% do not) rather than just that a handle exists.
%
% See also: SRNNModel2, srnn_param_preset

fprintf('=== Testing SRNNModel2 activation selection ===\n\n');
all_passed = true;
xg = linspace(-1, 2, 101);

%% Each name produces the documented function and derivative
m = SRNNModel2('activation', 'logistic', 'S_c', 0.4);
all_passed = check('logistic matches logisticSigmoid(x, S_c)', ...
    same(m.activation_function(xg), SRNNModel2.logisticSigmoid(xg, 0.4)) && ...
    same(m.activation_function_derivative(xg), SRNNModel2.logisticSigmoidDerivative(xg, 0.4))) ...
    && all_passed;

m = SRNNModel2('activation', 'piecewise', 'S_a', 0.9, 'S_c', 0.4);
all_passed = check('piecewise matches piecewiseSigmoid(x, S_a, S_c)', ...
    same(m.activation_function(xg), SRNNModel2.piecewiseSigmoid(xg, 0.9, 0.4)) && ...
    same(m.activation_function_derivative(xg), SRNNModel2.piecewiseSigmoidDerivative(xg, 0.9, 0.4))) ...
    && all_passed;

m = SRNNModel2('activation', 'tanh');
all_passed = check('tanh matches tanhActivation(x)', ...
    same(m.activation_function(xg), SRNNModel2.tanhActivation(xg)) && ...
    same(m.activation_function_derivative(xg), SRNNModel2.tanhActivationDerivative(xg))) ...
    && all_passed;

all_passed = check('the class default is logistic', ...
    strcmp(SRNNModel2().activation, 'logistic')) && all_passed;

%% The parameters mean what they say
m = SRNNModel2('activation', 'logistic', 'S_c', 0.4);
before = m.activation_function(0.4);
m.S_c = 0.0;
all_passed = check('S_c bites for logistic', ...
    abs(m.activation_function(0.4) - before) > 1e-6) && all_passed;
m.S_a = 0.1;
all_passed = check('S_a does NOT bite for logistic', ...
    same(m.activation_function(xg), SRNNModel2.logisticSigmoid(xg, 0.0))) && all_passed;

m = SRNNModel2('activation', 'piecewise', 'S_a', 0.9, 'S_c', 0.4);
before = m.activation_function(0.6);
m.S_a = 0.2;
all_passed = check('S_a bites for piecewise', ...
    abs(m.activation_function(0.6) - before) > 1e-6) && all_passed;
before = m.activation_function(0.6);
m.S_c = 0.0;
all_passed = check('S_c bites for piecewise', ...
    abs(m.activation_function(0.6) - before) > 1e-6) && all_passed;

m = SRNNModel2('activation', 'tanh');
before = m.activation_function(xg);
m.S_a = 0.1; m.S_c = 0.9;
all_passed = check('neither parameter bites for tanh', ...
    same(m.activation_function(xg), before)) && all_passed;

%% The handles are derived, never stale, and never settable
m = SRNNModel2('activation', 'piecewise', 'S_c', 0.4);
h_before = m.activation_function;
m.S_c = 0.1;
all_passed = check('a handle taken earlier does not go stale in the model', ...
    same(m.activation_function(xg), SRNNModel2.piecewiseSigmoid(xg, m.S_a, 0.1))) && all_passed;
all_passed = check('but a previously fetched handle keeps its own bound values', ...
    same(h_before(xg), SRNNModel2.piecewiseSigmoid(xg, m.S_a, 0.4))) && all_passed;

[threw, err] = capture(@() SRNNModel2('activation_function', @sin));
all_passed = check('setting activation_function errors, naming activation', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty') && ...
    contains(err.message, 'activation')) && all_passed;

[threw, err] = capture(@() SRNNModel2('activation_function_derivative', @cos));
all_passed = check('setting the derivative errors too', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty')) && all_passed;

%% An unknown name is rejected at build()
[threw, err] = capture(@() build_tiny(SRNNModel2('activation', 'sigmoidal')));
all_passed = check('unknown activation errors from build(), listing valid names', ...
    threw && contains(err.message, 'Unknown activation') && ...
    contains(err.message, 'piecewise')) && all_passed;

%% activation_custom overrides the named choice, and is validated
m = SRNNModel2('activation', 'tanh', ...
    'activation_custom', {@(x) 2*x, @(x) 2*ones(size(x))});
all_passed = check('activation_custom overrides the named activation', ...
    same(m.activation_function(xg), 2*xg) && ...
    same(m.activation_function_derivative(xg), 2*ones(size(xg)))) && all_passed;

[threw, ~] = capture(@() build_tiny(SRNNModel2('activation_custom', {@sin})));
all_passed = check('a 1-element activation_custom errors', threw) && all_passed;
[threw, ~] = capture(@() build_tiny(SRNNModel2('activation_custom', {1, 2})));
all_passed = check('a non-handle activation_custom errors', threw) && all_passed;

%% A built model runs end to end under each named activation
for nm = {'logistic', 'piecewise', 'tanh'}
    ok = true;
    try
        mm = SRNNModel2('activation', nm{1}, 'n', 12, 'indegree', 6, ...
            'T_range', [0 1], 'fs', 200, 'ode_solver', @ode_rk4, ...
            'lya_method', 'none', 'n_a_E', 1, 'n_b_E', 1);
        mm.build();
        mm.run();
    catch
        ok = false;
    end
    all_passed = check(sprintf('a %s model builds and runs', nm{1}), ok) && all_passed;
end

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function tf = same(a, b)
tf = isequal(size(a), size(b)) && max(abs(a(:) - b(:))) < 1e-15;
end

function build_tiny(m)
m.n = 6; m.indegree = 3; m.T_range = [0 0.1]; m.fs = 100; m.lya_method = 'none';
m.build();
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

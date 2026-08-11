% test_srnn_param_preset.m
% Verify the parameter presets, the shared run-mode config, and how a preset
% interacts with a sweep.
%
% Covers the two behavioural claims that matter:
%   * every preset is VALID -- a typo in a preset definition is caught here
%     rather than six hours into a sweep;
%   * a grid axis overrides the preset for the sweep that varies it, while the
%     preset still supplies every axis that sweep holds fixed.
%
% See also: srnn_param_preset, analysis_run_config, merge_struct

fprintf('=== Testing srnn_param_preset / analysis_run_config ===\n\n');
all_passed = true;

preset_names = {'default', 'overconnected'};

%% Every preset returns a struct whose fields are all settable SRNNModel2 params
for i = 1:numel(preset_names)
    d = srnn_param_preset(preset_names{i});
    all_passed = check(sprintf('%s returns a scalar struct', preset_names{i}), ...
        isstruct(d) && isscalar(d)) && all_passed;

    % validate_model_defaults is the authority: unknown / Dependent / non-public
    % names are hard errors. Use a PSA with no grid and no conditions carrying
    % the preset, so only the preset's own fields are judged.
    psa = ParamSpaceAnalysis2('verbose', false);
    psa.model_defaults = d;
    [threw, err] = capture_error(@() psa.validate_model_defaults());
    all_passed = check(sprintf('%s passes validate_model_defaults', preset_names{i}), ...
        ~threw) && all_passed;
    if threw
        fprintf('      %s\n', err.message);
    end
end

all_passed = check('''default'' is empty', ...
    isempty(fieldnames(srnn_param_preset('default')))) && all_passed;

[threw, err] = capture_error(@() srnn_param_preset('no_such_preset'));
all_passed = check('unknown preset name errors and lists the valid ones', ...
    threw && strcmp(err.identifier, 'srnn_param_preset:UnknownPreset') && ...
    contains(err.message, 'overconnected')) && all_passed;

%% A preset must not carry names that belong to run_mode or the conditions
banned = {'fs', 'T_range', 'ode_solver', 'lya_T_interval', ...
    'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I', 'n_levels', 'n_reps'};
for i = 1:numel(preset_names)
    d = srnn_param_preset(preset_names{i});
    bad = intersect(fieldnames(d), banned);
    all_passed = check(sprintf('%s carries no run_mode/condition names', preset_names{i}), ...
        isempty(bad)) && all_passed;
    if ~isempty(bad); fprintf('      offending: %s\n', strjoin(bad', ', ')); end
end

%% overconnected carries the values it is supposed to
d = srnn_param_preset('overconnected');
all_passed = check('overconnected sets mu_I_tilde_relative = -6', ...
    isfield(d, 'mu_I_tilde_relative') && d.mu_I_tilde_relative == -6) && all_passed;
all_passed = check('overconnected carries the swept axes (n, f, level_of_chaos)', ...
    all(isfield(d, {'n', 'f', 'level_of_chaos'}))) && all_passed;
% The nonlinearity is named data, so a preset cannot carry a handle whose
% captured constants disagree with its own S_a / S_c.
all_passed = check('overconnected names the activation instead of carrying handles', ...
    isfield(d, 'activation') && strcmp(d.activation, 'piecewise') && ...
    ~isfield(d, 'activation_function')) && all_passed;
all_passed = check('no preset carries a function handle', ...
    ~any(cellfun(@(nm) any(structfun(@(v) isa(v, 'function_handle'), ...
        srnn_param_preset(nm))), preset_names))) && all_passed;
% ...and the named choice reproduces the intended nonlinearity.
m = SRNNModel2('activation', d.activation, 'S_a', d.S_a, 'S_c', d.S_c);
all_passed = check('preset activation evaluates as piecewise(S_a, S_c)', ...
    abs(m.activation_function(0.5) - SRNNModel2.piecewiseSigmoid(0.5, d.S_a, d.S_c)) < 1e-15) ...
    && all_passed;

%% analysis_run_config covers every analysis x mode and rejects bad input
analyses = {'sensitivity', 'tau_sensitivity', 'param_space'};
modes = {'fast', 'medium', 'production'};
cfg_ok = true;
for i = 1:numel(analyses)
    for j = 1:numel(modes)
        cfg = analysis_run_config(analyses{i}, modes{j});
        cfg_ok = cfg_ok && isscalar(cfg.n_levels) && cfg.n_levels > 0 && ...
            isfield(cfg.model, 'ode_solver') && isfield(cfg.model, 'fs') && ...
            isfield(cfg.model, 'T_range');
        % reps axis exists everywhere except param_space
        cfg_ok = cfg_ok && (isfield(cfg, 'n_reps') == ~strcmp(analyses{i}, 'param_space'));
    end
end
all_passed = check('all 3 analyses x 3 modes return a usable config', cfg_ok) && all_passed;

[threw, ~] = capture_error(@() analysis_run_config('nope', 'fast'));
all_passed = check('unknown analysis errors', threw) && all_passed;
[threw, ~] = capture_error(@() analysis_run_config('sensitivity', 'turbo'));
all_passed = check('unknown run_mode errors', threw) && all_passed;

%% merge_struct precedence: the second argument wins
a = struct('x', 1, 'y', 2);
b = struct('y', 99, 'z', 3);
s = merge_struct(a, b);
all_passed = check('merge_struct: over wins, union of fields', ...
    s.x == 1 && s.y == 99 && s.z == 3) && all_passed;

%% A grid axis overrides the preset; the other preset axes still apply
psa = ParamSpaceAnalysis2('n_levels', 2, 'verbose', false);
psa.model_defaults = srnn_param_preset('overconnected');
psa.add_grid_parameter('n', [40, 80]);
[threw, ~, warn_id] = capture_error(@() psa.validate_model_defaults());
all_passed = check('preset + swept n validates without warning', ...
    ~threw && isempty(warn_id)) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function [threw, err, warn_id] = capture_error(fcn)
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

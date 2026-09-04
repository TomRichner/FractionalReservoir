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

% EVERY preset, not a hand-kept subset: the list is read from the function
% itself, so a preset added without a test cannot go unchecked. (This used to be
% {'default', 'overconnected'}, which left eight of eleven presets unexamined.)
preset_names = srnn_param_preset_names_for_test();

%% Every preset returns a struct whose fields are all settable model params
for i = 1:numel(preset_names)
    % The second output names the MODEL CLASS the overrides are written for --
    % the two classes have disjoint vocabularies, so validating a
    % SRNNCellTypePairs preset against SRNNModel2 would reject every per-type
    % field. Thread it through rather than assuming SRNNModel2.
    [d, model_class] = srnn_param_preset(preset_names{i});
    all_passed = check(sprintf('%s returns a scalar struct', preset_names{i}), ...
        isstruct(d) && isscalar(d)) && all_passed;

    % validate_model_defaults is the authority: unknown / Dependent / non-public
    % names are hard errors. Use a PSA with no grid and no conditions carrying
    % the preset, so only the preset's own fields are judged.
    psa = ParamSpaceAnalysis2('verbose', false);
    psa.model_class = model_class;
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
% lya_dt and lya_warmup join lya_T_interval here for the same reason: they are
% Lyapunov cost/fidelity knobs, so they belong to analysis_run_config (how much
% compute) rather than to a preset (which network).
banned = {'fs', 'T_range', 'ode_solver', 'lya_T_interval', 'lya_dt', ...
    'lya_warmup', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I', 'n_levels', 'n_reps'};
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
modes = run_mode_names();   % the canonical list; 'fast2' was removed 2026-09-03
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
all_passed = check('all 3 analyses x 4 modes return a usable config', cfg_ok) && all_passed;

%% The deterministic / stochastic solver pair
% sigma_u_noise is physics and lives in a preset, but sigma > 0 REQUIRES a
% stochastic integrator -- so the preset selects which of a cell's two solvers
% is used. Checking here keeps that contract pinned in one place.
noisy = struct('sigma_u_noise', 0.02);
quiet = struct('sigma_u_noise', 0);
det_ok = true; sto_ok = true; leak_ok = true; zero_ok = true;
for i = 1:numel(analyses)
    for j = 1:numel(modes)
        cd_ = analysis_run_config(analyses{i}, modes{j});
        cs  = analysis_run_config(analyses{i}, modes{j}, noisy);
        cz  = analysis_run_config(analyses{i}, modes{j}, quiet);
        % rk4 up to medium, ode45 only at production
        want_det = 'rk4';
        if strcmp(modes{j}, 'production'); want_det = 'ode45'; end
        det_ok  = det_ok  && strcmp(cd_.model.ode_solver, want_det) && ~cd_.is_stochastic;
        sto_ok  = sto_ok  && strcmp(cs.model.ode_solver, 'sra1')   &&  cs.is_stochastic;
        zero_ok = zero_ok && strcmp(cz.model.ode_solver, want_det) && ~cz.is_stochastic;
        % sde_solver must stay OUT of cfg.model, which becomes model_defaults
        leak_ok = leak_ok && ~isfield(cs.model, 'sde_solver') && ...
            ~isfield(cs.model, 'is_stochastic');
    end
end
all_passed = check('deterministic column is rk4 up to medium, ode45 at production', ...
    det_ok) && all_passed;
all_passed = check('a preset with sigma > 0 selects sra1 in every cell', ...
    sto_ok) && all_passed;
all_passed = check('a preset with sigma = 0 stays deterministic', zero_ok) && all_passed;
all_passed = check('sde_solver/is_stochastic never leak into cfg.model', ...
    leak_ok) && all_passed;

% The noise preset must actually reach that path end to end.
[dn, mcn] = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_4cond');
cfg_n = analysis_run_config('sensitivity', 'fast', dn);
all_passed = check('the noise preset selects sra1 at fast', ...
    strcmp(cfg_n.model.ode_solver, 'sra1') && strcmp(mcn, 'SRNNCellTypePairs')) && all_passed;
% ...and its merged model_defaults are what PSA will actually accept.
psa_n = ParamSpaceAnalysis2('verbose', false);
psa_n.model_class = mcn;
psa_n.model_defaults = merge_struct(dn, cfg_n.model);
[threw, err] = capture_error(@() psa_n.validate_model_defaults());
all_passed = check('noise preset + fast config validates as model_defaults', ...
    ~threw) && all_passed;
if threw; fprintf('      %s\n', err.message); end
[threw, ~] = capture_error(@() psa_n.validate_noise_settings());
all_passed = check('...and passes the sigma/solver pre-flight', ~threw) && all_passed;
% The same preset with the deterministic solver forced must be REJECTED.
psa_bad = ParamSpaceAnalysis2('verbose', false);
psa_bad.model_class = mcn;
psa_bad.model_defaults = merge_struct(dn, struct('ode_solver', 'rk4'));
[threw, err] = capture_error(@() psa_bad.validate_noise_settings());
all_passed = check('sigma > 0 forced onto rk4 is rejected', ...
    threw && contains(err.message, 'sra1')) && all_passed;

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
function names = srnn_param_preset_names_for_test()
% Every preset srnn_param_preset knows about.
%
% Read out of the function rather than kept as a second list here, so a preset
% added without touching this file is still covered. srnn_param_preset_names is
% local to srnn_param_preset.m and cannot be called directly, but the unknown-
% preset error enumerates them, and that message is asserted on a few lines
% above -- so the format is already pinned by this same test.
try
    srnn_param_preset('__definitely_not_a_preset__');
    error('test_srnn_param_preset:NoError', ...
        'srnn_param_preset accepted a bogus name; the list cannot be recovered.');
catch err
    if ~strcmp(err.identifier, 'srnn_param_preset:UnknownPreset')
        rethrow(err);
    end
end
marker = 'Valid presets:';
tail = extractAfter(err.message, marker);
names = strtrim(split(strtrim(erase(tail, '.')), ','))';
end

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

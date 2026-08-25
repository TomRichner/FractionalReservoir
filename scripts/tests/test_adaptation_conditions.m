% TEST_ADAPTATION_CONDITIONS Conditions state their SFA timescales explicitly.
%
% Covers the change that replaced srnn_adaptation_conditions' n_a_sfa COUNT with
% an sfa_timescales VECTOR, so a condition records which timescales it ran
% rather than leaving tau_a to the model's auto-fill.
%
% Two things this is really guarding:
%
%   1. logspace(a,b,1) returns 10^b. The auto-fill this replaces therefore gave
%      the SLOW 10 s end when asked for ONE timescale, silently. The test asserts
%      the trap still exists in MATLAB and that we no longer walk into it.
%   2. n_a and tau_a were two settable properties carrying one fact, and a
%      preset carrying tau_a of a different length from its conditions' n_a
%      failed outright -- which single_neuron_stf actually did.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: srnn_sfa_timescales, srnn_adaptation_conditions, srnn_param_preset

all_passed = true;

%% The ladder helper
fprintf('\n-- srnn_sfa_timescales --\n');
all_passed = check('K=3 is the standard ladder', ...
    isequal(srnn_sfa_timescales(3), logspace(log10(0.25), log10(10), 3))) && all_passed;
all_passed = check('K=1 is the FAST end (0.25)', ...
    isequal(srnn_sfa_timescales(1), 0.25)) && all_passed;
% The trap itself. If MATLAB ever changed this, the special case would be dead
% code and the reader should be told rather than left guessing.
all_passed = check('logspace(a,b,1) really does return the SLOW end', ...
    abs(logspace(log10(0.25), log10(10), 1) - 10) < 1e-12) && all_passed;
all_passed = check('K=0 is a 1x0 row, not []', ...
    isequal(size(srnn_sfa_timescales(0)), [1 0])) && all_passed;
all_passed = check('K=2 spans the same endpoints', ...
    isequal(srnn_sfa_timescales(2), [0.25 10])) && all_passed;

%% Pairs conditions carry tau_a, and it agrees with n_a
fprintf('\n-- SRNNCellTypePairs conditions --\n');
cond = srnn_adaptation_conditions('SRNNCellTypePairs');
all_passed = check('every condition carries tau_a', ...
    all(cellfun(@(c) isfield(c, 'tau_a') && iscell(c.tau_a), cond))) && all_passed;
% Conditions carry tau_a ALONE. n_a is Dependent on the model and read-only, so
% a condition setting it would now throw on every build rather than be ignored.
all_passed = check('no condition carries n_a', ...
    ~any(cellfun(@(c) isfield(c, 'n_a'), cond))) && all_passed;
all_passed = check('SFA conditions carry the ladder', ...
    isequal(cond{2}.tau_a{1}, srnn_sfa_timescales(3)) && ...
    isequal(cond{4}.tau_a{1}, srnn_sfa_timescales(3))) && all_passed;
all_passed = check('non-SFA conditions carry empty timescales', ...
    isempty(cond{1}.tau_a{1}) && isempty(cond{3}.tau_a{1})) && all_passed;
all_passed = check('SFA is on the FIRST type only', ...
    all(cellfun(@(c) all(cellfun(@isempty, c.tau_a(2:end))), cond))) && all_passed;

% A single-timescale request must reach the model as 0.25, not 10.
one = srnn_adaptation_conditions('SRNNCellTypePairs', [], srnn_sfa_timescales(1));
all_passed = check('one-timescale conditions carry 0.25', ...
    isequal(one{2}.tau_a{1}, 0.25) && numel(one{2}.tau_a{1}) == 1) && all_passed;

% C = 1 still produces 1-element rows/cells.
c1 = srnn_adaptation_conditions('SRNNCellTypePairs', [], srnn_sfa_timescales(3), 1);
all_passed = check('C=1 gives 1-element n_a and tau_a', ...
    all(cellfun(@(c) isscalar(c.tau_a), c1))) && all_passed;

%% SRNNModel2 still speaks counts, deliberately
fprintf('\n-- SRNNModel2 branch is untouched --\n');
m2 = srnn_adaptation_conditions('SRNNModel2');
all_passed = check('carries n_a_E / n_b_E, not tau_a', ...
    all(cellfun(@(c) isfield(c,'n_a_E') && isfield(c,'n_b_E') && ...
                     ~isfield(c,'tau_a'), m2))) && all_passed;
all_passed = check('SFA count derives from the timescale vector', ...
    isequal(m2{2}.n_a_E, 3) && isequal(m2{4}.n_a_E, 3)) && all_passed;
m2one = srnn_adaptation_conditions('SRNNModel2', [], srnn_sfa_timescales(1));
all_passed = check('a 1-element vector gives n_a_E = 1', ...
    isequal(m2one{2}.n_a_E, 1)) && all_passed;

%% No preset may carry tau_a: it is condition-owned
fprintf('\n-- presets do not carry tau_a --\n');
names = preset_names();
bad = {};
for k = 1:numel(names)
    d = srnn_param_preset(names{k});
    if isfield(d, 'tau_a'); bad{end+1} = names{k}; end %#ok<SAGROW>
end
all_passed = check(sprintf('no preset carries tau_a (%d checked)', numel(names)), ...
    isempty(bad)) && all_passed;
if ~isempty(bad); fprintf('     offenders: %s\n', strjoin(bad, ', ')); end

%% The regression this fixed: every preset can build every one of its conditions
% single_neuron_stf carried tau_a = {3} while its no_adaptation / std_only
% conditions set n_a = 0, so those two combinations threw
% "tau_a{1} must contain n_a(1) positive values" -- the figure only worked
% because it truncated tau_a by hand.
fprintf('\n-- every preset builds every one of its conditions --\n');
for k = 1:numel(names)
    [~, ~, cnd] = srnn_param_preset(names{k});
    ok = true; why = '';
    for j = 1:numel(cnd)
        try
            evalc(sprintf('build_one(''%s'', ''%s'');', names{k}, cnd{j}.name));
        catch ME
            ok = false; why = sprintf(' [%s: %s]', cnd{j}.name, ME.message);
            break
        end
    end
    all_passed = check(sprintf('%-24s builds all %d conditions%s', ...
        names{k}, numel(cnd), why), ok) && all_passed;
end

%% ------------------------------------------------------------------------
if all_passed
    fprintf('\ntest_adaptation_conditions: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_adaptation_conditions: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function build_one(preset_name, cond_name)
% Construct only -- build() is where the tau_a/n_a disagreement used to throw,
% and running would cost far more than the check is worth.
m = build_from_preset(preset_name, cond_name, 'lya_method', 'none', ...
    'T_range', [0 1], 'n', 8, 'indegree', 4);
% Only SRNNCellTypePairs has n_a / tau_a. Guard BEFORE dereferencing: MATLAB
% evaluates m.n_a eagerly, so a trailing `|| ~isprop(...)` does not save you.
if isprop(m, 'n_a') && isprop(m, 'tau_a')
    assert(isequal(m.n_a, cellfun(@numel, m.tau_a)), ...
        '%s/%s: n_a disagrees with numel(tau_a).', preset_name, cond_name);
end
end

function names = preset_names()
% The presets this refactor keeps working. Deliberately explicit rather than
% enumerated from srnn_param_preset, so a preset that quietly stops being
% covered shows up as an edit here.
names = {'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', 'bursting_pairs', ...
         'sompolinsky_pairs', 'single_neuron_stf', 'single_neuron_dualStd', ...
         'mc_esn', 'default', 'overconnected'};
end

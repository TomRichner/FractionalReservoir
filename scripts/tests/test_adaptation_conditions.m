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
one = srnn_adaptation_conditions('SRNNCellTypePairs', 'sfa_timescales', srnn_sfa_timescales(1));
all_passed = check('one-timescale conditions carry 0.25', ...
    isequal(one{2}.tau_a{1}, 0.25) && numel(one{2}.tau_a{1}) == 1) && all_passed;

% C = 1 still produces 1-element rows/cells.
c1 = srnn_adaptation_conditions('SRNNCellTypePairs', 'n_cell_types', 1);
all_passed = check('C=1 gives 1-element n_a and tau_a', ...
    all(cellfun(@(c) isscalar(c.tau_a), c1))) && all_passed;

%% The 'timescales' regime set
fprintf('\n-- the 7-regime set --\n');
% Explicit DUAL routes, because default_std_routes() is single-timescale and the
% _oneTS/multi contrast would then be invisible -- std_only and std_only_oneTS
% would be the same thing and the test would assert nothing.
dual = struct();
dual.E.E.std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
dual.I.I.std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
seven = srnn_adaptation_conditions('SRNNCellTypePairs', ...
    'synapse_config', dual, 'regimes', 'timescales');
want = {'no_adaptation','sfa_only_oneTS','sfa_only','std_only_oneTS', ...
        'std_only','sfa3_std1','sfa_and_std'};
got = cellfun(@(c) c.name, seven, 'UniformOutput', false);
% ORDER MATTERS: it is the column order of every sweep figure.
all_passed = check('7 regimes in the intended order', isequal(got, want)) && all_passed;

n_tau = cellfun(@(c) numel(c.tau_a{1}), seven);
all_passed = check('SFA timescale counts are [0 1 3 0 0 3 3]', ...
    isequal(n_tau, [0 1 3 0 0 3 3])) && all_passed;

n_std = cellfun(@(c) std_count(c.synapse_config), seven);
all_passed = check('STD timescale counts are [0 0 0 1 2 1 2]', ...
    isequal(n_std, [0 0 0 1 2 1 2])) && all_passed;

% The one-timescale routes are DERIVED, so they must be the first pair of the
% multi-timescale ones rather than an independently written value.
all_passed = check('_oneTS keeps the FIRST tau_rec/tau_rel pair', ...
    isequal(seven{4}.synapse_config.E.E.std.tau_rec, ...
            seven{5}.synapse_config.E.E.std.tau_rec(1)) && ...
    isequal(seven{4}.synapse_config.E.E.std.tau_rel, ...
            seven{5}.synapse_config.E.E.std.tau_rel(1))) && all_passed;

% The four shared names must mean exactly what they mean in the standard set,
% or results are not comparable between the 4- and 7-condition presets.
std4 = srnn_adaptation_conditions('SRNNCellTypePairs', 'synapse_config', dual);
shared = {'no_adaptation','sfa_only','std_only','sfa_and_std'};
same = true;
for s = shared
    a = seven{strcmp(got, s{1})};
    b = std4{strcmp(cellfun(@(c) c.name, std4, 'UniformOutput', false), s{1})};
    same = same && isequal(a, b);
end
all_passed = check('the 4 shared regimes are identical across both sets', same) && all_passed;

%% The 'single_multi' regime set
fprintf('\n-- the 3-regime set --\n');
% None / one timescale each / many timescales each. The whole point is that the
% only thing varying is HOW MANY timescales carry the adaptation, so the checks
% here are about counts and about the two shared regimes still matching.
three = srnn_adaptation_conditions('SRNNCellTypePairs', ...
    'synapse_config', dual, 'regimes', 'single_multi');
got3 = cellfun(@(c) c.name, three, 'UniformOutput', false);
all_passed = check('3 regimes in the intended order', ...
    isequal(got3, {'no_adaptation','sfa1_std1','sfa3_std2'})) && all_passed;

all_passed = check('SFA timescale counts are [0 1 3]', ...
    isequal(cellfun(@(c) numel(c.tau_a{1}), three), [0 1 3])) && all_passed;
all_passed = check('STD timescale counts are [0 1 2]', ...
    isequal(cellfun(@(c) std_count(c.synapse_config), three), [0 1 2])) && all_passed;

% BOTH mechanisms present together in both adapting regimes -- that is what
% separates this set from 'timescales', where the _only regimes isolate them.
all_passed = check('both adapting regimes carry SFA *and* STD', ...
    numel(three{2}.tau_a{1}) > 0 && std_count(three{2}.synapse_config) > 0 && ...
    numel(three{3}.tau_a{1}) > 0 && std_count(three{3}.synapse_config) > 0) && all_passed;

% sfa1_std1's single timescales must be DERIVED from the multi ones, not written
% out, so retuning the network cannot leave it describing an older one.
all_passed = check('sfa1_std1 takes the FIRST tau_a and the FIRST STD pair', ...
    isequal(three{2}.tau_a{1}, three{3}.tau_a{1}(1)) && ...
    isequal(three{2}.synapse_config.E.E.std.tau_rec, ...
            three{3}.synapse_config.E.E.std.tau_rec(1))) && all_passed;

% no_adaptation is the one name shared across all three sets, and it must still
% mean the same thing.
all_passed = check('no_adaptation is identical to the 7-regime set', ...
    isequal(three{strcmp(got3, 'no_adaptation')}, ...
            seven{strcmp(got, 'no_adaptation')})) && all_passed;

% sfa3_std2 is a second NAME for sfa_and_std, not a second regime. If these ever
% stop matching, one of the two switch cases has been edited in isolation and
% the "identical physics, different label" claim in both files is a lie.
a = three{strcmp(got3, 'sfa3_std2')};
b = seven{strcmp(got, 'sfa_and_std')};
all_passed = check('sfa3_std2 is byte-identical to sfa_and_std apart from the name', ...
    isequal(rmfield(a, 'name'), rmfield(b, 'name'))) && all_passed;
all_passed = check('...and they really do carry different names', ...
    ~strcmp(a.name, b.name)) && all_passed;

% Different titles is the whole reason the second name exists.
t3 = srnn_condition_titles();
all_passed = check('the two names get DIFFERENT display titles', ...
    ~strcmp(t3('sfa3_std2'), t3('sfa_and_std'))) && all_passed;
all_passed = check(sprintf('sfa1_std1 -> "%s"', t3('sfa1_std1')), ...
    strcmp(t3('sfa1_std1'), 'Single-Timescale Adaptation')) && all_passed;
all_passed = check(sprintf('sfa3_std2 -> "%s"', t3('sfa3_std2')), ...
    strcmp(t3('sfa3_std2'), 'Multiple-Timescale Adaptation')) && all_passed;

% Every regime needs a display title or a figure legend falls back to the raw
% snake_case name in one panel while reading properly in another.
titles3 = srnn_condition_titles();
all_passed = check('every 3-regime name has a display title', ...
    all(cellfun(@(n) isKey(titles3, n), got3))) && all_passed;

%% SRNNModel2 still speaks counts, deliberately
fprintf('\n-- SRNNModel2 branch is untouched --\n');
m2 = srnn_adaptation_conditions('SRNNModel2');
all_passed = check('carries n_a_E / n_b_E, not tau_a', ...
    all(cellfun(@(c) isfield(c,'n_a_E') && isfield(c,'n_b_E') && ...
                     ~isfield(c,'tau_a'), m2))) && all_passed;
all_passed = check('SFA count derives from the timescale vector', ...
    isequal(m2{2}.n_a_E, 3) && isequal(m2{4}.n_a_E, 3)) && all_passed;
m2one = srnn_adaptation_conditions('SRNNModel2', 'sfa_timescales', srnn_sfa_timescales(1));
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

function n = std_count(sc)
% How many depression timescales the E->E route of a synapse_config carries.
if isstruct(sc) && isfield(sc, 'E') && isfield(sc.E, 'E') && isfield(sc.E.E, 'std')
    n = numel(sc.E.E.std.tau_rec);
else
    n = 0;
end
end

function names = preset_names()
% The presets this refactor keeps working. Deliberately explicit rather than
% enumerated from srnn_param_preset, so a preset that quietly stops being
% covered shows up as an edit here.
names = {'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond', ...
         'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
         'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond', 'bursting_pairs', ...
         'sompolinsky_pairs', 'single_neuron_stf', 'single_neuron_dualStd', ...
         'mc_pairs_dualStd', 'default', 'overconnected'};
end

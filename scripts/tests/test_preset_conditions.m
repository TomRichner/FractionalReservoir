% TEST_PRESET_CONDITIONS Presets state their own conditions, and state them well.
%
% Was test_adaptation_conditions, which tested srnn_adaptation_conditions -- the
% shared regime-set builder, deleted when each preset took over stating its own
% conditions. What it was really guarding survives here, read off the PRESETS
% instead of off the builder's arguments:
%
%   * conditions carry tau_a explicitly, never n_a (which is Dependent and
%     read-only on SRNNCellTypePairs, so a condition setting it would throw on
%     every build);
%   * the one-timescale regimes are DERIVED from the multi-timescale ones, so
%     retuning a network cannot leave them describing an older one;
%   * regime ORDER, which is the column order of every sweep figure;
%   * every preset can BUILD every one of its own conditions.
%
% Plus what is new: validate_preset_conditions rejects each way a hand-written
% condition set can be malformed. That validator is the safety that moved when
% the shared builder went away -- the builder could not emit a bad set, whereas
% a preset can now be edited into one.
%
% NOT COVERED HERE, deliberately: the exact VALUES. test_preset_golden compares
% every preset's full output against a fixture captured before the refactor, so
% duplicating value checks here would only be a second place to update.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: validate_preset_conditions, test_preset_golden, srnn_param_preset,
%           srnn_condition_titles

all_passed = true;

[~, ~, four]  = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_4cond');
[~, ~, seven] = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_7cond');
[~, ~, three] = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_3cond');
name_of = @(cnd) cellfun(@(c) c.name, cnd, 'UniformOutput', false);

%% Pairs conditions carry tau_a, not n_a
fprintf('\n-- SRNNCellTypePairs conditions --\n');
all_passed = check('every condition carries tau_a as a cell', ...
    all(cellfun(@(c) isfield(c, 'tau_a') && iscell(c.tau_a), four))) && all_passed;
% n_a is Dependent on the model and read-only, so a condition setting it would
% throw on every build rather than be quietly ignored.
all_passed = check('no condition carries n_a', ...
    ~any(cellfun(@(c) isfield(c, 'n_a'), four))) && all_passed;
all_passed = check('SFA is on the FIRST cell type only', ...
    all(cellfun(@(c) all(cellfun(@isempty, c.tau_a(2:end))), four))) && all_passed;
% The empty entries must be 1x0, not 0x0 []. The class accepts both, but the
% golden fixture distinguishes them, so this is what keeps hand-written rows
% consistent with what was captured.
all_passed = check('empty tau_a entries are 1x0, not []', ...
    isequal(size(four{1}.tau_a{1}), [1 0])) && all_passed;

%% The 7-regime preset
fprintf('\n-- the 7-regime preset --\n');
got7 = name_of(seven);
% ORDER MATTERS: it is the column order of every sweep figure.
all_passed = check('7 regimes in the intended order', ...
    isequal(got7, {'no_adaptation','sfa1_std0','sfa3_std0','sfa0_std1', ...
                   'sfa0_std2','sfa3_std1','sfa3_std2'})) && all_passed;
all_passed = check('SFA timescale counts are [0 1 3 0 0 3 3]', ...
    isequal(cellfun(@(c) numel(c.tau_a{1}), seven), [0 1 3 0 0 3 3])) && all_passed;
all_passed = check('STD timescale counts are [0 0 0 1 2 1 2]', ...
    isequal(cellfun(@(c) std_count(c.synapse_config), seven), [0 0 0 1 2 1 2])) && all_passed;

% The one-timescale routes must be the FIRST pair of the multi-timescale ones.
% Written out independently they could drift; derived they cannot.
all_passed = check('_oneTS keeps the FIRST tau_rec/tau_rel pair', ...
    isequal(seven{4}.synapse_config.E.E.std.tau_rec, ...
            seven{5}.synapse_config.E.E.std.tau_rec(1)) && ...
    isequal(seven{4}.synapse_config.E.E.std.tau_rel, ...
            seven{5}.synapse_config.E.E.std.tau_rel(1))) && all_passed;
all_passed = check('sfa_only_oneTS takes the FIRST tau_a', ...
    isequal(seven{2}.tau_a{1}, seven{3}.tau_a{1}(1))) && all_passed;

%% Shared names mean the same thing across presets
% This is what makes a 4- and a 7-condition run comparable. Now that each preset
% writes its own conditions, nothing structural enforces it -- so it is asserted.
fprintf('\n-- shared regime names agree across presets --\n');
got4 = name_of(four);
same = true; differing = {};
for s = {'no_adaptation','sfa3_std0','sfa0_std2','sfa3_std2'}
    a = four{strcmp(got4, s{1})};
    b = seven{strcmp(got7, s{1})};
    if ~isequaln(a, b); same = false; differing{end+1} = s{1}; end %#ok<SAGROW>
end
all_passed = check('the 4 shared regimes are identical in the 4- and 7-cond presets', ...
    same) && all_passed;
if ~same; fprintf(2, '     differing: %s\n', strjoin(differing, ', ')); end

%% The 3-regime preset
fprintf('\n-- the 3-regime preset --\n');
got3 = name_of(three);
all_passed = check('3 regimes in the intended order', ...
    isequal(got3, {'no_adaptation','sfa1_std1','sfa3_std2'})) && all_passed;
all_passed = check('SFA timescale counts are [0 1 3]', ...
    isequal(cellfun(@(c) numel(c.tau_a{1}), three), [0 1 3])) && all_passed;
all_passed = check('STD timescale counts are [0 1 2]', ...
    isequal(cellfun(@(c) std_count(c.synapse_config), three), [0 1 2])) && all_passed;
% BOTH mechanisms present together in both adapting regimes -- that is what
% separates this set from the 7-regime one, where the _only regimes isolate them.
all_passed = check('both adapting regimes carry SFA *and* STD', ...
    numel(three{2}.tau_a{1}) > 0 && std_count(three{2}.synapse_config) > 0 && ...
    numel(three{3}.tau_a{1}) > 0 && std_count(three{3}.synapse_config) > 0) && all_passed;
all_passed = check('sfa1_std1 takes the FIRST tau_a and the FIRST STD pair', ...
    isequal(three{2}.tau_a{1}, three{3}.tau_a{1}(1)) && ...
    isequal(three{2}.synapse_config.E.E.std.tau_rec, ...
            three{3}.synapse_config.E.E.std.tau_rec(1))) && all_passed;

% ONE NAME, ONE REGIME, ACROSS PRESETS. sfa3_std2 used to be an alias invented so
% the 3-condition set could title its full regime differently from the 7-condition
% set's sfa_and_std. Now both presets simply call it sfa3_std2, and the two must be
% the same struct -- that is the invariant the whole naming scheme buys, and the
% one the old mechanism names could not express (sfa_and_std meant four different
% things across the ten presets).
all_passed = check('sfa3_std2 is identical in the 3- and 7-cond presets', ...
    isequaln(three{strcmp(got3, 'sfa3_std2')}, ...
             seven{strcmp(got7, 'sfa3_std2')})) && all_passed;
all_passed = check('no_adaptation is identical across the 3- and 7-cond presets', ...
    isequaln(three{strcmp(got3, 'no_adaptation')}, ...
             seven{strcmp(got7, 'no_adaptation')})) && all_passed;

% Titles are mostly structural, with these two deliberate exceptions: they are
% the 3-condition preset's adapting regimes, i.e. the paper's comparison, so they
% are titled in the manuscript's words rather than by their counts.
t = srnn_condition_titles();
all_passed = check(sprintf('sfa1_std1 -> "%s"', t('sfa1_std1')), ...
    strcmp(t('sfa1_std1'), 'Single-Timescale Adaptation')) && all_passed;
all_passed = check(sprintf('sfa3_std2 -> "%s"', t('sfa3_std2')), ...
    strcmp(t('sfa3_std2'), 'Multiple-Timescale Adaptation')) && all_passed;
all_passed = check('the two adapting regimes are titled differently', ...
    ~strcmp(t('sfa1_std1'), t('sfa3_std2'))) && all_passed;
% Legacy names stay in the map, or a figure regenerated from a run saved before
% the rename falls back to raw snake_case panel titles.
all_passed = check('legacy names still resolve to titles', ...
    all(cellfun(@(n) isKey(t, n), ...
        {'sfa_only','std_only','sfa_and_std','sfa_only_oneTS','std_only_oneTS'}))) ...
    && all_passed;

%% SRNNModel2 still speaks counts, deliberately
fprintf('\n-- SRNNModel2 presets speak counts --\n');
[~, m2_class, m2] = srnn_param_preset('default');
all_passed = check('''default'' is an SRNNModel2 preset', ...
    strcmp(m2_class, 'SRNNModel2')) && all_passed;
all_passed = check('carries n_a_E / n_b_E, not tau_a', ...
    all(cellfun(@(c) isfield(c,'n_a_E') && isfield(c,'n_b_E') && ...
                     ~isfield(c,'tau_a'), m2))) && all_passed;
all_passed = check('the SFA regimes ask for 3 timescales', ...
    isequal(m2{2}.n_a_E, 3) && isequal(m2{4}.n_a_E, 3)) && all_passed;

%% The validator rejects every malformed shape
% The safety that moved when the shared builder went away: it could not emit a
% bad set, whereas a preset can now be edited into one.
fprintf('\n-- validate_preset_conditions rejects --\n');
cls = 'SRNNCellTypePairs';
all_passed = rejects('a non-cell', ...
    @() validate_preset_conditions(three{1}, cls, 'p'), 'NotACell') && all_passed;
all_passed = rejects('an empty set', ...
    @() validate_preset_conditions({}, cls, 'p'), 'NoConditions') && all_passed;
% THE BRACE TRAP: struct('tau_a', <cell>) builds a 1xN struct ARRAY, one element
% per cell entry, rather than a scalar struct with a cell field. This is the one
% real hazard in writing conditions by hand, so it is checked explicitly.
bad = three; bad{2} = struct('name','x','tau_a',three{2}.tau_a);
all_passed = rejects('a struct ARRAY (the brace trap)', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'BadShape') && all_passed;
bad = three; bad{2} = rmfield(bad{2}, 'name');
all_passed = rejects('a condition with no name', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'BadName') && all_passed;
bad = three; bad{2}.name = bad{1}.name;
all_passed = rejects('a duplicate name', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'DuplicateName') && all_passed;
% The omission invariant: only the fields a condition ACTUALLY SETS are
% condition-owned, so one condition dropping a field silently inherits
% model_defaults instead of its siblings' value.
bad = three; bad{2} = rmfield(bad{2}, 'synapse_config');
all_passed = rejects('conditions that set different fields', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'FieldSetMismatch') && all_passed;
bad = three; for i = 1:numel(bad); bad{i}.synpase_config = struct(); end
all_passed = rejects('a typo''d field name', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'UnknownField') && all_passed;
bad = three; for i = 1:numel(bad); bad{i}.n_a = 1; end
all_passed = rejects('a Dependent, read-only property', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'UnknownField') && all_passed;
bad = three; bad{2}.name = 'not_a_titled_regime';
all_passed = rejects('a name with no display title', ...
    @() validate_preset_conditions(bad, cls, 'p'), 'NoTitle') && all_passed;

%% Names that claim timescale counts must tell the truth
% NOT "every condition must be named sfaX_stdY" -- conditions are general, and a
% preset comparing level_of_chaos or S_c would rightly use other names. This
% checks only names that OPT IN by matching the pattern, so the convention stays
% a convention. That is also why it lives here rather than in
% validate_preset_conditions, which stays free of any domain knowledge.
%
% The counts are what full_adaptation_condition sorts on, so a lying name would
% silently send the figures to the wrong regime.
fprintf('\n-- sfaX_stdY names match the actual timescale counts --\n');
liars = {};
for nm = preset_names()
    [~, cls, cnd] = srnn_param_preset(nm{1});
    for j = 1:numel(cnd)
        [claim, actual] = name_vs_counts(cnd{j}, cls);
        if isempty(claim); continue; end        % not a conforming name; skipped
        if ~isequal(claim, actual)
            liars{end+1} = sprintf('%s/%s claims [%s] but carries [%s]', ...
                nm{1}, cnd{j}.name, num2str(claim), num2str(actual)); %#ok<SAGROW>
        end
    end
end
all_passed = check('every sfaX_stdY name matches its condition', isempty(liars)) ...
    && all_passed;
if ~isempty(liars); fprintf(2, '     %s\n', strjoin(liars, '\n     ')); end

% The check must actually fire, or it is decoration. Take a real condition and
% rename it to claim counts it does not have.
[~, cls_p, cnd_p] = srnn_param_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_3cond');
mis = cnd_p{2};                       % genuinely sfa1_std1
mis.name = 'sfa3_std2';               % now a lie
[claim_m, actual_m] = name_vs_counts(mis, cls_p);
all_passed = check('a deliberately mis-named condition is caught', ...
    ~isequal(claim_m, actual_m)) && all_passed;

%% No preset may carry tau_a: it is condition-owned
fprintf('\n-- presets do not carry tau_a --\n');
names = preset_names();
bad_presets = {};
for k = 1:numel(names)
    d = srnn_param_preset(names{k});
    if isfield(d, 'tau_a'); bad_presets{end+1} = names{k}; end %#ok<SAGROW>
end
all_passed = check(sprintf('no preset carries tau_a (%d checked)', numel(names)), ...
    isempty(bad_presets)) && all_passed;
if ~isempty(bad_presets)
    fprintf(2, '     offenders: %s\n', strjoin(bad_presets, ', '));
end

%% Every preset builds every one of its conditions
% The regression this originally fixed: single_neuron_stf carried tau_a = {3}
% while its no_adaptation / std_only conditions set n_a = 0, so those two
% combinations threw "tau_a{1} must contain n_a(1) positive values" -- the figure
% only worked because it truncated tau_a by hand.
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
    fprintf('\ntest_preset_conditions: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_preset_conditions: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function ok = rejects(label, fn, want_suffix)
% The identifier must match, not merely "it threw": a validator that rejected
% everything for the wrong reason would otherwise pass every check here.
try
    fn();
    ok = false;
    fprintf(2, '  FAIL  %s -- accepted it\n', label);
catch ME
    ok = endsWith(ME.identifier, want_suffix);
    if ok
        fprintf('  PASS  %s -> %s\n', label, ME.identifier);
    else
        fprintf(2, '  FAIL  %s -> got %s, want *:%s\n', label, ME.identifier, want_suffix);
    end
end
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

function [claim, actual] = name_vs_counts(cond, model_class)
% What a sfaX_stdY[_stfZ] name CLAIMS, and what the condition actually carries,
% as [n_sfa n_std n_stf]. claim is [] for a name that does not match the pattern,
% which is how non-conforming names opt out of the check entirely.
% Two exact patterns, not one with an optional trailing group: MATLAB's regexp
% does not return a token for a trailing optional group even when it matches, so
% '_stf1' silently parsed as zero. This test caught that bug in
% full_adaptation_condition, which had the same defect.
claim = []; actual = [];
t3 = regexp(cond.name, '^sfa(\d+)_std(\d+)_stf(\d+)$', 'tokens', 'once');
t2 = regexp(cond.name, '^sfa(\d+)_std(\d+)$', 'tokens', 'once');
if ~isempty(t3)
    claim = cellfun(@str2double, t3);
elseif ~isempty(t2)
    claim = [cellfun(@str2double, t2), 0];
else
    return
end

if strcmp(model_class, 'SRNNModel2')
    % This class speaks COUNTS directly, and has no facilitation.
    actual = [cond.n_a_E, cond.n_b_E, 0];
    return
end

n_sfa = numel(cond.tau_a{1});
n_std = 0; n_stf = 0;
sc = cond.synapse_config;
if isstruct(sc) && ~isempty(fieldnames(sc))
    pre = fieldnames(sc);
    for a = 1:numel(pre)
        post = fieldnames(sc.(pre{a}));
        for b = 1:numel(post)
            r = sc.(pre{a}).(post{b});
            if isfield(r, 'std'); n_std = max(n_std, numel(r.std.tau_rec)); end
            if isfield(r, 'stf'); n_stf = max(n_stf, numel(r.stf.tau_fac)); end
        end
    end
end
actual = [n_sfa, n_std, n_stf];
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
% Deliberately explicit rather than enumerated from srnn_param_preset, so a
% preset that quietly stops being covered shows up as an edit here.
names = {'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond', ...
         'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
         'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond', 'bursting_pairs', ...
         'sompolinsky_pairs', 'single_neuron_stf', 'single_neuron_dualStd', ...
         'mc_pairs_dualStd', 'default', 'overconnected'};
end

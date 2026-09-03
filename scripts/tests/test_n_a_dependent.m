% TEST_N_A_DEPENDENT tau_a is authoritative; n_a is derived and read-only.
%
% n_a and tau_a used to be two settable properties carrying one fact, kept in
% agreement by a validate() check and bridged by an auto-fill that invented
% timescales from the count. n_a is now Dependent with private set access:
%
%   n_a = cellfun(@numel, tau_a)
%
% which is what the synapse side has always done (n_b = numel(tau_rec)).
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs, log_ladder, srnn_adaptation_conditions

all_passed = true;

%% n_a tracks tau_a
fprintf('\n-- n_a is derived from tau_a --\n');
m = pairs_model({log_ladder(0.25, 10, 3), zeros(1,0)});
all_passed = check('n_a = [3 0] from a 3-element ladder', ...
    isequal(m.n_a, [3 0])) && all_passed;

m2 = pairs_model({log_ladder(0.25, 10, 2), log_ladder(0.25, 10, 1)});
all_passed = check('n_a = [2 1] from [2-ladder, 1-ladder]', ...
    isequal(m2.n_a, [2 1])) && all_passed;

m0 = pairs_model({zeros(1,0), zeros(1,0)});
all_passed = check('n_a = [0 0] with no timescales', isequal(m0.n_a, [0 0])) && all_passed;

% Changing tau_a must move n_a with it, with no rebuild and no bookkeeping.
m0.tau_a = {[1 2 3 4], zeros(1,0)};
all_passed = check('n_a follows a later tau_a assignment', ...
    isequal(m0.n_a, [4 0])) && all_passed;

%% n_a cannot be set
fprintf('\n-- n_a is read-only --\n');
all_passed = check('assigning n_a throws SetProhibited', ...
    throws_id(@() set_n_a(m), 'MATLAB:class:SetProhibited')) && all_passed;
% Which also means it can never become a sweep axis, silently overriding every
% adaptation condition -- the hazard condition_field_names() only ever guarded
% for SRNNModel2's n_a_E.
all_passed = check('n_a is not a settable property of the class', ...
    ~ismember('n_a', settable_properties('SRNNCellTypePairs'))) && all_passed;
all_passed = check('tau_a IS settable', ...
    ismember('tau_a', settable_properties('SRNNCellTypePairs'))) && all_passed;

%% The auto-fill is gone
fprintf('\n-- no timescales are invented --\n');
% Previously an unset tau_a with n_a = 3 produced logspace(...), and n_a = 1
% produced the slow 10 s end. There is now nothing to infer from.
mn = pairs_model([]);          % tau_a unset entirely
all_passed = check('unset tau_a defaults to NO adaptation', ...
    isequal(mn.n_a, [0 0])) && all_passed;
all_passed = check('...and tau_a is a 1xC cell of 1x0 rows', ...
    iscell(mn.tau_a) && numel(mn.tau_a) == 2 && ...
    all(cellfun(@(t) isequal(size(t), [1 0]), mn.tau_a))) && all_passed;

%% Validation now reports tau_a, not an n_a disagreement
fprintf('\n-- validation --\n');
all_passed = check('a non-cell tau_a is rejected', ...
    throws_id(@() pairs_model(3), 'SRNNCellTypePairs:InvalidParams')) && all_passed;
all_passed = check('a wrong-length tau_a is rejected', ...
    throws_id(@() pairs_model({log_ladder(0.25, 10, 3)}), ...
              'SRNNCellTypePairs:InvalidParams')) && all_passed;
all_passed = check('a negative timescale is rejected', ...
    throws_id(@() pairs_model({[-1 2], zeros(1,0)}), ...
              'SRNNCellTypePairs:InvalidParams')) && all_passed;

%% Conditions no longer carry n_a at all
fprintf('\n-- conditions --\n');
cond = srnn_adaptation_conditions('SRNNCellTypePairs');
all_passed = check('no condition sets n_a', ...
    ~any(cellfun(@(c) isfield(c, 'n_a'), cond))) && all_passed;
all_passed = check('every condition sets tau_a', ...
    all(cellfun(@(c) isfield(c, 'tau_a'), cond))) && all_passed;
% Setting n_a would now throw, so a condition carrying it would break every
% build rather than being quietly ignored.
all_passed = check('a built model still reports n_a for reading', ...
    isequal(build_from_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
        'sfa_and_std', 'n', 8, 'indegree', 4, 'lya_method','none', ...
        'T_range',[0 1]).n_a, [3 0])) && all_passed;

%% C = 1 convenience: a bare numeric tau_a
fprintf('\n-- C = 1 --\n');
m1 = SRNNCellTypePairs('n_cellTypes', 1, 'cell_type_names', {'E'}, 'f', 1, ...
    'mu_tilde_relative', 0, 'sigma_tilde_relative', 1, 'n', 8, 'indegree', 4, ...
    'tau_a', 0.25, 'c', 0.5, 'lya_method', 'none', 'T_range', [0 1]);
m1.build();
all_passed = check('a numeric tau_a is wrapped at C = 1', ...
    iscell(m1.tau_a) && isequal(m1.tau_a{1}, 0.25) && isequal(m1.n_a, 1)) && all_passed;

%% ------------------------------------------------------------------------
if all_passed
    fprintf('\ntest_n_a_dependent: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_n_a_dependent: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function m = pairs_model(tau_a)
% A tiny two-type model. tau_a is passed straight through so the validation
% checks can hand it bad values.
args = {'n_cellTypes', 2, 'cell_type_names', {'E','I'}, 'f', [0.5 0.5], ...
        'mu_tilde_relative', [0 0; 0 0], 'sigma_tilde_relative', [1 1; 1 1], ...
        'n', 8, 'indegree', 4, 'c', [0.5 0], ...
        'lya_method', 'none', 'T_range', [0 1]};
if ~isempty(tau_a); args = [args, {'tau_a', tau_a}]; end
m = SRNNCellTypePairs(args{:});
evalc('m.build()');
end

function set_n_a(m)
m.n_a = [1 0];
end

function names = settable_properties(class_name)
% Public properties with public set access -- the same notion ParamSpaceAnalysis2
% uses to decide what may be a grid axis.
mc = meta.class.fromName(class_name);
p = mc.PropertyList;
keep = strcmp({p.GetAccess}, 'public') & strcmp({p.SetAccess}, 'public');
names = {p(keep).Name};
end

function tf = throws_id(f, id)
tf = false;
try
    f();
catch ME
    tf = strcmp(ME.identifier, id);
end
end

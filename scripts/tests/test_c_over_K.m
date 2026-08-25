% TEST_C_OVER_K Adaptation is scaled by c/K, and the change moved no numbers.
%
% c is the TOTAL adaptation budget; the model divides it by the number of SFA
% timescales in use. Every a_k relaxes to the rate, so sum_k a_k -> K*r at steady
% state whatever the taus are, and an unnormalized c*sum would scale linearly
% with K -- meaning a condition that changed the number of timescales would
% silently change adaptation STRENGTH too.
%
% THE FROZEN TABLE IS THE LOAD-BEARING PART. Rescaling every preset's c by its
% own K is supposed to leave the physics untouched, and "supposed to" is not
% evidence. These checksums were captured by RUNNING THE PRE-CHANGE CODE at
% commit cbcf637; they are not derived from anything in the current tree, so they
% cannot drift with it. If one moves, a preset's adaptation strength moved.
%
% The checksum covers x, r, a, b AND g. An x-only checksum is blind here: the
% single-neuron presets have W == 0, so x never sees r or b and every condition
% would report as identical.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs/get_params, SRNNModel2.effective_c

all_passed = true;

%% Frozen pre-change trajectory checksums
fprintf('\n-- physics unchanged (frozen at cbcf637) --\n');
% preset, n/indegree override, then one {condition, sum, norm} row per condition.
F = {
 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', {'n',60,'indegree',20}, {
    'no_adaptation',  5348.42091008, 141.460065715
    'sfa_only',       7714.98224774, 131.342709836
    'std_only',       42887.9484524, 162.324200553
    'sfa_and_std',    48209.9685713, 170.875882894 }
 'bursting_pairs', {}, {
    'no_adaptation',  56036.7024234, 429.243870676
    'sfa_only',       68236.1842781, 427.942567427
    'std_only',       20783.9501838, 127.513762718
    'sfa_and_std',    23035.5604914, 131.024526697 }
 'sompolinsky_pairs', {'n',60,'indegree',60}, {
    'no_adaptation',  2360.92676217, 129.267097157
    'sfa_only',       4205.09901870, 152.628453945 }
 'single_neuron_stf', {}, {
    'sfa_only',       93.9903825304, 5.55519680910
    'sfa_and_std',    647.781482404, 24.0779561777 }
 'single_neuron_dualStd', {}, {
    'no_adaptation',  175.026795357, 9.43683102318
    'sfa_only',       320.505680170, 10.9637823696
    'std_only',       424.368611267, 15.0522656270
    'sfa_and_std',    598.034634373, 16.6784274754 }
 'mc_esn', {'n',60,'indegree',20}, {
    'no_adaptation',  19107.4291486, 376.340171671
    'sfa_only',       18584.4201121, 336.323456565
    'std_only',       8176.46215537, 177.348510480
    'sfa_and_std',    10877.1677395, 169.120514596 }
 'overconnected', {'n',60,'indegree',40}, {
    'no_adaptation',  22173.2222384, 201.624718596
    'sfa_only',       21175.5140168, 170.518918800
    'sfa_and_std',    12513.0287070, 139.328370148 }
 'default', {'n',60,'indegree',20}, {
    'no_adaptation',  18377.7826083, 193.598106256
    'sfa_and_std',    14846.3717166, 132.179501549 }
};
for i = 1:size(F,1)
    rows = F{i,3};
    for k = 1:size(rows,1)
        [s, nrm] = traj(F{i,1}, rows{k,1}, F{i,2});
        ok = rel_eq(s, rows{k,2}) && rel_eq(nrm, rows{k,3});
        all_passed = check(sprintf('%-46s %-14s', F{i,1}, rows{k,1}), ok) && all_passed;
        if ~ok
            fprintf('      got  %.12g / %.12g\n      want %.12g / %.12g\n', ...
                s, nrm, rows{k,2}, rows{k,3});
        end
    end
end

%% The budget really is K-invariant
% The point of the change, stated as a property rather than a checksum: hold c
% fixed, vary K, and the steady-state adaptation must not move.
fprintf('\n-- total adaptation budget is independent of K --\n');
budget = nan(1,3);
for K = 1:3
    m = build_from_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
        'sfa_and_std', 'tau_a', {srnn_sfa_timescales(K), zeros(1,0)}, ...
'n', 8, 'indegree', 4, 'lya_method','none', 'T_range',[0 1]);
    p = m.get_params();
    budget(K) = p.c_eff(1) * K;      % (c/K)*sum(a) at steady state, where sum(a) = K*r
end
all_passed = check(sprintf('budget = %.4f / %.4f / %.4f at K = 1/2/3', budget), ...
    all(abs(budget - budget(1)) < 1e-12)) && all_passed;
all_passed = check('budget equals the preset''s c (0.5)', ...
    abs(budget(1) - 0.5) < 1e-12) && all_passed;

%% K = 0 must not divide by zero
fprintf('\n-- K = 0 --\n');
m0 = build_from_preset('celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
    'no_adaptation', 'n', 8, 'indegree', 4, 'lya_method','none', 'T_range',[0 1]);
p0 = m0.get_params();
all_passed = check('c_eff is finite at n_a = 0', all(isfinite(p0.c_eff))) && all_passed;
% SRNNModel2.effective_c is a private static, so exercise it through a built
% model instead of reaching into the class. A 0/0 there would surface as NaN in
% the trajectory rather than as an error.
[s0, n0] = traj('default', 'no_adaptation', {'n',8,'indegree',4});
all_passed = check('SRNNModel2 runs finite at n_a_E = 0', ...
    isfinite(s0) && isfinite(n0)) && all_passed;

%% Retired presets refuse rather than silently running at c/3
fprintf('\n-- retired presets --\n');
% NOTE the bare '..._dualStd' here, not '..._dualStd_4cond'. The retired parent
% and its live successor differ by a suffix, so a blanket rename over the repo
% will silently turn this check into "the live preset must error" -- which it did
% once already.
for nm = {'celltype_pairs', 'celltype_pairs_Sc0p2_noise0p025', ...
          'celltype_pairs_Sc0p2_noise0p025_dualStd'}
    all_passed = check(sprintf('%-46s errors', nm{1}), ...
        throws_id(@() srnn_param_preset(nm{1}), 'srnn_param_preset:RetiredPreset')) ...
        && all_passed;
end

%% ------------------------------------------------------------------------
if all_passed
    fprintf('\ntest_c_over_K: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_c_over_K: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function tf = rel_eq(got, want)
% The frozen values are printed to 12 significant figures, so compare relatively
% at that precision rather than demanding exact equality with a decimal literal.
tf = abs(got - want) <= 1e-9 * max(1, abs(want));
end

function [s, nrm] = traj(preset_name, cond_name, extra)
args = [{preset_name, cond_name}, extra, ...
        {'lya_method','none', 'T_range',[0 3], 'fs',100, 'rng_seeds',[42 42], ...
         'store_decimated_state',true, 'plot_freq',100}];
evalc('m = build_from_preset(args{:}); m.run();');
v = harvest(m.plot_data);
s = sum(v); nrm = norm(v);
end

function v = harvest(p)
v = [];
for f = {'x','r','a','b','g'}
    if isfield(p, f{1}); v = [v; flatten(p.(f{1}))]; end %#ok<AGROW>
end
end

function v = flatten(s)
v = [];
if isnumeric(s)
    v = s(:);
elseif isstruct(s)
    f = fieldnames(s);
    for j = 1:numel(f); v = [v; flatten(s.(f{j}))]; end %#ok<AGROW>
elseif iscell(s)
    for j = 1:numel(s); v = [v; flatten(s{j})]; end %#ok<AGROW>
end
end

function tf = throws_id(f, id)
tf = false;
try
    f();
catch ME
    tf = strcmp(ME.identifier, id);
end
end

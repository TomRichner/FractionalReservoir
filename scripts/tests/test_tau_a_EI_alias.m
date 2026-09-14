% test_tau_a_EI_alias.m - the tau_a_EI alias and the E+I tau sweep.
%
% Checks:
%   1. Setting tau_a_EI writes the ladder to BOTH types; tau_a_E still
%      touches type 1 only; reading tau_a_EI after tau_a_E moved E alone
%      errors SRNNCellTypePairs:TauAMismatch; a 3-type model refuses the
%      setter (TwoTypeAliasOnly).
%   2. The constructor accepts 'tau_a_EI' as a name-value (the sweep passes
%      it that way) and the model builds with equal ladders.
%   3. run_tau_sensitivity_analysis on the sfaEI scaled preset at 'fast' with
%      n = 40, 2 levels, 2 reps: chooses tau_a_EI (the condition's
%      ladders are identical on E and I) and writes tau_levels.mat / .md whose
%      E and I rows are equal at every level and whose last element is the
%      swept value; on the E-only 3cond_mu8p25 preset it chooses tau_a_E.
%
% ~1-2 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: run_tau_sensitivity_analysis, fig_sfa_EOC_allStd, SRNNCellTypePairs

fprintf('=== Testing tau_a_EI ===\n\n');
all_passed = true;
PEI = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdScaled_3cond_mu8p25';
PE  = 'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond_mu8p25';

%% 1. alias semantics
m = SRNNCellTypePairs('n', 10, 'indegree', 3, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [0.1 -0.1], 'sigma_tilde_relative', [0.01 0.01], ...
    'tau_a', {[0.25 1 4], [0.25 1 4]}, 'c', [0.5 0.5], 'T_range', [0 0.1], 'lya_method', 'none');
m.tau_a_EI = [0.25 2 9];
all_passed = check('tau_a_EI writes both ladders', isequal(m.tau_a{1}, [0.25 2 9]) && isequal(m.tau_a{2}, [0.25 2 9]) ...
    && isequal(m.tau_a_EI, [0.25 2 9])) && all_passed;
m.tau_a_E = [0.25 3 12];
all_passed = check('tau_a_E still moves type 1 only', isequal(m.tau_a{1}, [0.25 3 12]) && isequal(m.tau_a{2}, [0.25 2 9])) && all_passed;
all_passed = check('reading tau_a_EI with unequal ladders errors TauAMismatch', ...
    throws_id(@() m.tau_a_EI, 'SRNNCellTypePairs:TauAMismatch')) && all_passed;
m3 = SRNNCellTypePairs('n', 12, 'indegree', 3, 'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], 'mu_tilde_relative', [0.1 -0.2 -0.1], 'sigma_tilde_relative', [0.01 0.01 0.01], ...
    'T_range', [0 0.1], 'lya_method', 'none');
all_passed = check('three types refuse the setter', throws_id(@() set_ei(m3, [0.25 1]), 'SRNNCellTypePairs:TwoTypeAliasOnly')) && all_passed;

%% 2. constructor name-value
m2 = build_from_preset(PEI, 'sfa3_std2', 'n', 20, 'indegree', 5, 'T_range', [0 0.2], 'lya_method', 'none', 'tau_a_EI', [0.25 1.5 6]);
all_passed = check('constructor accepts tau_a_EI and both ladders carry it', ...
    isequal(m2.tau_a{1}, [0.25 1.5 6]) && isequal(m2.tau_a{2}, [0.25 1.5 6])) && all_passed;

%% 3. the sweep chooses the axis by the condition
out_root = fullfile(tempdir, 'tau_a_EI_test');
if exist(out_root, 'dir'); rmdir(out_root, 's'); end
picked = struct();
for P = {PEI, PE}
    ctx = resolve_run_context('tau_sensitivity', 'preset_name', P{1}, 'run_mode', 'fast', ...
        'output_dir', fullfile(out_root, P{1}(end-25:end)), 'save_figs', false, 'verbose', 'near-none');
    ctx.n_levels = 2; ctx.n_reps = 2;   % reps is a grid axis and needs >= 2 values
    ctx.model_defaults.n = 40; ctx.model_defaults.indegree = 8;
    ctx.model_defaults.T_range = [0 6]; ctx.model_defaults.lya_T_interval = [3 6]; ctx.model_defaults.lya_warmup = 1;
    ctx.model_defaults.lya_K = 2; ctx.model_defaults.lya_K_auto = false;
    ctx.output_dir = fullfile(out_root, P{1}(end-25:end));
    zz = evalc('d = run_tau_sensitivity_analysis(ctx);'); %#ok<NASGU>
    L = load(fullfile(d, 'tau_levels.mat'));
    picked.(matlab.lang.makeValidName(P{1}(end-25:end))) = L;
    fprintf('  %s -> axis %s, %d levels\n', P{1}, L.tau_axis, numel(L.levels));
end
LEI = picked.(matlab.lang.makeValidName(PEI(end-25:end)));
LE  = picked.(matlab.lang.makeValidName(PE(end-25:end)));
ok_ei = strcmp(LEI.tau_axis, 'tau_a_EI') && numel(LEI.levels) == 2 && ...
    all(arrayfun(@(l) isequal(l.tau_a_E, l.tau_a_I) && l.tau_a_E(end) == l.swept(end) && l.tau_a_E(1) == 0.25, LEI.levels));
ok_e = strcmp(LE.tau_axis, 'tau_a_E') && all(arrayfun(@(l) isempty(l.tau_a_I) && l.tau_a_E(end) == l.swept(end), LE.levels));
all_passed = check('sfaEI preset: tau_a_EI, E and I ladders equal at every level, slowest = swept', ok_ei) && all_passed;
all_passed = check('E-only preset: tau_a_E, I ladder empty', ok_e) && all_passed;
all_passed = check('tau_levels.md written', isfile(fullfile(fileparts(fullfile(out_root, PEI(end-25:end), 'x')), 'tau_levels.md')) || ...
    ~isempty(dir(fullfile(out_root, '**', 'tau_levels.md')))) && all_passed;
rmdir(out_root, 's');

fprintf('\n');
if all_passed
    fprintf('=== ALL tau_a_EI TESTS PASSED ===\n');
else
    fprintf('=== SOME tau_a_EI TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function set_ei(m, v)
m.tau_a_EI = v;
end

function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
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

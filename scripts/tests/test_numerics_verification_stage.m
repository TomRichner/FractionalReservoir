% test_numerics_verification_stage.m - run_numerics_verification end to end
% at the smallest counts, then the acceptance verdict and the figures.
%
% Checks:
%   1. A bad run mode errors :badMode (also in test_run_modes).
%   2. 'fast' with 1 reshoot trial and 2 LLE/QR trials into a temp run
%      directory: the .mat has results / cond_names / condition_titles /
%      settings / verdict; every condition carries reshoot, lle, qr, jacobian
%      (3 states, finite rel_fro_err, eig_err, max_abs_err, n_excluded >= 0,
%      N = the reduced network's state count) and a summary with
%      jac_rel_fro_max / jac_eig_err_max / jac_excluded_total;
%      settings.acceptance carries the pre-registered thresholds (fixed_on
%      2026-09-14); verdict has per_condition with the five checks A B L C J,
%      each with value / threshold / pass / text, and all_pass; the Jacobian
%      check itself passes (rel Frobenius <= 1e-6, eig <= 1e-6) -- it is the
%      one criterion the test asserts, since the others are finite-time
%      measurements with a scatter the fast mode cannot resolve;
%      numerics_verdict.md exists and names every check.
%   3. The figure variants 'solver', 'lya_method', 'ensemble' and 'jacobian'
%      write files; 'jacobian' errors NoJacobian on a .mat without the field.
%
% RUNTIME: the fast cost table still builds n = 500 networks for checks A/B
% (one seed, T_free 6 s at 400/800/1600 Hz plus an ode45 reference at 1e-10)
% and L (two seeds x two integrators, 10 s each), so expect ~5-10 min on a
% 12-worker pool, longer serial. Prints PASS/FAIL per check and a final
% banner. Assumes setup_paths.
%
% See also: run_numerics_verification, fig_numerics_verification,
%           test_numerics_probe, test_run_modes

fprintf('=== Testing run_numerics_verification + figures ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdScaled_3cond_mu8p25';

all_passed = check('bad run mode errors :badMode', ...
    throws_id(@() run_numerics_verification('run_mode', 'nope'), 'run_numerics_verification:badMode')) && all_passed;

run_dir = fullfile(tempdir, 'numerics_stage_test');
if exist(run_dir, 'dir'); rmdir(run_dir, 's'); end
mkdir(run_dir);
t0 = tic;
out = evalc(['mat_file = run_numerics_verification(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(run_dir, ''numerics_verification''), ''n_trials_reshoot'', 1, ''n_trials_lle'', 2);']);
fprintf('  (stage ran in %.0f s)\n', toc(t0));
all_passed = check('stage returns the .mat path and it exists', ischar(mat_file) && isfile(mat_file)) && all_passed;
D = load(mat_file);
all_passed = check('schema: results / cond_names / condition_titles / settings / verdict', ...
    all(isfield(D, {'results', 'cond_names', 'condition_titles', 'settings', 'verdict'}))) && all_passed;
R = D.results;
all_passed = check('every condition carries reshoot, lle, qr, jacobian, summary', ...
    numel(R) == 3 && all(isfield(R, {'reshoot', 'lle', 'qr', 'jacobian', 'summary'}))) && all_passed;
jac_ok = true; sum_ok = true;
for i = 1:numel(R)
    Jc = R(i).jacobian;
    jac_ok = jac_ok && numel(Jc) == 3 && all(isfield(Jc, {'rel_fro_err', 'eig_err', 'max_abs_err', 'n_excluded', 'N', 't', 'h'})) && ...
        all(isfinite([Jc.rel_fro_err])) && all(isfinite([Jc.eig_err])) && all([Jc.n_excluded] >= 0) && ...
        all([Jc.N] == Jc(1).N) && all(isfinite([Jc.t]));   % T_small starts negative, so t can be < 0
    S = R(i).summary;
    sum_ok = sum_ok && all(isfield(S, {'jac_rel_fro_max', 'jac_eig_err_max', 'jac_excluded_total'})) && ...
        S.jac_rel_fro_max == max([Jc.rel_fro_err]);
    fprintf('  %-14s Jacobian rel Frobenius %s, eig %s, excluded %s\n', R(i).name, ...
        mat2str([Jc.rel_fro_err], 2), mat2str([Jc.eig_err], 2), mat2str([Jc.n_excluded]));
end
all_passed = check('jacobian: 3 states per condition with finite errors', jac_ok) && all_passed;
all_passed = check('summary carries the Jacobian scalars', sum_ok) && all_passed;
acc = D.settings.acceptance;
all_passed = check('settings.acceptance carries the pre-registered thresholds', ...
    isstruct(acc) && strcmp(acc.fixed_on, '2026-09-14') && acc.J_rel_fro_max == 1e-6 && acc.L_abs_diff_max == 0.05) && all_passed;
V = D.verdict;
v_ok = isstruct(V) && isfield(V, 'per_condition') && isfield(V, 'all_pass') && numel(V.per_condition) == 3;
for i = 1:numel(V.per_condition)
    for c = {'A', 'B', 'L', 'C', 'J'}
        r = V.per_condition(i).(c{1});
        v_ok = v_ok && all(isfield(r, {'value', 'threshold', 'pass', 'text'})) && (isnan(r.pass) || r.pass == 0 || r.pass == 1);
    end
    fprintf('  %-14s A %s  B %s  L %s  C %s  J %s\n', V.per_condition(i).name, ptxt(V.per_condition(i).A.pass), ...
        ptxt(V.per_condition(i).B.pass), ptxt(V.per_condition(i).L.pass), ptxt(V.per_condition(i).C.pass), ptxt(V.per_condition(i).J.pass));
end
all_passed = check('verdict has the five checks per condition with value/threshold/pass/text', v_ok) && all_passed;
all_passed = check('the Jacobian criterion passes in every condition', ...
    all(arrayfun(@(c) c.J.pass == 1, V.per_condition))) && all_passed;
md_file = fullfile(run_dir, 'numerics_verification', 'numerics_verdict.md');
md = '';
if isfile(md_file); md = fileread(md_file); end
all_passed = check('numerics_verdict.md written and names every check and claim', ...
    ~isempty(md) && contains(md, '| J |') && contains(md, '| A |') && contains(md, 'Manuscript claims') && ...
    contains(md, 'finite differences')) && all_passed;

%% Figures
fig_ok = true; n_files = zeros(1, 4);
variants = {'solver', 'lya_method', 'ensemble', 'jacobian'};
for k = 1:numel(variants)
    try
        evalc(['res_k = fig_numerics_verification(''run_dir'', run_dir, ''variant'', variants{k}, ' ...
            '''out_dir'', fullfile(run_dir, [''fig_'' variants{k}]), ''visible'', false);']);
        n_files(k) = numel(res_k.files);
        close(res_k.figs);
    catch ME
        fprintf('  variant %s FAILED: %s\n', variants{k}, ME.message);
        fig_ok = false;
    end
end
all_passed = check(sprintf('all four figure variants write files (%s)', mat2str(n_files)), fig_ok && all(n_files >= 1)) && all_passed;
% a .mat without the jacobian field must be refused by the jacobian variant
D2 = D; D2.results = rmfield(D2.results, 'jacobian');
old_dir = fullfile(run_dir, 'old', 'numerics_verification'); mkdir(old_dir);
save(fullfile(old_dir, 'numerics_verification_data.mat'), '-struct', 'D2', '-v7.3');
all_passed = check('jacobian variant errors NoJacobian on a run without the field', ...
    throws_id(@() fig_numerics_verification('run_dir', fullfile(run_dir, 'old'), 'variant', 'jacobian', 'save', false, 'visible', false), ...
    'fig_numerics_verification:NoJacobian')) && all_passed;
rmdir(run_dir, 's');

fprintf('\n');
if all_passed
    fprintf('=== ALL numerics_verification stage TESTS PASSED ===\n');
else
    fprintf('=== SOME numerics_verification stage TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function s = ptxt(p)
if isnan(p); s = 'n/a '; elseif p; s = 'PASS'; else; s = 'FAIL'; end
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

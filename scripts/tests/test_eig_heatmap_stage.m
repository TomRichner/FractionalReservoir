% test_eig_heatmap_stage.m - run_eig_heatmap with the E:I imbalance examples,
% fig_eig_heatmap and fig_eig_heatmap_imbalance, end to end on a small network.
%
% Checks:
%   1. A bad run mode errors :badMode naming every mode (also in test_run_modes).
%   2. 'fast' with n_override 60 (so the examples run at 60 too), 5 states
%      per example, into a temp run directory: the .mat keeps every legacy
%      variable (evals_by_cond, lle_by_cond, lya_by_cond, num/spec abscissa,
%      J_times_by_cond, lle_window, lya_T_interval, settings, cond_names,
%      condition_titles) AND carries `examples`: 3 entries (inhibition-
%      dominant, reference, excitation-dominant) x 3 conditions with finite
%      lambda_1, mean rate in [0, 1] and B_E in [0, 1]; the reference entry's
%      eigenvalues equal the legacy variable; the two imbalanced examples
%      carry mu_EE_relative at 0.5x and 1.5x the preset's value and B_E
%      ordered inhibition < reference < excitation.
%   3. fig_eig_heatmap and fig_eig_heatmap_imbalance write files and the
%      latter its table; the imbalance figure errors :NoExamples on a .mat
%      without examples, and both error via resolve_data_file on a run
%      directory without the stage.
%
% ~2 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: run_eig_heatmap, fig_eig_heatmap, fig_eig_heatmap_imbalance

fprintf('=== Testing run_eig_heatmap + figures (imbalance examples) ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';

all_passed = check('bad run mode errors :badMode', ...
    throws_id(@() run_eig_heatmap('run_mode', 'nope'), 'run_eig_heatmap:badMode')) && all_passed;

run_dir = fullfile(tempdir, 'eig_heatmap_stage_test');
if exist(run_dir, 'dir'); rmdir(run_dir, 's'); end
mkdir(run_dir);
t0 = tic;
out = evalc(['mat_file = run_eig_heatmap(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(run_dir, ''eig_heatmap''), ''n_override'', 60, ''n_samples'', 5, ' ...
    '''n_samples_examples'', 5, ''use_parallel'', false);']);
fprintf('  (stage ran in %.0f s)\n', toc(t0));
all_passed = check('stage returns the .mat path and it exists', ischar(mat_file) && isfile(mat_file)) && all_passed;
D = load(mat_file);
legacy = {'evals_by_cond', 'lle_by_cond', 'lya_by_cond', 'num_abscissa_by_cond', ...
    'spec_abscissa_by_cond', 'J_times_by_cond', 'lle_window', 'lya_T_interval', ...
    'settings', 'cond_names', 'condition_titles'};
all_passed = check('legacy variables all present', all(isfield(D, legacy))) && all_passed;
all_passed = check('examples present: 3 entries, labels in order', isfield(D, 'examples') && ...
    numel(D.examples) == 3 && isequal({D.examples.label}, {'inhibition-dominant', 'reference', 'excitation-dominant'})) && all_passed;
E = D.examples;
n_cond = numel(D.cond_names);
fin = true; rng_ok = true; sz = true;
for e = 1:numel(E)
    fin = fin && all(isfinite(E(e).lle_by_cond)) && all(isfinite(E(e).mean_rate_by_cond)) && all(isfinite(E(e).B_E_by_cond));
    rng_ok = rng_ok && all(E(e).mean_rate_by_cond >= 0 & E(e).mean_rate_by_cond <= 1) && ...
        all(E(e).B_E_by_cond >= 0 & E(e).B_E_by_cond <= 1);
    sz = sz && numel(E(e).evals_by_cond) == n_cond && numel(E(e).lle_by_cond) == n_cond && ...
        all(cellfun(@(v) ~isempty(v), E(e).evals_by_cond)) && E(e).n == 60;
    fprintf('  %-20s mu_EE %6.3f  lambda_1 %s  <r> %s  B_E %s\n', E(e).label, E(e).mu_tilde_relative(1, 1), ...
        mat2str(E(e).lle_by_cond, 3), mat2str(E(e).mean_rate_by_cond, 3), mat2str(E(e).B_E_by_cond, 3));
end
all_passed = check('every example x condition: finite lambda_1, mean rate, B_E', fin) && all_passed;
all_passed = check('mean rate and B_E in [0, 1]; n = 60; eigenvalues pooled', rng_ok && sz) && all_passed;
i_ref = find(strcmp({E.label}, 'reference'));
all_passed = check('reference example equals the legacy variables', ...
    isequal(E(i_ref).evals_by_cond, D.evals_by_cond) && isequal(E(i_ref).lle_by_cond, D.lle_by_cond)) && all_passed;
d = srnn_param_preset(P);
mu_ref = d.mu_tilde_relative(1, 1);
i_inh = find(strcmp({E.label}, 'inhibition-dominant')); i_exc = find(strcmp({E.label}, 'excitation-dominant'));
all_passed = check('imbalanced examples carry mu_EE at 0.5x and 1.5x the preset', ...
    abs(E(i_inh).overrides.mu_EE_relative - 0.5 * mu_ref) < 1e-12 && abs(E(i_exc).overrides.mu_EE_relative - 1.5 * mu_ref) < 1e-12 && ...
    abs(E(i_inh).mu_tilde_relative(1, 1) - 0.5 * mu_ref) < 1e-9 && abs(E(i_exc).mu_tilde_relative(1, 1) - 1.5 * mu_ref) < 1e-9 && ...
    abs(E(i_ref).mu_tilde_relative(1, 1) - mu_ref) < 1e-9) && all_passed;
all_passed = check('B_E ordered inhibition < reference < excitation (every condition)', ...
    all(E(i_inh).B_E_by_cond < E(i_ref).B_E_by_cond) && all(E(i_ref).B_E_by_cond < E(i_exc).B_E_by_cond)) && all_passed;
all_passed = check('settings record the example labels and the override', ...
    isequal(D.settings.examples, {E.label}) && D.settings.n_samples_examples == 5 && D.settings.n_examples_override == 60) && all_passed;

%% Figures
evalc('res1 = fig_eig_heatmap(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig1''), ''visible'', false);');
all_passed = check(sprintf('fig_eig_heatmap writes files (%d)', numel(res1.files)), numel(res1.files) >= 1) && all_passed;
close(res1.figs);
evalc('res2 = fig_eig_heatmap_imbalance(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig2''), ''visible'', false);');
all_passed = check(sprintf('fig_eig_heatmap_imbalance writes files (%d) and the table', numel(res2.files)), ...
    numel(res2.files) >= 1 && isfile(fullfile(run_dir, 'fig2', 'Fig_Eig_Heatmap_Imbalance_table.md'))) && all_passed;
close(res2.figs);
old_dir = fullfile(tempdir, 'eig_heatmap_stage_test_old'); if ~exist(old_dir, 'dir'); mkdir(fullfile(old_dir, 'eig_heatmap')); end
Dold = rmfield(D, 'examples'); save(fullfile(old_dir, 'eig_heatmap', 'eig_heatmap_data.mat'), '-struct', 'Dold');
all_passed = check('imbalance figure errors :NoExamples on a pre-examples .mat', ...
    throws_id(@() fig_eig_heatmap_imbalance('run_dir', old_dir, 'save', false, 'visible', false), 'fig_eig_heatmap_imbalance:NoExamples')) && all_passed;
empty_dir = fullfile(tempdir, 'eig_heatmap_stage_test_empty'); if ~exist(empty_dir, 'dir'); mkdir(empty_dir); end
all_passed = check('both figures error on a run directory without the stage', ...
    throws_prefix(@() fig_eig_heatmap('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:') && ...
    throws_prefix(@() fig_eig_heatmap_imbalance('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:')) && all_passed;
rmdir(empty_dir, 's'); rmdir(old_dir, 's'); rmdir(run_dir, 's');

fprintf('\n');
if all_passed
    fprintf('=== ALL eig_heatmap stage TESTS PASSED ===\n');
else
    fprintf('=== SOME eig_heatmap stage TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
end

function ok = throws_prefix(fn, prefix)
ok = false;
try
    fn();
catch ME
    ok = startsWith(ME.identifier, prefix);
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

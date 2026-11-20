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
%   4. THE RITZ PASS (FD-033). Off by default: the .mat holds exactly the
%      variables it always did and `examples` has no ritz field. With
%      'ritz_K', 10 (reference example only, same seeds and states): every
%      legacy variable is identical to the run without it; ritz_by_cond
%      holds, per condition, H of size K x K x states x 2 scalings with
%      K = min(10, state dimension), every state valid and at the times of
%      J_times_by_cond, values = eig(H) with real parts summing to trace(H),
%      block fractions summing to 1, and the scales by kind of variable
%      (x alone without adaptation, where the scaled pass must then give the
%      raw Ritz values). A scaling name it does not know errors
%      :badRitzScaling.
%
% ~3 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
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
all_passed = check('no Ritz pass by default: the variables are the legacy ones and examples, no ritz field', ...
    isequal(sort(fieldnames(D)), sort([legacy, {'examples'}]')) && ~isfield(D.examples, 'ritz_by_cond') && ...
    ~isfield(D.settings, 'ritz_K')) && all_passed;
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

%% The Ritz pass (FD-033)
ritz_dir = fullfile(tempdir, 'eig_heatmap_stage_test_ritz');
if exist(ritz_dir, 'dir'); rmdir(ritz_dir, 's'); end
mkdir(ritz_dir);
ref_only = struct('label', 'reference', 'overrides', struct()); %#ok<NASGU>  used inside evalc
t0 = tic;
evalc(['mat_ritz = run_eig_heatmap(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(ritz_dir, ''eig_heatmap''), ''n_override'', 60, ''n_samples'', 5, ' ...
    '''examples'', ref_only, ''use_parallel'', false, ''ritz_K'', 10);']);
fprintf('  (stage with the Ritz pass, reference only, ran in %.0f s)\n', toc(t0));
DR = load(mat_ritz);
same = true;
for k = 1:numel(legacy)
    if ~ismember(legacy{k}, {'settings'}); same = same && isequaln(DR.(legacy{k}), D.(legacy{k})); end
end
all_passed = check('with the Ritz pass every legacy variable is identical to the run without it', same) && all_passed;
all_passed = check('settings record ritz_K and the scalings; the reference example carries the same structs', ...
    DR.settings.ritz_K == 10 && isequal(DR.settings.ritz_scalings, {'raw', 'type_sd'}) && ...
    isequaln(DR.examples(1).ritz_by_cond, DR.ritz_by_cond)) && all_passed;
shape = true; ident = true; kinds_ok = true; worst_trace = 0;
for i = 1:n_cond
    R = DR.ritz_by_cond{i};
    n_s = numel(DR.J_times_by_cond{i});
    shape = shape && R.K == min(10, R.N) && isequal(size(R.H), [R.K, R.K, n_s, 2]) && ...
        isequal(size(R.values), [R.K, n_s, 2]) && all(R.valid) && isequal(R.t(:), DR.J_times_by_cond{i}(:)) && ...
        all(isfinite(R.H), 'all') && all(isfinite(R.LE_spectrum), 'all');
    for sc = 1:2
        for k = 1:n_s
            worst_trace = max(worst_trace, abs(sum(real(R.values(:, k, sc))) - trace(R.H(:, :, k, sc))));
        end
    end
    ident = ident && max(abs(sum(R.block_fractions, 2) - 1), [], 'all') < 1e-12;
    if strcmp(DR.cond_names{i}, 'no_adaptation')
        d_raw = sort(real(R.values(:, :, 1))); d_sc = sort(real(R.values(:, :, 2)));
        kinds_ok = kinds_ok && isequal(R.scales.kinds, {'x'}) && max(abs(d_raw - d_sc), [], 'all') < 1e-8;
    else
        kinds_ok = kinds_ok && numel(R.scales.kinds) > 1 && all(R.scales.sd > 0);
    end
    fprintf('  %-14s N = %3d, K = %2d: lambda_1 %+.3f (raw) %+.3f (type_sd); x share of the subspace %.2f / %.2f; kinds %s\n', ...
        DR.cond_names{i}, R.N, R.K, R.LE_spectrum(1, 1), R.LE_spectrum(1, 2), ...
        mean(R.block_fractions(:, 1, 1)), mean(R.block_fractions(:, 1, 2)), strjoin(R.scales.kinds, ' '));
end
all_passed = check('per condition: H is K x K x states x 2, all states valid, at the sampled times', shape) && all_passed;
all_passed = check(sprintf('real parts of the Ritz values sum to trace(H) (worst %.1e)', worst_trace), worst_trace < 1e-9) && all_passed;
all_passed = check('block fractions sum to 1', ident) && all_passed;
all_passed = check('scales by kind: x alone without adaptation (scaled = raw there), several kinds with it', kinds_ok) && all_passed;
all_passed = check('an unknown scaling errors :badRitzScaling', ...
    throws_id(@() run_eig_heatmap('preset_name', P, 'run_mode', 'fast', 'ritz_K', 5, 'ritz_scalings', {'zscore'}), ...
    'run_eig_heatmap:badRitzScaling')) && all_passed;

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
rmdir(empty_dir, 's'); rmdir(old_dir, 's'); rmdir(run_dir, 's'); rmdir(ritz_dir, 's');

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

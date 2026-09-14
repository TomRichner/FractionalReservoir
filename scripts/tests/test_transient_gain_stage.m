% test_transient_gain_stage.m - run_transient_gain, fig_transient_gain and
% fig_transient_gain_excursions end to end on a small network.
%
% Checks:
%   1. A bad run mode errors :badMode naming every mode (also in test_run_modes).
%   2. 'fast' with n_override 60, TWO seeds (the assembly across seeds once
%      failed on the real network with one-seed tests green), 3 regular + up to 2 onset + 2
%      quiet samples, horizon 0.5 s, into a temp run directory: the .mat has
%      results / cond_names / condition_titles / settings; every condition
%      has trials, samples, N, n, n_onset_found, n_quiet_found; every sample
%      carries the schema (t, the five 3 x n_t readings, G_max, t_peak,
%      alignment, E fraction, participation, both v_opt, block fractions);
%      G(0) = 1 for every variant; worst >= every other reading at every t;
%      v_opt unit; the kind counts match the used counts in the trial record
%      and never exceed the found counts; regular samples lie in the sample
%      window; t_peak within the horizon.
%   3. Both figures write files and their tables; both error clearly on a
%      run directory without the stage.
%
% ~1 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: run_transient_gain, fig_transient_gain, fig_transient_gain_excursions,
%           test_transient_gain

fprintf('=== Testing run_transient_gain + figures ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';

all_passed = check('bad run mode errors :badMode', ...
    throws_id(@() run_transient_gain('run_mode', 'nope'), 'run_transient_gain:badMode')) && all_passed;

run_dir = fullfile(tempdir, 'transient_gain_stage_test');
if exist(run_dir, 'dir'); rmdir(run_dir, 's'); end
mkdir(run_dir);
t0 = tic;
out = evalc(['mat_file = run_transient_gain(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(run_dir, ''transient_gain''), ''n_override'', 60, ''n_seeds'', 2,' ...
    '''n_regular'', 3, ''n_excursion'', 2, ''horizon_s'', 0.5);']);
fprintf('  (stage ran in %.0f s)\n', toc(t0));
all_passed = check('stage returns the .mat path and it exists', ischar(mat_file) && isfile(mat_file)) && all_passed;
D = load(mat_file);
all_passed = check('schema: results / cond_names / condition_titles / settings', ...
    all(isfield(D, {'results', 'cond_names', 'condition_titles', 'settings'}))) && all_passed;
R = D.results;
all_passed = check('three conditions with trials, samples, counts, N, n', numel(R) == 3 && ...
    all(isfield(R, {'trials', 'samples', 'n_onset_found', 'n_quiet_found', 'N', 'n'}))) && all_passed;
need = {'t_sample', 'kind', 'seeds', 'local_rate_at_sample', 't', 'G_worst', 'G_noise', 'G_ei_diff', ...
    'G_ei_sum', 'G_lyap', 'G_max', 't_peak', 'align_opt_lyap', 'cos_opt_ei_diff', 'cos_opt_ei_sum', ...
    'frac_E_opt', 'participation_opt', 'v_opt_active', 'v_opt_frozen_x', 'lyap_block_fracs', 'seconds', ...
    'alpha_xx', 'omega_xx', 'alpha_full'};
schema_ok = true; g0 = true; ord = true; unit = true; counts_ok = true; win_ok = true; peak_ok = true;
win = D.settings.sample_window;
for i = 1:numel(R)
    smp = R(i).samples;
    schema_ok = schema_ok && all(isfield(smp, need));
    for s = 1:numel(smp)
        x = smp(s);
        nt = numel(x.t);
        schema_ok = schema_ok && isequal(size(x.G_worst), [3 nt]) && isequal(size(x.G_lyap), [3 nt]) && ...
            numel(x.G_max) == 3 && numel(x.v_opt_active) == R(i).n && numel(x.v_opt_frozen_x) == R(i).n && ...
            abs(sum(x.lyap_block_fracs) - 1) < 1e-9 && ismember(x.kind, {'regular', 'onset', 'quiet'});
        g0 = g0 && all(x.G_worst(:, 1) == 1) && all(x.G_noise(:, 1) == 1);
        ord = ord && all(x.G_worst(:) >= x.G_noise(:) - 1e-9) && all(x.G_worst(:) >= x.G_ei_diff(:) - 1e-9) && ...
            all(x.G_worst(:) >= x.G_ei_sum(:) - 1e-9) && all(x.G_worst(:) >= x.G_lyap(:) - 1e-9);
        schema_ok = schema_ok && isfinite(x.alpha_xx) && x.omega_xx >= x.alpha_xx - 1e-9 && (isnan(x.alpha_full) || isfinite(x.alpha_full));
        unit = unit && abs(norm(x.v_opt_active) - 1) < 1e-9 && abs(norm(x.v_opt_frozen_x) - 1) < 1e-9 && ...
            all(x.align_opt_lyap >= 0 & x.align_opt_lyap <= 1 + 1e-12);
        peak_ok = peak_ok && all(x.t_peak >= 0 & x.t_peak <= D.settings.horizon_s + 1e-9) && ...
            all(abs(x.G_max - max(x.G_worst, [], 2)') < 1e-9);
        if strcmp(x.kind, 'regular')
            win_ok = win_ok && x.t_sample >= win(1) - 1e-9 && x.t_sample <= win(2) + 1e-9;
        end
    end
    n_on = nnz(strcmp({smp.kind}, 'onset')); n_qu = nnz(strcmp({smp.kind}, 'quiet'));
    counts_ok = counts_ok && n_on == sum([R(i).trials.n_onset_used]) && n_qu == sum([R(i).trials.n_quiet_used]) && ...
        n_on <= R(i).n_onset_found && n_qu <= R(i).n_quiet_found && n_on <= 4 && n_qu <= 4 && ...
        nnz(strcmp({smp.kind}, 'regular')) == 6 && numel(R(i).trials) == 2 && isequal([R(i).trials.seeds], [1 2 2 3]);
    fprintf('  %-14s N %d, lambda_1 %+.3f, %d samples (%d onset, %d quiet; found %d / %d), G_max regular medians: %s\n', ...
        R(i).name, R(i).N, R(i).trials(1).LLE, numel(smp), n_on, n_qu, R(i).n_onset_found, R(i).n_quiet_found, ...
        mat2str(median(vertcat(smp(strcmp({smp.kind}, 'regular')).G_max), 1), 3));
end
all_passed = check('every sample carries the schema with 3 x n_t readings, n-vectors, fractions summing to 1', schema_ok) && all_passed;
all_passed = check('G(0) = 1 for every variant and sample', g0) && all_passed;
all_passed = check('worst >= noise, E/I diff, E/I sum, Lyapunov readings at every t', ord) && all_passed;
all_passed = check('v_opt unit; alignment in [0, 1]', unit) && all_passed;
all_passed = check('G_max = max of the curve, t_peak within the horizon', peak_ok) && all_passed;
all_passed = check('kind counts match the trial record and the caps; 3 regular per seed, 2 seeds', counts_ok) && all_passed;
all_passed = check('regular samples inside the sample window', win_ok) && all_passed;
all_passed = check('settings record T, horizon, counts, variants, directions, minutes', ...
    D.settings.T == 20 && D.settings.horizon_s == 0.5 && D.settings.n_regular == 3 && D.settings.n_excursion == 2 && ...
    numel(D.settings.variants) == 3 && numel(D.settings.directions) == 3 && isfinite(D.settings.minutes)) && all_passed;

%% Figures
evalc('res1 = fig_transient_gain(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig1''), ''visible'', false);');
all_passed = check(sprintf('fig_transient_gain writes files (%d) and the table', numel(res1.files)), ...
    numel(res1.files) >= 1 && isfile(fullfile(run_dir, 'fig1', 'Fig_Transient_Gain_table.md'))) && all_passed;
close(res1.figs);
evalc('res2 = fig_transient_gain_excursions(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig2''), ''visible'', false);');
all_passed = check(sprintf('fig_transient_gain_excursions writes files (%d) and the table', numel(res2.files)), ...
    numel(res2.files) >= 1 && isfile(fullfile(run_dir, 'fig2', 'Fig_Transient_Gain_Excursions_table.md'))) && all_passed;
close(res2.figs);
empty_dir = fullfile(tempdir, 'transient_gain_stage_test_empty'); if ~exist(empty_dir, 'dir'); mkdir(empty_dir); end
all_passed = check('both figures error on a run directory without the stage', ...
    throws_prefix(@() fig_transient_gain('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:') && ...
    throws_prefix(@() fig_transient_gain_excursions('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:')) && all_passed;
rmdir(empty_dir, 's'); rmdir(run_dir, 's');

fprintf('\n');
if all_passed
    fprintf('=== ALL transient_gain stage TESTS PASSED ===\n');
else
    fprintf('=== SOME transient_gain stage TESTS FAILED ===\n');
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

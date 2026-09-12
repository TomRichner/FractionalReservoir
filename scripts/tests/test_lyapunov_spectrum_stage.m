% test_lyapunov_spectrum_stage.m - run_lyapunov_spectrum and fig_lyapunov_spectrum
% end to end on a small network.
%
% Checks:
%   1. A bad run mode errors :badMode naming every mode (also in test_run_modes).
%   2. 'fast' with n_override 60, one seed, K 20, into a temp run directory:
%      the .mat has results / cond_names / condition_titles / settings; every
%      condition has noise_on AND noise_off runs; each run carries the
%      spectrum (K_used long, descending), n_positive, h_KS_bits, D_KY +
%      resolved, K_used in [20, 40], the finite-time curves, the block
%      fractions summing to 1, N and the seeds; noise off used rk4 and
%      sigma 0, noise on sra1 and the preset's sigma.
%   3. The stable regime's leading vector is in the SFA block; the
%      no-adaptation regime's is in x.
%   4. fig_lyapunov_spectrum on that run directory writes files and the
%      table; it errors clearly on a run directory without the stage.
%
% ~2 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: run_lyapunov_spectrum, fig_lyapunov_spectrum

fprintf('=== Testing run_lyapunov_spectrum + fig_lyapunov_spectrum ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';

all_passed = check('bad run mode errors :badMode', ...
    throws_id(@() run_lyapunov_spectrum('run_mode', 'nope'), 'run_lyapunov_spectrum:badMode')) && all_passed;

run_dir = fullfile(tempdir, 'lyapunov_spectrum_stage_test');
if exist(run_dir, 'dir'); rmdir(run_dir, 's'); end
mkdir(run_dir);
t0 = tic;
out = evalc(['mat_file = run_lyapunov_spectrum(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(run_dir, ''lyapunov_spectrum''), ''n_override'', 60, ''n_seeds'', 1, ''K'', 20);']);
fprintf('  (stage ran in %.0f s)\n', toc(t0));
all_passed = check('stage returns the .mat path and it exists', ischar(mat_file) && isfile(mat_file)) && all_passed;
D = load(mat_file);
all_passed = check('schema: results / cond_names / condition_titles / settings', ...
    all(isfield(D, {'results', 'cond_names', 'condition_titles', 'settings'}))) && all_passed;
R = D.results;
all_passed = check('three conditions, each with noise_on and noise_off', numel(R) == 3 && ...
    all(arrayfun(@(r) isfield(r, 'noise_on') && isfield(r, 'noise_off') && ~isempty(r.noise_on) && ~isempty(r.noise_off), R))) && all_passed;
ok = true; kok = true; frac_ok = true; solver_ok = true;
for i = 1:numel(R)
    for var = {'noise_on', 'noise_off'}
        r = R(i).(var{1});
        lam = r.LE_spectrum;
        ok = ok && numel(lam) == r.K_used && issorted(-lam) && isfinite(r.h_KS_bits) && ...
            isfinite(r.n_positive) && ismember(r.D_KY_resolved, [0 1]) && ...
            size(r.finite_LE_spectrum_t, 2) == min(10, r.K_used) && r.N > 0 && isequal(r.seeds, [1 2]);
        kok = kok && r.K_used >= 20 && r.K_used <= 40;
        frac_ok = frac_ok && abs(r.lead_frac_x + r.lead_frac_sfa + r.lead_frac_std + r.lead_frac_stf - 1) < 1e-9;
        if strcmp(var{1}, 'noise_off')
            solver_ok = solver_ok && strcmp(r.ode_solver, 'rk4') && r.sigma_u_noise == 0 && ~r.noise_on;
        else
            solver_ok = solver_ok && strcmp(r.ode_solver, 'sra1') && r.sigma_u_noise > 0 && r.noise_on;
        end
    end
    fprintf('  %-14s noise on : N %d, K_used %d, lambda_1 %+.3f, n_pos %d, h_KS %.2f, D_KY %.2f (resolved %d), sfa frac %.2f\n', ...
        R(i).name, R(i).noise_on.N, R(i).noise_on.K_used, R(i).noise_on.LLE, R(i).noise_on.n_positive, ...
        R(i).noise_on.h_KS_bits, R(i).noise_on.D_KY, R(i).noise_on.D_KY_resolved, R(i).noise_on.lead_frac_sfa);
end
all_passed = check('every run: spectrum of length K_used, descending; scalars finite; curves and seeds', ok) && all_passed;
all_passed = check('K_used within [K, 2K]', kok) && all_passed;
all_passed = check('block fractions sum to 1', frac_ok) && all_passed;
all_passed = check('noise off ran rk4 at sigma 0; noise on ran sra1 at the preset sigma', solver_ok) && all_passed;
all_passed = check('settings record T, K, K_max, n_seeds, variants, minutes', ...
    D.settings.T == 20 && D.settings.K == 20 && D.settings.K_max == 40 && D.settings.n_seeds == 1 && ...
    numel(D.settings.variants) == 2 && isfinite(D.settings.minutes)) && all_passed;
i_no = find(strcmp({R.name}, 'no_adaptation')); i_st = find(strcmp({R.name}, 'sfa3_std2'));
all_passed = check('no-adaptation leading vector in x; stable regime leading vector mostly SFA', ...
    R(i_no).noise_off.lead_frac_x == 1 && R(i_st).noise_off.lead_frac_sfa > 0.5) && all_passed;

%% Figure
fo = evalc('res = fig_lyapunov_spectrum(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig''), ''visible'', false);');
all_passed = check(sprintf('fig_lyapunov_spectrum writes files (%d)', numel(res.files)), numel(res.files) >= 1) && all_passed;
all_passed = check('...and the table markdown', isfile(fullfile(run_dir, 'fig', 'Fig_Lyapunov_Spectrum_table.md'))) && all_passed;
close(res.figs);
empty_dir = fullfile(tempdir, 'lyapunov_spectrum_stage_test_empty'); if ~exist(empty_dir, 'dir'); mkdir(empty_dir); end
all_passed = check('figure errors on a run directory without the stage', ...
    throws_prefix(@() fig_lyapunov_spectrum('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:')) && all_passed;
rmdir(empty_dir, 's'); rmdir(run_dir, 's');

fprintf('\n');
if all_passed
    fprintf('=== ALL lyapunov_spectrum stage TESTS PASSED ===\n');
else
    fprintf('=== SOME lyapunov_spectrum stage TESTS FAILED ===\n');
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

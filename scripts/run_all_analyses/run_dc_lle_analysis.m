% run_dc_lle_analysis.m
% Multi-seed DC-staircase Lyapunov analysis.
%
% Question: does tonic DC drive change the largest local Lyapunov exponent
% (LLE) of the bursting SRNN, and is the effect robust across networks?
%
% For each of n_seeds network/stimulus seeds we run ONE simulation in which a
% uniform DC is stepped through a staircase of levels (the same setup as
% scripts/presentations/Stability_Manuscript/fig_stim_engages_adaptation/
% bursting_SRNN_model_good_ex.m, single SFA+STD regime),
% compute the local LLE timeseries with Benettin's method, and bin it by DC
% level -- skipping the first psd_settle seconds after each step so SFA/STD have
% settled. Aggregating the per-level mean LLE ACROSS seeds gives a mean +/- std
% (seed variability) as a function of DC, plotted with confplot.
%
% Seed is the parallel axis (parfor); the DC levels live inside one staircase
% simulation, so DC is NOT a grid axis (it is a sub-field of input_config, which
% the grid-sweep framework cannot vary). This is why the analysis is a
% purpose-built parfor script rather than a ParamSpaceAnalysis2 run.
%
% Honors the run_all_analyses orchestrator contract: no clear/clc/close all, and
% it reads master_output_dir / master_save_figs / run_mode from the base
% workspace when present.
%
% See also: bursting_SRNN_model, dc_staircase_stimulus, confplot, SRNNModel2

% No clear/clc/close all on standalone runs, so base-workspace settings such as
% run_mode and save_figs (set in the console before running) survive into this
% script. run_all_analyses relies on the same: it never clears the sub-scripts.

%% Figure saving configuration
if exist('master_save_figs', 'var')
    if strcmp(master_save_figs, 'save_all_figs')
        save_figs = true;
    elseif strcmp(master_save_figs, 'save_no_figs')
        save_figs = false;
    end
end
if ~exist('save_figs', 'var')
    save_figs = false;
end

%% STD zero-floor configuration
% Matches the SRNNModel2 default (off) unless set in the base workspace.
if ~exist('std_zero_floor', 'var')
    std_zero_floor = false;
end

%% Setup paths
setup_paths();

%% Run mode: 'fast'/'medium'/'production'. Here it sets n_seeds and the solver
% (the staircase T_range is fixed by dc_levels x hold_dur, so the fs/T_range/LLE
% tuning used by the PSA sub-analyses does not apply). Set run_mode in the base
% workspace (or via run_all_analyses); defaults to 'production' when standalone.
if ~exist('run_mode', 'var'); run_mode = 'production'; end
switch run_mode
    case 'fast'
        n_seeds    = 25;
        dc_levels  = [0 0.0125 0.025 0.05 0.075 0.1 0.15 0.2 0.3];
        ode_solver = @ode_rk4;
        % ode_solver = @ode45;
    case 'medium'
        % Roughly halfway between fast and production: more seeds, still the
        % fast fixed-step solver to keep the (expensive) staircase runs quick.
        n_seeds    = 38;
        dc_levels  = [0 0.0125 0.025 0.05 0.075 0.1 0.15 0.2 0.3];
        ode_solver = @ode_rk4;
    case 'production'
        n_seeds    = 50;
        dc_levels  = [0 0.0125 0.025 0.05 0.075 0.1 0.15 0.2 0.3];
        ode_solver = @ode45;
    otherwise
        error('run_dc_lle_analysis:badMode', ...
            'Unknown run_mode ''%s'' (expected ''fast'', ''medium'', or ''production'').', run_mode);
end
nL = numel(dc_levels);
use_parallel = true;   % set false for serial debugging
fprintf('[run_dc_lle_analysis] run_mode=%s, n_seeds=%d, nL=%d, ode_solver=%s\n', ...
    run_mode, n_seeds, nL, func2str(ode_solver));

%% Network / model configuration (single SFA+STD "bursting" regime)
% Mirrors scripts/presentations/Stability_Manuscript/
% fig_stim_engages_adaptation/bursting_SRNN_model_good_ex.m.
n        = 50;      % total neurons
f        = 0.5;     % fraction excitatory
indegree = 4;       % expected in-degree

% Nonlinearity: piecewise sigmoid (both fn + derivative; derivative unused by
% Benettin but kept for parity with the bursting script / QR compatibility).
S_a = 0.9;
S_c = 0.5;
phi       = @(x) SRNNModel2.piecewiseSigmoid(x, S_a, S_c);
phi_deriv = @(x) SRNNModel2.piecewiseSigmoidDerivative(x, S_a, S_c);

% Stimulus timing
hold_dur        = 30;      % seconds each DC level is held
ramp_dur        = 10;      % ramp 0 -> dc_levels(1) over the first ramp_dur s
noise_intensity = 0.001;   % fs-invariant white-noise intensity
psd_settle      = 15;      % seconds skipped after each step before the LLE window
fs              = 400;
T_range         = [0, nL * hold_dur];

% Custom DC staircase generator (shared, on the path via setup_paths)
input_config = struct();
input_config.dc_levels       = dc_levels;
input_config.hold_dur        = hold_dur;
input_config.ramp_dur        = ramp_dur;
input_config.noise_intensity = noise_intensity;
input_config.intrinsic_drive = [];                      % required by the class
input_config.generator       = @dc_staircase_stimulus;

%% Output directory
project_root = fileparts(which('setup_paths'));
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
if exist('master_output_dir', 'var')
    base_dir = master_output_dir;
else
    base_dir = fullfile(project_root, 'data', 'param_space');
end
output_dir = fullfile(base_dir, sprintf('dc_lle_nSeeds_%d_%s', n_seeds, dt_str));
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
fprintf('[run_dc_lle_analysis] output_dir = %s\n', output_dir);

%% Run the seed sweep (parfor over seeds)
seed_results = cell(n_seeds, 1);
% Deterministic network seeds (offset 0) so every run draws the same networks.
% `seeds` records the actual RNG seeds used (offset + 1:n_seeds).
seed_offset = 0;
seeds = seed_offset + (1:n_seeds);

run_start = tic;
if use_parallel
    parfor s = 1:n_seeds
        seed_results{s} = run_one_seed(s, seed_offset, n, f, indegree, phi, phi_deriv, ...
            S_a, S_c, input_config, ode_solver, fs, T_range, ...
            hold_dur, psd_settle, nL, std_zero_floor);
    end
else
    for s = 1:n_seeds
        seed_results{s} = run_one_seed(s, seed_offset, n, f, indegree, phi, phi_deriv, ...
            S_a, S_c, input_config, ode_solver, fs, T_range, ...
            hold_dur, psd_settle, nL, std_zero_floor);
    end
end
fprintf('[run_dc_lle_analysis] %d seeds finished in %.1f s\n', n_seeds, toc(run_start));

%% Aggregate across seeds
% Per-seed per-level mean/std (n_seeds x nL)
per_seed_level_mean = cell2mat(cellfun(@(r) r.lvl_mean, seed_results, 'UniformOutput', false));
per_seed_level_std  = cell2mat(cellfun(@(r) r.lvl_std,  seed_results, 'UniformOutput', false));

% Across-seed mean / std of the per-level mean LLE (column vectors, length nL)
level_mean = mean(per_seed_level_mean, 1)';
level_std  = std(per_seed_level_mean, 0, 1)';

% Full local-LLE timeseries stacked across seeds (same t grid for all seeds).
t_lya = seed_results{1}.t_lya;
lens  = cellfun(@(r) numel(r.local_lya), seed_results);
if all(lens == lens(1))
    local_lya_all = cell2mat(cellfun(@(r) r.local_lya(:), seed_results', 'UniformOutput', false)); % nt x n_seeds
else
    warning('run_dc_lle_analysis:RaggedLLE', ...
        'local_lya lengths differ across seeds; storing as cell.');
    local_lya_all = cellfun(@(r) r.local_lya(:), seed_results, 'UniformOutput', false);
end

%% Save essential data
config = struct('n', n, 'f', f, 'indegree', indegree, 'S_a', S_a, 'S_c', S_c, ...
    'hold_dur', hold_dur, 'ramp_dur', ramp_dur, 'noise_intensity', noise_intensity, ...
    'psd_settle', psd_settle, 'fs', fs, 'T_range', T_range, ...
    'ode_solver', func2str(ode_solver), 'run_mode', run_mode);

dc_lle_results = struct( ...
    'dc_levels', dc_levels, 'hold_dur', hold_dur, 'psd_settle', psd_settle, ...
    'n_seeds', n_seeds, 'seeds', seeds, 'seed_offset', seed_offset, 't_lya', t_lya, ...
    'local_lya_all', local_lya_all, ...
    'per_seed_level_mean', per_seed_level_mean, 'per_seed_level_std', per_seed_level_std, ...
    'level_mean', level_mean, 'level_std', level_std, 'config', config);

save(fullfile(output_dir, 'dc_lle_results.mat'), 'dc_lle_results', '-v7.3');
fprintf('[run_dc_lle_analysis] saved dc_lle_results.mat\n');

% Copy this script to the output directory for reproducibility
copyfile([mfilename('fullpath') '.m'], output_dir);

%% Plot: mean local LLE vs DC level, +/- std across seeds
figure('Name', 'DC LLE: local Lyapunov vs DC level (across seeds)');
confplot(dc_levels, level_mean, level_std, level_std, [0 0 0.8; 0.7 0.8 1.0]);
yline(0, '--k', 'edge of chaos', 'HandleVisibility', 'off');
xlabel('DC level (input units)');
ylabel('mean local \lambda_1  (\pm std across seeds)');
title(sprintf('Largest local Lyapunov exponent vs DC stim level  (n_{seeds}=%d)', n_seeds));
grid off;

if save_figs
    fig_dir = fullfile(output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'dc_lle', [], {'fig', 'png'});
    fprintf('[run_dc_lle_analysis] figures saved to %s\n', fig_dir);
end

%% Summary
fprintf('\n=== DC LLE Analysis Summary ===\n');
fprintf('Output directory: %s\n', output_dir);
fprintf('DC levels: %s\n', mat2str(dc_levels));
for k = 1:nL
    fprintf('  DC = %6.4g :  lambda_1 = %+.4f  +/- %.4f (across %d seeds)\n', ...
        dc_levels(k), level_mean(k), level_std(k), n_seeds);
end
fprintf('Done! Results saved to: %s\n', output_dir);

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================
function res = run_one_seed(s, seed_offset, n, f, indegree, phi, phi_deriv, S_a, S_c, ...
        input_config, ode_solver, fs, T_range, hold_dur, psd_settle, nL, std_zero_floor)
    % Run one staircase simulation for seed s and return per-level local LLE.
    % seed_offset shifts the RNG seed so separate runs use distinct networks.
    rng_seed = s + seed_offset;
    model = SRNNModel2( ...
        'n', n, 'f', f, 'indegree', indegree, ...
        'mu_E_tilde_relative', 3, 'mu_I_tilde_relative', -4, ...
        'sigma_E_tilde_relative', 1, 'sigma_I_tilde_relative', 1, ...
        'level_of_chaos', 1.0, ...
        'n_a_E', 3, 'n_a_I', 0, 'c_E', 0.15/3, 'c_I', 0.15/3, ...
        'tau_a_E', logspace(log10(0.25), log10(10), 3), ...
        'n_b_E', 1, 'n_b_I', 0, 'tau_b_E_rec', 1, 'tau_b_E_rel', 0.25, ...
        'std_zero_floor', std_zero_floor, ...
        'tau_d', 0.1, 'S_a', S_a, 'S_c', S_c, ...
        'activation_function', phi, 'activation_function_derivative', phi_deriv, ...
        'input_config', input_config, 'u_ex_scale', 1.0, ...
        'fs', fs, 'T_range', T_range, 'ode_solver', ode_solver, ...
        'rng_seeds', [rng_seed rng_seed], 'lya_method', 'benettin', ...
        'store_full_state', false);
    model.build();
    model.run();

    ll = model.lya_results.local_lya;
    t  = model.lya_results.t_lya;

    lvl_mean = nan(1, nL);
    lvl_std  = nan(1, nL);
    for k = 1:nL
        lo  = (k-1)*hold_dur + psd_settle;
        hi  = k*hold_dur;
        sel = t > lo & t <= hi;
        lvl_mean(k) = mean(ll(sel));
        lvl_std(k)  = std(ll(sel));
    end

    res = struct('seed', s, 'local_lya', ll, 't_lya', t, ...
        'lvl_mean', lvl_mean, 'lvl_std', lvl_std, 'LLE', model.lya_results.LLE);
end

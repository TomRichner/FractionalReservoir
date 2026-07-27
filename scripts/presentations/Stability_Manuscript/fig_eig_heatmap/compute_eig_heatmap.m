% compute_eig_heatmap.m - Simulate + compute Jacobian-eigenvalue occupancy data
%
% Computing half of the eig-heatmap figure. Runs the four adaptation regimes
% (no adaptation / SFA / STD / SFA+STD) on a shared network, samples the
% instantaneous Jacobian at fixed intervals through each run, and pools the
% eigenvalues. The (expensive) eigenvalue computation is parallelized over the
% sampled time points via the Parallel Computing Toolbox. Results are saved to
% eig_heatmap_data.mat in this folder; plotting is done separately by
% fig_eig_heatmap.m, so the figure look can be iterated without re-simulating.
%
% Because the SRNN is nonlinear, the eigenvalues of the instantaneous Jacobian
% move around the complex plane as the state evolves; pooling them over many
% sampled times reveals how much time they spend in each region, and in
% particular to the right of the imaginary axis (Re > 0, locally unstable).
%
% All four conditions use the SAME rng_seeds, so W is identical and the
% conditions are directly comparable.
%
% Lifecycle per condition: SRNNModel2(...) -> build() -> run() -> sample.

clear; clc; close all;

setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% ---- Shared configuration (SRNNModel2 class defaults) ---------------------
% Same rng_seeds across all conditions => identical W (network seeded from
% rng_seeds(1)), so the adaptation regimes are compared on one network.
% level_of_chaos is exposed as a tuning knob: if the heatmaps look static, push
% it above 1 so the dynamics sit nearer/past the edge of chaos and the
% eigenvalues actually wander.
level_of_chaos = 3.0;

T_range    = [0 200];   % simulation window (s)
lle_window = 30;        % compute Benettin LLE over the LAST lle_window seconds
lya_T_interval = [T_range(2) - lle_window, T_range(2)];   % e.g. [170 200]

base_args = { ...
    'n', 300, 'f', 0.5, ...
    'level_of_chaos', level_of_chaos, ...
    'store_full_state', true, ...   % required so eigenvalue_time_series can read S_out
    'lya_method', 'benettin', ...   % finite-time LLE over the last lle_window seconds
    'lya_T_interval', lya_T_interval, ...
    'T_range', T_range, ...
    'rng_seeds', [1 2] };

% Four adaptation regimes (differ only in the adaptation counts). n_b_E is kept
% in {0,1} per the SRNNModel2 assumption.
condition_titles = {'No adaptation', 'SFA only', 'STD only', 'SFA + STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...     % no_adaptation
    {'n_a_E', 3, 'n_b_E', 0}, ...     % sfa_only
    {'n_a_E', 0, 'n_b_E', 1}, ...     % std_only
    {'n_a_E', 3, 'n_b_E', 1} };       % sfa_and_std
n_cond = numel(condition_titles);

%% ---- Sampling parameters --------------------------------------------------
n_samples = 300;    % Jacobian samples through the (post-warmup) run

%% ---- Run each condition and pool eigenvalues ------------------------------
% Conditions run serially; the eig loop inside eigenvalue_time_series uses
% parfor across the sampled time points (the expensive part).
evals_by_cond = cell(n_cond, 1);
lle_by_cond   = zeros(n_cond, 1);   % Benettin finite-time LLE over the last lle_window s
for i = 1:n_cond
    model = SRNNModel2(base_args{:}, condition_args{i}{:});
    model.build();
    model.run();   % run() also computes the Benettin LLE over lya_T_interval

    lle_by_cond(i) = model.lya_results.LLE;

    % Sample the Jacobian over the post-warmup window (t >= 0). Time vectors
    % typically start negative to allow transient settling before analysis.
    t_start = max(0, model.t_out(1));
    t_end   = model.t_out(end);
    J_times_sec = linspace(t_start, t_end, n_samples);

    evals_by_cond{i} = model.eigenvalue_time_series(J_times_sec, 'use_parallel', true);
    fprintf('Condition "%s": pooled %d eigenvalues, LLE = %.4f\n', ...
        condition_titles{i}, numel(evals_by_cond{i}), lle_by_cond(i));
end

%% ---- Report ---------------------------------------------------------------
fprintf('\nBenettin finite-time LLE over the last %g s ([%g %g]):\n', ...
    lle_window, lya_T_interval(1), lya_T_interval(2));
for i = 1:n_cond
    fprintf('  %-14s LLE = %+.4f\n', condition_titles{i}, lle_by_cond(i));
end

%% ---- Save for the plotting script -----------------------------------------
save_file = fullfile(this_dir, 'eig_heatmap_data.mat');
save(save_file, 'evals_by_cond', 'condition_titles', ...
    'base_args', 'condition_args', 'level_of_chaos', 'n_samples', ...
    'lle_by_cond', 'lle_window', 'lya_T_interval');
fprintf('\nSaved: %s\n', save_file);

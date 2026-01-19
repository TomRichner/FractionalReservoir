%% VAR_SRNN_comparison.m
% Simulate SRNN networks with different adaptation/depression configurations
% and export data for VAR analysis using the cscs_dynamics pipeline.
%
% This script runs 4 conditions:
%   1. None: No adaptation, no STD
%   2. SFA:  Spike-frequency adaptation only
%   3. STD:  Short-term depression only
%   4. Both: Both SFA and STD
%
% Each condition has baseline (no stim) and stim periods.
% Output is saved in cscs_dynamics-compatible format for VAR analysis.
%
% See also: export_SRNN_to_CSCS_format, SRNN_VAR_subject_structure

close all; clear; clc;
setup_paths();

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 1.25);
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTitleFontWeight', 'normal');

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Simulation parameters
sim_config = struct();

% Quick test settings: 500 Hz, 30 sec baseline, 30 sec stim
sim_config.fs = 100;               % Hz (no downsampling needed)
sim_config.T_baseline = 50;        % 30 seconds without stimulation
sim_config.T_stim = 50;            % 30 seconds with stimulation
sim_config.T_transient = 10;        % Seconds to discard from start of each block

% Development settings (comment out for production)
% sim_config.fs = 200;             % Hz
% sim_config.T_baseline = 30;      % seconds
% sim_config.T_stim = 30;          % seconds

% Neuron subset (empty = all, or specify number for random subset)
sim_config.n_neurons_subset = [];   % Not used when n is already small

% Network parameters
sim_config.n = 500;                  % Total neurons (reduced for speed)
sim_config.level_of_chaos = 1.7;    % Level of chaos
sim_config.u_ex_scale = 0.25;        % Stimulus scaling
sim_config.sigma_noise = 1e-6;       % Process noise std (set to 0 to disable)

% Random seeds for reproducibility
sim_config.rng_seeds = [105 25 33 44 55];
% sim_config.rng_seeds = [105 25 33 44 55]+1;

% Time config for plotting (negative time for settling before stim)
sim_config.T_settle = -20;          % Pre-simulation settling time

% Lyapunov computation
sim_config.compute_lyapunov = true;
sim_config.Lya_method = 'benettin';  % 'benettin', 'qr', or 'none'

% Output directory
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'VAR_SRNN');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Save figures and workspace
save_figs = true;
save_workspace = false;

%% ========================================================================
%  DEFINE CONDITIONS
%  ========================================================================

conditions = struct();

% Condition 1: No adaptation, no STD
conditions(1).name = 'SRNN_none';
conditions(1).n_a_E = 0;
conditions(1).n_b_E = 0;
conditions(1).description = 'No adaptation, no STD';

% Condition 2: SFA only
conditions(2).name = 'SRNN_SFA';
conditions(2).n_a_E = 3;
conditions(2).n_b_E = 0;
conditions(2).description = 'Spike-frequency adaptation only';

% Condition 3: STD only
conditions(3).name = 'SRNN_STD';
conditions(3).n_a_E = 0;
conditions(3).n_b_E = 1;
conditions(3).description = 'Short-term depression only';

% Condition 4: Both SFA and STD
conditions(4).name = 'SRNN_both';
conditions(4).n_a_E = 3;
conditions(4).n_b_E = 1;
conditions(4).description = 'Both SFA and STD';

n_conditions = length(conditions);

%% ========================================================================
%  RUN SIMULATIONS
%  ========================================================================

fprintf('=== VAR_SRNN_comparison ===\n');
fprintf('Running %d conditions for VAR analysis\n\n', n_conditions);
fprintf('Configuration:\n');
fprintf('  Sampling rate: %d Hz\n', sim_config.fs);
fprintf('  Baseline duration: %.1f sec\n', sim_config.T_baseline);
fprintf('  Stim duration: %.1f sec\n', sim_config.T_stim);
fprintf('  Compute Lyapunov: %s\n', mat2str(sim_config.compute_lyapunov));
fprintf('  Output directory: %s\n\n', output_dir);

% Store all results for plotting
all_results = cell(n_conditions, 1);
combined_runs = cell(n_conditions, 1);

for c = 1:n_conditions
    cond = conditions(c);

    fprintf('===========================================================\n');
    fprintf('Condition %d/%d: %s\n', c, n_conditions, cond.name);
    fprintf('  %s\n', cond.description);
    fprintf('  n_a_E = %d, n_b_E = %d\n', cond.n_a_E, cond.n_b_E);
    fprintf('===========================================================\n\n');

    % Clear persistent variables
    clear_SRNN_persistent();

    % Run simulation for this condition
    [result, lya_results, plot_data] = run_single_condition(cond, sim_config);

    % Store result
    all_results{c} = result;

    % Store for combined plotting
    run_data = struct();
    run_data.plot_data = plot_data;
    run_data.params = result.params;
    run_data.lya_results = lya_results;
    run_data.Lya_method = sim_config.Lya_method;
    combined_runs{c} = run_data;

    % Export to cscs_dynamics format
    export_SRNN_to_CSCS_format(result, cond.name, output_dir, sim_config);

    fprintf('\nCondition %s complete.\n\n', cond.name);
end

%% ========================================================================
%  PLOTTING
%  ========================================================================

fprintf('Generating plots...\n');

% Individual time series plots
for c = 1:n_conditions
    run_cell = {combined_runs{c}};
    fig_title = sprintf('%s (n_a_E=%d, n_b_E=%d)', ...
        conditions(c).name, conditions(c).n_a_E, conditions(c).n_b_E);

    figure('Position', [100 100 1200 800]);
    [~, ~] = plot_SRNN_combined_tseries(run_cell, 3, {'u_ex', 'x', 'br', 'a', 'b', 'lya'});
    sgtitle(fig_title, 'FontSize', 18);
end

% Combined comparison: none vs both
fprintf('Creating combined comparison plot (none vs both)...\n');
comparison_runs = {combined_runs{1}, combined_runs{4}};  % none vs both
figure('Position', [100 100 1400 900]);
[fig_combined, ~] = plot_SRNN_combined_tseries(comparison_runs, 5, {'u_ex', 'x', 'br', 'a', 'b', 'lya'});
sgtitle('SRNN Comparison: No Adaptation vs SFA+STD', 'FontSize', 18);

% Add letters to subplots
if exist('AddLetters2Plots', 'file')
    AddLetters2Plots(fig_combined, {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'}, ...
        'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.065, 'VShift', -0.025);
end

%% ========================================================================
%  SAVE FIGURES
%  ========================================================================

if save_figs
    figs_dir = fullfile(output_dir, 'figs');
    if ~exist(figs_dir, 'dir')
        mkdir(figs_dir);
    end

    % Save all open figures
    fig_handles = findobj('Type', 'figure');
    [~, idx] = sort([fig_handles.Number]);
    fig_handles = fig_handles(idx);

    save_some_figs_to_folder_2(figs_dir, 'VAR_SRNN', [fig_handles.Number], {'fig', 'png'});
    fprintf('Figures saved to: %s\n', figs_dir);
end

%% ========================================================================
%  GENERATE SUBJECT STRUCTURE
%  ========================================================================

fprintf('Generating subject structure...\n');
SRNN_VAR_subject_structure(conditions, sim_config, output_dir);

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n=== SIMULATION COMPLETE ===\n');
fprintf('Data saved to: %s\n', output_dir);
fprintf('Next steps:\n');
fprintf('  1. Run SRNN_VAR_load_simulated.m to preprocess data\n');
fprintf('  2. Run CSCS_run_ICA.m for ICA decomposition\n');
fprintf('  3. Run CSCS_run_VAR.m for VAR analysis\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [result, lya_results, plot_data] = run_single_condition(cond, sim_config)
%RUN_SINGLE_CONDITION Run SRNN simulation for a single condition
%
%   [result, lya_results, plot_data] = run_single_condition(cond, sim_config)
%
%   Simulates both baseline (no stim) and stim periods with optional
%   Lyapunov exponent computation.

% Initialize random number generator
rng(sim_config.rng_seeds(1));

% Set up base parameters
params = struct();
params.n = sim_config.n;
params.indegree = 20;
params.d = params.indegree / params.n;
params.f = 0.5;
params.G_stdev = 1 / sqrt(params.indegree);
params.mu_E = 1;
params.level_of_chaos = sim_config.level_of_chaos;
params.rng_seeds = sim_config.rng_seeds;

% Set adaptation parameters from condition
params.n_a_E = cond.n_a_E;
params.n_a_I = 0;
params.tau_a_E = logspace(log10(0.5), log10(10), max(1, params.n_a_E));
if params.n_a_E == 0
    params.tau_a_E = [];
end
params.tau_a_I = [];
params.c_E = 0.15 / 3;
params.c_I = 0.1;

% Set STD parameters from condition
params.n_b_E = cond.n_b_E;
params.n_b_I = 0;
params.tau_b_E_rel = 0.25;
params.tau_b_I_rel = 0.25;
params.tau_b_E_rec = 2;
params.tau_b_I_rec = 2;

% Activation function
S_a = 0.9;
S_c = 0.3;
params.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

% Dependent E/I parameters
params.n_E = round(params.f * params.n);
params.n_I = params.n - params.n_E;
params.N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + ...
    params.n_E * params.n_b_E + params.n_I * params.n_b_I + params.n;
params.E_indices = 1:params.n_E;
params.I_indices = params.n_E + 1:params.n;
imbalance = 1;
params.mu_I = -imbalance * params.f * params.mu_E / (1 - params.f);

% Create weight matrix
[W, ~, ~, ~] = create_W_matrix(params);
W_eigs = eig(W);
abscissa_0 = max(real(W_eigs));
gamma = 1 / abscissa_0;
params.W = params.level_of_chaos * gamma * W;
params.tau_d = 0.1;

% Time parameters
fs = sim_config.fs;
dt = 1 / fs;
T_total = sim_config.T_baseline + sim_config.T_stim;
T_settle = sim_config.T_settle;
T_range = [T_settle, T_total];

% Initial conditions
S0 = initialize_state(params);

% Generate external input with stim pattern
input_config.n_steps = 2;  % Baseline + Stim
input_config.step_density = 0.1;
input_config.amp = 0.5;
input_config.no_stim_pattern = [true, false];  % First half no stim, second half stim
input_config.intrinsic_drive = zeros(params.n, 1);

[u_stim, t_stim] = generate_external_input(params, T_total, fs, params.rng_seeds(2), input_config);

% Prepend settling period with zeros
if T_settle < 0
    t_pre = (T_settle:dt:-dt)';
    u_pre = zeros(params.n, length(t_pre));
    t_ex = [t_pre; t_stim];
    u_ex = [u_pre, u_stim];
else
    t_ex = t_stim;
    u_ex = u_stim;
end
u_ex = u_ex * sim_config.u_ex_scale;

% Generate process noise (pre-generated for reproducibility with ODE solver)
sigma_noise = sim_config.sigma_noise;
if sigma_noise > 0
    rng(sim_config.rng_seeds(3));  % Use seed 3 for noise
    noise_ex = sigma_noise * randn(params.n, length(t_ex));
else
    noise_ex = zeros(params.n, length(t_ex));
end

% Integrate
ode_solver = @ode45;
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params, noise_ex);
jac_wrapper = @(t, S) compute_Jacobian_fast(S, params);
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', dt, 'Jacobian', jac_wrapper);

fprintf('  Integrating equations...\n');
tic;
[t_out, S_out] = ode_solver(rhs, t_ex, S0, opts);
integration_time = toc;
fprintf('  Integration complete in %.1f sec\n', integration_time);

% Compute Lyapunov exponents
if sim_config.compute_lyapunov && ~strcmp(sim_config.Lya_method, 'none')
    fprintf('  Computing Lyapunov exponents (%s method)...\n', sim_config.Lya_method);
    T_interval = [max(0, T_range(1)+1), T_range(2)];
    lya_results = compute_lyapunov_exponents(sim_config.Lya_method, S_out, t_out, ...
        dt, fs, T_interval, params, opts, ode_solver, rhs, t_ex, u_ex);
    fprintf('  LLE = %.4f\n', lya_results.LLE);
else
    lya_results = struct('LLE', NaN, 't_lya', [], 'local_lya', []);
end

% Decimate for plotting
plot_deci = max(1, round(fs/20));
[t_plot, S_plot, plot_indices] = decimate_states(t_out, S_out, plot_deci);

% Unpack states for plotting
[x_plot, a_plot, b_plot, r_plot, br_plot] = unpack_and_compute_states(S_plot, params);

% Split external input for plotting
u_ex_plot = u_ex(:, plot_indices);
u_plot.E = u_ex_plot(params.E_indices, :);
u_plot.I = u_ex_plot(params.I_indices, :);

% Trim to positive time for plotting
keep_mask = t_plot >= 0;
t_plot = t_plot(keep_mask);
u_plot.E = u_plot.E(:, keep_mask);
u_plot.I = u_plot.I(:, keep_mask);
x_plot = trim_struct_data(x_plot, 2, keep_mask);
r_plot = trim_struct_data(r_plot, 2, keep_mask);
b_plot = trim_struct_data(b_plot, 2, keep_mask);
br_plot = trim_struct_data(br_plot, 2, keep_mask);
a_plot = trim_struct_data(a_plot, 3, keep_mask);

% Trim Lyapunov results
if ~isempty(lya_results.t_lya)
    keep_mask_lya = lya_results.t_lya >= 0;
    lya_results.t_lya = lya_results.t_lya(keep_mask_lya);
    if isfield(lya_results, 'local_lya')
        lya_results.local_lya = lya_results.local_lya(keep_mask_lya);
    end
    if isfield(lya_results, 'finite_lya')
        lya_results.finite_lya = lya_results.finite_lya(keep_mask_lya);
    end
end

% Create plot_data struct
plot_data = struct();
plot_data.t = t_plot;
plot_data.u = u_plot;
plot_data.x = x_plot;
plot_data.r = r_plot;
plot_data.a = a_plot;
plot_data.b = b_plot;
plot_data.br = br_plot;

% Unpack full states for VAR export (without decimation)
[x, a, b, r] = unpack_and_compute_states(S_out, params);

% Trim to positive time for export
export_mask = t_out >= 0;
t_export = t_out(export_mask);
x = trim_struct_data(x, 2, export_mask);
r = trim_struct_data(r, 2, export_mask);
a = trim_struct_data(a, 3, export_mask);
b = trim_struct_data(b, 2, export_mask);

% Separate baseline and stim periods
T_baseline = sim_config.T_baseline;
baseline_mask = t_export < T_baseline;
stim_mask = t_export >= T_baseline;

% Store results
result = struct();
result.params = params;
result.t = t_export;
result.S_out = S_out(export_mask, :);
result.x = x;
result.r = r;
result.a = a;
result.b = b;
result.u_ex = u_ex(:, export_mask);
result.t_ex = t_export;
result.fs = fs;
result.T_baseline = T_baseline;
result.baseline_mask = baseline_mask;
result.stim_mask = stim_mask;
result.lya_results = lya_results;
end

function s_out = trim_struct_data(s_in, dim, mask)
%TRIM_STRUCT_DATA Helper to trim fields of a struct along a dimension
s_out = s_in;
fields = fieldnames(s_in);
for i = 1:length(fields)
    val = s_in.(fields{i});
    if ~isempty(val)
        if dim == 2
            s_out.(fields{i}) = val(:, mask);
        elseif dim == 3
            s_out.(fields{i}) = val(:, :, mask);
        end
    end
end
end

function [t_out, S_out, params, lya_results, plot_data] = full_SRNN_run_SRNNModel(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config)
% FULL_SRNN_RUN_SRNNMODEL Runs the SRNN simulation using SRNNModel class.
%
% This is a refactored version that uses the SRNNModel class internally
% for cleaner code and improved W matrix parameterization via RMTMatrix.
%
% Inputs:
%   u_ex_scale    - Scaling factor for external input
%   n_a_E         - Number of adaptation time constants for E neurons
%   n_b_E         - Number of STD timescales for E neurons
%   level_of_chaos- Level of chaos parameter
%   rng_seeds     - Seeds for random number generation
%   save_dir      - Directory path to save outputs
%   save_figs     - Boolean, whether to save figures
%   save_workspace- Boolean, whether to save the workspace
%   note          - Optional string for naming figures (default: 'SRNN')
%   time_config   - Optional struct with fields:
%                   .T_range - [start, end] simulation time (default: [0, 50])
%                   .T_plot  - [start, end] plotting time (default: T_range)
%                   .J_periods - logical vector for which periods to compute
%                                J_eff/eigenspectrum (default: all true)

if nargin < 9 || isempty(note)
    note = 'SRNN';
end

if nargin < 10 || isempty(time_config)
    time_config = struct();
end
if ~isfield(time_config, 'T_range')
    time_config.T_range = [0, 50];
end
if ~isfield(time_config, 'T_plot')
    time_config.T_plot = time_config.T_range;
end
% J_periods default set after n_steps is known (see below)

%% Create and configure SRNNModel
model = SRNNModel();

% Network architecture (explicitly set for hand-tuning)
model.n = 300;
model.indegree = 100;
model.f = 0.5;
model.tau_d = 0.1;
model.fs = 100;
model.c_E = 0.25/3;

% RMT tilde-notation parameters (Harris 2023)
% Default val gives R=1: 1 / sqrt(n * alpha * (2 - alpha)) where alpha = indegree/n
alpha = model.indegree / model.n;
default_tilde_val = 1 / sqrt(model.n * alpha * (2 - alpha));
model.mu_E_tilde = 3.5*default_tilde_val;
model.mu_I_tilde = -3.5*default_tilde_val;
model.sigma_E_tilde = default_tilde_val/1;
model.sigma_I_tilde = default_tilde_val/1;
model.E_W = -0.5 / sqrt(model.n * alpha * (2 - alpha));
% model.zrs_mode = 'Partial_SZRS';  % ZRS mode: 'none', 'ZRS', 'SZRS', 'Partial_SZRS'
model.zrs_mode = 'none';


model.level_of_chaos = level_of_chaos;
model.rescale_by_abscissa = false;  % Scale W so abscissa matches level_of_chaos

% Adaptation parameters
model.n_a_E = n_a_E;
model.n_a_I = 0;
model.tau_a_E = logspace(log10(0.1), log10(10), n_a_E);

% STD parameters
model.n_b_E = n_b_E;
model.n_b_I = 0;
model.tau_b_E_rec = 1;  % Match old value
model.tau_b_E_rel = 0.5;

% Activation function
S_a = 0.9;
S_c = 0.40;
model.S_a = S_a;
model.S_c = S_c;
model.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
model.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

% Simulation settings
model.T_range = time_config.T_range;
model.T_plot = time_config.T_plot;
model.rng_seeds = rng_seeds;
model.u_ex_scale = u_ex_scale;

% Lyapunov settings
model.lya_method = 'benettin';
model.filter_local_lya = true;  % Filter before decimation to avoid edge effects
model.store_full_state = true;  % Need full state for Jacobian plots

% Input configuration (match old defaults)
model.input_config.n_steps = 3;
model.input_config.positive_only = true;
model.input_config.step_density_E = 0.15;
model.input_config.step_density_I = 0;
model.input_config.amp = 0.5;
model.input_config.no_stim_pattern = false(1, 3);
model.input_config.no_stim_pattern(1:2:end) = true;
model.input_config.intrinsic_drive = 0 * ones(model.n, 1);

%% Build and run model
model.build();
model.run();

%% Extract outputs for compatibility
t_out = model.t_out;
S_out = model.S_out;
params = model.get_params();
lya_results = model.lya_results;
plot_data = model.plot_data;

%% Plot time series using model method
model.plot();

%% Compute Jacobian eigenvalues at multiple time points
% Sample at the center of each period (both stim and non-stim)
T_stim = time_config.T_range(2);
n_steps = model.input_config.n_steps;
step_period = fix(T_stim / n_steps);

% Set default J_periods if not specified (all periods)
if ~isfield(time_config, 'J_periods')
    time_config.J_periods = true(1, n_steps);
end
J_periods = time_config.J_periods;

J_times_sec = [];
for k = 1:n_steps
    if J_periods(k)  % Only include periods where J_periods is true
        t_center = (k-1)*step_period + step_period/2;
        J_times_sec = [J_times_sec, t_center];
    end
end

% Convert times in seconds to indices
J_times = round((J_times_sec - t_out(1)) * model.fs) + 1;
J_times = unique(J_times);

fprintf('Computing Jacobian at %d time points...\n', length(J_times));
J_array = compute_Jacobian_at_indices(S_out, J_times, params);

% Compute eigenvalues for each Jacobian
eigenvalues_all = cell(length(J_times), 1);
for i = 1:length(J_times)
    eigenvalues_all{i} = eig(J_array(:,:,i));
end

%% Plot eigenvalue distributions on complex plane
n_plots = length(J_times);
if n_plots <= 4
    n_rows = 1;
    n_cols = n_plots;
else
    n_cols = ceil(sqrt(n_plots));
    n_rows = ceil(n_plots / n_cols);
end

figure('Position', [1312, 526, 600, 360]);
ax_handles = zeros(n_plots, 1);

% % Compute global axis limits across all eigenvalue sets
% all_real = [];
% all_imag = [];
% for i = 1:n_plots
%     evals = eigenvalues_all{i};
%     all_real = [all_real; real(evals)];
%     all_imag = [all_imag; imag(evals)];
% end
% global_xlim = [min(all_real), max(all_real)];
% global_ylim = [min(all_imag), max(all_imag)];

% % Add some padding (10% of range on each side)
% x_range = diff(global_xlim);
% y_range = diff(global_ylim);
% global_xlim = global_xlim + [-0.1, 0.1] * x_range;
% global_ylim = global_ylim + [-0.1, 0.1] * y_range;

% Hard-coded global axis limits for eigenvalue plots.  hard coded so they match between SFA and SFA+STD comparisons
global_xlim = [-28.0, 4.5];
global_ylim = [-17.0, 17.0];

for i = 1:n_plots
    ax_handles(i) = subplot(n_rows, n_cols, i);
    evals = eigenvalues_all{i};
    time_val = t_out(J_times(i));
    ax_handles(i) = plot_eigenvalues(evals, ax_handles(i), time_val, global_xlim, global_ylim);
    set(ax_handles(i), 'Color', 'none');
end

% Link axes of all subplots
linkaxes(ax_handles, 'xy');

%% Compute color limits from static W matrix (W sets the scale)
omit_diagonal_in_J_eff = false;
normalize_tau_d_for_J_eff = false;

% Prepare W matrix for plotting
W_plot = full(model.W);
if omit_diagonal_in_J_eff
    W_plot = W_plot - diag(diag(W_plot));
end

% Compute color limits from W (use 95th percentile of absolute non-zero values)
if omit_diagonal_in_J_eff
    W_values = W_plot(~eye(size(W_plot)));
else
    W_values = W_plot(:);
end
nonzero_abs_values = abs(W_values(W_values ~= 0));
max_abs = prctile(nonzero_abs_values, 95);
global_clim = [-max_abs, max_abs];

%% Plot static W matrix figure
figure('Position', [1312, 940, 600, 310]);

subplot(n_rows, n_cols, 1);
imagesc(W_plot);
colormap(redwhiteblue_colormap(256));
clim([-0.5, 0.5]);  % Hard-coded clim for W
axis square;
set(gca, 'XTick', [], 'YTick', []);
box off;
set(gca, 'Color', 'none');
set(gca, 'XColor', 'white', 'YColor', 'white', 'Layer', 'bottom');

%% Histogram of E and I synaptic weights from W
figure('Position', [100, 500, 400, 250]);

W_E = full(model.W(:, params.E_indices));
W_E_nonzero = W_E(W_E ~= 0);

W_I = full(model.W(:, params.I_indices));
W_I_nonzero = W_I(W_I ~= 0);

all_weights = [W_E_nonzero(:); W_I_nonzero(:)];
n_bins = 30;
bin_edges = linspace(min(all_weights), max(all_weights), n_bins + 1);

hold on;
histogram(W_E_nonzero, bin_edges, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
histogram(W_I_nonzero, bin_edges, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold off;

xlabel('Weight');
ylabel('Count');
legend('E', 'I', 'Location', 'northeast');
legend boxoff;
set(gca, 'Color', 'none');
box off;

%% Compute and plot J_eff at J_times (using W's color limits)
fprintf('Computing J_eff at %d time points...\n', length(J_times));

% Create params for J_eff plotting (optionally with normalized tau_d)
if normalize_tau_d_for_J_eff
    params_J_eff_plot = params;
    params_J_eff_plot.tau_d = 1;
else
    params_J_eff_plot = params;
end

% Compute J_eff for each time point
J_eff_array = zeros(params.n, params.n, length(J_times));
for i = 1:length(J_times)
    J_eff_array(:,:,i) = full(compute_J_eff(S_out(J_times(i),:)', params_J_eff_plot));

    if omit_diagonal_in_J_eff
        J_eff_array(:,:,i) = J_eff_array(:,:,i) - diag(diag(J_eff_array(:,:,i)));
    end
end

% Create J_eff figure with same layout as eigenvalue figure
figure('Position', [100, 100, 600, 310]);

for i = 1:n_plots
    subplot(n_rows, n_cols, i);
    imagesc(J_eff_array(:,:,i));
    colormap(redwhiteblue_colormap(256));
    clim([-5, 5]);  % Hard-coded clim for J_eff
    axis square;
    set(gca, 'XTick', [], 'YTick', []);
    box off;
    set(gca, 'Color', 'none');
    set(gca, 'XColor', 'white', 'YColor', 'white', 'Layer', 'bottom');
end

%% Create separate colorbar figure for W
figure('Position', [100, 346, 285, 154], 'Color', 'white');
ax = axes('Position', [0.3, 0.1, 0.3, 0.8]);
colormap(redwhiteblue_colormap(256));
cb = colorbar('Location', 'east');
clim([-0.5, 0.5]);
set(gca, 'Visible', 'off', 'Color', 'none');
set(cb, 'AxisLocation', 'out', 'Ticks', [-0.5, 0.5]);
ylabel(cb, 'W', 'Interpreter', 'tex', 'FontSize', 14);

%% Create separate colorbar figure for J_eff
figure('Position', [100, 500, 285, 154], 'Color', 'white');
ax = axes('Position', [0.3, 0.1, 0.3, 0.8]);
colormap(redwhiteblue_colormap(256));
cb = colorbar('Location', 'east');
clim([-5, 5]);
set(gca, 'Visible', 'off', 'Color', 'none');
set(cb, 'AxisLocation', 'out', 'Ticks', [-5, 5]);
ylabel(cb, 'J_{eff}', 'Interpreter', 'tex', 'FontSize', 14);

% %% Plot J_eff as directed graph
% figure('Position', [300, 400, 900, 400]);
%
% for i = 1:n_plots
%     subplot(n_rows, n_cols, i);
%     plot_J_eff_graph(J_eff_array(:,:,i), max_abs, global_clim);
%     set(gca, 'Color', 'none');
% end

%% eigenspectra plot of I-W
figure;
ax = gca;
eigs_diff = eig(-eye(params.n) + model.W);
plot_eigenvalues(eigs_diff, ax, 0);

%% Save results
if save_figs
    if ~isempty(save_dir)
        figs_folder = fullfile(save_dir, 'figs');
        if ~exist(figs_folder, 'dir')
            mkdir(figs_folder);
        end

        fig_handles = findobj('Type', 'figure');
        [~, idx] = sort([fig_handles.Number]);
        fig_handles = fig_handles(idx);
        fig_vec = [fig_handles.Number];

        save_some_figs_to_folder_2(figs_folder, note, fig_vec, {'fig', 'svg', 'png', 'pdf'});
    else
        warning('save_figs is true but save_dir is empty. Figures not saved.');
    end
end

if save_workspace
    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        save(fullfile(save_dir, 'SRNN_run.mat'), '-v7.3', '-nocompression');
    else
        warning('save_workspace is true but save_dir is empty. Workspace not saved.');
    end
end

end

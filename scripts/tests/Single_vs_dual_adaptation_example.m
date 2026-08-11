% close all; clear all; clc;  % Commented out for master script compatibility
% setup_paths();  % Commented out - called by master script

% Derive project root from script location for portable paths
% This script lives in scripts/tests/, so the project root is two directories
% up from its folder. (setup_paths is commented out above; resolve standalone.)
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figs_root = fullfile(project_root, 'figs');
output_folder_name = ['srnn_comparison_', datestr(now, 'yyyymmdd_HHMM')];

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextFontSize', 14);
set(groot, 'DefaultLineLineWidth', 0.75);
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTitleFontWeight', 'normal');

% NOTE: This script sets global MATLAB figure defaults that persist for the session.
% Run `reset(groot)` afterward to restore factory defaults if needed.

% Check for master override from run_all_figures.m
if exist('master_save_figs', 'var')
    if strcmp(master_save_figs, 'save_all_figs')
        save_figs = true;
    elseif strcmp(master_save_figs, 'save_no_figs')
        save_figs = false;
    end
end
if ~exist('save_figs', 'var')
    save_figs = false;  % Script default
end
save_workspace = false;

%% Shared simulation parameters
level_of_chaos = 1.0; % gamma from Sompolinsky

u_ex_scale = 1; % can change scale of stimulus.

rng_seeds = [42 42]; % seed 1 is for connection matrix, seed 2 is for stimulus

time_config.J_periods = [false true true];  % three periods: no-stim, stim, no-stim
time_config.T_range = [-15, 45]; % seconds
time_config.T_plot = [7.5, 45];  % seconds. Trim off half of first no-stim period

combined_runs = {}; % cell array to hold multiple simulations

%% ==================== PHASE 1: BUILD AND RUN ====================

%% Run 1: spike frequency adaptation (SFA) only
close all;
note = 'SFA_only';

n_a_E = 3; % three timeconstants of SFA
n_b_E = 0; % no short-term synaptic depression

save_dir = fullfile(figs_root, output_folder_name, note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);

% Create and configure SRNNModel
model1 = SRNNModel2();

% Network architecture (exactly matching full_SRNN_run_SRNNModel.m)
model1.n = 300;
model1.indegree = 100;
model1.f = 0.50;
model1.tau_d = 0.1;
model1.c_E = 0.25/3;

% RMT tilde-notation parameters (Harris 2023)
% F = 1/sqrt(N*alpha*(2-alpha)), the scaling factor yielding R=1 when
% all tilde parameters are equal (see parameter_table.md)
model1.mu_E_tilde_relative = 3.5;
model1.mu_I_tilde_relative = -3.5;
model1.sigma_E_tilde_relative = 1;
model1.sigma_I_tilde_relative = 1;
model1.E_W_relative = -0.5;
model1.zrs_mode = 'none';

model1.level_of_chaos = level_of_chaos;
model1.rescale_by_abscissa = false;

% Adaptation parameters
model1.n_a_E = n_a_E;
model1.n_a_I = 0;
model1.tau_a_E = logspace(log10(0.1), log10(10), n_a_E);

% STD parameters
model1.n_b_E = n_b_E;
model1.n_b_I = 0;
model1.tau_b_E_rec = 1;
model1.tau_b_E_rel = 0.5;

% Activation function
S_a = 0.9;
S_c = 0.40;
model1.S_a = S_a;
model1.S_c = S_c;
model1.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
model1.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

% Simulation settings
model1.T_range = time_config.T_range;
model1.T_plot = time_config.T_plot;
model1.rng_seeds = rng_seeds;
model1.u_ex_scale = u_ex_scale;

% Lyapunov settings
model1.lya_method = 'benettin';
model1.filter_local_lya = true;
model1.store_full_state = true;

% Input configuration
model1.input_config.n_steps = 3;
model1.input_config.positive_only = true;
model1.input_config.step_density_E = 0.15;
model1.input_config.step_density_I = 0;
model1.input_config.amp = 0.5;
model1.input_config.no_stim_pattern = false(1, 3);
model1.input_config.no_stim_pattern(1:2:end) = true;
model1.input_config.intrinsic_drive = 0 * ones(model1.n, 1);

% Build and run model
model1.build();
model1.run();

% Extract outputs
t_out_1 = model1.t_out;
S_out_1 = model1.S_out;
params_1 = model1.get_params();
lya_1 = model1.lya_results;
plot_data_1 = model1.plot_data;

run1.plot_data = plot_data_1;
run1.params = params_1;
run1.lya_results = lya_1;
run1.Lya_method = 'benettin';
combined_runs{1} = run1;

%% Run 2: SFA + STD
note = 'STD_and_SFA';
n_a_E = 3;
n_b_E = 1;

save_dir = fullfile(figs_root, output_folder_name, note);
fprintf('Running SRNN with u_ex_scale=%g, n_a_E=%d, n_b_E=%d, level_of_chaos=%g\n', u_ex_scale, n_a_E, n_b_E, level_of_chaos);

% Create and configure SRNNModel (same as Run 1, except n_b_E = 1)
model2 = SRNNModel2();

model2.n = 300;
model2.indegree = 100;
model2.f = 0.50;
model2.tau_d = 0.1;
model2.c_E = 0.25/3;

model2.mu_E_tilde_relative = 3.5;
model2.mu_I_tilde_relative = -3.5;
model2.sigma_E_tilde_relative = 1;
model2.sigma_I_tilde_relative = 1;
model2.E_W_relative = -0.5;
model2.zrs_mode = 'none';

model2.level_of_chaos = level_of_chaos;
model2.rescale_by_abscissa = false;

model2.n_a_E = n_a_E;
model2.n_a_I = 0;
model2.tau_a_E = logspace(log10(0.1), log10(10), n_a_E);

model2.n_b_E = n_b_E;  % KEY DIFFERENCE: 1 instead of 0
model2.n_b_I = 0;
model2.tau_b_E_rec = 1;
model2.tau_b_E_rel = 0.5;

S_a = 0.9;
S_c = 0.40;
model2.S_a = S_a;
model2.S_c = S_c;
model2.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
model2.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

model2.T_range = time_config.T_range;
model2.T_plot = time_config.T_plot;
model2.rng_seeds = rng_seeds;
model2.u_ex_scale = u_ex_scale;

model2.lya_method = 'benettin';
model2.filter_local_lya = true;
model2.store_full_state = true;

model2.input_config.n_steps = 3;
model2.input_config.positive_only = true;
model2.input_config.step_density_E = 0.15;
model2.input_config.step_density_I = 0;
model2.input_config.amp = 0.5;
model2.input_config.no_stim_pattern = false(1, 3);
model2.input_config.no_stim_pattern(1:2:end) = true;
model2.input_config.intrinsic_drive = 0 * ones(model2.n, 1);

model2.build();
model2.run();

t_out_4 = model2.t_out;
S_out_4 = model2.S_out;
params_4 = model2.get_params();
lya_4 = model2.lya_results;
plot_data_4 = model2.plot_data;

run4.plot_data = plot_data_4;
run4.params = params_4;
run4.lya_results = lya_4;
run4.Lya_method = 'benettin';
combined_runs{2} = run4;

%% ==================== PHASE 2: POST-PROCESSING ====================

% Compute eigenvalue sets for both runs
eig_data_1 = compute_eigenvalue_sets(model1, params_1, time_config);
eig_data_4 = compute_eigenvalue_sets(model2, params_4, time_config);

% Define which subplots share common axis limits
% Group -1: omitted from plotting (static W eigenvalues)
% Group 2: Jacobian eigenvalues (subplots 1-2) — shared limits across runs
common_lim_group = [-1, 2, 2];
padding_factor = 0.05;

% Compute axis limits across both runs
[xlims, ylims] = compute_grouped_axis_limits({eig_data_1, eig_data_4}, ...
    common_lim_group, padding_factor);

% Static W eigenvalue spectrum: plot only subplot 1, skip Jacobian subplots
common_lim_group_static = [1, -1, -1];
[xlims_static, ylims_static] = compute_grouped_axis_limits( ...
    {eig_data_1, eig_data_4}, common_lim_group_static, padding_factor);

%% ==================== PHASE 3: PLOTTING ====================

%% Plot Run 1: SFA only
model1.plot();

plot_jacobian_eigenvalues_fig(model1, params_1, time_config, eig_data_1, xlims, ylims, common_lim_group);
plot_jacobian_eigenvalues_fig(model1, params_1, time_config, eig_data_1, xlims_static, ylims_static, common_lim_group_static);
plot_W_matrix_fig(model1);
plot_weight_histogram_fig(model1, params_1);
plot_J_eff_matrices_fig(model1, S_out_1, t_out_1, params_1, time_config);
plot_W_colorbar_fig();
plot_J_eff_colorbar_fig();

%% Plot Run 2: SFA + STD
model2.plot();

plot_jacobian_eigenvalues_fig(model2, params_4, time_config, eig_data_4, xlims, ylims, common_lim_group);
plot_jacobian_eigenvalues_fig(model2, params_4, time_config, eig_data_4, xlims_static, ylims_static, common_lim_group_static);
plot_W_matrix_fig(model2);
plot_weight_histogram_fig(model2, params_4);
plot_J_eff_matrices_fig(model2, S_out_4, t_out_4, params_4, time_config);
plot_W_colorbar_fig();
plot_J_eff_colorbar_fig();

%% Plot Combined
[fig_handle, ~] = plot_SRNN_combined_tseries(combined_runs, 3, {'u_ex', 'x', 'br', 'a', 'b', 'lya'});

% Don't add letters to subplots for this manuscript.  They will be added in Affinity by hand for better alignment.
% AddLetters2Plots(fig_handle, {'(D)', '(E)', '(F)', '(G)', '(H)', '(I)'}, 'FontSize', 12, 'FontWeight', 'normal', 'HShift', -0.06, 'VShift', -0.02);

ylim([-1.9 1.9]) % y limits of the local lyapunov exponent

if save_figs
    save_dir_combined = fullfile(figs_root, output_folder_name);
    save_name_base = 'combined_comparison';

    % Use the existing helper function
    save_some_figs_to_folder_2(save_dir_combined, save_name_base, [], {'fig', 'svg', 'png', 'jp2'});
    fprintf('Combined plot saved to %s\n', save_dir_combined);
end


%% ==================== LOCAL FUNCTIONS ====================

function plot_jacobian_eigenvalues_fig(model, params, time_config, eig_data, xlims, ylims, common_lim_group)
%PLOT_JACOBIAN_EIGENVALUES_FIG Plot pre-computed Jacobian eigenvalues at sample times
%   Entries in common_lim_group with value -1 are omitted from the figure.

t_out = model.t_out;
J_times = eig_data.J_times;
eig_sets = eig_data.eig_sets;
circle_params = eig_data.circle_params;

% Determine which eig_sets to plot (skip group == -1)
plot_mask = (common_lim_group ~= -1);
plot_indices = find(plot_mask);
n_plot = length(plot_indices);

if n_plot <= 4
    n_rows = 1;
    n_cols = n_plot;
else
    n_cols = ceil(sqrt(n_plot));
    n_rows = ceil(n_plot / n_cols);
end

fh = figure;
ax_handles = zeros(n_plot, 1);

for k = 1:n_plot
    idx = plot_indices(k);
    ax_handles(k) = subplot(n_rows, n_cols, k);
    if idx == 1
        % Static eigenvalues of (1/tau_d)*(W - I) — pass circle_params
        ax_handles(k) = plot_eigenvalues(eig_sets{1}, ax_handles(k), 0, xlims(1,:), ylims(1,:), circle_params);
    else
        % Jacobian eigenvalues at sampled time point
        time_val = t_out(J_times(idx - 1));
        ax_handles(k) = plot_eigenvalues(eig_sets{idx}, ax_handles(k), time_val, xlims(idx,:), ylims(idx,:));
    end
    set(ax_handles(k), 'Color', 'none');
end

linkaxes(ax_handles, 'xy');
drawnow
set(fh,'Position',[100, 100, 360*n_plot, 360])
end

function plot_W_matrix_fig(model)
%PLOT_W_MATRIX_FIG Plot static W matrix with imagesc

W_plot = full(model.W);

figure('Position', [1312, 940, 600, 310]);
subplot(1, 1, 1);
imagesc(W_plot);
colormap(redwhiteblue_colormap(256));
clim([-0.5, 0.5]);  % Hard-coded clim for W
axis square;
set(gca, 'XTick', [], 'YTick', []);
box off;
set(gca, 'Color', 'none');
set(gca, 'XColor', 'white', 'YColor', 'white', 'Layer', 'bottom');
end

function plot_weight_histogram_fig(model, params)
%PLOT_WEIGHT_HISTOGRAM_FIG Histogram of E and I synaptic weights

W_E = full(model.W(:, params.E_indices));
W_E_nonzero = W_E(W_E ~= 0);

W_I = full(model.W(:, params.I_indices));
W_I_nonzero = W_I(W_I ~= 0);

all_weights = [W_E_nonzero(:); W_I_nonzero(:)];
n_bins = 50;
bin_edges = linspace(min(all_weights), max(all_weights), n_bins + 1);

figure('Position', [100, 500, 400, 250]);
hold on;
histogram(W_E_nonzero, bin_edges, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
histogram(W_I_nonzero, bin_edges, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold off;

xlabel('Weight');
ylabel('Count');
lg = legend('E', 'I', 'Location', 'northeast');
lg.Position(1) = lg.Position(1) + 0.08;
legend boxoff;
set(gca, 'Color', 'none');
box off;

% Add mu_tilde markers on x-axis
ax = gca;
y_bottom = ax.YLim(1);
mu_E_pos = model.mu_E_tilde + model.E_W;
mu_I_pos = model.mu_I_tilde + model.E_W;
text(mu_E_pos, y_bottom, '$\tilde{\mu}_E$', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);
text(mu_I_pos, y_bottom, '$\tilde{\mu}_I$', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);
end

function plot_J_eff_matrices_fig(model, S_out, t_out, params, time_config)
%PLOT_J_EFF_MATRICES_FIG Compute and plot J_eff at sample times

T_stim = time_config.T_range(2);
n_steps = model.input_config.n_steps;
step_period = fix(T_stim / n_steps);

if isfield(time_config, 'J_periods')
    J_periods = time_config.J_periods;
else
    J_periods = true(1, n_steps);
end

J_times_sec = [];
for k = 1:n_steps
    if J_periods(k)
        t_center = (k-1)*step_period + step_period/2;
        J_times_sec = [J_times_sec, t_center];
    end
end

J_times = round((J_times_sec - t_out(1)) * model.fs) + 1;
J_times = unique(J_times);

n_J_plots = length(J_times);
if n_J_plots <= 4
    n_rows = 1;
    n_cols = n_J_plots;
else
    n_cols = ceil(sqrt(n_J_plots));
    n_rows = ceil(n_J_plots / n_cols);
end

fprintf('Computing J_eff at %d time points...\n', length(J_times));

J_eff_array = zeros(params.n, params.n, length(J_times));
for idx = 1:length(J_times)
    J_eff_array(:,:,idx) = full(compute_J_eff(S_out(J_times(idx),:)', params));
end

figure('Position', [100, 100, 600, 310]);
for i_plot = 1:n_J_plots
    subplot(n_rows, n_cols, i_plot);
    imagesc(J_eff_array(:,:,i_plot));
    colormap(redwhiteblue_colormap(256));
    clim([-5, 5]);  % Hard-coded clim for J_eff
    axis square;
    set(gca, 'XTick', [], 'YTick', []);
    box off;
    set(gca, 'Color', 'none');
    set(gca, 'XColor', 'white', 'YColor', 'white', 'Layer', 'bottom');
end
end

function plot_W_colorbar_fig()
%PLOT_W_COLORBAR_FIG Create separate colorbar figure for W

figure('Position', [100, 346, 285, 154], 'Color', 'white');
ax = axes('Position', [0.3, 0.1, 0.3, 0.8]);
colormap(redwhiteblue_colormap(256));
cb = colorbar('Location', 'east');
clim([-0.5, 0.5]);
set(gca, 'Visible', 'off', 'Color', 'none');
set(cb, 'AxisLocation', 'out', 'Ticks', [-0.5, 0.5]);
ylabel(cb, 'W', 'Interpreter', 'tex', 'FontSize', 14);
end

function plot_J_eff_colorbar_fig()
%PLOT_J_EFF_COLORBAR_FIG Create separate colorbar figure for J_eff

figure('Position', [100, 500, 285, 154], 'Color', 'white');
ax = axes('Position', [0.3, 0.1, 0.3, 0.8]);
colormap(redwhiteblue_colormap(256));
cb = colorbar('Location', 'east');
clim([-5, 5]);
set(gca, 'Visible', 'off', 'Color', 'none');
set(cb, 'AxisLocation', 'out', 'Ticks', [-5, 5]);
ylabel(cb, 'J_{eff}', 'Interpreter', 'tex', 'FontSize', 14);
end

function eig_data = compute_eigenvalue_sets(model, params, time_config)
%COMPUTE_EIGENVALUE_SETS Compute eigenvalue sets for all eigenvalue subplots
%
%   eig_data = compute_eigenvalue_sets(model, params, time_config)
%
%   Returns a struct with:
%     eig_sets      - cell array of eigenvalue vectors (one per subplot)
%                     {1}: static eigenvalues of (1/tau_d)*(-I + W)
%                     {2}, {3}, ...: Jacobian eigenvalues at sampled times
%     J_times       - sample indices into t_out
%     circle_params - struct with .center and .radius for theoretical circle

t_out = model.t_out;
S_out = model.S_out;

% Determine time points for Jacobian evaluation
T_stim = time_config.T_range(2);
n_steps = model.input_config.n_steps;
step_period = fix(T_stim / n_steps);

if isfield(time_config, 'J_periods')
    J_periods = time_config.J_periods;
else
    J_periods = true(1, n_steps);
end

J_times_sec = [];
for k = 1:n_steps
    if J_periods(k)
        t_center = (k-1)*step_period + step_period/2;
        J_times_sec = [J_times_sec, t_center];
    end
end

J_times = round((J_times_sec - t_out(1)) * model.fs) + 1;
J_times = unique(J_times);

% Compute static eigenvalues (subplot 1)
eigs_static = 1/model.tau_d * eig(-eye(params.n) + model.W);

% Compute Jacobian eigenvalues at sampled times (subplots 2+)
fprintf('Computing Jacobian at %d time points...\n', length(J_times));
J_array = compute_Jacobian_at_indices(S_out, J_times, params);

eigenvalues_J = cell(length(J_times), 1);
for idx = 1:length(J_times)
    eigenvalues_J{idx} = eig(J_array(:,:,idx));
end

% Pack into eig_sets: {static, J1, J2, ...}
eig_sets = [{eigs_static}; eigenvalues_J];

% Theoretical circle for static eigenvalues
circle_center = -1 / model.tau_d;
circle_radius = model.R / model.tau_d;
circle_params = struct('center', circle_center, 'radius', circle_radius, 'outlier_threshold', 1.04);

% Return struct
eig_data.eig_sets = eig_sets;
eig_data.J_times = J_times;
eig_data.circle_params = circle_params;
end

function [xlims, ylims] = compute_grouped_axis_limits(eig_data_runs, common_lim_group, padding_factor)
%COMPUTE_GROUPED_AXIS_LIMITS Compute per-subplot axis limits from grouped eigenvalue data
%
%   [xlims, ylims] = compute_grouped_axis_limits(eig_data_runs, common_lim_group, padding_factor)
%
%   Inputs:
%     eig_data_runs    - cell array of eig_data structs (one per run)
%     common_lim_group - vector of group IDs (length = number of subplots per run)
%                        Subplots with the same group ID share axis limits
%                        across all runs.
%     padding_factor   - scalar padding fraction (e.g. 0.05 for 5% on each side)
%
%   Outputs:
%     xlims - n_subplots x 2 matrix of [xmin, xmax]
%     ylims - n_subplots x 2 matrix of [ymin, ymax]

n_subplots = length(common_lim_group);
unique_groups = unique(common_lim_group);

xlims = zeros(n_subplots, 2);
ylims = zeros(n_subplots, 2);

for g = 1:length(unique_groups)
    gid = unique_groups(g);
    if gid == -1
        continue;  % skip omitted subplots
    end
    subplot_mask = (common_lim_group == gid);

    % Collect all eigenvalues across all runs for subplots in this group
    all_eigs = [];
    for r = 1:length(eig_data_runs)
        eig_sets = eig_data_runs{r}.eig_sets;
        for s = 1:n_subplots
            if subplot_mask(s)
                all_eigs = [all_eigs; eig_sets{s}(:)];
            end
        end
    end

    % Compute extent from eigenvalue data
    xmin = min(real(all_eigs));
    xmax = max(real(all_eigs));
    ymin = min(imag(all_eigs));
    ymax = max(imag(all_eigs));

    % Apply padding
    x_range = xmax - xmin;
    y_range = ymax - ymin;
    xmin = xmin - padding_factor * x_range;
    xmax = xmax + padding_factor * x_range;
    ymin = ymin - padding_factor * y_range;
    ymax = ymax + padding_factor * y_range;

    % Ensure xmax >= 0 (preserving rule from plot_eigenvalues.m)
    if xmax < 0
        xmax = 0.05;
    end

    % Assign to all subplots in this group
    xlims(subplot_mask, :) = repmat([xmin, xmax], sum(subplot_mask), 1);
    ylims(subplot_mask, :) = repmat([ymin, ymax], sum(subplot_mask), 1);
end
end

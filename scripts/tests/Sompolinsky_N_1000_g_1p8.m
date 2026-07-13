%% Hopfield_Sompolinsky_chaos_demo.m
% Reproduce the Sompolinsky 1988 result showing that a random neural network
% with asymmetric random weights produces chaos when gamma > 1.
%
% Key parameters:
%   - tanh(x) nonlinearity
%   - No adaptation (n_a_E = n_a_I = n_b_E = n_b_I = 0)
%   - Random Gaussian weights with mean 0 and std 1/sqrt(N) (asymmetric)
%   - level_of_chaos (gamma) = 1.5 > 1 => chaotic regime
%   - Connection density = 1 (fully connected, alpha = 1)
%   - No external stimulus, random initial conditions
%
% Reference: Sompolinsky, H., Crisanti, A. & Sommers, H. J. Chaos in Random Neural Networks. Phys. Rev. Lett. 61, 259–262 (1988).

%% Setup
setup_paths;

% Figure saving configuration
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
% This script lives in scripts/tests/, so the project root is two directories
% up from its folder.
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figs_root = fullfile(project_root, 'figs');
output_folder_name = 'Sompolinsky1988';

%% Network parameters
N = 1000;               % Network size
gamma = 1.8;            % Level of chaos (Sompolinsky's g parameter)
% gamma > 1 => chaotic, gamma < 1 => stable fixed point

%% Create and configure SRNNModel
model = SRNNModel2();

% Network architecture
model.n = N;
model.f = 1;                       % All excitatory (E/I doesn't matter for Hopfield)
model.indegree = N;                % Fully connected (alpha = 1)
model.tau_d = 1;                 % Dendritic time constant (100 ms)

% RMT parameters for zero-mean, 1/sqrt(N) std Gaussian weights
% With alpha = 1 (full connectivity), the sparse parameters become dense ones
% Set mu_tilde = 0 and sigma_tilde = 1/sqrt(N) for both E and I
model.mu_E_tilde = 0;              % Zero mean excitatory
model.mu_I_tilde = 0;              % Zero mean inhibitory
model.sigma_E_tilde = 1/sqrt(N);   % Std = 1/sqrt(N)
model.sigma_I_tilde = 1/sqrt(N);   % Std = 1/sqrt(N)
model.E_W = 0;                     % No mean offset

model.level_of_chaos = gamma;      % This scales W by gamma
model.rescale_by_abscissa = false;
model.zrs_mode = 'none';           % No row-sum centering

% No adaptation (pure Hopfield network)
model.n_a_E = 0;
model.n_a_I = 0;
model.n_b_E = 0;
model.n_b_I = 0;

% Use tanh activation function
model.activation_function = @tanhActivation;
model.activation_function_derivative = @tanhActivationDerivative;

% Simulation settings
model.T_range = [-15, 50];           % 50 seconds simulation
model.fs = 200;                    % 400 Hz sampling
model.rng_seeds = [42, 1];         % Reproducible seeds

% No external input
model.input_config.n_steps = 1;
model.input_config.step_density = 0;
model.input_config.step_density_E = 0;
model.input_config.step_density_I = 0;
model.input_config.amp = 0;
model.input_config.no_stim_pattern = true;
model.input_config.intrinsic_drive = zeros(N, 1);
model.u_ex_scale = 0;

% Lyapunov settings
model.lya_method = 'benettin';
model.lya_T_interval = [15, model.T_range(2)];  % Skip first 15s transient
model.filter_local_lya = false;
model.store_full_state = false;    % Full state not needed
model.store_decimated_state = true;
model.plot_freq = 5;               % 5 Hz plot data

%% Build and run
model.build();
model.run();

%% Extract results
lya_results = model.lya_results;

fprintf('\n=== Sompolinsky Chaos Demo Results ===\n');
fprintf('Network size N = %d\n', N);
fprintf('Gamma (level_of_chaos) = %.2f\n', gamma);
fprintf('Theoretical spectral radius R = %.3f\n', model.R);
fprintf('Largest Lyapunov Exponent (LLE) = %.4f\n', lya_results.LLE);
if lya_results.LLE > 0
    fprintf('=> Network is CHAOTIC (LLE > 0)\n');
else
    fprintf('=> Network is NOT chaotic (LLE <= 0)\n');
end

%% Plot results
figure('Position', [100, 200, 1000, 600], 'Color', 'white');

% Plot 1: Time series of x (membrane potentials)
subplot(2, 1, 1);
t_plot = model.plot_data.t;
x_plot = model.plot_data.x.E;      % All neurons in .E since f=1

% Plot a subset of neurons for clarity
n_neurons_to_plot = min(1000, N);
neuron_indices = round(linspace(1, N, n_neurons_to_plot));

plot(t_plot, x_plot(neuron_indices, :)', 'LineWidth', 0.5);

xlabel('Time (s)');
ylabel('Membrane potential x');
title(sprintf('Hopfield Network Dynamics (\\gamma = %.1f, N = %d)', gamma, N));
xlim(model.T_range);
box off;
set(gca, 'Color', 'none');

% Plot 2: Local Lyapunov exponent over time
subplot(2, 1, 2);
t_lya = lya_results.t_lya;
local_lya = lya_results.local_lya;

plot(t_lya, local_lya, 'b', 'LineWidth', 1);
hold on;
yline(0, 'k--', 'LineWidth', 1);
yline(lya_results.LLE, 'r-', 'LineWidth', 1.5);
hold off;

xlabel('Time (s)');
ylabel('Local Lyapunov Exponent');
title(sprintf('Lyapunov Exponent (LLE = %.4f)', lya_results.LLE));
xlim([t_lya(1), t_lya(end)]);
legend('Local LE', '\lambda = 0', sprintf('LLE = %.3f', lya_results.LLE), ...
    'Location', 'best');
legend boxoff;
box off;
set(gca, 'Color', 'none');

sgtitle('Sompolinsky 1988: Chaos in Random Neural Networks', 'FontSize', 14, 'FontWeight', 'bold');

%% Optional: Plot eigenvalue spectrum of W
figure('Position', [200, 300, 500, 450], 'Color', 'white');

% Eigenvalues of the LTI Jacobian: (-I + W)/tau_d
J_lti = (-eye(N) + model.W) / model.tau_d;
eigs_J = eig(J_lti);

scatter(real(eigs_J), imag(eigs_J), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;

% Theoretical circle: center at -1/tau_d, radius R/tau_d
theta = linspace(0, 2*pi, 200);
center = -1 / model.tau_d;
radius = model.R / model.tau_d;
plot(center + radius*cos(theta), radius*sin(theta), 'r--', 'LineWidth', 1.5);

% Stability line
xline(0, 'k--', 'LineWidth', 1);

hold off;

xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title(sprintf('LTI Jacobian Eigenspectrum (R = %.2f)', model.R));
axis equal;
box off;
set(gca, 'Color', 'none');
legend('Eigenvalues', sprintf('Theoretical (R=%.2f)', model.R), 'Stability boundary', ...
    'Location', 'best');
legend boxoff;

%% Save figures
if save_figs
    save_dir = fullfile(figs_root, output_folder_name);

    % Save dynamics figure
    save_name_base = sprintf('Sompolinsky_N%d_dynamics', N);
    save_some_figs_to_folder_2(save_dir, save_name_base, 1, {'fig', 'svg', 'png', 'jp2'});
    fprintf('Dynamics plot saved to %s\n', save_dir);

    % Save eigenvalue figure
    save_name_base = sprintf('Sompolinsky_N%d_eigenspectrum', N);
    save_some_figs_to_folder_2(save_dir, save_name_base, 2, {'fig', 'svg', 'png', 'jp2'});
    fprintf('Eigenspectrum plot saved to %s\n', save_dir);
end

fprintf('\nDone. Two figures generated.\n');

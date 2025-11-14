%--- SRNN_reservoir_caller.m
% Simple calling script for SRNN_reservoir.m
% Demonstrates basic usage with spike-frequency adaptation

close all
clear all % must clear all due to use of persistent variables in SRNN_reservoir.m
clc

tic

%% Set random seed for reproducibility
seed = 42;
rng(seed, 'twister');

%% Network setup
n = 10; % number of neurons

% E/I composition
EI = 0.7; % fraction of excitatory neurons

% Connectivity parameters
mean_in_out_degree = 5; % desired mean number of connections in and out
density = mean_in_out_degree / (n - 1); % each neuron can make up to n-1 connections
sparsity = 1 - density;

% Weight structure
scale = 0.5 / 0.79782; % overall scaling factor of weights
w.EE = scale * 1;      % E to E connections
w.EI = scale * 1;      % E to I connections
w.IE = scale * 1;      % I to E connections
w.II = scale * 0.5;    % I to I connections
w.selfE = 0;           % self connections of E neurons
w.selfI = 0;           % self connections of I neurons

% Generate connectivity matrix
addpath('../reference_files'); % Add path to helper functions
[M, EI_vec] = generate_M_no_iso(n, w, sparsity, EI);
EI_vec = EI_vec(:); % make it a column vector
[E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

%% Time parameters
fs = 1000; % Sampling frequency (Hz)
dt = 1 / fs;
T = [-2 5]; % Time interval (seconds)

nt = round((T(2) - T(1)) * fs) + 1; % Number of time samples
t = linspace(T(1), T(2), nt)'; % Time vector (column)

%% External input (u_ex)
u_ex = zeros(n, nt);

% Add sine wave stimulus to first neuron
stim_b0 = 0.5; % baseline
amp = 0.5; % amplitude
dur = 2; % duration of stimulus (seconds)
f_sin = 1; % frequency (Hz)

% Apply stimulus starting at t = 1 second
stim_start = 1; % seconds
stim_start_idx = find(t >= stim_start, 1, 'first');
stim_end_idx = min(stim_start_idx + fix(fs * dur) - 1, nt);
stim_duration = stim_end_idx - stim_start_idx + 1;
t_stim = t(stim_start_idx:stim_end_idx) - t(stim_start_idx);
u_ex(1, stim_start_idx:stim_end_idx) = stim_b0 + amp * sin(2 * pi * f_sin * t_stim);

% Add DC offset with ramp-up to avoid initial transients
DC = 0.1;
ramp_duration = 1.5; % seconds (ramp up during negative time)
ramp_end_time = T(1) + ramp_duration;
ramp_logical = (t <= ramp_end_time);
ramp_profile = linspace(0, DC, sum(ramp_logical));
u_dc_profile = ones(1, nt) * DC;
u_dc_profile(ramp_logical) = ramp_profile;
u_ex = u_ex + u_dc_profile;

%% Adaptation parameters
tau_d = 0.025; % dendritic time constant (seconds)

% Number of adaptation time constants
n_a_E = 3; % number of SFA timescales for E neurons
n_a_I = 0; % number of SFA timescales for I neurons (typically 0)

% Define adaptation time constants
if n_a_E > 0
    tau_a_E = logspace(log10(0.3), log10(15), n_a_E); % seconds, 1 x n_a_E
else
    tau_a_E = [];
end
if n_a_I > 0
    tau_a_I = logspace(log10(0.3), log10(15), n_a_I); % seconds, 1 x n_a_I
else
    tau_a_I = [];
end

%% Package parameters
params.n = n;
params.n_E = n_E;
params.n_I = n_I;
params.E_indices = E_indices;
params.I_indices = I_indices;
params.n_a_E = n_a_E;
params.n_a_I = n_a_I;
params.W = M; % connection matrix
params.tau_d = tau_d;
params.tau_a_E = tau_a_E;
params.tau_a_I = tau_a_I;
params.activation_function = @(x) tanh(x);

%% Initial conditions
% Initialize adaptation variables for E neurons
a0_E = [];
if params.n_E > 0 && params.n_a_E > 0
    a0_E = zeros(params.n_E * params.n_a_E, 1);
end

% Initialize adaptation variables for I neurons
a0_I = [];
if params.n_I > 0 && params.n_a_I > 0
    a0_I = zeros(params.n_I * params.n_a_I, 1);
end

% Initialize dendritic states
x0 = zeros(n, 1);

% Combine into initial state vector
S0 = [a0_E; a0_I; x0];

N_sys_eqs = size(S0, 1); % Number of system equations / states

%% Integrate with ODE45
fprintf('Integrating SRNN with ode45...\n');

% Create wrapper function to pass extra parameters
SRNN_wrapper = @(tt, SS) SRNN_reservoir(tt, SS, t, u_ex, params);

% ODE options
ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Solve
[t_out, S_out] = ode45(SRNN_wrapper, t, S0, ode_options);

fprintf('Integration complete.\n');

%% Post-processing: Unpack state variables
% State organization in SRNN_reservoir: S = [a_E(:); a_I(:); x(:)]
% We need to unpack manually since the reference unpack_SRNN_state expects b_E and b_I

current_idx = 0;

% Adaptation states for E neurons (a_E)
len_a_E = n_E * n_a_E;
if len_a_E > 0
    a_E_ts = reshape(S_out(:, current_idx + (1:len_a_E))', n_E, n_a_E, length(t_out));
else
    a_E_ts = [];
end
current_idx = current_idx + len_a_E;

% Adaptation states for I neurons (a_I)
len_a_I = n_I * n_a_I;
if len_a_I > 0
    a_I_ts = reshape(S_out(:, current_idx + (1:len_a_I))', n_I, n_a_I, length(t_out));
else
    a_I_ts = [];
end
current_idx = current_idx + len_a_I;

% Dendritic states (x)
x_ts = S_out(:, current_idx + (1:n))'; % n x nt

%% Compute firing rates
% Since SRNN_reservoir doesn't have synaptic depression, we pass empty arrays for b_E and b_I
b_E_ts = [];
b_I_ts = [];

% The reference compute_dependent_variables expects u_d_ts, which is the dendritic state x
[r_ts, p_ts] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, x_ts, params);

%% Visualization
fprintf('Creating plots...\n');

figure('Position', [100, 100, 1200, 800]);

% Subplot 1: External input for selected neurons
subplot(3, 1, 1);
plot(t, u_ex(1, :), 'LineWidth', 1.5);
hold on;
if n >= 2
    plot(t, u_ex(2, :), 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('External Input');
title('External Input u_{ex}');
legend('Neuron 1', 'Neuron 2');
grid on;

% Subplot 2: Firing rates
subplot(3, 1, 2);
plot(t_out, r_ts', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Firing Rate');
title('Firing Rates r(t) for All Neurons');
grid on;
if n <= 10
    legend(arrayfun(@(i) sprintf('Neuron %d', i), 1:n, 'UniformOutput', false), ...
           'Location', 'eastoutside');
end

% Subplot 3: Adaptation variables for first E neuron
subplot(3, 1, 3);
if n_E > 0 && n_a_E > 0
    for k = 1:n_a_E
        plot(t_out, squeeze(a_E_ts(1, k, :)), 'LineWidth', 1.5);
        hold on;
    end
    xlabel('Time (s)');
    ylabel('Adaptation Variable');
    title(sprintf('Adaptation Variables a_{E} for Neuron %d (First E Neuron)', E_indices(1)));
    legend(arrayfun(@(k) sprintf('a_{%d} (\\tau=%.2fs)', k, tau_a_E(k)), ...
                    1:n_a_E, 'UniformOutput', false));
    grid on;
else
    text(0.5, 0.5, 'No adaptation variables for E neurons', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
end

sim_dur = toc;
fprintf('Simulation completed in %.2f seconds.\n', sim_dur);
fprintf('Simulation duration vs real time: %.2fx\n', sim_dur / (T(2) - T(1)));


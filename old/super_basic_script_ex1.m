close all
clear all
clc
rng(1)
% Setup parameters
params.n = 100;
f = 0.5; % fraction of neurons that are E
params.n_E = round(f*params.n);
params.n_I = params.n-params.n_E;
params.E_indices = 1:params.n_E;
params.I_indices = params.n_E+1:params.n;
params.n_a_E = 1;  % Two adaptation time constants for E neurons
params.n_a_I = 0;  % No adaptation for I neurons
mu_E = 1;
mu_I = -f*mu_E/(1-f);
M = [mu_E.*ones(params.n,params.n/2), mu_I.*ones(params.n,params.n/2)];
b_stdev = 1;
G = b_stdev*randn(params.n,params.n);
W = M+G;
d = 0.3; % density
Z = rand(params.n,params.n)>d;
W(Z) = 0;
W = W-mean(W,2);
params.W = W;
spectral_radius = b_stdev * sqrt(params.n * d);  % Random matrix theory: ρ ≈ σ√(Np) for sparse random matrix
level_of_chaos = 2;
params.tau_d = level_of_chaos/spectral_radius;  % 10 ms
% params.tau_a_E = logspace(log10(0.1), log10(10), params.n_a_E);  % Logarithmically spaced from 0.1 to 10
% params.tau_a_I = logspace(log10(0.1), log10(10), params.n_a_I);  % Logarithmically spaced from 0.1 to 10
params.tau_a_E = logspace(log10(2), log10(2), params.n_a_E);  % Logarithmically spaced from 0.1 to 10
params.tau_a_I = logspace(log10(5), log10(5), params.n_a_I);  % Logarithmically spaced from 0.1 to 10
params.activation_function = @(x) tanh(x);
params.activation_function_derivative = @(x) 1 - tanh(x).^2;
% params.activation_function = @(x) max(0,x);

% Initial conditions
N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n;
S0 = randn(N_sys_eqs, 1) * 0.01;  % Small random initial state

% External input
fs = 1000;  % Sampling frequency (Hz)
dt = 1/fs;
T = 60.0;    % Duration (s)
t_ex = (0:dt:T)';
nt = length(t_ex);

% Create sine and square wave stimulus similar to SRNN_basic_example.m
u_ex = zeros(params.n, nt);
stim_b0 = 0.5; 
amp = 0.5;
dur = 5; % duration of sine wave (shorter to fit in 1s window)
f_sin = 1.*ones(1,fs*dur);

% Square wave stimulus (starting at 0.1s)
start_idx_1 = fix(fs*nt/4);
stim_len_1 = min(fix(fs*dur), nt - start_idx_1);
u_ex(1, start_idx_1 + (1:stim_len_1)) = stim_b0 + amp.*sign(sin(2*pi*f_sin(1:stim_len_1).*t_ex(1:stim_len_1)'));

% Sine wave stimulus (starting at 0.5s)
start_idx_2 = fix(nt/2);
stim_len_2 = min(fix(fs*dur), nt - start_idx_2);
u_ex(1, start_idx_2 + (1:stim_len_2)) = stim_b0 + amp.*-cos(2*pi*f_sin(1:stim_len_2).*t_ex(1:stim_len_2)');

u_ex = u_ex(:,1:nt);

% Integrate
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', dt);
[t_out, S_out] = ode45(rhs, t_ex, S0, opts);

%% Unpack state vector and compute firing rates
% S_out is nt x N_sys_eqs, we need to unpack it
nt = size(S_out, 1);
current_idx = 0;

% Unpack adaptation states for E neurons (a_E)
len_a_E = params.n_E * params.n_a_E;
if len_a_E > 0
    a_E_ts = reshape(S_out(:, current_idx + (1:len_a_E))', params.n_E, params.n_a_E, nt);
else
    a_E_ts = [];
end
current_idx = current_idx + len_a_E;

% Unpack adaptation states for I neurons (a_I)
len_a_I = params.n_I * params.n_a_I;
if len_a_I > 0
    a_I_ts = reshape(S_out(:, current_idx + (1:len_a_I))', params.n_I, params.n_a_I, nt);
else
    a_I_ts = [];
end
current_idx = current_idx + len_a_I;

% Unpack dendritic states (x)
x_ts = S_out(:, current_idx + (1:params.n))';  % n x nt

% Compute firing rates
x_eff_ts = x_ts;  % n x nt

% Apply adaptation effect to E neurons
if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
    % sum(a_E_ts, 2) is n_E x 1 x nt, need to squeeze to n_E x nt
    sum_a_E = squeeze(sum(a_E_ts, 2));  % n_E x nt
    if size(sum_a_E, 1) ~= params.n_E  % Handle case where nt=1
        sum_a_E = sum_a_E';
    end
    x_eff_ts(params.E_indices, :) = x_eff_ts(params.E_indices, :) - sum_a_E;
end

% Apply adaptation effect to I neurons (though n_a_I = 0 in this example)
if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
    sum_a_I = squeeze(sum(a_I_ts, 2));  % n_I x nt
    if size(sum_a_I, 1) ~= params.n_I
        sum_a_I = sum_a_I';
    end
    x_eff_ts(params.I_indices, :) = x_eff_ts(params.I_indices, :) - sum_a_I;
end

r_ts = params.activation_function(x_eff_ts);  % n x nt

%% Plotting
figure('Position', [-1847         378        1200         800]);

% Plot all neurons
neuron_indices = 1:params.n;

% Subplot 1: External input
s(1) = subplot(4, 1, 1);
plot(t_out, u_ex(neuron_indices, :)');
xlabel('Time (s)');
ylabel('External Input');
title(sprintf('External Input u(t) - All %d Neurons', params.n));
grid on;

% Subplot 2: Dendritic states
s(2) = subplot(4, 1, 2);
plot(t_out, x_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Dendritic State x');
title(sprintf('Dendritic States x(t) - All %d Neurons', params.n));
grid on;

% Subplot 3: Adaptation variables (for all E neurons)
s(3) = subplot(4, 1, 3);
if ~isempty(a_E_ts)
    % Plot adaptation for all E neurons
    E_neurons_to_plot = neuron_indices(neuron_indices <= params.n_E);
    for n_idx = E_neurons_to_plot
        for k = 1:params.n_a_E
            plot(t_out, squeeze(a_E_ts(n_idx, k, :)));
            hold on;
        end
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Adaptation a');
    title(sprintf('Adaptation Variables for All %d E Neurons', params.n_E));
    grid on;
else
    text(0.5, 0.5, 'No adaptation variables', 'HorizontalAlignment', 'center');
    axis off;
end

% Subplot 4: Firing rates
s(4) = subplot(4, 1, 4);
plot(t_out, r_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Firing Rate r');
title(sprintf('Firing Rates r(t) - All %d Neurons', params.n));
grid on;

% Link x-axes of all time series plots
linkaxes(s,'x');

%% Compute Jacobian eigenvalues at multiple time points
% Select time indices every 0.05 seconds (50 ms)
dt_jacobian = fix(T/10); % time interval for Jacobian computation (seconds)
J_times = round(linspace(1, nt, floor(T/dt_jacobian) + 1));
J_times = unique(J_times); % ensure unique indices

fprintf('Computing Jacobian at %d time points...\n', length(J_times));
J_array = compute_Jacobian_at_indices(S_out, J_times, params);

% Compute eigenvalues for each Jacobian
eigenvalues_all = cell(length(J_times), 1);
for i = 1:length(J_times)
    eigenvalues_all{i} = eig(J_array(:,:,i));
end

%% Plot eigenvalue distributions on complex plane
n_plots = length(J_times);
n_cols = ceil(sqrt(n_plots));
n_rows = ceil(n_plots / n_cols);

figure('Position', [-1847, -500, 1400, 1000]);
ax_handles = zeros(n_plots, 1);

for i = 1:n_plots
    ax_handles(i) = subplot(n_rows, n_cols, i);
    
    % Get eigenvalues for this time point
    evals = eigenvalues_all{i};
    
    % Scatter plot on complex plane
    scatter(real(evals), imag(evals), 30, 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Add unit circle for reference (stability boundary for continuous systems at Re=0)
    hold on;
    plot([0 0], ylim, 'k--', 'LineWidth', 0.5); % imaginary axis
    plot(xlim, [0 0], 'k--', 'LineWidth', 0.5); % real axis
    hold off;
    
    % Labels and title
    xlabel('Real');
    ylabel('Imaginary');
    title(sprintf('t = %.3f s (idx %d)', t_out(J_times(i)), J_times(i)));
    grid on;
    axis equal;
end

% Link axes of all subplots
linkaxes(ax_handles, 'xy');

% Add overall title
sgtitle('Eigenvalue Distribution on Complex Plane');

fprintf('Jacobian analysis complete.\n');
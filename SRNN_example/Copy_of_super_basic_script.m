close all
clear all
clc

% Setup parameters
params.n = 100;
params.n_E = 80;
params.n_I = 20;
params.E_indices = 1:80;
params.I_indices = 81:100;
params.n_a_E = 2;  % Two adaptation time constants for E neurons
params.n_a_I = 0;  % No adaptation for I neurons
params.W = randn(100, 100) * 0.1;  % Example connectivity
params.tau_d = 0.01;  % 10 ms
params.tau_a_E = [0.1, 1.0];  % Fast and slow adaptation
params.tau_a_I = [];  % Empty since n_a_I = 0
% params.activation_function = @(x) tanh(x);
params.activation_function = @(x) max(0,x);

% Initial conditions
N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n;
S0 = randn(N_sys_eqs, 1) * 0.01;  % Small random initial state

% External input
fs = 1000;  % Sampling frequency (Hz)
T = 50;    % Duration (s)
t_ex = linspace(0, T, T*fs+1)';
nt = length(t_ex);

% Create sine and square wave stimulus similar to SRNN_basic_example.m
u_ex = zeros(params.n, nt);

amp = 0.5;
stim_b0 = amp; 
dur = 10; % duration of sine wave (shorter to fit in 1s window)
f_sin = 1.*ones(1,fs*dur);
dc = 0.1;

% Square wave stimulus (starting at 0.1s)
start_idx_1 = fix(fs*0.1);
stim_len_1 = min(fix(fs*dur), nt - start_idx_1);
u_ex(1, start_idx_1 + (1:stim_len_1)) = stim_b0 + amp.*sign(sin(2*pi*f_sin(1:stim_len_1).*t_ex(1:stim_len_1)'));

% Sine wave stimulus (starting at 0.5s)
start_idx_2 = fix(fs*0.5);
stim_len_2 = min(fix(fs*dur), nt - start_idx_2);
u_ex(1, start_idx_2 + (1:stim_len_2)) = stim_b0 + amp.*-cos(2*pi*f_sin(1:stim_len_2).*t_ex(1:stim_len_2)');

u_ex = u_ex(:,1:nt);

u_ex = u_ex + dc;

% Integrate
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
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
figure(1)
set(1,'Position', [-1853         530        1200         800]);

% Select neurons to plot
neurons_to_plot = min(5, params.n);  % Plot up to 5 neurons
neuron_indices = 1:neurons_to_plot;

% Subplot 1: External input
subplot(4, 1, 1);
plot(t_out, u_ex(neuron_indices, :)');
xlabel('Time (s)');
ylabel('External Input');
title('External Input u(t)');
legend(arrayfun(@(i) sprintf('Neuron %d', i), neuron_indices, 'UniformOutput', false), ...
       'Location', 'eastoutside');
grid on;

% Subplot 2: Dendritic states
subplot(4, 1, 2);
plot(t_out, x_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Dendritic State x');
title('Dendritic States x(t)');
legend(arrayfun(@(i) sprintf('Neuron %d', i), neuron_indices, 'UniformOutput', false), ...
       'Location', 'eastoutside');
grid on;

% Subplot 3: Adaptation variables (for first E neuron)
subplot(4, 1, 3);
if ~isempty(a_E_ts)
    for k = 1:params.n_a_E
        plot(t_out, squeeze(a_E_ts(1, k, :)), 'DisplayName', ...
             sprintf('\\tau_a = %.2f s', params.tau_a_E(k)));
        hold on;
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Adaptation a');
    title('Adaptation Variables for E Neuron 1');
    legend('Location', 'eastoutside');
    grid on;
else
    text(0.5, 0.5, 'No adaptation variables', 'HorizontalAlignment', 'center');
    axis off;
end

% Subplot 4: Firing rates
subplot(4, 1, 4);
plot(t_out, r_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Firing Rate r');
title('Firing Rates r(t)');
legend(arrayfun(@(i) sprintf('Neuron %d', i), neuron_indices, 'UniformOutput', false), ...
       'Location', 'eastoutside');
grid on;
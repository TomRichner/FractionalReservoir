close all
clear all
clc
rng(1)

% Lyapunov method selection
Lya_method = 'qr'; % 'benettin', 'qr', or 'none'

% Setup parameters
params.n = 100;
f = 0.5; % fraction of neurons that are E
params.n_E = round(f*params.n);
params.n_I = params.n-params.n_E;
params.E_indices = 1:params.n_E;
params.I_indices = params.n_E+1:params.n;
params.n_a_E = 0;  % Two adaptation time constants for E neurons
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
params.tau_a_E = logspace(log10(0.1), log10(10), params.n_a_E);  % Logarithmically spaced from 0.1 to 10
params.tau_a_I = logspace(log10(0.1), log10(10), params.n_a_I);  % Logarithmically spaced from 0.1 to 10
params.c_E = .2;  % Adaptation scaling for E neurons (scalar, typically 0-3)
params.c_I = .1;  % Adaptation scaling for I neurons (scalar, typically 0-3)
% params.activation_function = @(x) tanh(x);
% params.activation_function_derivative = @(x) 1 - tanh(x).^2;
% params.activation_function = @(x) min(max(0,x),1);
% params.activation_function_derivative = @(x) double(and(0<=x, x<=1));
params.activation_function = @(x) 1./(1 + exp(-4*x));
params.activation_function_derivative = @(x) 4*params.activation_function(x).*(1 - params.activation_function(x));

% Initial conditions
N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n;
% S0 = 0.1+randn(N_sys_eqs, 1) * 0.1;  % Small random initial state
S0 = 0.5*ones(N_sys_eqs,1);

% External input
rng(2);  % Fresh seed for u_ex to keep it independent of S0 size
fs = 100;  % Sampling frequency (Hz)
dt = 1/fs;
T = 100.0;    % Duration (s)
t_ex = (0:dt:T)';
nt = length(t_ex);

% Create random amplitude multidimensional sparse step function
% step_period = 20;  % Time period for each step (seconds)
step_density = 0.2;  % Fraction of neurons receiving input at each step
amp = 1;  % Amplitude scaling factor

% Calculate number of steps
% n_steps = ceil(T / step_period);
n_steps = 7;
step_period = fix(T/n_steps);

step_length = round(step_period * fs);  % Number of time points per step

% Precompute random sparse step matrix (vectorized)
% Generate all random amplitudes at once: n x n_steps
random_sparse_step = amp * randn(params.n, n_steps);

% Make it sparse based on step_density
sparse_mask = rand(params.n, n_steps) < step_density;
random_sparse_step = random_sparse_step .* sparse_mask;

% Option to zero out entire columns (steps with no stimulation)
no_stim_steps = false(1, n_steps);  % Logical vector: true = no stimulation for that step
no_stim_steps(1:2:end) = true;
% Example: no_stim_steps(1:10) = true;  % No stimulation for first 10 steps
random_sparse_step(:, no_stim_steps) = 0;

% Create u_ex using the precomputed random sparse steps
u_ex = zeros(params.n, nt);
for step_idx = 1:n_steps
    % Determine time indices for this step
    start_idx = (step_idx - 1) * step_length + 1;
    end_idx = min(step_idx * step_length, nt);
    
    if start_idx > nt
        break;
    end
    
    % Apply the step to all neurons (broadcasting the column vector)
    u_ex(:, start_idx:end_idx) = repmat(random_sparse_step(:, step_idx), 1, end_idx - start_idx + 1);
end

% add small intrinsic drive to neurons in u_ex
intrinsic_drive = -0.9+0.2*randn(params.n,1);
u_ex = u_ex+intrinsic_drive;

% Integrate
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', dt);
[t_out, S_out] = ode45(rhs, t_ex, S0, opts);

%% Compute Lyapunov exponent(s)
ode_solver = @ode45;
T_interval = [0, T];
lya_results = struct();

if ~strcmpi(Lya_method, 'none')
    % Adjust lya_dt based on method: QR needs longer interval
    if strcmpi(Lya_method, 'qr')
        lya_dt = 0.1;  % Longer interval for QR method
    else
        lya_dt = 5*dt;   % Standard interval for Benettin
    end
    lya_fs = 1/lya_dt;
    
    switch lower(Lya_method)
        case 'benettin'
            fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
            d0 = 1e-3;
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(S_out, t_out, dt, fs, d0, T_interval, lya_dt, params, opts, @SRNN_reservoir, t_ex, u_ex, ode_solver);
            fprintf('Largest Lyapunov Exponent: %.4f\n', LLE);
            lya_results.LLE = LLE; 
            lya_results.local_lya = local_lya; 
            lya_results.finite_lya = finite_lya; 
            lya_results.t_lya = t_lya;
            
        case 'qr'
            fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
            [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(S_out, t_out, lya_dt, params, ode_solver, opts, @SRNN_Jacobian_wrapper, T_interval, N_sys_eqs, fs);
            fprintf('Largest Lyapunov Exponent: %.4f\n', LE_spectrum(1));
            fprintf('Lyapunov Dimension: %.2f\n', compute_kaplan_yorke_dimension(LE_spectrum));
            lya_results.LE_spectrum = LE_spectrum; 
            lya_results.local_LE_spectrum_t = local_LE_spectrum_t; 
            lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t; 
            lya_results.t_lya = t_lya; 
            lya_results.N_sys_eqs = N_sys_eqs;
            
        otherwise
            error('Unknown Lyapunov method: %s. Use ''benettin'', ''qr'', or ''none''.', Lya_method);
    end
end

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

% Apply adaptation effect to E neurons (scaled by c_E)
if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
    % sum(a_E_ts, 2) is n_E x 1 x nt, need to squeeze to n_E x nt
    sum_a_E = squeeze(sum(a_E_ts, 2));  % n_E x nt
    if size(sum_a_E, 1) ~= params.n_E  % Handle case where nt=1
        sum_a_E = sum_a_E';
    end
    x_eff_ts(params.E_indices, :) = x_eff_ts(params.E_indices, :) - params.c_E * sum_a_E;
end

% Apply adaptation effect to I neurons (scaled by c_I)
if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
    sum_a_I = squeeze(sum(a_I_ts, 2));  % n_I x nt
    if size(sum_a_I, 1) ~= params.n_I
        sum_a_I = sum_a_I';
    end
    x_eff_ts(params.I_indices, :) = x_eff_ts(params.I_indices, :) - params.c_I * sum_a_I;
end

r_ts = params.activation_function(x_eff_ts);  % n x nt

%% Plotting
figure('Position', [-1847         378        1200         800]);

% Plot all neurons
neuron_indices = 1:params.n;

% Subplot 1: External input
s(1) = subplot(5, 1, 1);
plot(t_out, u_ex(neuron_indices, :)');
xlabel('Time (s)');
ylabel('External Input');

% Subplot 2: Dendritic states
s(2) = subplot(5, 1, 2);
plot(t_out, x_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Dendritic State x');

% Subplot 3: Adaptation variables (for all E neurons)
s(3) = subplot(5, 1, 3);
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

else
    text(0.5, 0.5, 'No adaptation variables', 'HorizontalAlignment', 'center');
    axis off;
end

% Subplot 4: Firing rates
s(4) = subplot(5, 1, 4);
plot(t_out, r_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Firing Rate r');

% Subplot 5: Lyapunov exponent(s)
if ~strcmpi(Lya_method, 'none')
    s(5) = subplot(5, 1, 5);
    
    if strcmpi(Lya_method, 'benettin')
        % Plot for Benettin's method (single exponent)
        [bL,aL] = butter(3,0.05/(lya_fs/2),'low');
        local_lya_filt = filtfilt(bL,aL,lya_results.local_lya);
        
        plot(lya_results.t_lya, lya_results.local_lya, 'LineWidth', 1.5)
        hold on
        plot(lya_results.t_lya, local_lya_filt,'LineWidth',1.5)
        yline(lya_results.LLE, '--r', 'LineWidth', 2)
        xlabel('Time (s)')
        ylabel('Local Lyapunov Exponent')
        title(sprintf('Local Lyapunov Exponent (LLE = %.4f)', lya_results.LLE))
        legend('Local', 'Filtered', 'Mean LLE', 'Location', 'best')
        
    elseif strcmpi(Lya_method, 'qr')
        % Plot for QR method (full spectrum)
        plot(lya_results.t_lya, lya_results.local_LE_spectrum_t, 'LineWidth', 1)
        hold on
        % Highlight the largest exponent
        plot(lya_results.t_lya, lya_results.local_LE_spectrum_t(:,1), 'r', 'LineWidth', 1)
        yline(0, '--k', 'LineWidth', 1)
        xlabel('Time (s)')
        ylabel('Lyapunov Exponents')
        
        % Add legend with final values
        legend_entries = cell(1, min(5, lya_results.N_sys_eqs));
        for i = 1:min(5, lya_results.N_sys_eqs)
            legend_entries{i} = sprintf('\\lambda_%d = %.3f', i, lya_results.LE_spectrum(i));
        end
        legend(legend_entries, 'Location', 'best')
    end
end

% Link x-axes of all time series plots
linkaxes(s,'x');

%% Compute Jacobian eigenvalues at multiple time points
% Sample at t=0 and 1 second before/after each step change
step_change_times = (0:step_period:(n_steps-1)*step_period);  % Times when steps change (0, 20, 40, ...)
J_times_sec = [0];  % Start with t=0

% Add 1 second before and after each step change (except t=0)
for i = 2:length(step_change_times)
    t_change = step_change_times(i);
    if t_change - 1 > 0  % Make sure we don't go below 0
        J_times_sec = [J_times_sec, t_change - 1];
    end
    if t_change + 1 <= T  % Make sure we don't exceed T
        J_times_sec = [J_times_sec, t_change + 1];
    end
end

% Convert times in seconds to indices (MATLAB is 1-indexed)
J_times = round(J_times_sec * fs) + 1;
J_times = unique(J_times);  % Ensure unique indices and sort

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
    title(sprintf('t = %.3f s', t_out(J_times(i))));
    axis equal;
end

% Link axes of all subplots
linkaxes(ax_handles, 'xy');

% Add overall title
sgtitle('Eigenvalue Distribution on Complex Plane');

fprintf('Jacobian analysis complete.\n');

%% Jacobian wrapper for lyapunov_spectrum_qr
function J_jac = SRNN_Jacobian_wrapper(tt, S, params)
    J_jac = compute_Jacobian(S, params);
end

%% Helper function to compute Kaplan-Yorke dimension
function D_KY = compute_kaplan_yorke_dimension(lambda)
    % Compute Kaplan-Yorke (Lyapunov) dimension from spectrum
    % lambda: sorted Lyapunov exponents (descending order)
    
    lambda = sort(lambda, 'descend');
    cumsum_lambda = cumsum(lambda);
    
    % Find largest j such that sum of first j exponents is non-negative
    j = find(cumsum_lambda >= 0, 1, 'last');
    
    if isempty(j)
        D_KY = 0;
    elseif j == length(lambda)
        D_KY = length(lambda);
    else
        D_KY = j + cumsum_lambda(j) / abs(lambda(j+1));
    end
end

close all
clear all
clc
rng(1)

% Lyapunov method selection
Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'

% Setup parameters
params.n = 200;
f = 0.5; % fraction of neurons that are E
params.n_E = round(f*params.n);
params.n_I = params.n-params.n_E;
params.E_indices = 1:params.n_E;
params.I_indices = params.n_E+1:params.n;

mu_E = 1;
mu_I = -f*mu_E/(1-f);
M = [mu_E.*ones(params.n,params.n/2), mu_I.*ones(params.n,params.n/2)];
b_stdev = 1;
G = b_stdev*randn(params.n,params.n);
W = M+G;
d = 15/params.n; % density
Z = rand(params.n,params.n)>d;
W(Z) = 0;
nonzero_mask = ~Z;
row_counts = sum(nonzero_mask, 2);
row_sums = sum(W, 2);
row_means = zeros(size(row_sums));
valid_rows = row_counts > 0;
row_means(valid_rows) = row_sums(valid_rows) ./ row_counts(valid_rows);
W = W - bsxfun(@times, row_means, nonzero_mask);
params.W = W;
spectral_radius = b_stdev * sqrt(params.n * d);  % Random matrix theory: ρ ≈ σ√(Np) for sparse random matrix
level_of_chaos = 1.7;
params.tau_d = level_of_chaos/spectral_radius;  % 10 ms

params.n_a_E = 0;  % Two adaptation time constants for E neurons
params.n_a_I = 0;  % No adaptation for I neurons
params.tau_a_E = logspace(log10(0.1), log10(10), params.n_a_E);  % Logarithmically spaced from 0.1 to 10
params.tau_a_I = logspace(log10(0.1), log10(10), params.n_a_I);  % Logarithmically spaced from 0.1 to 10
params.c_E = .2;  % Adaptation scaling for E neurons (scalar, typically 0-3)
params.c_I = .1;  % Adaptation scaling for I neurons (scalar, typically 0-3)

params.n_b_E = 0;  % Number of STD timescales for E neurons (0 or 1)
params.n_b_I = 0;  % Number of STD timescales for I neurons (0 or 1)
params.tau_b_E_rec = 2;  % STD recovery time constant for E neurons (s)
params.tau_b_E_rel = 0.5;  % STD release time constant for E neurons (s)
params.tau_b_I_rec = 0.5;  % STD recovery time constant for I neurons (s)
params.tau_b_I_rel = 2;  % STD release time constant for I neurons (s)
% params.activation_function = @(x) tanh(x);
% params.activation_function_derivative = @(x) 1 - tanh(x).^2;
% params.activation_function = @(x) min(max(0,x),1);
% params.activation_function_derivative = @(x) double(and(0<=x, x<=1));
params.activation_function = @(x) 1./(1 + exp(-4*x));
params.activation_function_derivative = @(x) 4*params.activation_function(x).*(1 - params.activation_function(x));

% Initial conditions
N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n_E * params.n_b_E + params.n_I * params.n_b_I + params.n;

% Initialize adaptation states
a0_E = [];
if params.n_a_E > 0
    a0_E = zeros(params.n_E * params.n_a_E, 1);
end

a0_I = [];
if params.n_a_I > 0
    a0_I = zeros(params.n_I * params.n_a_I, 1);
end

% Initialize STD states (b variables start at 1.0, no depression)
b0_E = [];
if params.n_b_E > 0
    b0_E = ones(params.n_E * params.n_b_E, 1);
end

b0_I = [];
if params.n_b_I > 0
    b0_I = ones(params.n_I * params.n_b_I, 1);
end

% Initialize dendritic states
x0 = 0.5*ones(params.n, 1);

% Pack state vector: [a_E; a_I; b_E; b_I; x]
S0 = [a0_E; a0_I; b0_E; b0_I; x0];

% Sanity check: compare dense and sparse Jacobians at initial state
J_ref = compute_Jacobian(S0, params);
J_fast = compute_Jacobian_fast(S0, params);
jac_diff_abs = norm(J_ref - full(J_fast), 'fro');
jac_diff_rel = jac_diff_abs / max(1, norm(J_ref, 'fro'));
fprintf('Jacobian check (fast vs. reference): abs=%.3e, rel=%.3e\n', jac_diff_abs, jac_diff_rel);
if jac_diff_rel > 1e-8
    error('Jacobian mismatch exceeds tolerance (rel=%.3e)', jac_diff_rel);
end

% External input
rng(2);  % Fresh seed for u_ex to keep it independent of S0 size
fs = 100;  % Sampling frequency (Hz)
dt = 1/fs;
T = 200.0;    % Duration (s)
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
            fprintf('Lyapunov Dimension: %.2f\n', compute_kaplan_yorke_dimension(LE_spectrum));
            lya_results.LE_spectrum = LE_spectrum; 
            lya_results.local_LE_spectrum_t = local_LE_spectrum_t; 
            lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t; 
            lya_results.t_lya = t_lya; 
            lya_results.N_sys_eqs = N_sys_eqs;

            % sort LE spectra by real component (descending) and keep index map
            [sorted_LE, sort_idx] = sort(real(lya_results.LE_spectrum), 'descend');
            lya_results.LE_spectrum = sorted_LE;
            lya_results.local_LE_spectrum_t = lya_results.local_LE_spectrum_t(:, sort_idx);
            lya_results.finite_LE_spectrum_t = lya_results.finite_LE_spectrum_t(:, sort_idx);
            lya_results.sort_idx = sort_idx;
            fprintf('Largest Lyapunov Exponent (sorted): %.4f\n', lya_results.LE_spectrum(1));
        otherwise
            error('Unknown Lyapunov method: %s. Use ''benettin'', ''qr'', or ''none''.', Lya_method);
    end
end

%% Unpack state vector and compute firing rates
% S_out is nt x N_sys_eqs, we need to unpack it
% State organization: [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
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

% Unpack b states for E neurons (b_E)
len_b_E = params.n_E * params.n_b_E;
if len_b_E > 0
    b_E_ts = S_out(:, current_idx + (1:len_b_E))';  % n_E x nt
else
    b_E_ts = [];
end
current_idx = current_idx + len_b_E;

% Unpack b states for I neurons (b_I)
len_b_I = params.n_I * params.n_b_I;
if len_b_I > 0
    b_I_ts = S_out(:, current_idx + (1:len_b_I))';  % n_I x nt
else
    b_I_ts = [];
end
current_idx = current_idx + len_b_I;

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

% Apply STD effect (b multiplicative factor)
b_ts = ones(params.n, nt);  % Initialize b = 1 for all neurons (no depression)
if params.n_b_E > 0 && ~isempty(b_E_ts)
    b_ts(params.E_indices, :) = b_E_ts;
end
if params.n_b_I > 0 && ~isempty(b_I_ts)
    b_ts(params.I_indices, :) = b_I_ts;
end

r_ts = params.activation_function(x_eff_ts);  % n x nt (without b modulation)

%% Plotting
figure('Position', [-1847         378        1200         800]);

% Plot all neurons
neuron_indices = 1:params.n;

% Subplot 1: External input
s(1) = subplot(6, 1, 1);
plot(t_out, u_ex(neuron_indices, :)');
xlabel('Time (s)');
ylabel('External Input');

% Subplot 2: Dendritic states
s(2) = subplot(6, 1, 2);
plot(t_out, x_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Dendritic State x');

% Subplot 3: Adaptation variables (for all E neurons)
s(3) = subplot(6, 1, 3);
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
s(4) = subplot(6, 1, 4);
plot(t_out, r_ts(neuron_indices, :)');
xlabel('Time (s)');
ylabel('Firing Rate r');

% Subplot 5: STD variables (b)
s(5) = subplot(6, 1, 5);
if params.n_b_E > 0 && ~isempty(b_E_ts)
    % Plot b for E neurons
    E_neurons_to_plot = neuron_indices(neuron_indices <= params.n_E);
    plot(t_out, b_E_ts(E_neurons_to_plot, :)');
    hold on;
end
if params.n_b_I > 0 && ~isempty(b_I_ts)
    % Plot b for I neurons
    I_neurons_to_plot = neuron_indices(neuron_indices > params.n_E) - params.n_E;
    plot(t_out, b_I_ts(I_neurons_to_plot, :)');
end
if (params.n_b_E == 0 || isempty(b_E_ts)) && (params.n_b_I == 0 || isempty(b_I_ts))
    text(0.5, 0.5, 'No STD variables', 'HorizontalAlignment', 'center');
    axis off;
else
    hold off;
    xlabel('Time (s)');
    ylabel('STD Variable b');
end

% Subplot 6: Lyapunov exponent(s)
if ~strcmpi(Lya_method, 'none')
    s(6) = subplot(6, 1, 6);
    
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
        plot_data = lya_results.local_LE_spectrum_t(:, end:-1:1);
        line_handles = plot(lya_results.t_lya, plot_data, 'LineWidth', 1);
        line_handles = line_handles(end:-1:1); % reorder handles to match sorted spectrum order
        hold on
        yline(0, '--k', 'LineWidth', 1)
        xlabel('Time (s)')
        ylabel('Lyapunov Exponents')
        
        % Add legend with final values (most positive exponents on top)
        legend_count = min(5, lya_results.N_sys_eqs);
        legend_entries = cell(1, legend_count);
        for i = 1:legend_count
            legend_entries{i} = sprintf('\\lambda_{%d} = %.3f', i, lya_results.LE_spectrum(i));
        end
        legend(line_handles(1:legend_count), legend_entries, 'Location', 'best')
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
    J_jac = compute_Jacobian_fast(S, params);
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

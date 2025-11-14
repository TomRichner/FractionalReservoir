close all
clear
clc

% Sompolinsky 1988 chaotic rate network integrated with ode45
%   dx_i/dt = -x_i + sum_j J_ij * tanh(x_j)
%   J_ij ~ N(0, g^2 / n), choose g slightly > 1 for weak chaos

% rng(2); % reproducibility

% Lyapunov method selection
Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'
plot_entropy_rate = false; % true: plot entropy production rate, false: plot full spectrum

% Parameters
n = 30;      % n
g = 1;              % gain (slightly chaotic when > 1)
fs = 250;          % sampling frequency (Hz)
dt = 1/fs;
tspan = [0, 200];       % seconds
t = linspace(tspan(1), tspan(2), (tspan(2)-tspan(1))*fs+1)'; % time vector

nt = length(t);
m = 1; % negative feedback
mu_E = 1.5;
mu_I = -1;
M = [mu_E.*ones(n,n/2), mu_I.*ones(n,n/2)];
% Connectivity: Gaussian with mean 0 and variance g^2/n
J = (g / sqrt(n)) * (randn(n, n)+M);
% J(rand(n,n)<0.5) = 0; 

% J=J-mean(J,2)./n; % row sum to zero

eigJ = eig(J);
figure(1)
plot(real(eigJ),imag(eigJ),'o')

% External input stimulus (n x time)
u = zeros(n, nt); % User can customize this as needed

stim_change_interval = 200;
% Calculate number of time steps per stimulus and number of stimuli needed
steps_per_stim = round(stim_change_interval * fs);
num_stimuli = ceil(nt / steps_per_stim);

% Pre-generate all stimulus vectors
amp = 0;
stim_vectors = amp .* randn(n, num_stimuli);

% Fill u with the pre-generated stimulus vectors using a for loop
for i = 1:num_stimuli
    start_idx = (i-1) * steps_per_stim + 1;
    end_idx = min(i * steps_per_stim, nt);
    u(:, start_idx:end_idx) = repmat(stim_vectors(:, i), 1, end_idx - start_idx + 1);
end

% Initial condition (small random perturbation)
x0 = 0.05+0.1 * randn(n, 1);

% ODE right-hand side (now with external stimulus)
rhs = @(t_current, x) rhs_somp(t_current, x, t, u, J, m);

% Integrate
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'MaxStep', dt);
[t_out, X] = ode45(rhs, t, x0, opts);

%% Compute Lyapunov exponent(s)
params.n = n;
params.J = J;
params.m = m;
N_sys_eqs = n;
T = tspan;
ode_solver = @ode45;

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
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t_out, dt, fs, d0, T, lya_dt, params, opts, @rhs_somp_wrapper, t, u, ode_solver);
            fprintf('Largest Lyapunov Exponent: %.4f\n', LLE);
            lya_results.LLE = LLE; 
            lya_results.local_lya = local_lya; 
            lya_results.finite_lya = finite_lya; 
            lya_results.t_lya = t_lya;
            
        case 'qr'
            fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
            [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(X, t_out, lya_dt, params, ode_solver, opts, @rhs_somp_Jacobian, T, N_sys_eqs, fs);
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

%% Plots
% Plot: activity over time (all neurons)
figure(2);
set(2,'Color', 'w')
set(gcf,'Position',[-2302         726        1574         543])
imagesc(t_out, 1:n, X.');
axis xy tight
colormap turbo
colorbar
xlabel('Time (s)')
ylabel('Neuron index')
% title(sprintf('Sompolinsky network dynamics (n = %d, g = %.2f)', n, g))

% Also show trajectories of a subset of neurons and Lyapunov exponent(s)
if ~strcmpi(Lya_method, 'none')
    figure(3);
    clf;
    set(gcf,'Position',[-2244          82        1340         1000])

    % Subplot 1: External input stimulus
    subplot(3,1,1)
    plot(t,u')
    xlabel('Time (s)')
    ylabel('Stimulus u(t)')
    ylim(1.2*ylim);

    % Subplot 2: Sample neuron trajectories
    subplot(3,1,2)
    numToShow = n
    plot(t_out, X(:, 1:numToShow))
    xlabel('Time (s)')
    ylabel('Dendritic activity, x(t)')

    % Subplot 3: Lyapunov exponent(s)
    subplot(3,1,3)
    
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
        % Plot for QR method (full spectrum or entropy production rate)
        if plot_entropy_rate
            % Compute local entropy production rate (sum of positive Lyapunov exponents)
            % Convert from nats/s to bits/s by multiplying by 1/ln(2)
            nats_to_bits = 1/log(2);
            local_entropy_rate = sum(max(lya_results.local_LE_spectrum_t, 0), 2) * nats_to_bits;
            mean_entropy_rate = sum(max(lya_results.LE_spectrum, 0)) * nats_to_bits;
            
            % Plot entropy production rate
            plot(lya_results.t_lya, local_entropy_rate, 'b', 'LineWidth', 1.5)
            hold on
            yline(mean_entropy_rate, '--r', 'LineWidth', 2)
            yline(0, '--k', 'LineWidth', 1)
            xlabel('Time (s)')
            ylabel('Entropy Production (bits/s)')
            legend('Local Entropy Rate', sprintf('Mean = %.3f bits/s', mean_entropy_rate), 'Location', 'best')
        else
            % Plot full spectrum (original behavior)
            plot(lya_results.t_lya, lya_results.local_LE_spectrum_t, 'LineWidth', 1)
            hold on
            % Highlight the largest exponent
            plot(lya_results.t_lya, lya_results.local_LE_spectrum_t(:,1), 'r', 'LineWidth', 1)
            yline(0, '--k', 'LineWidth', 1)
            xlabel('Time (s)')
            ylabel('Lyapunov Exponents')

            % title(sprintf('Lyapunov Spectrum (LLE = %.4f)', lya_results.LE_spectrum(1)))
            legend_entries = cell(1, min(5, lya_results.N_sys_eqs));
            for i = 1:min(5, lya_results.N_sys_eqs)
                legend_entries{i} = sprintf('\\lambda_%d = %.3f', i, lya_results.LE_spectrum(i));
            end
            legend(legend_entries, 'Location', 'best')
        end
    end
end

%% Plot: Eigenvalue distribution of Jacobian at different times
figure(4);
clf;
set(gcf, 'Color', 'w')
eig_J_times = 50:50:tspan(2);
num_times = length(eig_J_times);
set(gcf, 'Position', [-2302, 100, 400*num_times, 400])

for i = 1:num_times
    % Find closest time index in t_out
    [~, time_idx] = min(abs(t_out - eig_J_times(i)));
    actual_time = t_out(time_idx);
    
    % Extract state at this time
    x_t = X(time_idx, :)';
    
    % Compute Jacobian at this state
    J_jac = rhs_somp_Jacobian(actual_time, x_t, params);
    
    % Compute eigenvalues
    eig_J_t = eig(J_jac);
    
    % Create subplot
    subplot(1, num_times, i)
    hold on
    
    % Plot unit circle
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1)
    
    % Plot axis lines at 0
    xline(0, 'k-', 'LineWidth', 0.5)
    yline(0, 'k-', 'LineWidth', 0.5)
    
    % Plot eigenvalues
    plot(real(eig_J_t), imag(eig_J_t), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1)
    
    % Format subplot
    axis equal

    xlabel('Real')
    ylabel('Imaginary')
    title(sprintf('t = %.0f s', actual_time))
    
    % Set axis limits to show unit circle nicely
    max_val = max(abs(eig_J_t));
    lim = max(1.2, max_val * 1.1);
    xlim([-lim, lim])
    ylim([-lim, lim])
end

sgtitle('Eigenvalue Distribution of Jacobian')

%% Subfunction: ODE right-hand side with external stimulus
function dx_dt = rhs_somp(t_current, x, t_vec, u, J, m)
    % Persistent variables for efficient interpolation
    persistent u_interpolant t_vec_last;
    
    % Create or update griddedInterpolant if time vector has changed
    if isempty(u_interpolant) || isempty(t_vec_last) || ...
       numel(t_vec_last) ~= numel(t_vec) || t_vec_last(1) ~= t_vec(1) || t_vec_last(end) ~= t_vec(end)
        
        % Use 'linear' interpolation and 'none' for extrapolation
        % This returns NaN for out-of-bounds queries to catch errors
        u_interpolant = griddedInterpolant(t_vec, u', 'linear', 'none');
        t_vec_last = t_vec;
    end
    
    % Interpolate u at current time (u_interpolant returns 1 x n, transpose to n x 1)
    u_t = u_interpolant(t_current)';
    
    % Compute derivative: dx/dt = -x + J * relu(x) + u(t)
    % relu with saturation at 5
    % dx_dt = 0.5-m.*x + J * min(max(x, 0), 1) + u_t;
    % dx_dt = -m.*x + J * tanh(max(x, 0)) + u_t;
    dx_dt = -m.*x + J * tanh(x) + u_t;
end

%% Wrapper function for benettin_algorithm
function dx_dt = rhs_somp_wrapper(tt, XX, t_ex, u_ex, params)
    dx_dt = rhs_somp(tt, XX, t_ex, u_ex, params.J, params.m);
end

%% Jacobian function for lyapunov_spectrum_qr
function J_jac = rhs_somp_Jacobian(tt, x, params)
    % Jacobian of dx/dt = -m*x + J*tanh(x) + u(t)
    % d/dx[tanh(x)] = sech(x)^2 = 1 - tanh(x)^2
    % J_jac = -m*I + J * diag(1 - tanh(x).^2)
    n = params.n;
    J_jac = -params.m * eye(n) + params.J * diag(1 - tanh(x).^2);
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

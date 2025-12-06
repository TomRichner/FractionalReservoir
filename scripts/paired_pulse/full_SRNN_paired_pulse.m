function [t_out, S_out, params, lya_results] = full_SRNN_paired_pulse(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note)
% FULL_SRNN_PAIRED_PULSE Runs the SRNN simulation with paired-pulse paradigm.
%
% Inputs:
%   u_ex_scale    - Scaling factor for external input
%   n_a_E         - Number of adaptation time constants for E neurons
%   n_b_E         - Number of STD timescales for E neurons
%   level_of_chaos- Level of chaos parameter
%   save_dir      - Directory path to save outputs
%   save_figs     - Boolean, whether to save figures
%   save_workspace- Boolean, whether to save the workspace
%   note          - Optional string for naming figures (default: 'SRNN_PP')

if nargin < 8 || isempty(note)
    note = 'SRNN_PP';
end

% background of all figs should be white
set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1);
set(groot, 'DefaultAxesLineWidth', 1);
set(groot, 'DefaultAxesTitleFontWeight', 'normal');

% params.rng_seeds = [5 3 3 4 5];
params.rng_seeds = [5 3 3 4 5];
rng(params.rng_seeds(1))

% Lyapunov method selection
Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'

%% time 
fs = 2000;  % Higher sampling rate for short pulses
dt = 1/fs;
T = 10;    % Duration (s)
t_ex = (0:dt:T)';
nt = length(t_ex);

%% input
% u_ex_scale passed as argument

% set simulation parameters
params.n = 50;
params.indegree = 10; % expected indegree
params.d = params.indegree/params.n; % expected density
params.f = 0.5; % fraction of neurons that are E
params.G_stdev = 1/sqrt(params.indegree); 
params.mu_E = 1;
params.level_of_chaos = level_of_chaos; % Passed as argument

%% set adaptation params
params.n_a_E = n_a_E;  % Passed as argument
params.n_a_I = 0;  % 0 to 3 adaptation for I neurons
params.tau_a_E = logspace(log10(0.05), log10(25), params.n_a_E);  % Logarithmically spaced from 0.1 to 10
params.tau_a_I = logspace(log10(0.1), log10(25), params.n_a_I);  % Logarithmically spaced from 0.1 to 10
params.c_E = 0.25*1/3;  % Adaptation scaling for E neurons (scalar, typically 0-3)
params.c_I = .1;  % Adaptation scaling for I neurons (scalar, typically 0-3)

params.n_b_E = n_b_E;  % Passed as argument
params.n_b_I = 0;  % Number of STD timescales for I neurons (0 or 1)
params.tau_b_E_rel = 0.05;  % STD release time constant for E neurons (s)
params.tau_b_I_rel = 0.1;  % STD release time constant for I neurons (s)
params.tau_b_E_rec = 1;  % STD recovery time constant for E neurons (s)
params.tau_b_I_rec = 1;  % STD recovery time constant for I neurons (s)

%% set activaiton function
S_a = 0.9;
S_c = 0.45;
params.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

%% dependent E vs I params
params.n_E = round(params.f*params.n);
params.n_I = params.n-params.n_E;
params.N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n_E * params.n_b_E + params.n_I * params.n_b_I + params.n;
params.E_indices = 1:params.n_E;
params.I_indices = params.n_E+1:params.n;
params.mu_I = -params.f*params.mu_E/(1-params.f);


%% make W connection matrix
[W, M, G, Z] = create_W_matrix(params); % should use params.d not recompute it
% Scale W based on abscissa to control level of chaos
W_eigs = eig(W);
radius_0 = max(abs(W_eigs))
abscissa_0 = max(real(W_eigs))

gamma = 1 / abscissa_0;  % Scaling to reach edge of chaos
params.W = params.level_of_chaos * gamma * W;

% verify new abscissa and spectral radius
W_eigs_1 = eig(params.W);
radius_1 = max(abs(W_eigs_1))
abscissa_1 = max(real(W_eigs_1))

params.tau_d = 0.1;  % dentdritic time constant, seconds

%% ICs
S0 = initialize_state(params);

%% External input
% Configure paired pulse input
pulse_config.amp = 3;
pulse_config.width_ms = 0.5;
pulse_config.neuron_frac = 0.5; % Arbitrary choice, or could be passed
% Generate external input
[u_ex, t_ex] = generate_paired_pulse_input(params, T, fs, params.rng_seeds(2), pulse_config);
u_ex = u_ex.*u_ex_scale; % way to turn it off easily

% Integrate
ode_solver = @ode45;
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
jac_wrapper = @(t, S) compute_Jacobian_fast(S, params);
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', 1*dt, 'Jacobian', jac_wrapper);

disp('Integrating Equations...')
tic
[t_out, S_out] = ode_solver(rhs, t_ex, S0, opts);
toc

%% Compute Lyapunov exponent(s)
T_interval = [1, T];
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
lya_results = compute_lyapunov_exponents(Lya_method, S_out, t_out, dt, fs, T_interval, params, opts, ode_solver, rhs, t_ex, u_ex);

%% Decimate state vector for plotting
% Need less decimation to see pulses
% plot_deci = round(fs/500); % Approx 500 Hz plotting rate? Or just keep all points?
% Actually, let's ensure we catch the pulses. 
plot_deci = 1; % No decimation for paired pulse to ensure we see the 2ms gap
[t_plot, S_plot, plot_indices] = decimate_states(t_out, S_out, plot_deci);

%% Unpack state vector and compute firing rates
[x_plot, a_plot, b_plot, r_plot] = unpack_and_compute_states(S_plot, params);

% Split external input into E and I
u_ex_plot = u_ex(:, plot_indices);
u_plot.E = u_ex_plot(params.E_indices, :);
u_plot.I = u_ex_plot(params.I_indices, :);

%% Plotting
params.plot_mean_dendrite = true;
[~, ~] = plot_SRNN_paired_pulse_tseries(t_plot, u_plot, x_plot, r_plot, a_plot, b_plot, params, lya_results, Lya_method);

% Jacobian Analysis for Paired Pulse - DISABLED
% The original analysis was based on long steps. 
% For now, we skip the Jacobian/Eigenvalue plots as they require specific time points 
% that might need to be manually selected relative to the pulses (e.g., before/during/after).

%% Save results
if save_figs
    if ~isempty(save_dir)
        figs_folder = fullfile(save_dir, 'figs');
        if ~exist(figs_folder, 'dir')
            mkdir(figs_folder);
        end
        
        % Find all open figures
        fig_handles = findobj('Type', 'figure');
        % Sort by figure number
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


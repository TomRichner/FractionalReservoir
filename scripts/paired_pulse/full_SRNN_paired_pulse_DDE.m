function [t_out, S_out, params] = full_SRNN_paired_pulse_DDE(u_ex_scale, n_a_E, n_b_E, level_of_chaos, save_dir, save_figs, save_workspace, note)
% FULL_SRNN_PAIRED_PULSE_DDE Runs the SRNN simulation with paired-pulse paradigm using DDEs.
%
% Implements a paired E-I architecture with delayed inhibition and paired pulse input.

if nargin < 8 || isempty(note)
    note = 'SRNN_PP_DDE';
end

% background of all figs should be white
set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1);
set(groot, 'DefaultAxesLineWidth', 1);
set(groot, 'DefaultAxesTitleFontWeight', 'normal');

params.rng_seeds = [5 3 3 4 5];
rng(params.rng_seeds(1))

%% time 
fs = 2000;  % Higher sampling rate for short pulses
dt = 1/fs;
T = 10;    % Duration (s)
t_ex = (0:dt:T)';
nt = length(t_ex);

%% Simulation Parameters (DDE Paired Setup)
params.n = 100; 
% Ensure even n for pairing
if mod(params.n, 2) ~= 0
    params.n = params.n + 1;
    warning('Adjusting n to %d to allow equal E-I pairing.', params.n);
end
params.n_E = params.n / 2;
params.n_I = params.n / 2;
params.f = 0.5;
params.E_indices = 1:params.n_E;
params.I_indices = (params.n_E + 1):params.n;

params.indegree = 15; 
params.d = params.indegree/params.n;
params.mu_E = 1;
params.mu_I = -1; % Explicitly set for the paired matrix logic
params.G_stdev = 1/sqrt(params.indegree); 
params.level_of_chaos = level_of_chaos; 

%% Adaptation Params
params.n_a_E = n_a_E;
params.n_a_I = 0; 
params.tau_a_E = logspace(log10(0.25), log10(25), params.n_a_E); 
params.tau_a_I = logspace(log10(0.25), log10(25), params.n_a_I);
params.c_E = 0.25*1/3;  
params.c_I = .1;  

%% STD Params
params.n_b_E = n_b_E; 
params.n_b_I = 0; 
params.tau_b_E_rel = 0.5; 
params.tau_b_I_rel = 0.5; 
params.tau_b_E_rec = 2; 
params.tau_b_I_rec = 2; 

%% Activation Function
S_a = 0.9;
S_c = 0.4;
params.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

%% Create Paired Connectivity Matrix
[W, M, G, Z] = create_W_matrix(params);

% Calculate Dependent Params
params.N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + ...
                   params.n_E * params.n_b_E + params.n_I * params.n_b_I + params.n;

% Scale W for Chaos
W_eigs = eig(W);
abscissa_0 = max(real(W_eigs));
radius_0 = max(abs(W_eigs));

if abscissa_0 > 0
    gamma = 1 / abscissa_0; 
    params.W = params.level_of_chaos * gamma * W;
else
    % Fallback if abscissa is negative (stable)
    params.W = params.level_of_chaos * W; 
end

% Verify
W_eigs_1 = eig(params.W);
radius_1 = max(abs(W_eigs_1));
abscissa_1 = max(real(W_eigs_1));
fprintf('W scaled: Radius=%.2f, Abscissa=%.2f\n', radius_1, abscissa_1);

params.tau_d = 0.1; 

%% Setup DDE Components
% We assume E outputs are instant (Lag 0) and I outputs are delayed (Lag tau)
tau_syn_delay = 0.004; % 2 ms delay (kept same as DDE run, 0.004 seems to be 4ms?)
% In full_SRNN_run_DDE it was 0.004. Note comment said "2 ms delay" but value is 0.004. 
% I will trust the code value 0.004.
params.lags = [tau_syn_delay];

W_inst = zeros(params.n);
W_inst(:, params.E_indices) = params.W(:, params.E_indices);

W_delayed = zeros(params.n);
W_delayed(:, params.I_indices) = params.W(:, params.I_indices);

params.W_components = {W_inst, W_delayed};

%% ICs
S0 = initialize_state(params);

%% External Input (Paired Pulse)
pulse_config.amp = 3;
pulse_config.width_ms = 0.5;
pulse_config.neuron_frac = 0.5; 
% Generate external input
[u_ex, t_ex] = generate_paired_pulse_input(params, T, fs, params.rng_seeds(2), pulse_config);
u_ex = u_ex.*u_ex_scale;

%% DDE Solver
disp('Integrating DDEs...');
tic
history = S0; % Constant history for t < 0
options = ddeset('RelTol', 1e-5, 'AbsTol', 1e-5); 

sol = dde23(@(t,y,S_delay) SRNN_reservoir_DDE(t,y,S_delay, t_ex, u_ex, params), params.lags, history, [0 T], options);
toc

% Resample to uniform time grid matching t_ex
t_out = t_ex;
S_out = deval(sol, t_out)'; % deval returns cols as time points, we want rows

%% Decimate and Plot
plot_deci = 1; % No decimation for paired pulse to ensure we see the gap
[t_plot, S_plot, plot_indices] = decimate_states(t_out, S_out, plot_deci);

[x_plot, a_plot, b_plot, r_plot] = unpack_and_compute_states(S_plot, params);

u_ex_plot = u_ex(:, plot_indices);
u_plot.E = u_ex_plot(params.E_indices, :);
u_plot.I = u_ex_plot(params.I_indices, :);

% Placeholder for Lyapunov results
lya_results = [];
Lya_method = 'none';

%% Plotting
params.plot_mean_dendrite = true;
[~, ~] = plot_SRNN_paired_pulse_tseries(t_plot, u_plot, x_plot, r_plot, a_plot, b_plot, params, lya_results, Lya_method);

%% Save Results
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
        save(fullfile(save_dir, 'SRNN_run_DDE.mat'), '-v7.3', '-nocompression');
    else
        warning('save_workspace is true but save_dir is empty. Workspace not saved.');
    end
end

end


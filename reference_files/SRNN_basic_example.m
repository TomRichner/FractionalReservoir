%--- SRNN_caller.m

close all
clear all % must clear all due to use of persistant variables in SRNN.m
clc

tic

%% 
seed = 42;
rng(seed,'twister');

%% Network
n = 10; % number of neurons

Lya_method = 'benettin'; % 'benettin', 'qr', 'svd', or 'none'
use_Jacobian = false;

mean_in_out_degree = 5; % desired mean number of connections in and out
density = mean_in_out_degree/(n-1); % each neuron can make up to n-1 connections with other neurons
sparsity = 1-density;

EI = 0.7;
scale = 0.5/0.79782; % overall scaling factor of weights
w.EE = scale*1; % E to E. Change to scale*2 for bursting
w.EI = scale*1; % E to I connections
w.IE = scale*1; % I to E
w.II = scale*.5; % I to I
w.selfE = 0;    % self connections of E neurons
w.selfI = 0;    % self connections of I neurons

[M, EI_vec] = generate_M_no_iso(n,w,sparsity, EI);
EI_vec = EI_vec(:); % make it a column
[E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec);

%% Time
fs = 1000; %Plotting sample frequency
dt = 1/fs;
T = [-10 10];

T_lya_1 = T(1);

nt = round((T(2)-T(1))*fs)+1; % Number of plotting samples
t = linspace(T(1), T(2), nt)'; % Plotting time vector

%% u_ex, external input, stimulation

u_ex = zeros(n, nt);
% sine and square wave stim
stim_b0 = 0.5; amp = 0.5;
dur = 3; % duration of sine
f_sin = 1.*ones(1,fs*dur);
% f_sin = logspace(log10(0.5),log10(3),fs*5);
u_ex(1,-t(1)*fs+fix(fs*6)+(1:fix(fs*dur))) = stim_b0+amp.*sign(sin(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))'));
u_ex(1,-t(1)*fs+fix(fs*1)+(1:fix(fs*dur))) = stim_b0+amp.*-cos(2*pi*f_sin(1:fix(fs*dur)).*t(1:fix(fs*dur))');
u_ex = u_ex*1;
u_ex = u_ex(:,1:nt);
DC = 0.1;

% Ramp up to DC over the first 5 seconds to avoid big IC transient
ramp_duration = 5; % seconds
ramp_logical = t <= (T(1) + ramp_duration);
ramp_profile = linspace(0, DC, sum(ramp_logical));
u_dc_profile = ones(1, nt) * DC;
u_dc_profile(ramp_logical) = ramp_profile;
u_ex = u_ex + u_dc_profile;

%% parameters

tau_STD = 0.5; % scalar, time constant of synaptic depression

% Define number of timescales for E and I neurons separately
n_a_E = 3; % typically 3, number of SFA timescales for E neurons
n_a_I = 0; % typically 0, number of SFA timescales for I neurons (typically 0)
n_b_E = 1; % typically 1 or 2, number of STD timescales for E neurons
n_b_I = 0; % typically 0, number of STD timescales for I neurons (typically 0)

% Define tau_a and tau_b for E and I neurons
% Ensure these are empty if the corresponding n_a_X or n_b_X is 0
if n_a_E > 0
    tau_a_E = logspace(log10(0.3), log10(15), n_a_E); % s, 1 x n_a_E
else
    tau_a_E = [];
end
if n_a_I > 0
    tau_a_I = logspace(log10(0.3), log10(15), n_a_I); % s, 1 x n_a_I 
else
    tau_a_I = [];
end

if n_b_E > 0
    tau_b_E = logspace(log10(0.6), log10(9), n_b_E);  % s, 1 x n_b_E
    if n_b_E == 1 % Specific condition from original code
        tau_b_E = 4*tau_STD;
    end
else
    tau_b_E = [];
end
if n_b_I > 0
    tau_b_I = logspace(log10(0.6), log10(9), n_b_I); % s, 1 x n_b_I
    if n_b_I == 1 % Retain similar logic if ever used
        tau_b_I = 4*tau_STD;
    end
else
    tau_b_I = [];
end


tau_d = 0.025; % s, scalar

if n_a_E > 0
    c_SFA = (1/n_a_E) * double(EI_vec == 1); % n x 1, Example: SFA only for E neurons
else
    c_SFA = zeros(n, 1);
end
% c_SFA(I_indices) = 0; % Explicitly set to 0 for I if desired, or rely on n_a_I = 0
F_STD = 1 * double(EI_vec == 1); % n x 1, Example: STD only for E neurons
% F_STD(I_indices) = 0; % Explicitly set to 0 for I if desired, or rely on n_b_I = 0

params = package_params(n_E, n_I, E_indices, I_indices, ...
                        n_a_E, n_a_I, n_b_E, n_b_I, ...
                        tau_a_E, tau_a_I, tau_b_E, tau_b_I, ...
                        tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec);


%% Initial Conditions
a0_E = [];
if params.n_E > 0 && params.n_a_E > 0
    a0_E = zeros(params.n_E * params.n_a_E, 1);
end

a0_I = [];
if params.n_I > 0 && params.n_a_I > 0
    a0_I = zeros(params.n_I * params.n_a_I, 1);
end

b0_E = [];
if params.n_E > 0 && params.n_b_E > 0
    b0_E = ones(params.n_E * params.n_b_E, 1);
end

b0_I = [];
if params.n_I > 0 && params.n_b_I > 0
    b0_I = ones(params.n_I * params.n_b_I, 1);
end

u_d0 = zeros(n, 1);

X_0 = [a0_E; a0_I; b0_E; b0_I; u_d0];

N_sys_eqs = size(X_0,1); % Number of system equations / states

%% Integrate with ODE solver

SRNN_wrapper = @(tt,XX) SRNN(tt,XX,t,u_ex,params); % inline wrapper function to add t, u_ex, and params to SRNN

% wrap ode_RKn to limit the exposure of extra parameters for usage to match builtin integrators
solver_method = 6; % 5 is classic RK4
deci = 1; % deci > 1 does not work for benettin's method.  Need to fix this

% these ode_options only apply if using ode15s
ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e8, 'MaxStep',dt, 'InitialStep', 0.1*dt); % RelTol must be less than perturbation d0, which is 1e-3

ode_solver = @ode15s; % stiff ode solver

[~, X] = ode_solver(SRNN_wrapper, t, X_0, ode_options);

% This block computes LLEs for Phase 2, and decides whether to keep them or revert to Phase 1 results
lya_results = struct();

if ~strcmpi(Lya_method, 'none')
    if strcmpi(Lya_method,'qr')
        lya_dt = 4*tau_d;
    else
        lya_dt = 0.5*tau_d;
    end
    lya_calc_start_idx = find(t >= T_lya_1, 1, 'first');
    if isempty(lya_calc_start_idx)
        error('Could not find T_lya_1 in time vector t. Check T and T_lya_1 values.');
    end
    X_for_lya = X(lya_calc_start_idx:end, :);
    t_for_lya = t(lya_calc_start_idx:end);
    
    switch lower(Lya_method)
        case 'svd'
            error('svd not verified to be working yet')
            fprintf('Computing full Lyapunov spectrum using SVD method...\n');
            [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_svd(X_for_lya, t_for_lya, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);
        case 'qr'
            fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
            [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(X_for_lya, t_for_lya, lya_dt, params, ode_solver, ode_options, @SRNN_Jacobian, T, N_sys_eqs, fs);
        case 'benettin'
            fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
            d0 = 1e-3;
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X_for_lya, t_for_lya, dt, fs, d0, T, lya_dt, params, ode_options, @SRNN, t, u_ex, ode_solver);
    end

    if strcmpi(Lya_method, 'benettin')
        lya_results.LLE = LLE; lya_results.local_lya = local_lya; lya_results.finite_lya = finite_lya; lya_results.t_lya = t_lya;
    else % qr or svd
        lya_results.LE_spectrum = LE_spectrum; lya_results.local_LE_spectrum_t = local_LE_spectrum_t; lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t; lya_results.t_lya = t_lya; lya_results.N_sys_eqs = N_sys_eqs;
    end

end

%% Convert state X to named variables and compute dependent variables for plotting and comparisons to analytic method
[a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts] = unpack_SRNN_state(X, params);
[r, p] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params);

%% Make plots using the plotting function

% Call the plotting function
if ~strcmpi(Lya_method, 'none') && ~isempty(fieldnames(lya_results))
    SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method, lya_results);
else
    SRNN_tseries_plot(t, u_ex, r, a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params, T, Lya_method);
end

sim_dur = toc

simulation_duration_vs_real_time = sim_dur./(T(2)-T(1))
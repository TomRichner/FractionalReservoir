% test_SRNN2_defaults.m - Run SRNNModel2 with all defaults and plot
%
% Constructs SRNNModel2 with n_a_E=3 (SFA) and n_b_E=1 (STD), all other
% parameters at defaults (n=300, T_range=[0,50], lya_method='benettin').

close all; clear; clc;

%% Create model
rng_seeds = [1 2] + 21;
model = SRNNModel2('tau_d', 0.1, 'level_of_chaos',1,'f',2/3,'n_a_E', 3, 'n_b_E', 1, 'tau_b_E_rec', [1], 'S_c', 0.35, 'ode_solver', @ode_rk4, 'fs', 200, 'check_connectivity', true, 'T_range', [0 40], 'rng_seeds',rng_seeds);

F = model.default_val;
model.mu_E_tilde =  3*F;   % default
model.mu_I_tilde = -6*F;   % default is -4*F; doubled since there are half as many I neurons


%% Build, run, and plot
model.build();
model.run();
model.plot();

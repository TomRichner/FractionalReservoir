% test_SRNN2_defaults.m - Run SRNNModel2 with all defaults and plot
%
% Constructs SRNNModel2 with n_a_E=3 (SFA) and n_b_E=1 (STD), all other
% parameters at defaults (n=300, T_range=[0,50], lya_method='benettin').
% Activation is overridden to piecewiseSigmoid(x, S_a, S_c).

close all; clear; clc;

rng_seeds = [19 20]+4;

%% Create model
model = SRNNModel2('n',300,'indegree',100,'check_connectivity',true,'level_of_chaos', 1,'S_c',0.0,'tau_d',1,'tau_b_E_rec',1,'rng_seeds',rng_seeds,'c_E', 1/3,'n_a_E', 3, 'n_b_E',1, 'f', 2/3, 'T_range', [0 40], ...
    'ode_solver', @ode_rk4, 'fs', 200);

%% E:I ratio of 2:1 (f=2/3), with per-synapse inhibition doubled to compensate
% The _relative properties are in multiples of the RMT normalization factor
% F = 1/sqrt(n*alpha*(2-alpha)) (Harris 2023), which depends only on n and
% alpha, so it is unchanged by f.
model.mu_E_tilde_relative =  3;   % default
model.mu_I_tilde_relative = -6;   % default is -4; doubled since there are half as many I neurons

%% Use the piecewise sigmoid instead of the default logistic sigmoid
model.activation_function            = @(x) SRNNModel2.piecewiseSigmoid(x, model.S_a, model.S_c);
model.activation_function_derivative = @(x) SRNNModel2.piecewiseSigmoidDerivative(x, model.S_a, model.S_c);

%% Build, run, and plot
model.build();
model.run();
model.plot();

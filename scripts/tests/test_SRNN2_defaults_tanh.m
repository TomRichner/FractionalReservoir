% test_SRNN2_defaults.m - Run SRNNModel2 with all defaults and plot
%
% Constructs SRNNModel2 with n_a_E=3 (SFA) and n_b_E=1 (STD), all other
% parameters at defaults (n=300, T_range=[0,50], lya_method='benettin').

close all; clear; clc;

%% Create model
model = SRNNModel2('n_a_E', 3, 'n_b_E', 1, 'f', 0.5, 'level_of_chaos',1.5);

%% Use tanh activation instead of the default piecewise sigmoid
model.activation_function = @SRNNModel2.tanhActivation;
model.activation_function_derivative = @SRNNModel2.tanhActivationDerivative;

%% Build, run, and plot
model.build();
model.run();
model.plot();

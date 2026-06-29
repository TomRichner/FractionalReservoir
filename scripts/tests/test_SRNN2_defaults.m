% test_SRNN2_defaults.m - Run SRNNModel2 with all defaults and plot
%
% Constructs SRNNModel2 with n_a_E=3 (SFA) and n_b_E=1 (STD), all other
% parameters at defaults (n=300, T_range=[0,50], lya_method='benettin').

close all; clear; clc;

% Add paths
setup_paths();

%% Create model
model = SRNNModel2('n_a_E', 3, 'n_b_E', 1);

%% Build, run, and plot
model.build();
model.run();
model.plot();

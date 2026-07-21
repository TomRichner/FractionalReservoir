% test_SRNN2_defaults.m - Run SRNNModel2 with all defaults and plot
%
% Constructs SRNNModel2 with n_a_E=3 (SFA) and n_b_E=1 (STD), all other
% parameters at defaults (n=300, T_range=[0,50], lya_method='benettin').

close all; clear; clc;

rng_seeds = [1 3] + 4;

% Add paths
setup_paths();

%% Create model
model = SRNNModel2('tau_d',0.1,'rng_seeds',rng_seeds,'tau_b_E_rel',0.25,'std_zero_floor',false,'c_E', 0.15/3,'S_c', 0.35,'n_a_E', 3, 'n_b_E',1, 'tau_a_E', logspace(log10(0.25), log10(10), 3),  'tau_b_E_rec', [2], 'f', 0.5, 'level_of_chaos', 10, 'T_range', [0 30]);

%% Build, run, and plot
model.build();
model.run();
model.plot();

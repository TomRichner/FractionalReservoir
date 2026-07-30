% test_SRNN2_balanced_mu.m - SRNNModel2 with balanced E/I means and plot
%
% Same as test_SRNN2_defaults.m (n_a_E=3 SFA, n_b_E=1 STD, f=0.5, all other
% parameters at defaults) but with mu_E_tilde = 3*F and mu_I_tilde = -3*F
% instead of the default 3*F / -4*F, i.e. E and I means balanced in magnitude.
%
% F = 1/sqrt(n*alpha*(2-alpha)) is the RMT normalization factor (Harris 2023),
% exposed as the dependent property default_val.

close all; clear; clc;

%% Create model
model = SRNNModel2('tau_a_E',logspace(log10(0.25), log10(10), 3),'n_a_E', 3, 'n_b_E', 1, 'f', 2/3,'indegree',100,'T_range',[0 100],'S_c',0,'rng_seeds',[2 3],'lya_method','none','c_E',0.5/3);
% model = SRNNModel2('tau_a_E',logspace(log10(0.25), log10(10), 3),'n_a_E', 3, 'n_b_E', 0, 'f', 2/3,'indegree',100,'T_range',[0 100],'S_c',0,'rng_seeds',[2 3],'lya_method','none','c_E',0.5/3);

% Override the E/I population means (must be set before build())
F = model.default_val;
model.mu_E_tilde =  10 * F;
model.mu_I_tilde = -20 * F;

% stim steps
n_steps = 3;
model.input_config.n_steps = n_steps;
model.input_config.no_stim_pattern = true(1, n_steps);        % start all silent
model.input_config.no_stim_pattern(2:2:end) = false;      % drive blocks 2, 4

%% Build, run, and plot
model.build();
model.run();
model.plot();

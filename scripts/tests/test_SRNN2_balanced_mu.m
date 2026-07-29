% test_SRNN2_balanced_mu.m - SRNNModel2 with balanced E/I means and plot
%
% Same as test_SRNN2_defaults.m (n_a_E=3 SFA, n_b_E=1 STD, f=0.5, all other
% parameters at defaults) but with mu_E_tilde = 3*F and mu_I_tilde = -3*F
% instead of the default 3*F / -4*F, i.e. E and I means balanced in magnitude.
%
% F = 1/sqrt(n*alpha*(2-alpha)) is the RMT normalization factor (Harris 2023),
% exposed as the dependent property default_val.

close all; clear; clc;

% Add paths
setup_paths();

%% Create model
model = SRNNModel2('n_a_E', 3, 'n_b_E', 1, 'f', 0.5);

% Override the E/I population means (must be set before build())
F = model.default_val;
model.mu_E_tilde =  3 * F;
model.mu_I_tilde = -3 * F;

%% Build, run, and plot
model.build();
model.run();
model.plot();

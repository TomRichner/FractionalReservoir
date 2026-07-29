% test_SRNN2_multi_std.m - Run SRNNModel2 with multi-timescale STD and plot
%
% Constructs SRNNModel2 with n_a_E=3 (SFA) and n_b_E=2 (two STD timescales,
% tau_b_E_rec = [1, 3] s), all other parameters at defaults (n=300,
% T_range=[0,50], lya_method='benettin'). The synaptic depression factor
% is the product of the two per-timescale resources, prod_m b_{i,m}.

close all; clearvars -except rng_seeds; clc;



if exist('rng_seeds', 'var')
    rng_seeds = rng_seeds + 1
else
    rng_seeds = [1 2]
end

%% Create model
tau_d = 0.1;
n_a_E = 3;
tau_a_E = logspace(log10(0.5), log10(15), n_a_E);
n_b_E = 1;
tau_b_E_rec = [1];
c_E = 1/n_a_E;
level_of_chaos = 1;
std_zero_floor = false;
T_range = [-5 15];
model = SRNNModel2('T_range', T_range, 'tau_d',tau_d, 'rng_seeds', rng_seeds, 'n_a_E', n_a_E, 'n_b_E', n_b_E, 'level_of_chaos',level_of_chaos,'tau_a_E',tau_a_E, 'tau_b_E_rec', tau_b_E_rec, 'c_E', c_E, 'std_zero_floor', std_zero_floor, 'store_full_state', true);

%% Build, run, and plot
model.build();
model.plot_W_spectrum()
model.plot_weight_histogram()
model.run();
model.plot_eigenvalues(mean(model.T_range));   % Jacobian spectrum at simulation midpoint (25 s)
model.plot();

% test_SRNNCellTypePairs_defaults.m - Run SRNNCellTypePairs at the SRNNModel2
% default operating point and plot.
%
% The Pairs analogue of test_SRNN2_defaults.m, which is just
%
%   model = SRNNModel2('n_a_E', 3, 'n_b_E', 1, 'f', 0.5);
%
% SRNNCellTypePairs needs a little more because it is general over cell types:
% the type list and the RMT block statistics are required rather than defaulted.
% Everything else below is the class default, chosen to match SRNNModel2's:
%
%   n = 300, indegree = 100          (class defaults in both)
%   mu_tilde_relative    = [3 -4]    matches mu_E/I_tilde_relative
%   sigma_tilde_relative = [1 1]     matches sigma_E/I_tilde_relative
%   tau_a = logspace(0.25, 10, 3)    filled in by both when n_a = 3
%   c = 0.15/3                       SFA scaling, same default
%   activation 'logistic', S_c = 0.4, tau_d = 0.1, T_range [0 50], fs 400
%
% The one structural difference worth understanding: SRNNModel2's n_b_E = 1 puts
% depression on every OUTGOING excitatory synapse. Here that is spelled out as
% two routes, E->E and E->I -- which is exactly the freedom this class buys.
% Delete the E->I line and depression applies to E->E connections only.
%
% See also: test_SRNN2_defaults, SRNNCellTypePairs, test_SRNNCellTypePairs_blocks

close all; clear; clc;

%% Short-term depression on all outgoing E synapses (SRNNModel2's n_b_E = 1)
% tau_rec / tau_rel match SRNNModel2's tau_b_E_rec = 1, tau_b_E_rel = 0.25.
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
% synapse_config.I.E.stf = struct('tau_dec', 5, 'tau_fac', 0.75, 'G',1.5);
synapse_config.I.I.std = struct('tau_rec', 4, 'tau_rel', 1);
% synapse_config.E.I.std = struct('tau_rec', 1, 'tau_rel', 0.25);
% synapse_config.I.E.std = struct('tau_rec', 1, 'tau_rel', 0.25);

%% Create model
Tend = 25;
model = SRNNCellTypePairs( ...
    'rng_seeds',[19 20]+7, ...
    'activation', 'piecewise', ...
    'S_c', 0.0, ...
    'S_a', 0.8, ...
    'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], ...
    'mu_tilde_relative',    [4 -3; 3 -3], ...   % multiples of F, (post <- pre)
    'sigma_tilde_relative', [1 1; 1 1], ...
    'level_of_chaos', 1.3, ...
    'n_a', [3 0], ...                     % SFA on E only, 3 timescales
    'c', [0.5/3, 0], ...
    ... % Per-neuron setpoints: one mean and one SD per cell type. Uncomment to
    ... % give E and I different excitability plus cell-to-cell spread. Left
    ... % out, every neuron shares the scalar S_c above.
    ... 'mu_S_c',    [0.0 -0.1], ...
    ... 'sigma_S_c', [0.05 0.05], ...
    'T_range', [-5 Tend],...
    'ode_solver', 'rk4', ...
    'fs', 200, ...
    'T_plot', [0 Tend], ...
    'lya_T_interval', [max(0,Tend-10), Tend], ...
    'synapse_config', synapse_config);

model.input_config.intrinsic_drive = 0.1 * ones(model.n, 1);

% plot_eigenvalues reads the Jacobian off the full trajectory, so the full state
% has to be kept. At n=300 over [0 20] at fs=200 that is ~4001 x N_sys_eqs
% doubles -- tens of MB, fine here, but drop it if you scale n up.
model.store_full_state = true;

%% Build and run
model.build();
model.run();

%% Plots
model.plot();                       % compact summary: one panel per quantity

% Per-cell-type view: one COLUMN per type, every individual neuron trace, with
% b and g collapsed across routes as prod(b) and coloured by target type.
model.plot_celltypes();

% Effective-Jacobian spectrum at three times: early (still settling), middle,
% and late. Watch whether adaptation pulls eigenvalues back inside the unit
% circle as the network engages.
model.plot_eigenvalues([2 10 18]);

% W itself, imaged with a diverging colormap so zero is white and the sign of
% each synapse reads at a glance. Cell-type boundaries mark the (post, pre)
% blocks, so a raised mu_EE shows as one quadrant going redder. This is the
% SCALED W -- what was simulated -- not the generator's unscaled output.
model.plot_W();

% Eigenvalues of W itself. This is where a raised mu_EE shows up directly --
% the outlier leaving the bulk disk of radius R, rather than only as the
% lambda_O number printed below.
model.plot_W_spectrum();

% Weight distribution per cell type. The negative tail on the E columns is the
% Gaussian sampler, not a bug: Dale's law here is statistical, so a fraction
% ~normcdf(-mu_tilde/sigma_tilde) of E synapses come out negative.
model.plot_weight_histogram();

%% Report the RMT predictions the block generator provides
fprintf('\n=== SRNNCellTypePairs defaults ===\n');
fprintf('F (normalization)      = %.6g\n', model.default_val);
fprintf('mu_tilde (post<-pre)   = %s\n', mat2str(model.mu_tilde, 4));
fprintf('bulk radius R          = %.4f\n', model.R);
fprintf('outlier eigenvalues    = %s\n', mat2str(round(model.lambda_O(:)', 4)));
fprintf('state dimension        = %d\n', model.N_sys_eqs);
fprintf('LLE                    = %.6f\n', model.lya_results.LLE);

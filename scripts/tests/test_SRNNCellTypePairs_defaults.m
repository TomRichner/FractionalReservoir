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
synapse_config.E.E.std = struct('tau_rec', 1, 'tau_rel', 0.25);
synapse_config.E.I.std = struct('tau_rec', 1, 'tau_rel', 0.25);
% synapse_config.I.E.std = struct('tau_rec', 1, 'tau_rel', 0.5);

%% Create model
model = SRNNCellTypePairs( ...
    'rng_seeds',[19 20], ...
    'activation', 'piecewise', ...
    'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], ...
    'mu_tilde_relative',    [3 -4; 3 -4], ...   % multiples of F, (post <- pre)
    'sigma_tilde_relative', [1 1; 1 1], ...
    'level_of_chaos', 1.2, ...
    'n_a', [3 0], ...                     % SFA on E only, 3 timescales
    'T_range', [0 20],...
    'synapse_config', synapse_config);

%% go faster
model.ode_solver = @ode_rk4;
model.fs = 200;

%% Build, run, and plot
model.build();
model.run();
model.plot();

%% Report the RMT predictions the block generator provides
fprintf('\n=== SRNNCellTypePairs defaults ===\n');
fprintf('F (normalization)      = %.6g\n', model.default_val);
fprintf('mu_tilde (post<-pre)   = %s\n', mat2str(model.mu_tilde, 4));
fprintf('bulk radius R          = %.4f\n', model.R);
fprintf('outlier eigenvalues    = %s\n', mat2str(round(model.lambda_O(:)', 4)));
fprintf('state dimension        = %d\n', model.N_sys_eqs);
fprintf('LLE                    = %.6f\n', model.lya_results.LLE);

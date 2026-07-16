 % Fig_example_cellTypeNet.m - Example cell-typed reservoir for the STF presentation
%
% Presentation-figure version of test_celltypes_defaults.m. Builds one
% SRNNModelCellTypes network with SFA + STD + STF on, parameterized from the
% committed Campagnola matrices (4 cell types: Pyr/Pvalb/Sst/Vip -- block
% connectivity, per-type STD/STF release-model params, fitted SFA tau_a, per-type
% n_a with non-adapting Pvalb). The by-cell-type grid figure is saved to this
% folder.
%
% Cell-type fractions are overridden so the network is 50% excitatory (Pyr): the
% remaining 50% is split among the inhibitory types (Pvalb/Sst/Vip) keeping their
% Campagnola relative proportions -- i.e. the default inhibitory ratios
% [0.08 0.07 0.05] rescaled by 0.50/0.20 = 2.5 -> [0.20 0.175 0.125].

close all; clear; clc;
setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% Cell-type fractions: 50% Pyr, inhibitory rescaled to keep relative proportions
% default: [0.80 0.08 0.07 0.05] (pyr pvalb sst vip); inhibitory sums to 0.20.
% Target pyr = 0.50 -> inhibitory must sum to 0.50; rescale inh by 0.50/0.20 = 2.5.
frac_E = 0.8;
inh_default   = [0.08, 0.07, 0.05];                 % pvalb, sst, vip (default)
inh_rescaled  = inh_default * ((1-frac_E) / sum(inh_default));
cell_frac    = [frac_E, inh_rescaled]               % -> [0.50 0.20 0.175 0.125]

%% Create model: 500 neurons, all 4 cell types from Campagnola data
% Piecewise sigmoid nonlinearity, S_a=0.9 (shape), S_c=0.35 (center), r_max=1.
model = SRNNModelCellTypes('fs',500, ...
    'T_range',[0 25], ...
    'n',300, 'lya_method', 'benettin', ...
    'type_fractions', cell_frac, ...
    'c_gain', 0.0, ...
    'S_a', 0.9, 'S_c', 0.35, ...
    'activation_function', @(x) SRNNModelBase.piecewiseSigmoid(x, 0.9, 0.35, 10), ...
    'activation_function_derivative', @(x) SRNNModelBase.piecewiseSigmoidDerivative(x, 0.9, 0.35, 10));

% Stronger external stimulus drive (default input_config.amp = 0.5)
model.input_config.amp = 8.0;

%% Build, run, and plot
model.build();

% Report realized per-type neuron counts and the data-derived per-type SFA
fprintf('\nPer-type composition and SFA (Campagnola data):\n');
for T = 1:model.n_types
    j = find(model.type_of == T, 1);
    fprintf('  %-6s n=%3d (%.1f%%)  adaptation_index=%.4f  adapting=%d  tau_a=%.3f s  tau_d=%.4f s\n', ...
        model.type_names{T}, nnz(model.type_of == T), 100*nnz(model.type_of == T)/model.n, ...
        model.adapt_index(T), ismember(j, model.ad_idx), ...
        model.tau_a_type(T), model.tau_d_type(T));
end

model.run();
[fig, ~] = model.plot_by_celltype();

% Weight matrix colored by presynaptic cell type (+ separate gradient legend)
[fig_W, fig_W_legend] = model.plot_W_bytype();

%% Save
% save_some_figs_to_folder_2(this_dir, 'Fig_example_cellTypeNet', fig.Number, {'fig', 'png', 'svg'});

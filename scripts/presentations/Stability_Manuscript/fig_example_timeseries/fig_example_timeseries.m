% fig_example_timeseries.m - Example SRNNModel2 time-series figure
%
% Essentially scripts/tests/test_SRNN2_defaults.m (SRNNModel2 at defaults with
% n_a_E=3 SFA and n_b_E=1 STD), but using the PIECEWISE sigmoid activation
% instead of the class-default logistic sigmoid. The piecewise variant is
% evaluated at the class-default parameters (S_a=0.9, S_c=0.4). model.plot()
% renders the stimulus, dendritic state x, firing rate r, synaptic output
% b.*r, and the adaptation/depression variables.

close all; clear; clc;

% Add paths
setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% Piecewise sigmoid activation (class-default S_a, S_c)
S_a = 0.9;   % piecewiseSigmoid slope param (class default)
S_c = 0.4;   % piecewiseSigmoid center param (class default)

%% Create model
model = SRNNModel2('n_a_E', 3, 'n_b_E', 1, 'f', 0.5, ...
    'activation', 'logistic', 'S_a', S_a, 'S_c', S_c);

%% Build, run, and plot
model.build();
model.run();
[fig_handle, ~] = model.plot();

%% Save the figure
save_some_figs_to_folder_2(this_dir, 'fig_example_timeseries', fig_handle.Number, {'fig', 'png', 'svg'});

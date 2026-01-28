% test_zrs_modes.m
% Test script comparing SRNN behavior with different ZRS modes
%
% This script builds and runs 3 networks with identical parameters except
% for the ZRS mode, then plots their activity time series and eigenvalue
% spectra for comparison.

clear; close all; clc;

%% Setup
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

%% Define ZRS modes to test
zrs_modes = {'none', 'SZRS', 'Partial_SZRS'};
n_modes = length(zrs_modes);

%% Run simulations and plot
for i = 1:n_modes
    fprintf('\n=== Testing ZRS mode: %s ===\n', zrs_modes{i});

    % Create model with this ZRS mode (all other params use defaults)
    model = SRNNModel('zrs_mode', zrs_modes{i}, 'level_of_chaos', 1.15, 'S_c', 0, 'n', 500, 'rng_seeds', [1 2]);

    % Build and run
    model.build();

    % Plot eigenvalue spectra (W and LTI Jacobian)
    [fig_spec, ~] = model.plot_W_spectrum();
    sgtitle(fig_spec, sprintf('Eigenvalue Spectra - ZRS mode: %s', zrs_modes{i}), ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    % Run simulation
    model.run();

    % Use the built-in plot method for time series
    [fig, ~] = model.plot();

    % Update figure title to include ZRS mode
    sgtitle(fig, sprintf('Time Series - ZRS mode: %s', zrs_modes{i}), 'FontWeight', 'bold', 'Interpreter', 'none');
end

fprintf('\n=== All simulations complete ===\n');

% run_paired_pulse_mi.m
% Example script demonstrating how to use the PairedPulseMIAnalysis class
%
% This script performs paired-pulse stimulation experiments to measure how
% well the neural response during a probe pulse (pulse 2) encodes the
% amplitude of a preceding conditioning pulse (pulse 1).
%
% The analysis compares mutual information (MI) across the four adaptation
% conditions, using the SAME network connectivity for fair comparison.
%
% Key metrics:
%   - MI vs delay curve: how MI evolves during pulse 2
%   - Peak MI: maximum information about pulse 1 amplitude
%   - LLE (Lyapunov) and mean rate distributions
%
% See also: PairedPulseMIAnalysis, SRNNModel

clear;
clc;
close all;

%% Setup paths
setup_paths();

%% Create PairedPulseMIAnalysis object
pp = PairedPulseMIAnalysis(...
    'n_networks', 5, ...           % Number of different network realizations
    'n_levels', 5, ...              % Number of levels per grid parameter
    'note', 'bigger_mi', ...             % Optional note for folder naming
    'verbose', true ...             % Print progress during execution
    );

% Copy this script to the output directory
copyfile([mfilename('fullpath') '.m'], pp.output_dir);

%% Configure model defaults
% Benettin LLE is enabled by default for the new plotting methods
pp.model_defaults.n = 20;                   % Network size
pp.model_defaults.fs = 200;                 % Sampling frequency
pp.model_defaults.T_range = [0, 200];       % Simulation time
pp.model_defaults.level_of_chaos = 1.5;     % Abscissa scaling
% pp.model_defaults.lya_method = 'benettin'; % Already default

%% Optional: Add grid parameters for parameter space exploration
% Uncomment to test across parameter combinations:
%
pp.add_grid_parameter('level_of_chaos', [0.5, 3]);
% pp.add_grid_parameter('n', [50, 200]);

%% Debug: Visualize time series for parameter tuning
% Run simulations to see stimulus and response before full analysis
fprintf('\n=== Debug: Visualizing time series for parameter tuning ===\n');
pp.plot_debug(1, 'no_adaptation');
pp.plot_debug(1, 'sfa_only');
pp.plot_debug(1, 'std_only');
pp.plot_debug(1, 'sfa_and_std');
fprintf('Press any key to continue with full analysis (or Ctrl+C to abort and tune parameters)...\n');
pause;

%% Run the analysis
% Loops over networks and conditions, computing LLE, mean rate, and MI vs delay
pp.run();

%% Create all plots
% MI vs delay comparison
pp.plot();

% LLE and mean rate distributions per condition
pp.plot_histograms();

% LLE vs delay heatmap (mean MI)
pp.plot_mi_imagesc();

% LLE vs MI slice with bootstrap CIs
pp.plot_mi_vs_lle('delay_window', [0.5, 1.5],'lle_range', [-5, 3]);

%% Display summary
fprintf('\n=== Paired-Pulse MI Analysis Summary ===\n');
fprintf('Output directory: %s\n', pp.output_dir);
fprintf('Number of networks: %d\n', pp.n_networks);
fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, pp.conditions, 'UniformOutput', false), ', '));

%% Optional: Load and replot results from previous run
% pp_loaded = PairedPulseMIAnalysis();
% pp_loaded.load_results('/path/to/paired_pulse_MI_...');
% pp_loaded.plot_histograms();
% pp_loaded.plot_mi_imagesc();
% pp_loaded.plot_mi_vs_lle();

fprintf('\nDone! Results saved to: %s\n', pp.output_dir);


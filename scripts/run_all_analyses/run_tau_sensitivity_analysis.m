% run_tau_sensitivity_analysis.m
% Sensitivity analysis for tau_a_E(end) and tau_b_E_rec parameters
%
% Analyzes the effect of the maximum SFA time constant (tau_a_E_max)
% and STD recovery time constant (tau_b_E_rec) on SRNN dynamics.
%
% Uses ParamSpaceAnalysis2 with add_vector_parameter for tau_a_E
% and add_grid_parameter for tau_b_E_rec.
%
% See also: ParamSpaceAnalysis2, SRNNModel2

if ~exist('master_output_dir', 'var')
    clear;
    clc;
    close all;
end

%% Figure saving configuration
if exist('master_save_figs', 'var')
    if strcmp(master_save_figs, 'save_all_figs')
        save_figs = true;
    elseif strcmp(master_save_figs, 'save_no_figs')
        save_figs = false;
    end
end
if ~exist('save_figs', 'var')
    save_figs = false;
end

%% Setup paths
setup_paths();

%% Analysis Configuration
% Run mode: 'fast' for quick checks, 'production' for full-size runs.
% Set run_mode in the base workspace (or via run_all_analyses); defaults to
% 'production' when this script is run standalone.
if ~exist('run_mode', 'var'); run_mode = 'production'; end
switch run_mode
    case 'fast',       n_levels = 5;  n_reps = 5;
    case 'production', n_levels = 25; n_reps = 50;
    otherwise, error('run_tau_sensitivity_analysis:badMode', ...
        'Unknown run_mode ''%s'' (expected ''fast'' or ''production'').', run_mode);
end
fprintf('[run_tau_sensitivity_analysis] run_mode=%s, n_levels=%d, n_reps=%d\n', run_mode, n_levels, n_reps);
note = 'tau_timescales';

% Condition: SFA + STD (n_a_E=3, n_b_E=1)
condition = {struct('name', 'sfa_and_std', 'n_a_E', 3, 'n_b_E', 1)};

%% 1. tau_a_E(end) sweep — vector parameter
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity: tau_a_E(end) [5, 60] ===\n');
fprintf('========================================\n');

psa_tau_a = ParamSpaceAnalysis2(...
    'n_levels', n_levels, ...
    'batch_size', 25, ...
    'note', sprintf('%s_tau_a_E_max', note), ...
    'randomize_order', false, ...
    'verbose', true ...
    );
psa_tau_a.folder_prefix = 'tau_sensitivity';
if exist('master_output_dir', 'var')
    psa_tau_a.output_dir = master_output_dir;
end

psa_tau_a.set_conditions(condition);

% tau_a_E is a vector of length n_a_E=3, logspaced from 0.25 to max
% We sweep the max (last element) from 5 to 60
psa_tau_a.add_vector_parameter('tau_a_E', ...
    'vary_element', 'last', ...
    'fixed_value', 0.25, ...
    'vary_range', [5, 60], ...
    'n_elements', 3, ...
    'spacing', 'log', ...
    'level_spacing', 'linear');

psa_tau_a.add_grid_parameter('reps', 1:n_reps);

psa_tau_a.run();

% Copy script for reproducibility
copyfile([mfilename('fullpath') '.m'], psa_tau_a.output_dir);

% Plot
psa_tau_a.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1]);
psa_tau_a.plot_sensitivity('metric', 'mean_rate');

save(fullfile(psa_tau_a.output_dir, 'psa_object.mat'), 'psa_tau_a');

if save_figs
    fig_dir = fullfile(psa_tau_a.output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'tau_sensitivity_tau_a', [], {'fig', 'png'});
    fprintf('Figures saved to %s\n', fig_dir);
end
close all;

%% 2. tau_b_E_rec sweep — scalar parameter
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity: tau_b_E_rec [5, 60] ===\n');
fprintf('========================================\n');

psa_tau_b = ParamSpaceAnalysis2(...
    'n_levels', n_levels, ...
    'batch_size', 25, ...
    'note', sprintf('%s_tau_b_E_rec', note), ...
    'randomize_order', false, ...
    'verbose', true ...
    );
psa_tau_b.folder_prefix = 'tau_sensitivity';
if exist('master_output_dir', 'var')
    psa_tau_b.output_dir = master_output_dir;
end

psa_tau_b.set_conditions(condition);

psa_tau_b.add_grid_parameter('tau_b_E_rec', [5, 60]);
psa_tau_b.add_grid_parameter('reps', 1:n_reps);

psa_tau_b.run();

copyfile([mfilename('fullpath') '.m'], psa_tau_b.output_dir);

psa_tau_b.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1]);
psa_tau_b.plot_sensitivity('metric', 'mean_rate');

save(fullfile(psa_tau_b.output_dir, 'psa_object.mat'), 'psa_tau_b');

if save_figs
    fig_dir = fullfile(psa_tau_b.output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'tau_sensitivity_tau_b', [], {'fig', 'png'});
    fprintf('Figures saved to %s\n', fig_dir);
end
close all;

%% Summary
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity Analysis Complete ===\n');
fprintf('tau_a_E results: %s\n', psa_tau_a.output_dir);
fprintf('tau_b_E_rec results: %s\n', psa_tau_b.output_dir);
fprintf('========================================\n');

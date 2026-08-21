% Fig_2_fraction_excitatory_analysis.m
% Parameter space analysis: fraction excitatory (f) sweep
%
% Sweeps f=[0.4, 0.6] with 5 reps across 4 adaptation conditions.
% Generates unit histograms, paired beeswarm plots, and combined figure.
%
% Ported from ConnectivityAdaptation to use ParamSpaceAnalysis2 / SRNNModel2.
%
% See also: ParamSpaceAnalysis2, SRNNModel2

%% Configuration
% Edit these two. They used to be read out of the base workspace
% (`master_save_figs`, `master_output_dir`) because run_all_analyses once called
% this script; that call is commented out and the protocol is gone. No script in
% this repo passes settings to another through workspace variables -- see
% CLAUDE.md. If this analysis is ever wired back into a pipeline, make it a
% function taking these as arguments.
save_figs  = false;
output_dir = '';    % '' = let ParamSpaceAnalysis2 create its own dated folder

%% Setup paths
setup_paths();

% Derive project_root from setup_paths.m (at the repo root) so this tolerates living
% in a subdirectory such as scripts/EI_balance/.
project_root = fileparts(which('setup_paths'));
figs_root = fullfile(project_root, 'figs');

%% Create ParamSpaceAnalysis2 object
N = 300;
indegree = 100;

psa = ParamSpaceAnalysis2(...
    'n_levels', 5, ...
    'batch_size', 25, ...
    'note', 'frac_exc', ...
    'verbose', true ...
    );
psa.folder_prefix = 'fig2';
if ~isempty(output_dir)
    psa.output_dir = output_dir;
end

%% Add parameters to the grid
psa.add_grid_parameter('f', [0.4, 0.6]);
psa.add_grid_parameter('reps', 1:3);

%% Configure model defaults
psa.model_defaults.n = N;
psa.model_defaults.indegree = indegree;

% Timing
psa.model_defaults.T_range = [-15, 45];
psa.model_defaults.tau_d = 0.1;

% RMT tilde-parameters (Harris 2023), in multiples of F = default_val
psa.model_defaults.mu_E_tilde_relative = 3.5;
psa.model_defaults.mu_I_tilde_relative = -3.5;
psa.model_defaults.sigma_E_tilde_relative = 1;
psa.model_defaults.sigma_I_tilde_relative = 1;
psa.model_defaults.E_W_relative = -0.5;
psa.model_defaults.zrs_mode = 'none';
psa.model_defaults.level_of_chaos = 1.0;
psa.model_defaults.rescale_by_abscissa = false;

% Adaptation parameters
psa.model_defaults.c_E = 0.25/3;
psa.model_defaults.tau_a_E = logspace(log10(0.1), log10(10), 3);
psa.model_defaults.tau_b_E_rec = 1;
psa.model_defaults.tau_b_E_rel = 0.5;

% Activation function (use SRNNModel2 static methods)
S_a = 0.9;
S_c = 0.40;
psa.model_defaults.activation = 'piecewise';
psa.model_defaults.S_a = S_a;
psa.model_defaults.S_c = S_c;

% Input configuration
psa.model_defaults.u_ex_scale = 1.0;
psa.model_defaults.input_config.n_steps = 3;
psa.model_defaults.input_config.positive_only = true;
psa.model_defaults.input_config.step_density = 0.15;
psa.model_defaults.input_config.step_density_E = 0.15;
psa.model_defaults.input_config.step_density_I = 0;
psa.model_defaults.input_config.amp = 0.5;
psa.model_defaults.input_config.no_stim_pattern = [true false true];
psa.model_defaults.input_config.intrinsic_drive = zeros(N, 1);

% Lyapunov settings
psa.model_defaults.lya_method = 'benettin';
psa.store_local_lya = true;
psa.store_local_lya_dt = 0.1;

%% Run
psa.run();

% Save PSA object
save_file = fullfile(psa.output_dir, 'psa_object.mat');
save(save_file, 'psa');
fprintf('PSA object saved to: %s\n', save_file);

% Copy script for reproducibility
copyfile([mfilename('fullpath') '.m'], psa.output_dir);

%% Plot results
[~, figs_hist] = load_and_make_unit_histograms(psa.output_dir, 'Metrics', {'br','lle'});
fig_paired_swarm = load_and_plot_lle_by_stim_period(psa.output_dir, 'transient_skip', 3, 'periods_to_plot', [0 1 1]);

% Combine into 3x4 layout
fig_combined = concatenate_figs([figs_hist, fig_paired_swarm], 'vertical', 'HideTitlesAfterFirstRow', true);
drawnow
try
    AddLetters2Plots(fig_combined, ...
        {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)'}, ...
        'FontSize', 14, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.04);
catch ME
    warning('AddLetters2Plots failed (non-critical): %s', ME.message);
end

%% Save figures
if save_figs
    save_dir = fullfile(figs_root, 'fraction_excitatory_analysis');
    save_some_figs_to_folder_2(save_dir, 'fraction_excitatory', [], {'fig', 'png'});
    fprintf('Figures saved to %s\n', save_dir);

    data_source_file = fullfile(save_dir, 'data_source.txt');
    fid = fopen(data_source_file, 'w');
    fprintf(fid, 'Figures generated from data in:\n%s\n', psa.output_dir);
    fclose(fid);
end

%% Summary
fprintf('\n=== Parameter Space Analysis Summary ===\n');
fprintf('Output directory: %s\n', psa.output_dir);
fprintf('Grid parameters: %s\n', strjoin(psa.grid_params, ', '));
fprintf('Levels per parameter: %d\n', psa.n_levels);
fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false), ', '));
fprintf('\nDone! Results saved to: %s\n', psa.output_dir);

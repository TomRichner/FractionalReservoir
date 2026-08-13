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

% No clear/clc/close all on standalone runs, so base-workspace settings such as
% run_mode and save_figs (set in the console before running) survive into this
% script. run_all_analyses relies on the same: it never clears the sub-scripts.

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
cfg = analysis_run_config('tau_sensitivity', run_mode);
n_levels = cfg.n_levels;
n_reps = cfg.n_reps;
fprintf('[run_tau_sensitivity_analysis] run_mode=%s, n_levels=%d, n_reps=%d, ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    run_mode, n_levels, n_reps, cfg.model.ode_solver, cfg.model.fs, ...
    cfg.model.T_range(1), cfg.model.T_range(2));

% Which network to simulate, and with which model class. run_all_analyses sets
% both from srnn_param_preset; standalone runs get SRNNModel2 class defaults.
if exist('master_model_overrides', 'var')
    preset_defaults = master_model_overrides;
else
    preset_defaults = srnn_param_preset('default');
end
if exist('master_model_class', 'var')
    model_class = master_model_class;
else
    model_class = 'SRNNModel2';
end

note = 'tau_timescales';

% Condition: SFA + STD, taken from the shared four-regime set rather than
% respelled here, so this sweep uses exactly the regime the other analyses call
% 'sfa_and_std'. Whichever class is in play, it gives E three SFA timescales --
% which is what n_elements = 3 below is coupled to.
if exist('master_conditions', 'var')
    all_conditions = master_conditions;
else
    all_conditions = srnn_adaptation_conditions(model_class);
end
is_sfa_and_std = cellfun(@(c) strcmp(c.name, 'sfa_and_std'), all_conditions);
condition = all_conditions(is_sfa_and_std);

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
psa_tau_a.model_class = model_class;
psa_tau_a.integer_params = {'n', 'indegree'};
if exist('master_output_dir', 'var')
    psa_tau_a.output_dir = master_output_dir;
end
% Parameter preset first, run_mode timings second, so run_mode keeps final say
% over ode_solver / fs / T_range / lya_T_interval.
psa_tau_a.model_defaults = merge_struct(preset_defaults, cfg.model);

psa_tau_a.set_conditions(condition);

% tau_a_E is a vector of length 3, logspaced from 0.25 to max, and we sweep the
% max (last element) from 5 to 60. It is a real property on SRNNModel2 and a
% scalar-row alias onto tau_a{1} on SRNNCellTypePairs, so the same axis name
% works for both. n_elements stays coupled BY HAND to the condition's three SFA
% timescales above.
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

% NOTE: tau_b_E_rec sweep skipped — section 2 is commented out below.
%{
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
psa_tau_b.model_defaults = merge_struct(preset_defaults, cfg.model);

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
%}
% end of skipped tau_b_E_rec sweep

%% Summary
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity Analysis Complete ===\n');
fprintf('tau_a_E results: %s\n', psa_tau_a.output_dir);
% fprintf('tau_b_E_rec results: %s\n', psa_tau_b.output_dir);  % skipped
fprintf('========================================\n');

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
switch run_mode
    case 'fast'
        % Fast iteration: fewer levels/reps, half the sample rate, short time
        % range. fs=200 keeps Benettin's lya_dt/dt guard satisfied (4>=3);
        % T_range=[0,20] with an explicit 10 s LLE window [10,20]. NOTE: this
        % window is far shorter than the swept tau (up to 60 s), so long-tau
        % LLE reflects a transient — accepted for fast qualitative iteration.
        n_levels = 7;  n_reps = 7;  ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 20];  lya_T_interval_mode = [10, 20];
    case 'medium'
        % Medium: roughly halfway between fast and production. ode45 at fs=200,
        % 11 levels x 25 reps. T_range=[0,50] with the auto LLE window
        % ([] -> [15,50]) so the analysis window is longer relative to the swept
        % tau (up to 60 s) than the [10,20] used by the other medium sub-scripts.
        n_levels = 11; n_reps = 25; ode_solver_mode = @ode45;
        fs_mode = 200;  T_range_mode = [0, 50];  lya_T_interval_mode = [];
    case 'production'
        n_levels = 25; n_reps = 50; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 50];  lya_T_interval_mode = [];
    otherwise, error('run_tau_sensitivity_analysis:badMode', ...
        'Unknown run_mode ''%s'' (expected ''fast'', ''medium'', or ''production'').', run_mode);
end
fprintf('[run_tau_sensitivity_analysis] run_mode=%s, n_levels=%d, n_reps=%d, ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    run_mode, n_levels, n_reps, func2str(ode_solver_mode), fs_mode, T_range_mode(1), T_range_mode(2));
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
psa_tau_a.model_defaults.ode_solver = ode_solver_mode;  % fast=ode_rk4, production=ode45
psa_tau_a.model_defaults.fs = fs_mode;                  % fast=200 (default 400)
psa_tau_a.model_defaults.T_range = T_range_mode;        % fast=[0,20] (default [0,50])
if ~isempty(lya_T_interval_mode)
    psa_tau_a.model_defaults.lya_T_interval = lya_T_interval_mode;  % fast: LLE window [10,20]
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
psa_tau_b.model_defaults.ode_solver = ode_solver_mode;  % fast=ode_rk4, production=ode45
psa_tau_b.model_defaults.fs = fs_mode;                  % fast=200 (default 400)
psa_tau_b.model_defaults.T_range = T_range_mode;        % fast=[0,20] (default [0,50])
if ~isempty(lya_T_interval_mode)
    psa_tau_b.model_defaults.lya_T_interval = lya_T_interval_mode;  % fast: LLE window [10,20]
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
%}
% end of skipped tau_b_E_rec sweep

%% Summary
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity Analysis Complete ===\n');
fprintf('tau_a_E results: %s\n', psa_tau_a.output_dir);
% fprintf('tau_b_E_rec results: %s\n', psa_tau_b.output_dir);  % skipped
fprintf('========================================\n');

%% replot_tau_sensitivity.m
% Reload saved ParamSpaceAnalysis2 objects and regenerate tau sensitivity
% figures with wider y-axis limits (lambda_1: [-1.5, 1.5]).
%
% No simulations are re-run — only plotting from saved results.
% Figures are saved to a NEW timestamped subfolder to avoid overwriting
% the originals.

%% Setup
setup_paths();
close all;

%% Configuration
data_root = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'data', 'param_space', 'run_all_mar_02_26_15_30');

% New y-axis (histogram) range for LLE plots
lle_hist_range = [-1.5, 1.5];

% Create a new output folder for the replots (never overwrites originals)
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
replot_dir = fullfile(data_root, sprintf('replot_%s', dt_str));
if ~exist(replot_dir, 'dir')
    mkdir(replot_dir);
end
fprintf('Replot output directory:\n  %s\n\n', replot_dir);

%% 1. Reload and replot tau_a_E sensitivity
tau_a_dir = fullfile(data_root, ...
    'tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_7_mar_02_26_15_30');

fprintf('Loading tau_a_E results from:\n  %s\n', tau_a_dir);
S = load(fullfile(tau_a_dir, 'psa_object.mat'));
psa_tau_a = S.psa_tau_a;

% Redirect output_dir so internal saves go to the new folder
psa_tau_a.output_dir = replot_dir;

% Plot with new range
psa_tau_a.plot_sensitivity('metric', 'LLE', 'hist_range', lle_hist_range);
psa_tau_a.plot_sensitivity('metric', 'mean_rate');

% Save figures using the project's standard export utility
fig_dir = fullfile(replot_dir, 'figures');
save_some_figs_to_folder_2(fig_dir, 'tau_sensitivity_tau_a', [], {'fig', 'png'});

fprintf('\nDone. Figures saved to:\n  %s\n', fig_dir);

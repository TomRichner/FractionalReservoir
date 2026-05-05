%% replot_sensitivity.m
% Reload saved ParamSpaceAnalysis2 objects from a 1D sensitivity run and
% regenerate figures with a wider LLE y-axis range. mean_rate is replotted
% with its default range.
%
% No simulations are re-run -- only plotting from saved results.
% Figures are saved to a NEW timestamped subfolder so originals are not
% overwritten.

%% Setup
setup_paths();
close all;

%% Configuration
data_root = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'data', 'param_space', 'run_all_mar_02_26_17_12');

% Subfolders produced by run_sensitivity_analysis.m (one per swept param)
sensitivity_subdirs = {
    '1D_sensitivity_sensitivity_n_nLevs_25_mar_02_26_20_03',   'n';
    '1D_sensitivity_sensitivity_f_nLevs_25_mar_02_26_23_32',   'f';
    '1D_sensitivity_sensitivity_S_c_nLevs_25_mar_03_26_04_32', 'S_c';
};

% New y-axis (histogram) range for LLE plots
lle_hist_range = [-1, 1];

% Create a new output folder for the replots (never overwrites originals)
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
replot_dir = fullfile(data_root, sprintf('replot_sensitivity_%s', dt_str));
if ~exist(replot_dir, 'dir')
    mkdir(replot_dir);
end
fprintf('Replot output directory:\n  %s\n\n', replot_dir);

%% Loop over sensitivity runs
for k = 1:size(sensitivity_subdirs, 1)
    sub = sensitivity_subdirs{k, 1};
    param_name = sensitivity_subdirs{k, 2};
    src_dir = fullfile(data_root, sub);

    fprintf('Loading %s sensitivity from:\n  %s\n', param_name, src_dir);
    S = load(fullfile(src_dir, 'psa_object.mat'));
    psa = S.psa;

    % Redirect output_dir so any internal saves go to the new folder
    psa.output_dir = replot_dir;

    % Replot with new LLE range; mean_rate uses default range
    psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_hist_range);
    psa.plot_sensitivity('metric', 'mean_rate');

    % Save figures
    fig_dir = fullfile(replot_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, ...
        sprintf('sensitivity_%s', param_name), [], {'fig', 'png'});

    close all;
    fprintf('  -> figures saved to %s\n\n', fig_dir);
end

fprintf('Done. All replots in:\n  %s\n', replot_dir);

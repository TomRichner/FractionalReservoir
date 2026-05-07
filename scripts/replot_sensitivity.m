%% replot_sensitivity.m
% Reload saved ParamSpaceAnalysis2 objects from a 1D sensitivity run and
% regenerate figures with a wider LLE y-axis range. mean_rate is replotted
% with its default range.
%
% Master-driveable:
%   master_replot_data_root - override the data_root below
%   master_lle_hist_range   - override the lle_hist_range below
%
% No simulations are re-run -- only plotting from saved results.
% Figures are saved to a NEW timestamped subfolder so originals are not
% overwritten.

%% Setup
if ~exist('master_replot_data_root', 'var')
    close all;
end
setup_paths();

%% Configuration
if exist('master_replot_data_root', 'var')
    data_root = master_replot_data_root;
else
    data_root = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'data', 'param_space', 'run_all_mar_02_26_17_12');
end

% New y-axis (histogram) range for LLE plots
if exist('master_lle_hist_range', 'var')
    lle_hist_range = master_lle_hist_range;
else
    lle_hist_range = [-1, 1];
end

% Locate 1D sensitivity subfolders inside data_root
sens_listing = dir(fullfile(data_root, '1D_sensitivity_*'));
sens_listing = sens_listing([sens_listing.isdir]);
if isempty(sens_listing)
    error('replot_sensitivity:NotFound', ...
        'No 1D_sensitivity_* subfolder found in:\n  %s', data_root);
end

% Create a new output folder for the replots (never overwrites originals)
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
replot_dir = fullfile(data_root, sprintf('replot_sensitivity_%s', dt_str));
if ~exist(replot_dir, 'dir')
    mkdir(replot_dir);
end
fprintf('Replot output directory:\n  %s\n\n', replot_dir);

%% Loop over sensitivity runs
for k = 1:length(sens_listing)
    src_dir = fullfile(sens_listing(k).folder, sens_listing(k).name);

    psa_file = fullfile(src_dir, 'psa_object.mat');
    if ~exist(psa_file, 'file')
        warning('replot_sensitivity:NoPsaFile', ...
            'Skipping (no psa_object.mat): %s', src_dir);
        continue;
    end

    fprintf('Loading sensitivity from:\n  %s\n', src_dir);
    S = load(psa_file);
    psa = S.psa;

    % Infer swept parameter name from the loaded PSA
    swept = setdiff(psa.grid_params, {'reps'}, 'stable');
    if isempty(swept)
        warning('replot_sensitivity:NoSweptParam', ...
            'Skipping (no swept param found): %s', src_dir);
        continue;
    end
    param_name = swept{1};

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

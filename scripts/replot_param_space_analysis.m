%% replot_param_space_analysis.m
% Reload a saved ParamSpaceAnalysis2 object from a parameter-space run and
% regenerate figures without re-running simulations.
%
% Master-driveable: if `master_replot_data_root` is set in the workspace,
% that root is used; otherwise the hardcoded `data_root` below is used.
%
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

% Locate the param_space subfolder inside data_root.
% run_param_space_analysis2 produces folders named param_space_<note>_nLevs_*_<dt>
ps_listing = dir(fullfile(data_root, 'param_space_*'));
ps_listing = ps_listing([ps_listing.isdir]);
if isempty(ps_listing)
    error('replot_param_space_analysis:NotFound', ...
        'No param_space_* subfolder found in:\n  %s', data_root);
end

% Create a new output folder for the replots (never overwrites originals)
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
replot_dir = fullfile(data_root, sprintf('replot_param_space_%s', dt_str));
if ~exist(replot_dir, 'dir')
    mkdir(replot_dir);
end
fprintf('Replot output directory:\n  %s\n\n', replot_dir);

%% Loop over param_space runs (typically just one per run_all)
for k = 1:length(ps_listing)
    src_dir = fullfile(ps_listing(k).folder, ps_listing(k).name);

    fprintf('Loading param_space results from:\n  %s\n', src_dir);
    psa_file = fullfile(src_dir, 'psa_object.mat');
    if ~exist(psa_file, 'file')
        warning('replot_param_space_analysis:NoPsaFile', ...
            'Skipping (no psa_object.mat): %s', src_dir);
        continue;
    end
    S = load(psa_file);
    psa = S.psa;

    % Redirect output_dir so any internal saves go to the new folder
    psa.output_dir = replot_dir;

    % Replot
    psa.plot('metric', 'LLE');
    psa.plot('metric', 'mean_rate');

    % Save figures
    fig_dir = fullfile(replot_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'param_space', [], {'fig', 'png'});

    close all;
    fprintf('  -> figures saved to %s\n\n', fig_dir);
end

fprintf('Done. All replots in:\n  %s\n', replot_dir);

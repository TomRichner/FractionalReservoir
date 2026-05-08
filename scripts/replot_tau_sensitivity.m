function replot_dir = replot_tau_sensitivity(data_root, lle_hist_range)
% REPLOT_TAU_SENSITIVITY Reload tau sensitivity PSA objects and regenerate figures.
%
%   replot_dir = REPLOT_TAU_SENSITIVITY(data_root)
%   replot_dir = REPLOT_TAU_SENSITIVITY(data_root, lle_hist_range)
%
% Inputs:
%   data_root      : path to a run_all_<dt> folder containing tau_sensitivity_*
%                    subdirectories produced by run_tau_sensitivity_analysis.m
%   lle_hist_range : optional 2-element [ymin ymax] for the LLE histogram
%                    (default [-1.5, 1.5])
%
% Output:
%   replot_dir     : path to the new replot_tau_<dt> folder under data_root
%
% No simulations are re-run -- only plotting from saved results.

    if nargin < 2 || isempty(lle_hist_range)
        lle_hist_range = [-1.5, 1.5];
    end

    setup_paths();

    tau_listing = dir(fullfile(data_root, 'tau_sensitivity_*'));
    tau_listing = tau_listing([tau_listing.isdir]);
    if isempty(tau_listing)
        error('replot_tau_sensitivity:NotFound', ...
            'No tau_sensitivity_* subfolder found in:\n  %s', data_root);
    end

    dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
    replot_dir = fullfile(data_root, sprintf('replot_tau_%s', dt_str));
    if ~exist(replot_dir, 'dir')
        mkdir(replot_dir);
    end
    fprintf('Replot output directory:\n  %s\n\n', replot_dir);

    for k = 1:length(tau_listing)
        src_dir = fullfile(tau_listing(k).folder, tau_listing(k).name);

        psa_file = fullfile(src_dir, 'psa_object.mat');
        if ~exist(psa_file, 'file')
            warning('replot_tau_sensitivity:NoPsaFile', ...
                'Skipping (no psa_object.mat): %s', src_dir);
            continue;
        end

        fprintf('Loading tau sensitivity from:\n  %s\n', src_dir);
        S = load(psa_file);
        fns = fieldnames(S);
        psa = S.(fns{1});  % saved as psa_tau_a / psa_tau_b

        swept = setdiff(psa.grid_params, {'reps'}, 'stable');
        if isempty(swept)
            param_name = 'tau';
        else
            param_name = swept{1};
        end

        psa.output_dir = replot_dir;

        psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_hist_range);
        psa.plot_sensitivity('metric', 'mean_rate');

        fig_dir = fullfile(replot_dir, 'figures');
        save_some_figs_to_folder_2(fig_dir, ...
            sprintf('tau_sensitivity_%s', param_name), [], {'fig', 'png'});

        close all;
        fprintf('  -> figures saved to %s\n\n', fig_dir);
    end

    fprintf('Done. All replots in:\n  %s\n', replot_dir);
end

function replot_dir = replot_param_space_analysis(data_root)
% REPLOT_PARAM_SPACE_ANALYSIS Reload param-space PSA object and regenerate figures.
%
%   replot_dir = REPLOT_PARAM_SPACE_ANALYSIS(data_root)
%
% Inputs:
%   data_root  : path to a run_all_<dt> folder containing a param_space_*
%                subdirectory produced by run_param_space_analysis.m
%
% Output:
%   replot_dir : path to the new replot_param_space_<dt> folder under data_root
%
% No simulations are re-run -- only plotting from saved results.

    setup_paths();

    ps_listing = dir(fullfile(data_root, 'param_space_*'));
    ps_listing = ps_listing([ps_listing.isdir]);
    if isempty(ps_listing)
        error('replot_param_space_analysis:NotFound', ...
            'No param_space_* subfolder found in:\n  %s', data_root);
    end

    dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
    replot_dir = fullfile(data_root, sprintf('replot_param_space_%s', dt_str));
    if ~exist(replot_dir, 'dir')
        mkdir(replot_dir);
    end
    fprintf('Replot output directory:\n  %s\n\n', replot_dir);

    for k = 1:length(ps_listing)
        src_dir = fullfile(ps_listing(k).folder, ps_listing(k).name);

        psa_file = fullfile(src_dir, 'psa_object.mat');
        if ~exist(psa_file, 'file')
            warning('replot_param_space_analysis:NoPsaFile', ...
                'Skipping (no psa_object.mat): %s', src_dir);
            continue;
        end

        fprintf('Loading param_space results from:\n  %s\n', src_dir);
        psa = ParamSpaceAnalysis2.from_dir(src_dir);

        psa.output_dir = replot_dir;

        psa.plot('metric', 'LLE');
        psa.plot('metric', 'mean_rate');

        fig_dir = fullfile(replot_dir, 'figures');
        save_some_figs_to_folder_2(fig_dir, 'param_space', [], {'fig', 'png'});

        close all;
        fprintf('  -> figures saved to %s\n\n', fig_dir);
    end

    fprintf('Done. All replots in:\n  %s\n', replot_dir);
end

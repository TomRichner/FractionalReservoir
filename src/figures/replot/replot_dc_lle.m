function replot_dir = replot_dc_lle(data_root)
% REPLOT_DC_LLE Reload dc_lle_results.mat and regenerate the DC-vs-LLE figure.
%
%   replot_dir = REPLOT_DC_LLE(data_root)
%
% Inputs:
%   data_root : path to a run_all_<dt> folder containing one or more
%               dc_lle_nSeeds_* subdirectories produced by run_dc_lle_analysis.m
%
% Output:
%   replot_dir : path to the new replot_dc_lle_<dt> folder under data_root
%
% No simulations are re-run -- the confplot figure is rebuilt from the saved
% dc_lle_results.mat alone (mean local LLE vs DC level, +/- std across seeds).

    setup_paths();

    listing = dir(fullfile(data_root, 'dc_lle_nSeeds_*'));
    listing = listing([listing.isdir]);
    if isempty(listing)
        error('replot_dc_lle:NotFound', ...
            'No dc_lle_nSeeds_* subfolder found in:\n  %s', data_root);
    end

    dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
    replot_dir = fullfile(data_root, sprintf('replot_dc_lle_%s', dt_str));
    if ~exist(replot_dir, 'dir')
        mkdir(replot_dir);
    end
    fprintf('Replot output directory:\n  %s\n\n', replot_dir);

    fig_dir = fullfile(replot_dir, 'figures');

    for k = 1:length(listing)
        src_dir = fullfile(listing(k).folder, listing(k).name);

        res_file = fullfile(src_dir, 'dc_lle_results.mat');
        if ~exist(res_file, 'file')
            warning('replot_dc_lle:NoResults', ...
                'Skipping (no dc_lle_results.mat): %s', src_dir);
            continue;
        end

        fprintf('Loading DC LLE results from:\n  %s\n', src_dir);
        S = load(res_file);
        r = S.dc_lle_results;

        figure('Name', 'DC LLE: local Lyapunov vs DC level (across seeds)');
        confplot(r.dc_levels, r.level_mean, r.level_std, r.level_std, ...
            [0 0 0.8; 0.7 0.8 1.0]);
        yline(0, '--k', 'edge of chaos', 'HandleVisibility', 'off');
        xlabel('DC level (input units)');
        ylabel('mean local \lambda_1  (\pm std across seeds)');
        title(sprintf('Largest local Lyapunov exponent vs DC stim level  (n_{seeds}=%d)', ...
            r.n_seeds));
        grid off;

        save_some_figs_to_folder_2(fig_dir, ...
            sprintf('dc_lle_nSeeds_%d', r.n_seeds), [], {'fig', 'png'});

        close all;
        fprintf('  -> figures saved to %s\n\n', fig_dir);
    end

    fprintf('Done. All replots in:\n  %s\n', replot_dir);
end

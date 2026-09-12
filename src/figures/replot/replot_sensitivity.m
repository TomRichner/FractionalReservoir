function replot_dir = replot_sensitivity(data_root, lle_hist_range, n_bins, rate_n_bins)
% REPLOT_SENSITIVITY Reload 1D sensitivity PSA objects and regenerate figures.
%
%   replot_dir = REPLOT_SENSITIVITY(data_root)
%   replot_dir = REPLOT_SENSITIVITY(data_root, lle_hist_range)
%   replot_dir = REPLOT_SENSITIVITY(data_root, lle_hist_range, n_bins)
%   replot_dir = REPLOT_SENSITIVITY(data_root, lle_hist_range, n_bins, rate_n_bins)
%
% Inputs:
%   data_root      : path to a run_all_<dt> folder containing 1D_sensitivity_*
%                    subdirectories produced by run_sensitivity_analysis.m
%   lle_hist_range : optional 2-element [ymin ymax] for the LLE histogram
%                    (default [-1, 1]); mean_rate uses its plot_sensitivity
%                    default range of [0, 1]
%   n_bins         : optional number of LLE histogram bins. When omitted, uses
%                    plot_sensitivity's default (35).
%   rate_n_bins    : optional number of mean_rate histogram bins. When omitted,
%                    uses plot_sensitivity's default (35). Pass the same value
%                    as n_bins to give both metrics the same vertical
%                    resolution.
%
% Note: n_bins counts linspace EDGES, so the plotted image has n_bins+1 rows
% (n_bins-1 interior bins plus the two -inf/+inf overflow bins).
%
% Output:
%   replot_dir     : path to the new replot_sensitivity_<dt> folder under data_root
%
% No simulations are re-run -- only plotting from saved results.

    if nargin < 2 || isempty(lle_hist_range)
        lle_hist_range = [-2,2];
    end
    if nargin < 3
        n_bins = [];
    end
    if nargin < 4
        rate_n_bins = [];
    end

    setup_paths();

    sens_listing = dir(fullfile(data_root, '1D_sensitivity_*'));
    sens_listing = sens_listing([sens_listing.isdir]);
    if isempty(sens_listing)
        error('replot_sensitivity:NotFound', ...
            'No 1D_sensitivity_* subfolder found in:\n  %s', data_root);
    end

    dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
    replot_dir = fullfile(data_root, sprintf('replot_sensitivity_%s', dt_str));
    if ~exist(replot_dir, 'dir')
        mkdir(replot_dir);
    end
    fprintf('Replot output directory:\n  %s\n\n', replot_dir);

    % Per-metric plot_sensitivity args (constant across params); include n_bins
    % only when the caller provided it, else plot_sensitivity uses its default.
    % One plot_sensitivity call per registry measure with in_sheets set
    % (sweep_metrics). LLE keeps the caller's range; every other measure
    % takes the registry's sens_range. n_bins applies to LLE, rate_n_bins to
    % the rest (historical names; both default to plot_sensitivity's own).
    specs = sweep_metrics();
    specs = specs([specs.in_sheets]);
    metric_args = cell(1, numel(specs));
    for mi = 1:numel(specs)
        args = {'metric', specs(mi).field};
        if strcmp(specs(mi).field, 'LLE')
            args = [args, {'hist_range', lle_hist_range}]; %#ok<AGROW>
            if ~isempty(n_bins); args = [args, {'n_bins', n_bins}]; end %#ok<AGROW>
        elseif ~isempty(rate_n_bins)
            args = [args, {'n_bins', rate_n_bins}]; %#ok<AGROW>
        end
        metric_args{mi} = args;
    end

    for k = 1:length(sens_listing)
        src_dir = fullfile(sens_listing(k).folder, sens_listing(k).name);

        psa_file = fullfile(src_dir, 'psa_object.mat');
        if ~exist(psa_file, 'file')
            warning('replot_sensitivity:NoPsaFile', ...
                'Skipping (no psa_object.mat): %s', src_dir);
            continue;
        end

        fprintf('Loading sensitivity from:\n  %s\n', src_dir);
        psa = ParamSpaceAnalysis2.from_dir(src_dir);

        swept = setdiff(psa.grid_params, {'reps'}, 'stable');
        if isempty(swept)
            warning('replot_sensitivity:NoSweptParam', ...
                'Skipping (no swept param found): %s', src_dir);
            continue;
        end
        param_name = swept{1};

        psa.output_dir = replot_dir;

        for mi = 1:numel(metric_args)
            psa.plot_sensitivity(metric_args{mi}{:});
        end

        fig_dir = fullfile(replot_dir, 'figures');
        save_some_figs_to_folder_2(fig_dir, ...
            sprintf('sensitivity_%s', param_name), [], {'fig', 'png'});

        close all;
        fprintf('  -> figures saved to %s\n\n', fig_dir);
    end

    fprintf('Done. All replots in:\n  %s\n', replot_dir);
end

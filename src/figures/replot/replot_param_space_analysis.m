function replot_dir = replot_param_space_analysis(data_root, opts)
% REPLOT_PARAM_SPACE_ANALYSIS Reload param-space PSA object and regenerate figures.
%
%   replot_dir = REPLOT_PARAM_SPACE_ANALYSIS(data_root)
%   replot_dir = REPLOT_PARAM_SPACE_ANALYSIS(data_root, 'conditions', names, 'verbose', level)
%
% Inputs:
%   data_root    : path to a run_all_<dt> folder containing a param_space_*
%                  subdirectory produced by run_param_space_analysis.m
%   'conditions' : cell of condition names to draw as columns, LEFT TO RIGHT
%                  IN THE ORDER GIVEN (FD-006). Empty (the default) draws every
%                  condition the sample holds, in its own order, as before the
%                  option existed. Resolved by select_conditions, so a name the
%                  run does not hold is an error. Applied to the LOADED PSA
%                  object's conditions list only (psa.plot draws
%                  psa.conditions in order; results are keyed by name).
%   'verbose'    : 'verbose' | 'minimal' (default) | 'near-none'
%
% Output:
%   replot_dir : path to the new replot_param_space_<dt> folder under data_root
%
% No simulations are re-run -- only plotting from saved results. Only the
% figures this function draws are saved and closed (it used to save every open
% figure and then 'close all', which in a batch destroyed the figures of the
% entries run before it).

arguments
    data_root (1,:) char
    opts.conditions = {}
    opts.verbose (1,:) char = 'minimal'
end

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
    vprintf(opts.verbose, 'verbose', 'Replot output directory:\n  %s\n\n', replot_dir);

    for k = 1:length(ps_listing)
        src_dir = fullfile(ps_listing(k).folder, ps_listing(k).name);

        psa_file = fullfile(src_dir, 'psa_object.mat');
        if ~exist(psa_file, 'file')
            warning('replot_param_space_analysis:NoPsaFile', ...
                'Skipping (no psa_object.mat): %s', src_dir);
            continue;
        end

        vprintf(opts.verbose, 'verbose', 'Loading param_space results from:\n  %s\n', src_dir);
        psa = ParamSpaceAnalysis2.from_dir(src_dir);

        % The columns, in the requested order (empty = all, unchanged).
        held = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
        [~, idx] = select_conditions(held, opts.conditions);
        psa.conditions = psa.conditions(idx);

        psa.output_dir = replot_dir;

        before = findall(groot, 'Type', 'figure');
        specs = sweep_metrics();
        for spec = specs([specs.in_sheets])
            psa.plot('metric', spec.field);
        end
        after = findall(groot, 'Type', 'figure');
        made = after(~ismember(after, before));
        if isempty(made)
            % Never pass an empty list on: save_some_figs_to_folder_2 reads []
            % as "every open figure".
            warning('replot_param_space_analysis:NoFigures', 'psa.plot drew nothing for %s.', src_dir);
            continue;
        end

        fig_dir = fullfile(replot_dir, 'figures');
        save_some_figs_to_folder_2(fig_dir, 'param_space', sort([made.Number]), {'fig', 'png'});

        close(made);
        vprintf(opts.verbose, 'verbose', '  -> figures saved to %s\n\n', fig_dir);
    end

    vprintf(opts.verbose, 'verbose', 'Done. All replots in:\n  %s\n', replot_dir);
end

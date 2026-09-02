function [fig1, fig2] = replot_memory_capacity(mat_file, out_dir)
% REPLOT_MEMORY_CAPACITY Reload a saved MC results .mat and regenerate figures.
%
%   replot_memory_capacity(mat_file)
%   replot_memory_capacity(mat_file, out_dir)
%
% Convenience wrapper: loads the `results_all` struct saved by
% run_memory_capacity.m and calls plot_memory_capacity -- no simulation is
% re-run, so you can iterate on figure styling from a finished run.
%
% Inputs:
%   mat_file : path to a <run_tag>_results.mat produced by run_memory_capacity.
%   out_dir  : optional save folder for the figures (default: the mat_file's
%              folder). Pass '' to display only, without saving.

    % Make plot_memory_capacity resolvable regardless of cwd.
    setup_paths();

    if ~isfile(mat_file)
        error('replot_memory_capacity:NotFound', 'File not found:\n  %s', mat_file);
    end
    if nargin < 2; out_dir = fileparts(mat_file); end

    S = load(mat_file);
    if ~isfield(S, 'results_all')
        error('replot_memory_capacity:BadFile', ...
            '%s does not contain ''results_all''.', mat_file);
    end

    [fig1, fig2] = plot_memory_capacity(S.results_all, out_dir);
end

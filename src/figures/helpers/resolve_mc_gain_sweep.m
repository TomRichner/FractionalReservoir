function D = resolve_mc_gain_sweep(src)
% RESOLVE_MC_GAIN_SWEEP The gain-sweep result a block C figure was asked to draw.
%
%   D = RESOLVE_MC_GAIN_SWEEP('data', cfg.data, 'data_file', cfg.data_file, 'run_dir', cfg.run_dir)
%   D = RESOLVE_MC_GAIN_SWEEP('run_dir', run_dir, 'grid', 'selfcoupling')
%   D = RESOLVE_MC_GAIN_SWEEP('run_dir', run_dir, 'grid', 'bias', 'variant', 'EI3to2')
%
% The three ways every block C figure takes its result, in this order:
%   data       a result in memory: read_mc_gain_sweep's output, or anything
%              it reads (a table, a struct with T)
%   data_file  a results .mat, named explicitly
%   run_dir    a run directory. WHICH RESULT in it is named by two tags, as
%              the bundle's mc_gain_sweep stage names its runs
%              (mc_gain_sweep_stage_plan):
%                grid     the grid's tag ('bias', 'selfcoupling'); empty
%                         for a run directory with one unnamed grid, and
%                         for run_mc_gain_sweep's own output directory
%                variant  a variant network of that grid ('EI3to2'); empty
%                         for the grid itself
%              The result is <run_dir>/<folder>/mc_gain_sweep/
%              mc_gain_sweep_results.mat, with the folder formed by
%              mc_gain_sweep_folder, the function the stage writes by.
% Nothing outside the run directory is read, and a run that does not hold
% the result asked for is an error, not a search elsewhere: not in another
% grid's folder, and not outside the run (resolve_data_file).
%
% See also: read_mc_gain_sweep, resolve_data_file, mc_gain_sweep_folder

arguments
    src.data                 = []
    src.data_file (1,:) char = ''
    src.run_dir   (1,:) char = ''
    src.grid      (1,:) char = ''
    src.variant   (1,:) char = ''
end

if ~isempty(src.data)
    D = read_mc_gain_sweep(src.data);
    return
end
stage_dir = fullfile(src.run_dir, mc_gain_sweep_folder(src.grid, src.variant), 'mc_gain_sweep');
file = resolve_data_file(src.data_file, src.run_dir, {stage_dir}, 'mc_gain_sweep_results.mat', ...
    'Run run_mc_gain_sweep first, or pass data_file');
D = read_mc_gain_sweep(file);
end

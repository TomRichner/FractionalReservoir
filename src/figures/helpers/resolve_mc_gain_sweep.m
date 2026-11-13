function D = resolve_mc_gain_sweep(src)
% RESOLVE_MC_GAIN_SWEEP The gain-by-bias result a block C figure was asked to draw.
%
%   D = RESOLVE_MC_GAIN_SWEEP('data', cfg.data, 'data_file', cfg.data_file, 'run_dir', cfg.run_dir)
%   D = RESOLVE_MC_GAIN_SWEEP('run_dir', run_dir, 'variant', 'EI3to2')
%
% The three ways every block C figure takes its result, in this order:
%   data       a result in memory: read_mc_gain_sweep's output, or anything
%              it reads (a table, a struct with T)
%   data_file  a results .mat, named explicitly
%   run_dir    a run directory: <run_dir>/mc_gain_sweep/
%              mc_gain_sweep_results.mat, where run_mc_gain_sweep writes it.
%              With variant, the result of that variant network instead,
%              <run_dir>/mc_gain_sweep_<variant>/mc_gain_sweep/ (a bundle's
%              mc_gain_sweep stage; mc_gain_sweep_variant_folder).
% Nothing outside the run directory is read, and a run that holds no result
% is an error, not a search elsewhere (resolve_data_file).
%
% See also: read_mc_gain_sweep, resolve_data_file, mc_gain_sweep_variant_folder

arguments
    src.data                 = []
    src.data_file (1,:) char = ''
    src.run_dir   (1,:) char = ''
    src.variant   (1,:) char = ''
end

if ~isempty(src.data)
    D = read_mc_gain_sweep(src.data);
    return
end
stage_dir = src.run_dir;
if ~isempty(src.variant)
    stage_dir = fullfile(src.run_dir, mc_gain_sweep_variant_folder(src.variant));
end
file = resolve_data_file(src.data_file, src.run_dir, ...
    {fullfile(stage_dir, 'mc_gain_sweep')}, 'mc_gain_sweep_results.mat', ...
    'Run run_mc_gain_sweep first, or pass data_file');
D = read_mc_gain_sweep(file);
end

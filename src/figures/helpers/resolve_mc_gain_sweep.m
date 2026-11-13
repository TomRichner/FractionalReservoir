function D = resolve_mc_gain_sweep(data_file, run_dir)
% RESOLVE_MC_GAIN_SWEEP Find and read the block C result for a figure.
%
%   D = RESOLVE_MC_GAIN_SWEEP(cfg.data_file, cfg.run_dir)
%
% An explicit data_file wins; otherwise the newest mc_gain_sweep*.mat in
% <run_dir>/mc_gain_sweep, and nothing outside the run directory
% (resolve_data_file). The table is read by read_mc_gain_sweep, the one
% reader. The folder and file pattern are this file's guess at
% run_mc_gain_sweep's layout (Plan C section 2), written once, here.
%
% See also: read_mc_gain_sweep, resolve_data_file, fig_mc_gain_bias_maps

f = resolve_data_file(char(data_file), char(run_dir), ...
    {fullfile(char(run_dir), 'mc_gain_sweep')}, 'mc_gain_sweep*.mat', ...
    'Run run_mc_gain_sweep first, or pass data_file');
D = read_mc_gain_sweep(f);
end

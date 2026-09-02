%% combine_two_runs.m
% Pool two run_all_<dt> runs into combined sensitivity/tau/param-space (and DC)
% figures via combine_runs.
%
% Runs pooled here:
%   run_all_jul_05_26_22_52  (medium; no run_manifest -> seed offset 0)
%   run_all_jul_06_26_12_47  (medium; has run_manifest + DC)
%
% NOTE ON SEED OFFSETS: combine_runs refuses to pool runs that used the SAME
% network seed offset (same config + same offset => identical networks, i.e.
% duplicate data). jul_05 predates the offset feature (offset 0). If jul_06 also
% used offset 0, you'll get a 'runs share a network_seed_offset' error -- in that
% case generate an independent run first:  run_index = 1; run_all_analyses  (->
% offset 1e6), then pool THAT with jul_05. combine_runs prints each run's offset.
%
% DC won't pool: only jul_06 has a dc_lle_* folder, so DC plots from it alone.

clear; clc; close all;

%% Paths
setup_paths();
project_root = fileparts(which('setup_paths'));

data_root = fullfile(project_root, 'data', 'param_space');
run_dirs = { ...
    fullfile(data_root, 'run_all_jul_05_26_22_52'), ...
    fullfile(data_root, 'run_all_jul_06_26_12_47') };

for i = 1:numel(run_dirs)
    if ~isfolder(run_dirs{i})
        error('combine_two_runs:MissingRun', 'Run folder not found:\n  %s', run_dirs{i});
    end
end

%% Combine (writes to data/param_space/combined_runs_<dt>/figures/)
combined_dir = combine_runs(run_dirs);

fprintf('\nCombined figures written to:\n  %s\n', combined_dir);

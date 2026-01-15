function export_SRNN_to_CSCS_format(result, subject_name, output_dir, sim_config)
%EXPORT_SRNN_TO_CSCS_FORMAT Export SRNN simulation to cscs_dynamics-compatible format
%
%   export_SRNN_to_CSCS_format(result, subject_name, output_dir, sim_config)
%
%   Converts SRNN simulation output to the format expected by the
%   cscs_dynamics VAR analysis pipeline.
%
%   Inputs:
%       result       - Struct from run_single_condition with fields:
%                      .t, .x, .fs, .params, .baseline_mask, .stim_mask
%       subject_name - Subject identifier (e.g., 'SRNN_none')
%       output_dir   - Directory to save output files
%       sim_config   - Simulation configuration struct
%
%   Output format matches CSCSreview_load_revived_parallel.m:
%       data{block}  - [time x channels] matrix
%       Fs_final     - Sampling rate
%       Files        - Cell array of channel names
%       header       - Metadata struct
%
%   See also: VAR_SRNN_comparison, SRNN_VAR_load_simulated

fprintf('  Exporting %s to CSCS format...\n', subject_name);

% Extract dendritic states (x) as the observed variable
% x is a struct with .E and .I fields, each [n_neurons x n_timepoints]
x_E = result.x.E;  % [n_E x nt]
x_I = result.x.I;  % [n_I x nt]

% Combine E and I neurons
x_all = [x_E; x_I];  % [n x nt]

% Apply neuron subset if specified
n_total = size(x_all, 1);
if ~isempty(sim_config.n_neurons_subset) && sim_config.n_neurons_subset < n_total
    rng(sim_config.rng_seeds(3));  % Reproducible subset selection
    subset_idx = randperm(n_total, sim_config.n_neurons_subset);
    subset_idx = sort(subset_idx);  % Keep order consistent
    x_all = x_all(subset_idx, :);
    n_channels = sim_config.n_neurons_subset;
    fprintf('    Using random subset of %d neurons\n', n_channels);
else
    subset_idx = 1:n_total;
    n_channels = n_total;
end

% Transpose to [time x channels] format (matching SEEG convention)
x_data = x_all';  % [nt x n_channels]

% Split into baseline and stim blocks
baseline_data = x_data(result.baseline_mask, :);
stim_data = x_data(result.stim_mask, :);

% Create data cell array (matching cscs_dynamics format)
% Block 1: baseline, Block 2: stim
data = cell(1, 2);
data{1} = baseline_data;
data{2} = stim_data;

% Sampling rate
Fs_final = result.fs;
Fs_initial = result.fs;  % Same for simulation

% Create channel names (neuron IDs)
Files = cell(n_channels, 1);
for i = 1:n_channels
    neuron_idx = subset_idx(i);
    if neuron_idx <= result.params.n_E
        Files{i} = sprintf('E_%03d.sim', neuron_idx);
    else
        Files{i} = sprintf('I_%03d.sim', neuron_idx - result.params.n_E);
    end
end

% Create header struct with metadata
header = struct();
header.simulation = true;
header.condition = subject_name;
header.n_a_E = result.params.n_a_E;
header.n_b_E = result.params.n_b_E;
header.level_of_chaos = result.params.level_of_chaos;
header.n_neurons = result.params.n;
header.n_E = result.params.n_E;
header.n_I = result.params.n_I;
header.T_baseline = result.T_baseline;
header.rng_seeds = result.params.rng_seeds;
header.subset_idx = subset_idx;
header.creation_date = datestr(now);

% Create Label array (all neurons labeled as "nSOZ" equivalent = 3)
Label = 3 * ones(1, n_channels);

% Save to file
output_file = fullfile(output_dir, sprintf('%s.mat', subject_name));
save(output_file, 'data', 'Fs_final', 'Fs_initial', 'Files', 'header', 'Label', '-v7.3');

fprintf('    Saved: %s\n', output_file);
fprintf('    Block 1 (baseline): %d samples x %d channels\n', size(data{1}, 1), size(data{1}, 2));
fprintf('    Block 2 (stim): %d samples x %d channels\n', size(data{2}, 1), size(data{2}, 2));
end

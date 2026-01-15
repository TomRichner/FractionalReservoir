function SRNN_VAR_subject_structure(conditions, sim_config, output_dir)
%SRNN_VAR_SUBJECT_STRUCTURE Generate subject structure for SRNN simulations
%
%   SRNN_VAR_subject_structure(conditions, sim_config, output_dir)
%
%   Creates a subject structure matching CSCS_subject_structure.m format.
%   Each SRNN condition becomes a pseudo-subject.
%
%   Inputs:
%       conditions - Struct array with fields: .name, .n_a_E, .n_b_E, .description
%       sim_config - Simulation configuration struct
%       output_dir - Directory containing simulation data and for saving structure
%
%   Output:
%       Saves SRNN_subjects.mat to output_dir
%
%   See also: VAR_SRNN_comparison, CSCS_subject_structure

fprintf('  Creating SRNN subject structure...\n');

n_conditions = length(conditions);

% Common parameters
deci_mode = 'none';  % Simulation already at desired rate
password = '';       % Not used for simulation
BlockLength = 16;    % minutes (placeholder)
BlockLengthFinal = 15;

% Determine number of channels
if ~isempty(sim_config.n_neurons_subset)
    n_channels = sim_config.n_neurons_subset;
else
    n_channels = sim_config.n;
end

% Initialize subjects struct
subjects = struct([]);

for idx = 1:n_conditions
    cond = conditions(idx);

    subjects(idx).SaveName = cond.name;
    subjects(idx).exclude_subject = false;
    subjects(idx).FileDir = output_dir;  % Directory with simulated data

    % Generate channel names
    Files = cell(n_channels, 1);
    for i = 1:n_channels
        if i <= sim_config.n * 0.5  % Assuming f=0.5 for E/I split
            Files{i} = sprintf('E_%03d.sim', i);
        else
            Files{i} = sprintf('I_%03d.sim', i - round(sim_config.n * 0.5));
        end
    end
    subjects(idx).Files = Files;

    % Labels: all neurons treated as "nSOZ" equivalent (Label = 3)
    subjects(idx).Label = 3 * ones(1, n_channels);

    % No bad channels in simulation
    subjects(idx).bad_Files = {};

    % Block start times not used for simulation (use empty)
    subjects(idx).BlockStartTimes = [];

    % Sampling rate
    subjects(idx).Fs_final = sim_config.fs;

    % Block indices
    % Block 1 = baseline (no stim), Block 2 = stim
    subjects(idx).baseline_block_indices = [1];
    subjects(idx).stim_block_indices = [2];

    % Polarity and time offset (not applicable)
    subjects(idx).polrev = 0;
    subjects(idx).TimeOffSet = 0;

    % Sleep state (not applicable)
    subjects(idx).isAsleep = [false, false];

    % Simulation-specific fields
    subjects(idx).is_simulation = true;
    subjects(idx).n_a_E = cond.n_a_E;
    subjects(idx).n_b_E = cond.n_b_E;
    subjects(idx).description = cond.description;
end

% Display summary
fprintf('\n=== SRNN Subject Structure ===\n');
fprintf('Total subjects: %d\n', length(subjects));
for i = 1:length(subjects)
    fprintf('  %d. %s (n_a_E=%d, n_b_E=%d)\n', i, subjects(i).SaveName, ...
        subjects(i).n_a_E, subjects(i).n_b_E);
end

% Save structure
save_filename = fullfile(output_dir, 'SRNN_subjects.mat');
save(save_filename, 'subjects', 'deci_mode', 'password', 'BlockLength', ...
    'BlockLengthFinal', 'output_dir', '-v7.3');
fprintf('\nSaved subject structure to: %s\n', save_filename);
end

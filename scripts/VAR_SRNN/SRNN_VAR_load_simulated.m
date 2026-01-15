%% SRNN_VAR_load_simulated.m
% Load and preprocess simulated SRNN data for VAR analysis
%
% This script is analogous to CSCSreview_load_revived_parallel.m but
% loads data from simulation .mat files instead of .mef files.
%
% Preprocessing steps:
%   1. Load simulated data from .mat files
%   2. Apply optional Butterworth filtering
%   3. Downsample if needed (for production runs at 5000 Hz)
%   4. Save in format compatible with ICA/VAR pipeline
%
% Prerequisites:
%   1. Run VAR_SRNN_comparison.m to generate simulation data
%
% See also: VAR_SRNN_comparison, CSCS_run_ICA, CSCS_run_VAR

close all
clc

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Input/output directories
input_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'VAR_SRNN');
output_dir = input_dir;  % Save preprocessed data in same directory

% Subject file
subject_file = fullfile(input_dir, 'SRNN_subjects.mat');

% Filtering options
%   'none'   - No filtering (keep simulation as-is)
%   'butter' - Butterworth low-pass filter (anti-aliasing before downsample)
filter_mode = 'none';

% Downsampling
% Set to target Fs if production data is at 5000 Hz
% Set to 0 or same as input Fs to skip downsampling
Fs_target = 500;  % Target sampling rate (Hz)

% Skip subjects that are already processed
skip_completed = false;

%% ========================================================================
%  LOAD SUBJECT STRUCTURE
%  ========================================================================

fprintf('=== SRNN_VAR_load_simulated ===\n');
fprintf('Loading and preprocessing simulated SRNN data\n\n');

if ~exist(subject_file, 'file')
    error('Subject file not found: %s\nRun VAR_SRNN_comparison.m first.', subject_file);
end

load(subject_file, 'subjects', 'output_dir');
fprintf('Loaded %d subjects from %s\n\n', length(subjects), subject_file);

%% ========================================================================
%  PROCESS EACH SUBJECT
%  ========================================================================

n_subjects = length(subjects);
success = false(n_subjects, 1);
errors = cell(n_subjects, 1);

for i = 1:n_subjects
    s = subjects(i);

    % Skip excluded subjects
    if isfield(s, 'exclude_subject') && s.exclude_subject
        fprintf('Skipping %s: excluded\n', s.SaveName);
        success(i) = true;
        continue;
    end

    fprintf('Processing %d/%d: %s\n', i, n_subjects, s.SaveName);

    try
        % Load simulation data
        input_file = fullfile(input_dir, sprintf('%s.mat', s.SaveName));
        if ~exist(input_file, 'file')
            error('Input file not found: %s', input_file);
        end

        loaded = load(input_file);
        data = loaded.data;
        Fs_original = loaded.Fs_final;
        Files = loaded.Files;
        header = loaded.header;
        Label = loaded.Label;

        n_blocks = length(data);
        n_channels = size(data{1}, 2);

        fprintf('  Loaded: %d blocks, %d channels, Fs = %d Hz\n', ...
            n_blocks, n_channels, Fs_original);

        % Apply filtering if requested
        if ~strcmp(filter_mode, 'none')
            fprintf('  Applying %s filter...\n', filter_mode);
            data = apply_filter(data, Fs_original, filter_mode, Fs_target);
        end

        % Downsample if needed
        if Fs_target > 0 && Fs_target < Fs_original
            fprintf('  Downsampling: %d Hz -> %d Hz\n', Fs_original, Fs_target);
            data = downsample_data(data, Fs_original, Fs_target);
            Fs_final = Fs_target;
        else
            Fs_final = Fs_original;
        end

        % Save preprocessed data
        output_file = fullfile(output_dir, sprintf('%s.mat', s.SaveName));
        Fs_initial = Fs_original;
        save(output_file, 'data', 'Fs_final', 'Fs_initial', 'Files', 'header', 'Label', '-v7.3');

        fprintf('  Saved: %s\n', output_file);
        for b = 1:n_blocks
            fprintf('    Block %d: %d samples (%.1f min)\n', b, ...
                size(data{b}, 1), size(data{b}, 1) / Fs_final / 60);
        end

        success(i) = true;

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        errors{i} = ME.message;
    end
end

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n=== Processing Complete ===\n');
fprintf('Successful: %d/%d\n', sum(success), n_subjects);

if any(~success)
    fprintf('\nFailed subjects:\n');
    for i = 1:n_subjects
        if ~success(i)
            fprintf('  - %s: %s\n', subjects(i).SaveName, errors{i});
        end
    end
end

fprintf('\nNext steps:\n');
fprintf('  1. Run CSCS_run_ICA.m (or skip with skip_ica flag)\n');
fprintf('  2. Run CSCS_run_VAR.m for VAR analysis\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function data_out = apply_filter(data_in, Fs, mode, Fs_target)
%APPLY_FILTER Apply filtering to data blocks

n_blocks = length(data_in);
data_out = cell(size(data_in));

switch mode
    case 'butter'
        % Butterworth low-pass for anti-aliasing
        if Fs_target > 0 && Fs_target < Fs
            cutoff = 0.45 * Fs_target;
        else
            cutoff = 0.45 * Fs;
        end
        [b, a] = butter(2, cutoff / (Fs/2), 'low');
        for i = 1:n_blocks
            data_out{i} = filtfilt(b, a, data_in{i});
        end
    otherwise
        data_out = data_in;
end
end

function data_out = downsample_data(data_in, Fs_original, Fs_target)
%DOWNSAMPLE_DATA Downsample data blocks

n_blocks = length(data_in);
data_out = cell(size(data_in));

factor = round(Fs_original / Fs_target);

for i = 1:n_blocks
    data_out{i} = data_in{i}(1:factor:end, :);
end
end

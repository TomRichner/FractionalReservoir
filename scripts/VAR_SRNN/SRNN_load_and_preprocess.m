%% SRNN_load_and_preprocess.m
% Load and preprocess simulated SRNN data to match CSCS output format exactly.
%
% This script is analogous to CSCSreview_load_revived_parallel.m but for
% simulation data. It applies the same filtering and decimation as the SEEG
% preprocessing pipeline.
%
% Processing steps:
%   1. Load raw simulation .mat files
%   2. Apply Butterworth lowpass filter (anti-aliasing)
%   3. Downsample from Fs_initial to Fs_final
%   4. Detrend each block
%   5. Save with all CSCS-compatible fields
%
% Prerequisites:
%   1. Run VAR_SRNN_comparison.m to generate raw simulation data
%
% See also: VAR_SRNN_comparison, CSCSreview_loadfunction_revived

close all
clc

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Input/output directories
script_dir = fileparts(mfilename('fullpath'));
input_dir = fullfile(script_dir, '..', '..', 'data', 'VAR_SRNN');
output_dir = fullfile(script_dir, '..', '..', 'data', 'VAR_SRNN_processed');

% Subject structure file
subject_file = fullfile(input_dir, 'SRNN_subjects.mat');

% Decimation settings
% Options: 'no_deci', 'butter', 'use_decimate', 'box_car_deci'
deci_mode = 'no_deci';  % No decimation - simulation already at target Fs

% Target sampling rate (Hz)
% Production: 500 Hz with 'butter' mode (requires simulation at >= 500 Hz)
% Development: use same as simulation Fs (typically 200 Hz) with no_deci
Fs_final_target = 100;

% Plotting settings
plot_comparison = true;         % Generate before/after comparison plots
plot_time_window = [0, 120];     % Time window to plot (seconds)
save_figs = true;               % Save comparison figures

% Detrend settings
% Options: 'none', 'dc', 'linear'
%   'none'   - No detrending
%   'dc'     - Remove DC offset (mean) only
%   'linear' - Remove DC and linear slope (MATLAB's detrend default)
detrend_mode = 'dc';

% Create output directory if needed
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

%% ========================================================================
%  LOAD SUBJECT STRUCTURE
%  ========================================================================

fprintf('=== SRNN_load_and_preprocess ===\n');
fprintf('Preprocessing simulated SRNN data for CSCS pipeline\n\n');

if ~exist(subject_file, 'file')
    error('Subject file not found: %s\nRun VAR_SRNN_comparison.m first.', subject_file);
end

load(subject_file, 'subjects');
fprintf('Loaded %d subjects from %s\n\n', length(subjects), subject_file);

fprintf('Configuration:\n');
fprintf('  Input dir: %s\n', input_dir);
fprintf('  Output dir: %s\n', output_dir);
fprintf('  Decimation mode: %s\n', deci_mode);
fprintf('  Target Fs_final: %d Hz\n\n', Fs_final_target);

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
        % Load raw simulation data
        input_file = fullfile(input_dir, sprintf('%s.mat', s.SaveName));
        if ~exist(input_file, 'file')
            error('Input file not found: %s', input_file);
        end

        loaded = load(input_file);
        data = loaded.data;
        Fs_initial = loaded.Fs_initial;
        Files = loaded.Files;
        header = loaded.header;
        Label = loaded.Label;
        channel_labels = loaded.channel_labels;

        n_blocks = length(data);
        n_channels = size(data{1}, 2);

        fprintf('  Loaded: %d blocks, %d channels, Fs_initial = %d Hz\n', ...
            n_blocks, n_channels, Fs_initial);

        % Store pre-decimation data for comparison plot
        if plot_comparison
            data_orig = data;
        end

        % Validate: cannot upsample
        if Fs_initial < Fs_final_target
            error('Cannot upsample: Fs_initial (%d) < Fs_final_target (%d)', ...
                Fs_initial, Fs_final_target);
        end

        % Validate: Fs_final must divide evenly into Fs_initial
        if strcmp(deci_mode, 'no_deci')
            Fs_final = Fs_initial;
            deci_factor = 1;
        else
            if mod(Fs_initial, Fs_final_target) ~= 0
                error('Fs_final (%d Hz) does not divide evenly into Fs_initial (%d Hz)', ...
                    Fs_final_target, Fs_initial);
            end
            Fs_final = Fs_final_target;
            deci_factor = Fs_initial / Fs_final;
        end

        fprintf('  Decimation factor: %d (Fs: %d -> %d Hz)\n', ...
            deci_factor, Fs_initial, Fs_final);

        % Process each block
        for b = 1:n_blocks
            block_data = data{b};
            n_samples_orig = size(block_data, 1);

            % Process each channel
            for ch = 1:n_channels
                tmp = block_data(:, ch);

                % Apply decimation based on mode
                switch deci_mode
                    case 'no_deci'
                        % No decimation - keep original data

                    case 'use_decimate'
                        % MATLAB's decimate function (includes anti-aliasing filter)
                        tmp = decimate(tmp, deci_factor);

                    case 'box_car_deci'
                        % Box-car decimation: average every deci_factor samples
                        n_samples = floor(length(tmp) / deci_factor) * deci_factor;
                        tmp = tmp(1:n_samples);
                        tmp = mean(reshape(tmp, deci_factor, []), 1)';

                    case 'butter'
                        % Butterworth lowpass filter then downsample
                        % Remove DC before filtering to avoid filtfilt edge effects
                        dc_offset = mean(tmp);
                        tmp = tmp - dc_offset;
                        % Cutoff at 0.45*Fs_final normalized to Nyquist (Fs_initial/2)
                        Wn = 0.45 * Fs_final / (Fs_initial / 2);
                        [b_filt, a_filt] = butter(2, Wn, 'low');
                        tmp = filtfilt(b_filt, a_filt, tmp);
                        % Add DC back (will be removed later if detrending)
                        tmp = tmp + dc_offset;
                        % Downsample by taking every deci_factor-th sample
                        tmp = tmp(1:deci_factor:end);

                    otherwise
                        error('Unknown decimation mode: %s', deci_mode);
                end

                % Detrend based on mode
                switch detrend_mode
                    case 'none'
                        % No detrending
                    case 'dc'
                        % Remove DC offset (mean) only
                        tmp = tmp - mean(tmp);
                    case 'linear'
                        % Remove DC and linear slope
                        tmp = detrend(tmp);
                    otherwise
                        error('Unknown detrend mode: %s', detrend_mode);
                end

                % Store back (resize on first channel)
                if ch == 1
                    n_samples_final = length(tmp);
                    processed_block = zeros(n_samples_final, n_channels);
                end
                processed_block(:, ch) = tmp;
            end

            data{b} = processed_block;
            fprintf('    Block %d: %d -> %d samples\n', b, n_samples_orig, n_samples_final);
        end

        % Set Fs = Fs_final for CSCS compatibility
        Fs = Fs_final;

        % Save processed data with all CSCS-compatible fields
        output_file = fullfile(output_dir, sprintf('%s.mat', s.SaveName));
        save(output_file, 'data', 'Fs', 'Fs_final', 'Fs_initial', 'Files', 'Label', ...
            'header', 'deci_mode', 'channel_labels', '-v7.3');

        fprintf('  Saved: %s\n', output_file);

        % Generate comparison plot
        if plot_comparison
            fprintf('  Generating comparison plot...\n');

            % Use block 1 (baseline) for comparison
            x_orig = data_orig{1};
            x_deci = data{1};

            % Create time vectors
            n_samples_orig = size(x_orig, 1);
            n_samples_deci = size(x_deci, 1);
            t_orig = (0:n_samples_orig-1)' / Fs_initial;
            t_deci = (0:n_samples_deci-1)' / Fs_final;

            % Set up plotting parameters
            plot_params = struct();
            plot_params.time_window = plot_time_window;
            plot_params.n_E = header.n_E;
            plot_params.linewidth_orig = 1.5;
            plot_params.linewidth_deci = 0.5;
            plot_params.detrend_mode = detrend_mode;

            % Create comparison figure
            fig = plot_decimation_comparison(t_orig, x_orig, t_deci, x_deci, ...
                plot_params, s.SaveName);

            % Save figure
            if save_figs
                fig_file = fullfile(output_dir, sprintf('%s_decimation.png', s.SaveName));
                saveas(fig, fig_file);
                fprintf('  Saved figure: %s\n', fig_file);
            end
        end

        success(i) = true;

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        errors{i} = ME.message;
    end
end

%% ========================================================================
%  UPDATE SUBJECT STRUCTURE WITH NEW OUTPUT DIR
%  ========================================================================

fprintf('\nUpdating subject structure...\n');

% Update FileDir in subjects to point to processed data
for i = 1:n_subjects
    subjects(i).FileDir = output_dir;
    subjects(i).Fs_final = Fs_final_target;
end

% Save updated subjects structure
processed_subject_file = fullfile(output_dir, 'SRNN_subjects.mat');
BlockLength = 16;  % placeholder
BlockLengthFinal = 15;
password = '';
save(processed_subject_file, 'subjects', 'deci_mode', 'password', ...
    'BlockLength', 'BlockLengthFinal', 'output_dir', '-v7.3');
fprintf('Saved updated subjects to: %s\n', processed_subject_file);

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

fprintf('\nOutput directory: %s\n', output_dir);
fprintf('\nNext steps:\n');
fprintf('  1. Run CSCS_check_for_noise.m (optional)\n');
fprintf('  2. Run CSCS_run_ICA.m\n');
fprintf('  3. Run CSCS_pick_ICA_components.m\n');
fprintf('  4. Run CSCS_run_VAR.m\n');

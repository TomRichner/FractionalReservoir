classdef test_VAR_SRNN_export < matlab.unittest.TestCase
    %TEST_VAR_SRNN_EXPORT Unit tests for VAR SRNN export functionality
    %
    %   Tests the data format compatibility between SRNN simulation output
    %   and the cscs_dynamics VAR analysis pipeline.
    %
    %   Run tests:
    %       results = runtests('test_VAR_SRNN_export');
    %
    %   See also: export_SRNN_to_CSCS_format, SRNN_VAR_subject_structure

    properties
        TestDataDir
        TempDir
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            % Add paths
            script_dir = fileparts(mfilename('fullpath'));
            addpath(fullfile(script_dir, '..', '..', 'src'));
            addpath(genpath(fullfile(script_dir, '..', '..', 'src')));

            % Create temp directory for test outputs
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
        end
    end

    methods (TestClassTeardown)
        function cleanupTempDir(testCase)
            % Remove temp directory
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)
        function test_export_format_structure(testCase)
            %TEST_EXPORT_FORMAT_STRUCTURE Verify output has required fields

            % Create mock result
            result = create_mock_result();
            sim_config = create_mock_sim_config();

            % Export
            export_SRNN_to_CSCS_format(result, 'test_subject', testCase.TempDir, sim_config);

            % Load and verify
            output_file = fullfile(testCase.TempDir, 'test_subject.mat');
            testCase.verifyTrue(exist(output_file, 'file') == 2, ...
                'Output file should exist');

            loaded = load(output_file);

            % Check required fields
            testCase.verifyTrue(isfield(loaded, 'data'), 'Should have data field');
            testCase.verifyTrue(isfield(loaded, 'Fs_final'), 'Should have Fs_final field');
            testCase.verifyTrue(isfield(loaded, 'Files'), 'Should have Files field');
            testCase.verifyTrue(isfield(loaded, 'header'), 'Should have header field');
            testCase.verifyTrue(isfield(loaded, 'Label'), 'Should have Label field');
        end

        function test_data_dimensions(testCase)
            %TEST_DATA_DIMENSIONS Verify data has correct [time x channels] format

            result = create_mock_result();
            sim_config = create_mock_sim_config();

            export_SRNN_to_CSCS_format(result, 'test_dims', testCase.TempDir, sim_config);

            loaded = load(fullfile(testCase.TempDir, 'test_dims.mat'));

            % Check data is cell array
            testCase.verifyTrue(iscell(loaded.data), 'data should be cell array');

            % Check each block has [time x channels] format
            for b = 1:length(loaded.data)
                block_data = loaded.data{b};
                n_samples = size(block_data, 1);
                n_channels = size(block_data, 2);

                % Time should be first dimension (larger than channels for typical data)
                testCase.verifyGreaterThan(n_samples, 0, 'Should have samples');
                testCase.verifyGreaterThan(n_channels, 0, 'Should have channels');
            end
        end

        function test_channel_count_matches_files(testCase)
            %TEST_CHANNEL_COUNT_MATCHES_FILES Verify Files matches data columns

            result = create_mock_result();
            sim_config = create_mock_sim_config();

            export_SRNN_to_CSCS_format(result, 'test_channels', testCase.TempDir, sim_config);

            loaded = load(fullfile(testCase.TempDir, 'test_channels.mat'));

            n_files = length(loaded.Files);
            n_data_channels = size(loaded.data{1}, 2);

            testCase.verifyEqual(n_files, n_data_channels, ...
                'Number of Files should match number of data channels');
        end

        function test_sampling_rate(testCase)
            %TEST_SAMPLING_RATE Verify Fs_final is set correctly

            result = create_mock_result();
            sim_config = create_mock_sim_config();
            expected_fs = result.fs;

            export_SRNN_to_CSCS_format(result, 'test_fs', testCase.TempDir, sim_config);

            loaded = load(fullfile(testCase.TempDir, 'test_fs.mat'));

            testCase.verifyEqual(loaded.Fs_final, expected_fs, ...
                'Fs_final should match simulation sampling rate');
        end

        function test_subject_structure_format(testCase)
            %TEST_SUBJECT_STRUCTURE_FORMAT Verify subject struct has required fields

            % Create mock conditions and config
            conditions = struct('name', {'SRNN_test1', 'SRNN_test2'}, ...
                'n_a_E', {0, 3}, ...
                'n_b_E', {0, 1}, ...
                'description', {'Test 1', 'Test 2'});
            sim_config = create_mock_sim_config();

            % Generate subject structure
            SRNN_VAR_subject_structure(conditions, sim_config, testCase.TempDir);

            % Load and verify
            loaded = load(fullfile(testCase.TempDir, 'SRNN_subjects.mat'));

            testCase.verifyTrue(isfield(loaded, 'subjects'), ...
                'Should have subjects field');

            subjects = loaded.subjects;

            % Check required fields
            required_fields = {'SaveName', 'Fs_final', 'baseline_block_indices', ...
                'stim_block_indices', 'exclude_subject', 'Files', 'Label'};
            for f = 1:length(required_fields)
                testCase.verifyTrue(isfield(subjects, required_fields{f}), ...
                    sprintf('Subjects should have %s field', required_fields{f}));
            end
        end

        function test_neuron_subset_option(testCase)
            %TEST_NEURON_SUBSET_OPTION Verify random subset selection works

            result = create_mock_result();
            sim_config = create_mock_sim_config();
            sim_config.n_neurons_subset = 20;  % Use subset of 20

            export_SRNN_to_CSCS_format(result, 'test_subset', testCase.TempDir, sim_config);

            loaded = load(fullfile(testCase.TempDir, 'test_subset.mat'));

            n_channels = size(loaded.data{1}, 2);
            testCase.verifyEqual(n_channels, 20, ...
                'Should have 20 channels when subset specified');
        end

        function test_baseline_stim_split(testCase)
            %TEST_BASELINE_STIM_SPLIT Verify data is split into baseline and stim blocks

            result = create_mock_result();
            sim_config = create_mock_sim_config();

            export_SRNN_to_CSCS_format(result, 'test_split', testCase.TempDir, sim_config);

            loaded = load(fullfile(testCase.TempDir, 'test_split.mat'));

            % Should have exactly 2 blocks
            testCase.verifyEqual(length(loaded.data), 2, ...
                'Should have 2 blocks (baseline + stim)');

            % Both blocks should have data
            testCase.verifyGreaterThan(size(loaded.data{1}, 1), 0, ...
                'Baseline block should have samples');
            testCase.verifyGreaterThan(size(loaded.data{2}, 1), 0, ...
                'Stim block should have samples');
        end
    end
end

%% Helper Functions

function result = create_mock_result()
%CREATE_MOCK_RESULT Create a mock result structure for testing

% Parameters
n = 100;
n_E = 50;
n_I = 50;
fs = 200;
T = 10;  % seconds
nt = T * fs;

% Create params struct
params = struct();
params.n = n;
params.n_E = n_E;
params.n_I = n_I;
params.n_a_E = 0;
params.n_b_E = 0;
params.level_of_chaos = 1.7;
params.rng_seeds = [1 2 3 4 5];

% Create mock dendritic states
x.E = randn(n_E, nt);
x.I = randn(n_I, nt);

r.E = abs(x.E);
r.I = abs(x.I);

a.E = [];
a.I = [];

b.E = ones(n_E, nt);
b.I = ones(n_I, nt);

% Time vector
t = (0:nt-1)' / fs;

% Create masks for baseline/stim split
T_baseline = T / 2;
baseline_mask = t < T_baseline;
stim_mask = t >= T_baseline;

% Assemble result
result = struct();
result.params = params;
result.t = t;
result.x = x;
result.r = r;
result.a = a;
result.b = b;
result.fs = fs;
result.T_baseline = T_baseline;
result.baseline_mask = baseline_mask;
result.stim_mask = stim_mask;
result.u_ex = randn(n, nt);
result.t_ex = t;
result.S_out = randn(nt, n);
end

function sim_config = create_mock_sim_config()
%CREATE_MOCK_SIM_CONFIG Create a mock simulation config for testing

sim_config = struct();
sim_config.fs = 200;
sim_config.T_baseline = 5;
sim_config.T_stim = 5;
sim_config.n_neurons_subset = [];  % Use all
sim_config.n = 100;
sim_config.rng_seeds = [1 2 3 4 5];
end

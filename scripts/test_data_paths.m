% test_data_paths.m
% Verifies that analysis classes default to the correct data directories

clear; clc;
addpath(genpath('../src'));
fprintf('=== Checking Analysis Class Defaults ===\n');

project_root = fileparts(fileparts(mfilename('fullpath')));
data_dir = fullfile(project_root, 'data');

%% 1. SensitivityAnalysis
fprintf('\n--- SensitivityAnalysis ---\n');
sa = SensitivityAnalysis();
expected = fullfile(data_dir, 'sensitivity');
if contains(sa.output_dir, expected)
    fprintf('PASS: Output dir is %s\n', sa.output_dir);
else
    fprintf('FAIL: Expected %s, got %s\n', expected, sa.output_dir);
end

%% 2. PairedPulseMIAnalysis
fprintf('\n--- PairedPulseMIAnalysis ---\n');
pp = PairedPulseMIAnalysis();
expected = fullfile(data_dir, 'paired_pulse');
if contains(pp.output_dir, expected)
    fprintf('PASS: Output dir is %s\n', pp.output_dir);
else
    fprintf('FAIL: Expected %s, got %s\n', expected, pp.output_dir);
end

%% 3. ParamSpaceAnalysis
fprintf('\n--- ParamSpaceAnalysis ---\n');
psa = ParamSpaceAnalysis();
expected = fullfile(data_dir, 'param_space');
if contains(psa.output_dir, expected)
    fprintf('PASS: Output dir is %s\n', psa.output_dir);
else
    fprintf('FAIL: Expected %s, got %s\n', expected, psa.output_dir);
end

fprintf('\n=== Done ===\n');

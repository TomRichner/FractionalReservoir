%% test_paired_beeswarm.m
% Test script demonstrating paired_beeswarm functionality on synthetic data

%% Setup - add beeswarm to path
addpath('/Users/richner.thomas/Desktop/local_code/beeswarm');

%% Test 1: Basic usage with 2 conditions
figure('Name', 'Test 1: Basic 2 Conditions', 'Position', [100 100 400 400]);
N = 30;
data = [randn(N,1), randn(N,1) + 1];  % Condition 2 shifted up by 1
paired_beeswarm(data);
title('Basic Usage: 2 Conditions');
ylabel('Value');

%% Test 2: Three conditions with per-subject colors
figure('Name', 'Test 2: Per-Subject Colors', 'Position', [520 100 400 400]);
N = 25;
data = [randn(N,1), randn(N,1) + 0.5, randn(N,1) + 1.5];
colors = parula(N);  % One color per subject
paired_beeswarm(data, 'Colors', colors, 'Labels', {'Baseline', 'Treatment', 'Follow-up'});
title('Per-Subject Colors: 3 Conditions');
ylabel('Value');

%% Test 3: Single color for all
figure('Name', 'Test 3: Single Color', 'Position', [940 100 400 400]);
N = 40;
data = [randn(N,1)*0.5 + 2, randn(N,1)*0.8 + 1];  % Different spread and means
paired_beeswarm(data, 'Colors', [0.2 0.4 0.8], 'Alpha', 0.7, 'LineWidth', 1);
title('Single Color for All');
ylabel('Measurement');

%% Test 4: Different sort styles
figure('Name', 'Test 4: Sort Styles', 'Position', [100 520 800 400]);
N = 35;
data_for_sort = [randn(N,1), randn(N,1) + 1];

sortStyles = {'nosort', 'up', 'hex', 'square'};
for i = 1:length(sortStyles)
    subplot(1, 4, i);
    paired_beeswarm(data_for_sort, 'SortStyle', sortStyles{i}, ...
        'Colors', lines(N), 'Alpha', 0.6);
    title(['SortStyle: ' sortStyles{i}]);
    ylabel('Value');
end

%% Test 5: Large N with many conditions
figure('Name', 'Test 5: Large N, 4 Conditions', 'Position', [520 520 500 400]);
N = 50;
data = randn(N, 4);
data(:,2) = data(:,1) + 0.3 + randn(N,1)*0.2;  % Correlated with condition 1
data(:,3) = data(:,1) - 0.5 + randn(N,1)*0.3;  % Anti-correlated
data(:,4) = randn(N,1) + 2;  % Independent

colors = cool(N);
paired_beeswarm(data, 'Colors', colors, ...
    'Labels', {'A', 'B', 'C', 'D'}, ...
    'MarkerSize', 0.8, 'LineWidth', 0.3, 'Alpha', 0.4);
title('50 Subjects, 4 Conditions');
ylabel('Score');

%% Test 6: Extreme paired differences (clear treatment effect)
figure('Name', 'Test 6: Clear Treatment Effect', 'Position', [1040 520 400 400]);
N = 20;
baseline = randn(N,1) * 0.5;
treatment = baseline + 2 + randn(N,1)*0.3;  % Strong consistent increase
data = [baseline, treatment];

% Color by baseline value
colors = winter(N);
[~, sortIdx] = sort(baseline);
reorderedColors = zeros(N, 3);
reorderedColors(sortIdx, :) = colors;

paired_beeswarm(data, 'Colors', reorderedColors, ...
    'Labels', {'Pre', 'Post'}, 'LineWidth', 1.5, 'Alpha', 0.8);
title('Clear Treatment Effect');
ylabel('Response');

%% Test 7: Handling NaN values
figure('Name', 'Test 7: Missing Data (NaN)', 'Position', [100 50 400 400]);
N = 25;
data = [randn(N,1), randn(N,1) + 0.5];
% Introduce some NaN values
data(5, 2) = NaN;
data(10, 1) = NaN;
data(15, :) = NaN;

paired_beeswarm(data, 'Colors', copper(N), 'Labels', {'Cond1', 'Cond2'});
title('Handling NaN Values');
ylabel('Value');

%% Summary
fprintf('\n=== test_paired_beeswarm completed ===\n');
fprintf('Generated 7 test figures demonstrating:\n');
fprintf('  1. Basic 2-condition usage\n');
fprintf('  2. Per-subject coloring with 3 conditions\n');
fprintf('  3. Single color for all subjects\n');
fprintf('  4. Different beeswarm sort styles\n');
fprintf('  5. Large N with 4 conditions\n');
fprintf('  6. Clear treatment effect visualization\n');
fprintf('  7. Handling of missing data (NaN)\n');

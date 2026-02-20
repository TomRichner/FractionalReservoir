%% Memory Capacity (MC) experiment — paper-ready, replicated, paired design
% Compares 3 SRNN conditions:
%   1) Baseline (no adaptation)
%   2) SFA only
%   3) SFA + STD
%
% Key "fairness" choices:
%   - Paired trials: same (W, W_in, input u(t)) across conditions within each trial
%   - Primary input: i.i.d. white noise (prevents inflated MC from autocorrelated inputs)
%   - Many trials (different seeds) to support statistical inference
%   - Paper-ready plots (distributions + mean curves w/ CI)
%
% Requirements:
%   - SRNN_ESN_reservoir class must support input_type 'white'
%   - verify_shared_build(esn, differing_props, shared_fields) must exist

clear; clc; close all;

%% Add paths
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

%% -------------------- Global experiment settings --------------------
% Network / dynamics
n = 300;                    % Number of neurons
level_of_chaos = 1.0;       % Chaos level

% Sampling
fs = 200;                   % Hz

% MC protocol (seconds -> samples)
T_wash_sec  = 20;
T_train_sec = 50;
T_test_sec  = 50;

T_wash  = T_wash_sec  * fs;
T_train = T_train_sec * fs;
T_test  = T_test_sec  * fs;

d_max_sec = 3.0;            % seconds
d_max = round(d_max_sec * fs);

% Input (PRIMARY: white / i.i.d.)
input_type = 'white';       % 'white' | 'bandlimited' | 'one_over_f'
u_f_cutoff = 5;             % only used if bandlimited
u_alpha    = 1;             % only used if one_over_f (1=pink, 2=red)

% Trials / seeds
n_trials = 50;              % >=30 recommended; 50-100 typical
seed_net_base  = 1000;      % deterministic seed schedule
seed_stim_base = 2000;

% Analysis knobs
R2_threshold_for_horizon = 0.10;
n_boot = 2000;              % bootstrap resamples for CI
n_perm = 10000;             % permutation sign-flip count

% Output directory
script_dir   = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
out_dir = fullfile(project_root, 'data', 'memory_capacity', 'paper_ready');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_tag = sprintf('MC_%s_%s_trials%d', input_type, timestamp, n_trials);

%% -------------------- Conditions --------------------
condition_names = {'Baseline', 'SFA', 'SFA+STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...   % Baseline
    {'n_a_E', 3, 'n_b_E', 0}, ...   % SFA only
    {'n_a_E', 3, 'n_b_E', 1}, ...   % SFA + STD
};
n_cond = numel(condition_names);

%% -------------------- Shared base config template --------------------
% NOTE: rng_seeds will be set per trial: [seed_net, seed_stim]
base_args_template = { ...
    'n', n, ...
    'fs', fs, ...
    'level_of_chaos', level_of_chaos, ...
    'tau_d', 0.1, ...              % Dendritic time constant (s)
    'S_c', 0.4, ...                % Nonlinearity bias (center)
    'S_a', 0.9, ...                % Fraction of nonlinearity with slope 1
    'n_a_I', 0, ...                % no SFA for I neurons (all conditions)
    'n_b_I', 0, ...                % no STD for I neurons (all conditions)
    'c_E', 0.15/3, ...             % adaptation strength for E neurons
    'tau_a_E', [0.1, 1.0, 10], ... % SFA time constants (s)
    'tau_b_E_rec', 1.0, ...        % STD recovery (s)
    'tau_b_E_rel', 0.25, ...       % STD release (s)
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'u_alpha', u_alpha, ...
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max ...
};

%% -------------------- Preallocate logs --------------------
MC_trials  = nan(n_trials, n_cond);
H_trials   = nan(n_trials, n_cond);          % horizon delay (seconds) at R2>threshold
R2_trials  = nan(n_trials, n_cond, d_max);   % per-delay R2

seed_net   = seed_net_base  + (1:n_trials);
seed_stim  = seed_stim_base + (1:n_trials);

% (Optional) keep per-trial structs (can be large; enable if you need internals)
store_internal_results = false;
internal_results = cell(n_trials, n_cond);

%% -------------------- Main loop: paired trials --------------------
fprintf('\n==== Running %d paired trials (%s input) ====\n', n_trials, input_type);

for k = 1:n_trials
    fprintf('\n--- Trial %d / %d | seeds: net=%d, stim=%d ---\n', ...
        k, n_trials, seed_net(k), seed_stim(k));

    % Per-trial base args includes trial seeds (ensures paired W, W_in, u(t))
    base_args = [{'rng_seeds', [seed_net(k), seed_stim(k)]}, base_args_template];

    % Build all 3 conditions
    esn = cell(1, n_cond);
    for i = 1:n_cond
        esn{i} = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
        esn{i}.build();
    end

    % Verify identical structure + identical input across conditions (fairness check)
    verify_shared_build(esn, {'n_a_E','n_b_E'}, {'W','W_in','u_scalar','u_ex','t_ex'});

    % Run MC
    for i = 1:n_cond
        [mc_i, r2_i, res_i] = esn{i}.run_memory_capacity();

        MC_trials(k,i) = mc_i;
        r2_i = r2_i(:);
        if numel(r2_i) ~= d_max
            error('Trial %d condition %s: expected d_max=%d, got %d', ...
                k, condition_names{i}, d_max, numel(r2_i));
        end
        R2_trials(k,i,:) = r2_i;

        % horizon (in seconds): max delay with R2 > threshold
        idx = find(r2_i > R2_threshold_for_horizon, 1, 'last');
        if isempty(idx)
            H_trials(k,i) = 0;
        else
            H_trials(k,i) = idx / fs;
        end

        if store_internal_results
            internal_results{k,i} = res_i;
        end
    end
end

fprintf('\n==== Done. Computing summary + figures... ====\n');

%% -------------------- Summary statistics --------------------
MC_mean = mean(MC_trials, 1, 'omitnan');
MC_sem  = std(MC_trials, 0, 1, 'omitnan') / sqrt(n_trials);

H_mean  = mean(H_trials, 1, 'omitnan');
H_sem   = std(H_trials, 0, 1, 'omitnan') / sqrt(n_trials);

% Bootstrap 95% CI for mean MC and horizon
MC_ci = bootstrap_mean_ci(MC_trials, n_boot, 0.05);
H_ci  = bootstrap_mean_ci(H_trials,  n_boot, 0.05);

% Mean R2(d) and CI across trials
R2_mean = squeeze(mean(R2_trials, 1, 'omitnan')); % [cond x d]
R2_ci   = bootstrap_mean_ci_3d(R2_trials, n_boot, 0.05); % struct with lo/hi [cond x d]

% Cumulative MC curves (mean across trials)
R2_cum_mean = cumsum(R2_mean, 2);
R2_cum_ci_lo = cumsum(R2_ci.lo, 2);
R2_cum_ci_hi = cumsum(R2_ci.hi, 2);

%% -------------------- Paired statistical tests (Total MC) --------------------
% Pairwise comparisons with sign-flip permutation test + paired effect size
pairs = [1 2; 1 3; 2 3]; % (Baseline vs SFA), (Baseline vs SFA+STD), (SFA vs SFA+STD)
pair_labels = { ...
    'Baseline vs SFA', ...
    'Baseline vs SFA+STD', ...
    'SFA vs SFA+STD' ...
};

stats = struct();
for p = 1:size(pairs,1)
    a = pairs(p,1); b = pairs(p,2);
    x = MC_trials(:,a);
    y = MC_trials(:,b);
    [pval, diff_mean] = paired_signflip_permutation_pvalue(x, y, n_perm);
    d_z = paired_cohens_dz(x, y);
    stats(p).pair = pair_labels{p};
    stats(p).mean_diff = diff_mean;
    stats(p).p_perm = pval;
    stats(p).cohens_dz = d_z;
end

%% -------------------- Paper-ready figures --------------------
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultLegendInterpreter','none');

delay_s = (1:d_max) / fs;

% Colors (simple, print-friendly)
colors = [0.45 0.45 0.45;   % Baseline
          0.20 0.50 0.80;   % SFA
          0.80 0.30 0.25];  % SFA+STD

% Figure 1: Total MC distribution (paired) + mean/CI
fig1 = figure('Color','w','Position',[100 100 1050 420]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% 1A: Paired scatter with lines
nexttile; hold on; grid on; box on;
xpos = 1:n_cond;
for k = 1:n_trials
    plot(xpos, MC_trials(k,:), '-', 'Color', [0 0 0 0.15], 'LineWidth', 1);
end
for i = 1:n_cond
    scatter(i*ones(n_trials,1), MC_trials(:,i), 20, ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.65);
end
% Mean + 95% CI
for i = 1:n_cond
    plot(i, MC_mean(i), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot([i i], [MC_ci.lo(i) MC_ci.hi(i)], '-', 'Color', 'k', 'LineWidth', 2);
end
xlim([0.5 n_cond+0.5]);
set(gca,'XTick',xpos,'XTickLabel',condition_names);
ylabel('Total Memory Capacity (sum R^2)');
title('Total MC (paired trials)');

% 1B: Horizon distribution (optional but helpful)
nexttile; hold on; grid on; box on;
for k = 1:n_trials
    plot(xpos, H_trials(k,:), '-', 'Color', [0 0 0 0.15], 'LineWidth', 1);
end
for i = 1:n_cond
    scatter(i*ones(n_trials,1), H_trials(:,i), 20, ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.65);
end
for i = 1:n_cond
    plot(i, H_mean(i), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot([i i], [H_ci.lo(i) H_ci.hi(i)], '-', 'Color', 'k', 'LineWidth', 2);
end
xlim([0.5 n_cond+0.5]);
set(gca,'XTick',xpos,'XTickLabel',condition_names);
ylabel(sprintf('Memory horizon (s), R^2 > %.2f', R2_threshold_for_horizon));
title('Horizon (paired trials)');

sgtitle(sprintf('Memory Capacity Comparison (%s input, n=%d, trials=%d)', input_type, n, n_trials));

saveas(fig1, fullfile(out_dir, [run_tag '_Fig1_MC_Distributions.png']));
print(fig1, fullfile(out_dir, [run_tag '_Fig1_MC_Distributions.pdf']), '-dpdf', '-painters');

% Figure 2: R^2(d) mean curve with 95% CI shading
fig2 = figure('Color','w','Position',[100 100 1050 420]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% 2A: R^2(d)
nexttile; hold on; grid on; box on;
for i = 1:n_cond
    shaded_ci(delay_s, R2_ci.lo(i,:), R2_ci.hi(i,:), colors(i,:), 0.18);
    plot(delay_s, R2_mean(i,:), '-', 'Color', colors(i,:), 'LineWidth', 2);
end
xlabel('Delay (s)');
ylabel('R^2(d)');
title('Per-delay memory (mean ± 95% CI)');
legend(condition_names,'Location','northeast');

% 2B: Cumulative MC vs delay
nexttile; hold on; grid on; box on;
for i = 1:n_cond
    shaded_ci(delay_s, R2_cum_ci_lo(i,:), R2_cum_ci_hi(i,:), colors(i,:), 0.18);
    plot(delay_s, R2_cum_mean(i,:), '-', 'Color', colors(i,:), 'LineWidth', 2);
end
xlabel('Delay (s)');
ylabel('Cumulative MC (sum_{j<=d} R^2(j))');
title('Cumulative memory (mean ± 95% CI)');
legend(condition_names,'Location','southeast');

sgtitle(sprintf('Delay Profile of Memory (%s input, n=%d, trials=%d)', input_type, n, n_trials));

saveas(fig2, fullfile(out_dir, [run_tag '_Fig2_R2_Curves.png']));
print(fig2, fullfile(out_dir, [run_tag '_Fig2_R2_Curves.pdf']), '-dpdf', '-painters');

%% -------------------- Save results --------------------
results_all = struct();
results_all.run_tag = run_tag;
results_all.timestamp = timestamp;

results_all.settings = struct();
results_all.settings.n = n;
results_all.settings.fs = fs;
results_all.settings.level_of_chaos = level_of_chaos;
results_all.settings.T_wash = T_wash;
results_all.settings.T_train = T_train;
results_all.settings.T_test = T_test;
results_all.settings.d_max = d_max;
results_all.settings.input_type = input_type;
results_all.settings.u_f_cutoff = u_f_cutoff;
results_all.settings.u_alpha = u_alpha;
results_all.settings.n_trials = n_trials;
results_all.settings.R2_threshold_for_horizon = R2_threshold_for_horizon;
results_all.settings.n_boot = n_boot;
results_all.settings.n_perm = n_perm;

results_all.conditions = condition_names;
results_all.seeds = table(seed_net(:), seed_stim(:), 'VariableNames', {'seed_net','seed_stim'});

results_all.MC_trials = MC_trials;
results_all.H_trials  = H_trials;
results_all.R2_trials = R2_trials;

results_all.summary = struct();
results_all.summary.MC_mean = MC_mean;
results_all.summary.MC_sem  = MC_sem;
results_all.summary.MC_ci95 = MC_ci;
results_all.summary.H_mean  = H_mean;
results_all.summary.H_sem   = H_sem;
results_all.summary.H_ci95  = H_ci;
results_all.summary.R2_mean = R2_mean;
results_all.summary.R2_ci95 = R2_ci;
results_all.summary.stats   = stats;

if store_internal_results
    results_all.internal_results = internal_results;
end

mat_file = fullfile(out_dir, [run_tag '_results.mat']);
save(mat_file, 'results_all', '-v7.3');

% Also save a compact CSV for total MC and horizon
T = array2table([MC_trials H_trials], ...
    'VariableNames', [strcat('MC_',condition_names), strcat('H_',condition_names)]);
T.seed_net  = seed_net(:);
T.seed_stim = seed_stim(:);
T = movevars(T, {'seed_net','seed_stim'}, 'Before', 1);
writetable(T, fullfile(out_dir, [run_tag '_MC_Horizon.csv']));

% Save a readable text summary
txt_file = fullfile(out_dir, [run_tag '_summary.txt']);
fid = fopen(txt_file,'w');
fprintf(fid, 'Run: %s\n', run_tag);
fprintf(fid, 'Input: %s | n=%d | fs=%d Hz | trials=%d | d_max=%.2f s\n\n', ...
    input_type, n, fs, n_trials, d_max/fs);

fprintf(fid, 'Total MC (mean ± SEM) and 95%% CI:\n');
for i=1:n_cond
    fprintf(fid, '  %s: %.3f ± %.3f | CI [%.3f, %.3f]\n', ...
        condition_names{i}, MC_mean(i), MC_sem(i), MC_ci.lo(i), MC_ci.hi(i));
end
fprintf(fid, '\nHorizon (R^2>%.2f) (mean ± SEM) and 95%% CI:\n', R2_threshold_for_horizon);
for i=1:n_cond
    fprintf(fid, '  %s: %.3f ± %.3f s | CI [%.3f, %.3f]\n', ...
        condition_names{i}, H_mean(i), H_sem(i), H_ci.lo(i), H_ci.hi(i));
end

fprintf(fid, '\nPaired permutation tests (sign-flip) on Total MC:\n');
for p=1:numel(stats)
    fprintf(fid, '  %s: mean diff = %.3f | p_perm = %.4g | Cohen''s dz = %.3f\n', ...
        stats(p).pair, stats(p).mean_diff, stats(p).p_perm, stats(p).cohens_dz);
end
fclose(fid);

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', mat_file, txt_file, ...
    fullfile(out_dir, [run_tag '_Fig1_MC_Distributions.pdf']), ...
    fullfile(out_dir, [run_tag '_Fig2_R2_Curves.pdf']));

%% ==================== Local helper functions ====================

function ci = bootstrap_mean_ci(X, n_boot, alpha)
% X: [N x C]
% Returns ci.lo, ci.hi for each column (mean CI)
    if nargin < 3; alpha = 0.05; end
    [N, C] = size(X);
    boot_means = nan(n_boot, C);
    for b = 1:n_boot
        idx = randi(N, N, 1);
        boot_means(b,:) = mean(X(idx,:), 1, 'omitnan');
    end
    ci.lo = quantile(boot_means, alpha/2, 1);
    ci.hi = quantile(boot_means, 1-alpha/2, 1);
end

function ci = bootstrap_mean_ci_3d(R2_trials, n_boot, alpha)
% R2_trials: [N x C x D]
% Returns ci.lo and ci.hi as [C x D] for mean across N
    if nargin < 3; alpha = 0.05; end
    [N, C, D] = size(R2_trials);
    boot_means = nan(n_boot, C, D);
    for b = 1:n_boot
        idx = randi(N, N, 1);
        boot_means(b,:,:) = squeeze(mean(R2_trials(idx,:,:), 1, 'omitnan')); % [C x D]
    end
    ci.lo = squeeze(quantile(boot_means, alpha/2, 1));     % [C x D]
    ci.hi = squeeze(quantile(boot_means, 1-alpha/2, 1));   % [C x D]
end

function [pval, mean_diff] = paired_signflip_permutation_pvalue(x, y, n_perm)
% Robust paired permutation test via sign flips on differences.
% H0: mean(x-y) = 0
    x = x(:); y = y(:);
    d = x - y;
    d = d(~isnan(d));
    mean_diff = mean(d);
    N = numel(d);
    if N == 0
        pval = NaN; return;
    end
    % Permute by random sign flips
    perm_stats = zeros(n_perm,1);
    for r = 1:n_perm
        s = 2*(rand(N,1) > 0.5) - 1; % +/-1
        perm_stats(r) = mean(s .* d);
    end
    % two-sided p-value
    pval = mean(abs(perm_stats) >= abs(mean_diff));
end

function dz = paired_cohens_dz(x, y)
% Cohen's dz for paired designs: mean(diff)/std(diff)
    d = x(:) - y(:);
    d = d(~isnan(d));
    dz = mean(d) / std(d, 0);
end

function shaded_ci(x, lo, hi, rgb, alpha_fill)
% Draw shaded confidence interval [lo, hi] around x
    x = x(:)'; lo = lo(:)'; hi = hi(:)';
    X = [x, fliplr(x)];
    Y = [hi, fliplr(lo)];
    fill(X, Y, rgb, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
end
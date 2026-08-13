%% Memory Capacity (MC) experiment — paper-ready, replicated, paired design
% Compares 4 SRNN conditions:
%   1) Baseline (no adaptation)
%   2) SFA only
%   3) STD only
%   4) SFA + STD
%
% Key "fairness" choices:
%   - Paired trials: same (W, W_in, input u(t)) across conditions within each trial
%   - Input: sample-and-hold i.i.d. values (fair for a low-pass reservoir; MC is
%     measured in hold units, i.i.d. at the hold rate -> no autocorrelation inflation)
%   - Many trials (different seeds) to support statistical inference
%   - Paper-ready plots (distributions + mean curves w/ CI)
%
% Requirements:
%   - SRNN_ESN_reservoir class must support input_type 'sample_hold'
%   - verify_shared_build(esn, differing_props, shared_fields) must exist

clear; clc; close all;

%% Add paths
setup_paths();

%% -------------------- Global experiment settings --------------------
% Network / dynamics
n = 300;                    % Number of neurons (matches SRNNModel2 default)
f = 0.6;                    % Fraction excitatory (matches SRNNModel2 default)
level_of_chaos = 2.0;       % Above the edge of chaos (logistic mean slope < 1 raises it); SRNNModel2 default is 1.0

% Sampling
fs = 200;                   % Hz (matches SRNNModel2 default)

% MC protocol (seconds -> samples). N_train = T_train_sec/T_hold hold-samples
% must exceed n features for the readout to be well posed (here ~667 > 300).
% T_wash=10 s is short vs tau_a_E=10 s (matches the example dial-in); raise it
% for the final figure if the SFA conditions need more settling.
T_wash_sec  = 10;
T_train_sec = 600;
T_test_sec  = 150;

T_wash  = T_wash_sec  * fs;
T_train = T_train_sec * fs;
T_test  = T_test_sec  * fs;

d_max_sec = 15.0;           % seconds
d_max = round(d_max_sec * fs);

% Input: sample-and-hold i.i.d. -- MC is measured in hold units (fair for a
% low-pass reservoir, free of autocorrelation inflation). See SRNN_ESN_reservoir.
input_type = 'sample_hold'; % 'white' | 'bandlimited' | 'one_over_f' | 'sample_hold'
u_f_cutoff = 5;             % only used if bandlimited
u_alpha    = 1;             % only used if one_over_f (1=pink, 2=red)
tau_d      = 0.1;           % Dendritic time constant (s)
T_hold     = 0.3;           % sample_hold: hold each i.i.d. value this long (s); sets the MC delay increment (matches compute_memory_capacity_example.m)

% Readout signal for the MC regression: 'rate' reads out r = phi(x_eff) (STD's b
% NOT exposed); 'synaptic' reads out br = b.*r (exposes the STD state to the
% linear readout; unchanged for the no-STD conditions since br == r there).
readout_signal = 'synaptic'; % 'rate' | 'synaptic'

% Trials / seeds
n_trials = 30;     % 10 for fast, 50 for production   % VALIDATION pass (fast); restore to 50-100 for real runs
seed_net_base  = 3000;      % deterministic seed schedule
seed_stim_base = 4000;

% Analysis knobs
R2_threshold_for_horizon = 0.10;
n_boot = 2000;              % bootstrap resamples for CI
n_perm = 10000;             % permutation sign-flip count

% Output directory
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out_dir = fullfile(project_root, 'data', 'memory_capacity', 'paper_ready');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_tag = sprintf('MC_%s_%s_trials%d', input_type, timestamp, n_trials);

%% -------------------- Conditions --------------------
condition_names = {'Baseline', 'SFA', 'STD', 'SFA+STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...   % Baseline
    {'n_a_E', 3, 'n_b_E', 0}, ...   % SFA only
    {'n_a_E', 0, 'n_b_E', 1}, ...   % STD only
    {'n_a_E', 3, 'n_b_E', 1}, ...   % SFA + STD
    };
n_cond = numel(condition_names);

%% -------------------- Shared base config template --------------------
% NOTE: rng_seeds will be set per trial: [seed_net, seed_stim]
base_args_template = { ...
    'n', n, ...
    'f', f, ...                    % fraction excitatory (matches SRNNModel2 default)
    'fs', fs, ...
    'level_of_chaos', level_of_chaos, ...
    'ode_solver', 'rk4', ...    % fixed-step solver (fast); SRNNModel2 default is 'ode45'
    'tau_d', tau_d, ...            % Dendritic time constant (s; matches SRNNModel2 default)
    'S_c', 0.35, ...                % Nonlinearity bias (center); matches SRNNModel2 default
    'S_a', 0.9, ...                % Fraction of nonlinearity with slope 1 (matches SRNNModel2 default; unused by the logistic)
    'n_a_I', 0, ...                % no SFA for I neurons (all conditions; matches SRNNModel2 default)
    'n_b_I', 0, ...                % no STD for I neurons (all conditions; matches SRNNModel2 default)
    'c_E', 0.5/3, ...             % adaptation strength for E neurons (matches SRNNModel2 default)
    'tau_b_E_rec', 1.0, ...        % STD recovery (s; matches SRNNModel2 default)
    'tau_b_E_rel', 0.25, ...       % STD release (s; matches SRNNModel2 default)
    'std_zero_floor', false, ...   % matches SRNNModel2 default (no zero-floor rescale)
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'u_alpha', u_alpha, ...
    'T_hold', T_hold, ...          % hold duration for input_type='sample_hold'
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max ...
    };

%% -------------------- Delay grid + preallocate logs --------------------
% Size the delay grid up front so R2_trials is fully allocated before the parfor
% (a parfor sliced output cannot be lazily grown). Mirrors the hold-unit
% decimation in SRNN_ESN_reservoir/run_memory_capacity.
if strcmpi(input_type, 'sample_hold')
    hold_len = round(T_hold * fs);
else
    hold_len = 1;
end
d_max_eff = max(1, floor(d_max / hold_len));
delay_s   = (1:d_max_eff) * (hold_len / fs);   % delay axis (seconds)

MC_trials  = nan(n_trials, n_cond);
H_trials   = nan(n_trials, n_cond);            % horizon delay (seconds) at R2>threshold
R2_trials  = nan(n_trials, n_cond, d_max_eff); % per-delay R2

seed_net   = seed_net_base  + (1:n_trials);
seed_stim  = seed_stim_base + (1:n_trials);

% (Optional) keep per-trial structs (can be large; enable if you need internals)
store_internal_results = false;
internal_results = cell(n_trials, n_cond);

%% -------------------- One-time fairness check --------------------
% verify_shared_build needs all conditions built at once (memory-heavy). The
% paired invariant (same seeds + only n_a_E/n_b_E differ => identical
% W/W_in/u_scalar/u_ex/t_ex) is STRUCTURAL, so verify it once on the first seed
% pair; the trial loop then builds one condition at a time to keep memory low.
fprintf('\nVerifying shared build (fairness check, trial 1)...\n');
chk_args = [{'rng_seeds', [seed_net(1), seed_stim(1)]}, base_args_template];
esn_chk = cell(1, n_cond);
for i = 1:n_cond
    esn_chk{i} = SRNN_ESN_reservoir(chk_args{:}, condition_args{i}{:});
    esn_chk{i}.build();
end
% tau_a_E now uses the SRNNModel2 default (auto-filled per n_a_E), so it legitimately differs across conditions.
verify_shared_build(esn_chk, {'n_a_E','n_b_E','tau_a_E'}, {'W','W_in','u_scalar','u_ex','t_ex'});
clear esn_chk chk_args;

%% -------------------- Main loop: paired trials (parallel) --------------------
% Trials are independent (each uses its own seed_net(k)/seed_stim(k) and the ESN
% seeds rng() explicitly at build), so parfor is order-independent and matches
% the serial result. Degrades to serial if no parallel pool is available.
% Each condition is built/run/released one at a time, with a lean MC-only run
% (no stored state time series, no LLE), to keep per-worker memory low.
fprintf('\n==== Running %d paired trials (%s input) ====\n', n_trials, input_type);

parfor k = 1:n_trials
    fprintf('--- Trial %d / %d | seeds: net=%d, stim=%d ---\n', ...
        k, n_trials, seed_net(k), seed_stim(k));
    
    % Per-trial base args include the trial seeds (ensures paired W, W_in, u(t))
    base_args = [{'rng_seeds', [seed_net(k), seed_stim(k)]}, base_args_template];
    
    % Per-trial accumulators (assigned to the sliced outputs once, below)
    mc_row  = nan(1, n_cond);
    h_row   = nan(1, n_cond);
    r2_row  = nan(n_cond, d_max_eff);
    int_row = cell(1, n_cond);
    
    for i = 1:n_cond
        esn_i = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
        esn_i.build();
        [mc_i, r2_i, res_i] = esn_i.run_memory_capacity('store_timeseries', false, 'verbose', false, ...
            'readout_signal', readout_signal);
        r2_i = r2_i(:)';
        if numel(r2_i) ~= d_max_eff
            error('Trial %d condition %s: expected %d delays, got %d', ...
                k, condition_names{i}, d_max_eff, numel(r2_i));
        end
        mc_row(i)   = mc_i;
        r2_row(i,:) = r2_i;
        
        % horizon (in seconds): largest delay with R2 > threshold
        idx = find(r2_i > R2_threshold_for_horizon, 1, 'last');
        if isempty(idx); h_row(i) = 0; else; h_row(i) = delay_s(idx); end
        
        if store_internal_results; int_row{i} = res_i; end
        esn_i = [];   % release this condition's memory before building the next
    end
    
    % Sliced assignments -- each output written exactly once per trial k
    MC_trials(k,:)        = mc_row;
    H_trials(k,:)         = h_row;
    R2_trials(k,:,:)      = reshape(r2_row, [1, n_cond, d_max_eff]);
    internal_results(k,:) = int_row;
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
% (Cumulative curves are recomputed inside plot_memory_capacity from these.)

%% -------------------- Paired statistical tests (Total MC) --------------------
% Pairwise comparisons with sign-flip permutation test + paired effect size
% All pairwise comparisons of the 4 conditions (Baseline, SFA, STD, SFA+STD).
pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
pair_labels = { ...
    'Baseline vs SFA', ...
    'Baseline vs STD', ...
    'Baseline vs SFA+STD', ...
    'SFA vs STD', ...
    'SFA vs SFA+STD', ...
    'STD vs SFA+STD' ...
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

%% -------------------- Figures --------------------
% Plotting is factored into plot_memory_capacity so it can be iterated without
% recomputing. It is called after the results are saved, below. To re-plot a
% finished run later:  replot_memory_capacity('<run_tag>_results.mat')

%% -------------------- Save results --------------------
results_all = struct();
results_all.run_tag = run_tag;
results_all.timestamp = timestamp;

results_all.settings = struct();
results_all.settings.n = n;
results_all.settings.f = f;
results_all.settings.fs = fs;
results_all.settings.level_of_chaos = level_of_chaos;
results_all.settings.S_c = 0.35;
results_all.settings.std_zero_floor = true;
results_all.settings.tau_d = tau_d;
results_all.settings.T_wash = T_wash;
results_all.settings.T_train = T_train;
results_all.settings.T_test = T_test;
results_all.settings.d_max = d_max;
results_all.settings.d_max_eff = d_max_eff;      % delays actually scored (hold units)
results_all.settings.delay_s = delay_s;          % delay axis (seconds)
results_all.settings.input_type = input_type;
results_all.settings.T_hold = T_hold;
results_all.settings.hold_len = round(T_hold * fs);
results_all.settings.u_f_cutoff = u_f_cutoff;
results_all.settings.u_alpha = u_alpha;
results_all.settings.ode_solver = 'ode_rk4';
% Nonlinearity fingerprint from a bare construct (esn is a parfor temporary; no
% build needed -- the activation is set by the constructor).
probe_esn = SRNN_ESN_reservoir(base_args_template{:});
results_all.settings.activation = func2str(probe_esn.activation_function);
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

%% -------------------- Plot (from the saved results struct) --------------------
plot_memory_capacity(results_all, out_dir);

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
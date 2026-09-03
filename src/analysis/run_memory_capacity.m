function mat_file = run_memory_capacity(opts)
% RUN_MEMORY_CAPACITY Paired-trial memory-capacity experiment across 4 regimes.
%
%   mat_file = RUN_MEMORY_CAPACITY()
%   mat_file = RUN_MEMORY_CAPACITY('run_mode', 'fast')
%   mat_file = RUN_MEMORY_CAPACITY('run_mode', 'production', 'output_dir', d)
%
% Compares Baseline / SFA / STD / SFA+STD on their ability to reconstruct a
% delayed input from a linear readout, and returns the path of the saved
% <run_tag>_results.mat.
%
% FAIRNESS. Within a trial all four conditions share the same (W, W_in, u(t)) --
% same seeds, and only the adaptation counts differ -- so the comparison is
% paired. verify_shared_build asserts that invariant once, structurally, on the
% first seed pair. Trials differ only by seed, so parfor is order-independent.
%
% TWO KNOBS, KEPT APART, exactly as the sweep pipeline does it:
%   preset_name  WHICH NETWORK. srnn_param_preset('mc_pairs_dualStd') by default.
%   run_mode     HOW MUCH COMPUTE: trials, training/test duration, delay
%                horizon, bootstrap/permutation counts, and fs. See the table
%                in mc_run_config below.
%
% The MC protocol settings live in run_mode rather than in the preset because
% they size the EXPERIMENT, not the network -- which is what lets a 'fast' and a
% 'production' run be the same network measured with different effort.
%
% WAS A SCRIPT (looped_memory_capacity.m) with a 60-line config block, no
% arguments, and no way to be called from a master script.
%
% SRNNCellTypePairs ONLY, since 2026-09-02: SRNN_ESN_reservoir was re-parented
% off SRNNModel2, so memory capacity is no longer the one part of the paper on a
% different model class.
%
% See also: SRNN_ESN_reservoir, plot_memory_capacity, replot_memory_capacity,
%           srnn_param_preset, verify_shared_build

arguments
    opts.preset_name            (1,:) char    = 'mc_pairs_dualStd'
    opts.run_mode               (1,:) char    = 'production'
    opts.output_dir             (1,:) char    = ''
    opts.save_figs              (1,1) logical = true
    opts.verbose                (1,1) logical = true
    opts.store_internal_results (1,1) logical = false
    % Integrator override. Empty means "decide from the preset": a deterministic
    % scheme at sigma_u_noise = 0, the stochastic one above it, the same rule
    % analysis_run_config applies to the sweeps. Name one explicitly to force it
    % -- 'sra1' at sigma = 0 is legal (sde_fixed_step treats an absent noise
    % tensor as zero) and is how an MC run is made to use the same integrator as
    % the rest of the analyses even with noise off.
    opts.ode_solver             (1,:) char    = ''
end

setup_paths();

%% Resolve the two knobs
[preset, model_class, conditions] = srnn_param_preset(opts.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_memory_capacity:BadModelClass', ...
        ['Preset ''%s'' is written for %s, but SRNN_ESN_reservoir now ' ...
         'subclasses SRNNCellTypePairs and cannot take its parameters.'], ...
        opts.preset_name, model_class);
end
cfg = mc_run_config(opts.run_mode, preset, opts.ode_solver);

fs = cfg.fs;
T_wash  = round(cfg.T_wash_sec  * fs);
T_train = round(cfg.T_train_sec * fs);
T_test  = round(cfg.T_test_sec  * fs);
d_max   = round(cfg.d_max_sec   * fs);

fprintf('[memory_capacity] preset=%s  run_mode=%s\n', opts.preset_name, opts.run_mode);
fprintf('[memory_capacity] trials=%d  fs=%d  T_train=%gs  T_test=%gs  d_max=%gs\n', ...
    cfg.n_trials, fs, cfg.T_train_sec, cfg.T_test_sec, cfg.d_max_sec);

%% Output directory
project_root = fileparts(which('setup_paths'));
if isempty(opts.output_dir)
    out_dir = fullfile(project_root, 'data', 'memory_capacity');
else
    out_dir = fullfile(opts.output_dir, 'memory_capacity');
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

timestamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
run_tag = sprintf('MC_%s_%s_trials%d', cfg.input_type, timestamp, cfg.n_trials);

%% Conditions
% From the preset, which states its conditions itself, rather than hardcoded here.
% They used to be four literal {'n_a_E', 0, 'n_b_E', 0} pairs, which was a
% second definition of the four regimes and could drift from the sweeps'.
%
% Names are now the project's snake_case KEYS (no_adaptation, sfa_only,
% std_only, sfa_and_std), not the MC figure's old display strings
% ('Baseline', 'SFA', ...). The keys are what gets saved, so an MC run
% directory names its conditions the same way every other run does;
% srnn_condition_titles supplies the display text at plot time.
condition_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
condition_args  = cellfun(@(c) namevalue(rmfield(c, 'name')), conditions, ...
    'UniformOutput', false);
n_cond = numel(condition_names);

%% Base args: preset (the network) + protocol (the experiment)
% Nothing is stripped from the preset here. On SRNNModel2 tau_a_E had to be
% removed, because n_a_E was the settable property and tau_a_E auto-filled from
% it, so passing both would fight. SRNNCellTypePairs inverts that -- tau_a IS
% the property and n_a is Dependent on it, read-only -- so the conditions carry
% tau_a directly and there is nothing to strip.
base_args_template = [ ...
    namevalue(preset), ...
    {'fs',             fs, ...
     'ode_solver',     cfg.ode_solver, ...
     'input_type',     cfg.input_type, ...
     'u_f_cutoff',     cfg.u_f_cutoff, ...
     'u_alpha',        cfg.u_alpha, ...
     'T_hold',         cfg.T_hold, ...
     'T_wash',         T_wash, ...
     'T_train',        T_train, ...
     'T_test',         T_test, ...
     'd_max',          d_max}];

%% Delay grid + preallocation
% Sized up front so R2_trials is fully allocated before the parfor (a sliced
% output cannot be lazily grown). Mirrors the hold-unit decimation inside
% SRNN_ESN_reservoir/run_memory_capacity.
if strcmpi(cfg.input_type, 'sample_hold')
    hold_len = round(cfg.T_hold * fs);
else
    hold_len = 1;
end
d_max_eff = max(1, floor(d_max / hold_len));
delay_s   = (1:d_max_eff) * (hold_len / fs);

n_trials  = cfg.n_trials;
MC_trials = nan(n_trials, n_cond);
H_trials  = nan(n_trials, n_cond);
R2_trials = nan(n_trials, n_cond, d_max_eff);

seed_net  = cfg.seed_net_base  + (1:n_trials);
seed_stim = cfg.seed_stim_base + (1:n_trials);

internal_results = cell(n_trials, n_cond);
store_internal   = opts.store_internal_results;
readout_signal   = cfg.readout_signal;
R2_thresh        = cfg.R2_threshold_for_horizon;

%% One-time fairness check
fprintf('\nVerifying shared build (fairness check, trial 1)...\n');
chk_args = [{'rng_seeds', [seed_net(1), seed_stim(1)]}, base_args_template];
esn_chk = cell(1, n_cond);
for i = 1:n_cond
    esn_chk{i} = SRNN_ESN_reservoir(chk_args{:}, condition_args{i}{:});
    esn_chk{i}.build();
end
% tau_a and synapse_config are exactly what the conditions vary, so they are
% EXPECTED to differ; everything else -- the network, the input weights, the
% stimulus -- must be shared, which is what makes the trials paired.
verify_shared_build(esn_chk, {'tau_a','synapse_config'}, ...
    {'W','W_in','u_scalar','u_ex','t_ex'});
clear esn_chk chk_args;

%% Main loop: paired trials
fprintf('\n==== Running %d paired trials (%s input) ====\n', n_trials, cfg.input_type);

parfor k = 1:n_trials
    fprintf('--- Trial %d / %d | seeds: net=%d, stim=%d ---\n', ...
        k, n_trials, seed_net(k), seed_stim(k));

    base_args = [{'rng_seeds', [seed_net(k), seed_stim(k)]}, base_args_template];

    mc_row  = nan(1, n_cond);
    h_row   = nan(1, n_cond);
    r2_row  = nan(n_cond, d_max_eff);
    int_row = cell(1, n_cond);

    for i = 1:n_cond
        esn_i = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
        esn_i.build();
        [mc_i, r2_i, res_i] = esn_i.run_memory_capacity( ...
            'store_timeseries', false, 'verbose', false, ...
            'readout_signal', readout_signal);
        r2_i = r2_i(:)';
        if numel(r2_i) ~= d_max_eff
            error('Trial %d condition %s: expected %d delays, got %d', ...
                k, condition_names{i}, d_max_eff, numel(r2_i));
        end
        mc_row(i)   = mc_i;
        r2_row(i,:) = r2_i;

        idx = find(r2_i > R2_thresh, 1, 'last');
        if isempty(idx); h_row(i) = 0; else; h_row(i) = delay_s(idx); end

        if store_internal; int_row{i} = res_i; end
        esn_i = []; %#ok<NASGU>  release before building the next condition
    end

    MC_trials(k,:)        = mc_row;
    H_trials(k,:)         = h_row;
    R2_trials(k,:,:)      = reshape(r2_row, [1, n_cond, d_max_eff]);
    internal_results(k,:) = int_row;
end

fprintf('\n==== Done. Computing summary + figures... ====\n');

%% Summary statistics
MC_mean = mean(MC_trials, 1, 'omitnan');
MC_sem  = std(MC_trials, 0, 1, 'omitnan') / sqrt(n_trials);
H_mean  = mean(H_trials, 1, 'omitnan');
H_sem   = std(H_trials, 0, 1, 'omitnan') / sqrt(n_trials);

MC_ci = bootstrap_mean_ci(MC_trials, cfg.n_boot, 0.05);
H_ci  = bootstrap_mean_ci(H_trials,  cfg.n_boot, 0.05);

R2_mean = squeeze(mean(R2_trials, 1, 'omitnan'));            % [cond x d]
R2_ci   = bootstrap_mean_ci_3d(R2_trials, cfg.n_boot, 0.05); % lo/hi [cond x d]

%% Paired statistical tests (total MC)
pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
pair_labels = {'Baseline vs SFA', 'Baseline vs STD', 'Baseline vs SFA+STD', ...
               'SFA vs STD', 'SFA vs SFA+STD', 'STD vs SFA+STD'};
stats = struct();
for p = 1:size(pairs,1)
    x = MC_trials(:, pairs(p,1));
    y = MC_trials(:, pairs(p,2));
    [pval, diff_mean] = paired_signflip_permutation_pvalue(x, y, cfg.n_perm);
    stats(p).pair       = pair_labels{p};
    stats(p).mean_diff  = diff_mean;
    stats(p).p_perm     = pval;
    stats(p).cohens_dz  = paired_cohens_dz(x, y);
end

%% Assemble and save
% The settings struct is what plot_memory_capacity, plot_memory_capacity_combined
% and the figure READMEs read, and what makes a saved run self-describing. Every
% value below is taken from what ACTUALLY RAN. The old script hardcoded three of
% them and two were wrong: it recorded std_zero_floor = true when base_args set
% false, and ode_solver = 'ode_rk4' (a handle-era name) when the property was
% 'rk4'. It also stored the activation as a func2str'd handle; the nonlinearity
% is named data now, so the name is the honest record.
probe = SRNN_ESN_reservoir(base_args_template{:});

results_all = struct();
results_all.run_tag   = run_tag;
results_all.timestamp = timestamp;

s = struct();
s.preset_name = opts.preset_name;
s.run_mode    = opts.run_mode;
s.model_class = 'SRNN_ESN_reservoir';
s.n = probe.n;  s.f = probe.f;  s.fs = probe.fs;
s.level_of_chaos = probe.level_of_chaos;
s.tau_d = probe.tau_d;
s.S_c = probe.S_c;  s.S_a = probe.S_a;
s.activation = probe.activation;
s.std_zero_floor = probe.std_zero_floor;
s.ode_solver = probe.ode_solver;
s.sigma_u_noise = probe.sigma_u_noise;
% Per-cell-type now, not the E/I scalars: c is a 1 x C row, cell_type_names
% says what its entries mean, and the STD time constants live per ROUTE inside
% synapse_config rather than as tau_b_E_rec/tau_b_E_rel.
s.cell_type_names = probe.cell_type_names;
s.c = probe.c;
s.tau_a = probe.tau_a;
s.synapse_config = probe.synapse_config;
s.T_wash = T_wash;  s.T_train = T_train;  s.T_test = T_test;
s.T_wash_sec = cfg.T_wash_sec;  s.T_train_sec = cfg.T_train_sec;
s.T_test_sec = cfg.T_test_sec;
s.d_max = d_max;  s.d_max_sec = cfg.d_max_sec;
s.d_max_eff = d_max_eff;      % delays actually scored (hold units)
s.delay_s = delay_s;          % delay axis (seconds)
s.input_type = cfg.input_type;
s.T_hold = cfg.T_hold;  s.hold_len = hold_len;
s.u_f_cutoff = cfg.u_f_cutoff;  s.u_alpha = cfg.u_alpha;
s.readout_signal = readout_signal;
s.n_trials = n_trials;
s.R2_threshold_for_horizon = R2_thresh;
s.n_boot = cfg.n_boot;  s.n_perm = cfg.n_perm;
s.seed_net_base = cfg.seed_net_base;  s.seed_stim_base = cfg.seed_stim_base;
results_all.settings = s;

results_all.conditions = condition_names;
results_all.seeds = table(seed_net(:), seed_stim(:), ...
    'VariableNames', {'seed_net','seed_stim'});
results_all.MC_trials = MC_trials;
results_all.H_trials  = H_trials;
results_all.R2_trials = R2_trials;

results_all.summary = struct( ...
    'MC_mean', MC_mean, 'MC_sem', MC_sem, 'MC_ci95', MC_ci, ...
    'H_mean',  H_mean,  'H_sem',  H_sem,  'H_ci95',  H_ci, ...
    'R2_mean', R2_mean, 'R2_ci95', R2_ci, 'stats', stats);

if store_internal
    results_all.internal_results = internal_results;
end

mat_file = fullfile(out_dir, [run_tag '_results.mat']);
save(mat_file, 'results_all', '-v7.3');

% Compact CSV for total MC and horizon
T = array2table([MC_trials H_trials], 'VariableNames', ...
    [strcat('MC_', matlab.lang.makeValidName(condition_names)), ...
     strcat('H_',  matlab.lang.makeValidName(condition_names))]);
T.seed_net  = seed_net(:);
T.seed_stim = seed_stim(:);
T = movevars(T, {'seed_net','seed_stim'}, 'Before', 1);
writetable(T, fullfile(out_dir, [run_tag '_MC_Horizon.csv']));

write_summary_txt(fullfile(out_dir, [run_tag '_summary.txt']), run_tag, s, ...
    condition_names, MC_mean, MC_sem, MC_ci, H_mean, H_sem, H_ci, stats);

%% Plot
if opts.save_figs
    plot_memory_capacity(results_all, out_dir);
end

fprintf('\nSaved:\n  %s\n', mat_file);
end

%% ======================================================================
%  RUN-MODE TABLE
%  ======================================================================
function cfg = mc_run_config(run_mode, preset_defaults, ode_solver_override)
% Cost/fidelity settings for one MC run. The analogue of analysis_run_config,
% kept separate because the MC protocol has no n_levels/n_reps and its cost is
% driven by trials x training duration rather than by a grid.
%
% Trials buy statistical power on the paired tests; T_train buys a well-posed
% readout (N_train = T_train_sec/T_hold hold-samples must exceed n features);
% d_max_sec sets how far the horizon can reach.
%
% preset_defaults is read ONLY to choose the integrator -- see select_solver
% below. It is the same third-argument arrangement analysis_run_config uses, and
% for the same reason: sigma_u_noise lives in the preset while the integrator is
% a run-mode setting, so the one place they meet has to see both.
if nargin < 2; preset_defaults = struct(); end
if nargin < 3; ode_solver_override = ''; end
switch run_mode
    case 'fast'
        % Smoke test: the plumbing, not the numbers. N_train = 60/0.3 = 200
        % hold-samples against n = 300 features, so the readout is
        % UNDER-DETERMINED and the MC values are not meaningful. Ridge
        % regularization keeps it from blowing up. Do not read results from it.
        cfg = pack(4, 10, 60, 30, 5, 200, 500, 200);
    case 'medium'
        % N_train = 300/0.3 = 1000 > 300 features: well posed, usable numbers.
        cfg = pack(15, 10, 300, 90, 10, 1000, 2000, 200);
    case 'production'
        % The paper's settings. N_train = 600/0.3 = 2000 hold-samples.
        cfg = pack(30, 10, 600, 150, 15, 2000, 10000, 200);
    otherwise
        error('run_memory_capacity:badMode', ...
            'Unknown run_mode ''%s'' (expected fast, medium or production).', run_mode);
end

cfg.ode_solver = select_solver(cfg, preset_defaults, ode_solver_override);
end

function solver = select_solver(cfg, preset_defaults, override)
% Pick the integrator, mirroring analysis_run_config's rule.
%
% This used to be a hardcoded cfg.ode_solver = 'rk4', which meant a preset with
% sigma_u_noise > 0 could not run MC at all: the model rejects a deterministic
% scheme above sigma = 0, so the run died at the first build.
%
% Auto (override empty): the deterministic scheme at sigma = 0, the stochastic
% one above it. Explicit: whatever is named, still validated against sigma.
%
% WHY AN OVERRIDE AND NOT JUST AUTO-SELECTION. 'sra1' is legal at sigma = 0 --
% sde_fixed_step reads has_noise = ~isempty(noise) && noise.sigma ~= 0, so with
% no noise tensor it runs as a deterministic two-stage scheme -- and naming it
% is how an MC run uses the SAME integrator as the sweeps even with noise off.
% Auto-selection alone cannot express that, because at sigma = 0 it has no
% reason to prefer sra1 over rk4.
if ~isempty(override)
    solver = override;
else
    is_stochastic = isstruct(preset_defaults) && ...
        isfield(preset_defaults, 'sigma_u_noise') && ...
        any(preset_defaults.sigma_u_noise(:) > 0);
    if is_stochastic
        solver = cfg.sde_solver;
    else
        solver = cfg.det_solver;
    end
end

% Fail here rather than at the first build, and through the check both model
% classes and ParamSpaceAnalysis2's pre-flight already share.
sigma = 0;
if isstruct(preset_defaults) && isfield(preset_defaults, 'sigma_u_noise')
    sigma = preset_defaults.sigma_u_noise;
end
check_noise_settings(sigma, solver, 'run_memory_capacity');
end

function cfg = pack(n_trials, T_wash_sec, T_train_sec, T_test_sec, d_max_sec, ...
        n_boot, n_perm, fs)
cfg = struct();
cfg.n_trials    = n_trials;
cfg.T_wash_sec  = T_wash_sec;
cfg.T_train_sec = T_train_sec;
cfg.T_test_sec  = T_test_sec;
cfg.d_max_sec   = d_max_sec;
cfg.n_boot      = n_boot;
cfg.n_perm      = n_perm;
cfg.fs          = fs;

% Protocol constants -- identical in every mode, so they are not knobs.
% sample_hold: i.i.d. values held for T_hold, so MC is measured in HOLD UNITS.
% That is the fair choice for a low-pass reservoir and is free of the
% autocorrelation inflation that white input produces.
cfg.input_type     = 'sample_hold';
cfg.T_hold         = 0.3;      % s; sets the MC delay increment
cfg.u_f_cutoff     = 5;        % only used by 'bandlimited'
cfg.u_alpha        = 1;        % only used by 'one_over_f'
% The two integrators this run mode would use, deterministic and stochastic.
% select_solver picks between them from the preset's sigma_u_noise, or takes an
% explicit override. Named here rather than chosen here so the choice sits in
% one place -- the same split analysis_run_config uses.
cfg.det_solver     = 'rk4';    % fixed-step; fast and adequate at sigma = 0
cfg.sde_solver     = 'sra1';   % Rossler SRA1: same cost as heun, strong order 1.5
% 'synaptic' reads out br = b.*r, exposing the STD state to the linear readout.
% Unchanged for the no-STD conditions, where br == r.
cfg.readout_signal = 'synaptic';
cfg.R2_threshold_for_horizon = 0.10;
cfg.seed_net_base  = 3000;
cfg.seed_stim_base = 4000;
end

%% ======================================================================
%  HELPERS
%  ======================================================================
function nv = namevalue(s)
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f); nv{2*i-1} = f{i}; nv{2*i} = s.(f{i}); end
end

function write_summary_txt(path, run_tag, s, condition_names, MC_mean, MC_sem, ...
        MC_ci, H_mean, H_sem, H_ci, stats)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));   %#ok<*NASGU> closes fid on any exit path
fprintf(fid, 'Run: %s\n', run_tag);
fprintf(fid, 'Preset: %s | run_mode: %s\n', s.preset_name, s.run_mode);
fprintf(fid, 'Input: %s | n=%d | fs=%d Hz | trials=%d | d_max=%.2f s\n\n', ...
    s.input_type, s.n, s.fs, s.n_trials, s.d_max_sec);
fprintf(fid, 'Total MC (mean +/- SEM) and 95%% CI:\n');
for i = 1:numel(condition_names)
    fprintf(fid, '  %s: %.3f +/- %.3f | CI [%.3f, %.3f]\n', ...
        condition_names{i}, MC_mean(i), MC_sem(i), MC_ci.lo(i), MC_ci.hi(i));
end
fprintf(fid, '\nHorizon (R^2>%.2f) (mean +/- SEM) and 95%% CI:\n', ...
    s.R2_threshold_for_horizon);
for i = 1:numel(condition_names)
    fprintf(fid, '  %s: %.3f +/- %.3f s | CI [%.3f, %.3f]\n', ...
        condition_names{i}, H_mean(i), H_sem(i), H_ci.lo(i), H_ci.hi(i));
end
fprintf(fid, '\nPaired permutation tests (sign-flip) on Total MC:\n');
for p = 1:numel(stats)
    fprintf(fid, '  %s: mean diff = %.3f | p_perm = %.4g | Cohen''s dz = %.3f\n', ...
        stats(p).pair, stats(p).mean_diff, stats(p).p_perm, stats(p).cohens_dz);
end
end

function ci = bootstrap_mean_ci(X, n_boot, alpha)
% X: [N x C]. Returns ci.lo, ci.hi per column (bootstrap CI of the mean).
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
% R2_trials: [N x C x D]. Returns ci.lo/ci.hi as [C x D].
[N, C, D] = size(R2_trials);
boot_means = nan(n_boot, C, D);
for b = 1:n_boot
    idx = randi(N, N, 1);
    boot_means(b,:,:) = reshape(mean(R2_trials(idx,:,:), 1, 'omitnan'), [1, C, D]);
end
ci.lo = reshape(quantile(boot_means, alpha/2,   1), [C, D]);
ci.hi = reshape(quantile(boot_means, 1-alpha/2, 1), [C, D]);
end

function [pval, mean_diff] = paired_signflip_permutation_pvalue(x, y, n_perm)
% Paired permutation test via sign flips on the differences. H0: mean(x-y) = 0.
d = x(:) - y(:);
d = d(~isnan(d));
mean_diff = mean(d);
N = numel(d);
if N == 0; pval = NaN; return; end
perm_stats = zeros(n_perm, 1);
for r = 1:n_perm
    sgn = 2*(rand(N,1) > 0.5) - 1;
    perm_stats(r) = mean(sgn .* d);
end
pval = mean(abs(perm_stats) >= abs(mean_diff));   % two-sided
end

function dz = paired_cohens_dz(x, y)
% Cohen's dz for paired designs: mean(diff)/std(diff).
d = x(:) - y(:);
d = d(~isnan(d));
dz = mean(d) / std(d, 0);
end

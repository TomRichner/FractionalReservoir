%% test_mc_gain_sweep.m
% run_mc_gain_sweep: memory capacity and lambda_1 over a gain-by-bias grid.
%
% A 50-neuron version of the mu7chain memory-capacity preset (model_args),
% a protocol shortened THROUGH ARGUMENTS for this test only (2 s washout,
% 12 s training, 6 s test, delays to 1.5 s; a 6 s Lyapunov run), a 2 x 2
% grid, one seed, two conditions, serial in the client.
%
% Checks:
%   1. The stage runs: the consolidated file, the CSV, one job file per
%      job, provenance; the table has the documented columns, one row per
%      (gain, bias, seed, condition), finite values, nothing failed; the
%      settings carry what the figures read.
%   2. The stage's memory capacity IS run_memory_capacity's. One grid
%      point, both conditions: a reservoir built from the same arguments
%      and scored by the existing full-trajectory path gives the same R^2
%      curve and total, isequal. Tolerance zero: the lean path repeats the
%      existing path's steps bit for bit (test_mc_lean_path says why), and
%      the scoring is the same method.
%   3. A rerun into the same directory runs nothing (job files untouched)
%      and gives the same table; every run_mode name is accepted.
%   4. A larger grid in the same directory runs only the new job.
%   5. A job file from a different protocol is refused, not mixed in.
%   6. A failing job (gain NaN) is recorded as failed with its message and
%      the other job finishes; retry_failed = false then leaves it alone.
%   7. An unknown run_mode, condition or protocol field errors before any
%      compute.
%   8. THE SECOND AXIS IS THE CALLER'S. The bias run's job files carry no
%      trace of the option. A run with second_axis = 'self_feedback' (2
%      gains x 2 k, one seed): the column, the job files and the settings
%      are named by the property; each condition's row is bit-identical to
%      a reservoir built by hand with self_feedback set, on the same drawn
%      W in both conditions (paired); both readers, mc_sweep_quantity, the
%      statistics and every block C figure take the axis and its wording
%      ("self-feedback k (R_W)") from the result; options that do not fit
%      the axis are refused before any compute.
%
% Prints [PASS]/[FAIL] per check and errors at the end if anything failed.
% About 80 s. Assumes setup_paths has run.

n_fail = 0;
preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
model_args = {'n', 50, 'indegree', 25};
protocol = struct('T_wash_sec', 2, 'T_train_sec', 12, 'T_test_sec', 6, 'd_max_sec', 1.5);
conds = {'no_adaptation', 'sfa3_stdc2'};
gains = [0.8 1.6];
biases = [0.1 0.3];
seed = 3001;
common = {'preset_name', preset_name, 'model_args', model_args, 'protocol', protocol, ...
    'conditions', conds, 'seeds', seed, 'lya_T_sec', 6, 'lya_acc_sec', 3, ...
    'use_parallel', false, 'verbose', 'near-none'};

tmp = tempname;
mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));

%% ---- 1. the stage runs ----
fprintf('=== 1. the stage runs ===\n');
[out_dir, results_file] = run_mc_gain_sweep(common{:}, 'gain_levels', gains, ...
    'bias_levels', biases, 'output_dir', tmp);
L = load(results_file);
T = L.T;  s = L.settings;
job_files = dir(fullfile(out_dir, 'jobs', 'job_*.mat'));
n_fail = n_fail + report(strcmp(out_dir, fullfile(tmp, 'mc_gain_sweep')) ...
    && exist(results_file, 'file') == 2 ...
    && exist(fullfile(out_dir, 'mc_gain_sweep_results.csv'), 'file') == 2 ...
    && exist(fullfile(out_dir, 'git_provenance.txt'), 'file') == 2 ...
    && numel(job_files) == 4, ...
    'results .mat, .csv, provenance and 4 job files in <output_dir>/mc_gain_sweep');
want_cols = {'gain', 'bias', 'seed', 'seed_stim', 'condition', 'mc_total', 'horizon', ...
    'R2_delay', 'lle_input', 'lle_noinput', 'rate_E', 'rate_I', 'frac_pinned', ...
    'frac_silent', 'seconds', 'failed', 'message'};
n_fail = n_fail + report(isequal(T.Properties.VariableNames, want_cols), ...
    'the table has the documented columns, in order');
keys = unique(T(:, {'gain', 'bias', 'seed', 'condition'}), 'rows');
n_fail = n_fail + report(height(T) == 8 && height(keys) == 8 ...
    && isequal(unique(T.gain)', gains) && isequal(unique(T.bias)', biases) ...
    && all(T.seed == seed) && all(T.seed_stim == seed + 1000), ...
    'one row per (gain, bias, seed, condition): 2 x 2 x 1 x 2 = 8; stimulus seed = seed + 1000');
num = [T.mc_total, T.horizon, T.lle_input, T.lle_noinput, T.rate_E, T.rate_I, ...
    T.frac_pinned, T.frac_silent, T.seconds];
n_fail = n_fail + report(~any(T.failed) && all(isfinite(num(:))) && all(isfinite(T.R2_delay(:))) ...
    && all(T.mc_total > 0) && all(abs(T.mc_total - sum(T.R2_delay, 2)) < 1e-12) ...
    && all(T.frac_pinned >= 0 & T.frac_pinned <= 1) && all(T.rate_E >= 0 & T.rate_E <= 1), ...
    sprintf('nothing failed, every value finite, MC = sum of the R2 curve (MC %.2f to %.2f)', ...
    min(T.mc_total), max(T.mc_total)));
n_fail = n_fail + report(any(T.lle_input ~= T.lle_noinput), ...
    sprintf('lambda_1 with and without the input differ (%+.3f vs %+.3f in row 1)', ...
    T.lle_input(1), T.lle_noinput(1)));
n_fail = n_fail + report(isequal(s.conditions, conds) && size(T.R2_delay, 2) == numel(s.delay_s) ...
    && numel(s.delay_s) == 5 && abs(s.delay_s(1) - 0.3) < 1e-12 ...
    && strcmp(s.preset_name, preset_name) && s.preset_gain == 1 && s.preset_bias == 0.2 ...
    && s.n == 50 && s.n_jobs == 4 && s.n_failed == 0, ...
    'settings: conditions, delay_s (5 delays of 0.3 s), preset_name, preset_gain 1, preset_bias 0.2');

%% ---- 2. equal to the existing path ----
fprintf('\n=== 2. the stage''s memory capacity is run_memory_capacity''s ===\n');
[d, ~, conditions] = srnn_param_preset(preset_name);
P = mc_protocol();
fs = P.fs;
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
for c = 1:numel(conds)
    cond = conditions{strcmp(names, conds{c})};
    cond_args = struct2namevalue(rmfield(cond, 'name'));
    preset_args = struct2namevalue(d);
    esn = SRNN_ESN_reservoir('rng_seeds', [seed, seed + 1000], preset_args{:}, model_args{:}, ...
        'fs', fs, 'ode_solver', 'sra1', 'input_type', P.input_type, 'T_hold', P.T_hold, ...
        'T_wash', protocol.T_wash_sec * fs, 'T_train', protocol.T_train_sec * fs, ...
        'T_test', protocol.T_test_sec * fs, 'd_max', protocol.d_max_sec * fs, ...
        cond_args{:}, 'level_of_chaos', gains(2), 'S_c', biases(1), 'verbose', 'near-none');
    esn.build();
    [mc_ref, r2_ref] = esn.run_memory_capacity('readout_signal', P.readout_signal, ...
        'store_timeseries', false, 'verbose', false);
    row = T(T.gain == gains(2) & T.bias == biases(1) & strcmp(T.condition, conds{c}), :);
    n_fail = n_fail + report(height(row) == 1 && isequal(row.R2_delay, r2_ref(:)') ...
        && row.mc_total == mc_ref, ...
        sprintf('%s at gain %.1f, bias %.1f: R2 curve and MC bit-identical (MC = %.4f, max |dR2| = %.1e)', ...
        conds{c}, gains(2), biases(1), mc_ref, max(abs(row.R2_delay - r2_ref(:)'))));
end

%% ---- 3. rerun skips; every run mode accepted ----
fprintf('\n=== 3. rerun ===\n');
stamps = sort([job_files.datenum]);
modes = run_mode_names();
same = true;
for k = 1:numel(modes)
    [~, f2] = run_mc_gain_sweep(common{:}, 'gain_levels', gains, 'bias_levels', biases, ...
        'output_dir', tmp, 'run_mode', modes{k});
    L2 = load(f2);
    same = same && isequaln(L2.T, T) && L2.settings.n_run == 0;
end
again = dir(fullfile(out_dir, 'jobs', 'job_*.mat'));
n_fail = n_fail + report(same && isequal(sort([again.datenum]), stamps), ...
    sprintf('a rerun under each of the %d run modes runs no job, leaves the job files alone and gives the same table', ...
    numel(modes)));

%% ---- 4. a larger grid runs only the new job ----
fprintf('\n=== 4. extending the grid ===\n');
[~, f4] = run_mc_gain_sweep(common{:}, 'gain_levels', [gains 1.2], 'bias_levels', biases(1), ...
    'output_dir', tmp);
L4 = load(f4);
old_rows = L4.T(L4.T.gain ~= 1.2, :);
n_fail = n_fail + report(L4.settings.n_run == 1 && L4.settings.n_jobs == 3 && height(L4.T) == 6 ...
    && isequaln(old_rows, T(T.bias == biases(1), :)), ...
    'three gains at one bias: one new job run, the two finished ones read back unchanged');

%% ---- 5. a different experiment is refused ----
fprintf('\n=== 5. a different protocol in the same directory ===\n');
other = protocol;
other.T_train_sec = 10;
args5 = common;
args5{find(strcmp(common, 'protocol')) + 1} = other;
n_fail = n_fail + report(strcmp(error_id(@() run_mc_gain_sweep(args5{:}, 'gain_levels', gains, ...
    'bias_levels', biases, 'output_dir', tmp)), 'run_mc_gain_sweep:DifferentRun'), ...
    'job files made with another protocol raise run_mc_gain_sweep:DifferentRun');

%% ---- 6. a failing job ----
fprintf('\n=== 6. a failing job ===\n');
% A NaN gain builds (NaN weights) and the trajectory is not finite; the
% Lyapunov estimator warns and the stage records the condition as failed.
% The stage's own failure line and those warnings are captured, not
% printed: run_all_tests reads a printed FAIL as a failed test.
tmp6 = fullfile(tmp, 'failing');
args6 = [common, {'gain_levels', [NaN gains(1)], 'bias_levels', biases(1), 'output_dir', tmp6}];
warn_state = warning('off', 'all');
said = evalc('[out6, f6] = run_mc_gain_sweep(args6{:});');
warning(warn_state);
n_fail = n_fail + report(contains(said, 'NonFiniteLyapunov'), ...
    'the failure is announced on the console at any verbosity');
L6 = load(f6);
bad = isnan(L6.T.gain);
n_fail = n_fail + report(nnz(bad) == 2 && all(L6.T.failed(bad)) && all(isnan(L6.T.mc_total(bad))) ...
    && all(~cellfun(@isempty, L6.T.message(bad))) && L6.settings.n_failed == 2, ...
    sprintf('the NaN-gain job is recorded as failed, with its message (%s)', ...
    strtok(L6.T.message{find(bad, 1)}, ' ')));
good = L6.T(~bad, :);
n_fail = n_fail + report(~any(good.failed) && all(isfinite(good.mc_total)) ...
    && isequaln(good(:, {'mc_total', 'R2_delay', 'lle_input'}), ...
    T(T.gain == gains(1) & T.bias == biases(1), {'mc_total', 'R2_delay', 'lle_input'})), ...
    'the other job finished, with the values it has in the first run');
before = dir(fullfile(out6, 'jobs', 'job_*.mat'));
[~, f6b] = run_mc_gain_sweep(common{:}, 'gain_levels', [NaN gains(1)], 'bias_levels', biases(1), ...
    'output_dir', tmp6, 'retry_failed', false);
L6b = load(f6b);
after = dir(fullfile(out6, 'jobs', 'job_*.mat'));
n_fail = n_fail + report(L6b.settings.n_run == 0 && isequal([before.datenum], [after.datenum]), ...
    'retry_failed = false leaves the failed job alone');

%% ---- 7. bad arguments fail before any compute ----
fprintf('\n=== 7. bad arguments ===\n');
n_fail = n_fail + report(strcmp(error_id(@() run_mc_gain_sweep(common{:}, 'run_mode', 'no_such_mode')), ...
    'run_mc_gain_sweep:badMode'), 'an unknown run_mode raises run_mc_gain_sweep:badMode');
args7 = common;
args7{find(strcmp(common, 'conditions')) + 1} = {'no_such_condition'};
n_fail = n_fail + report(strcmp(error_id(@() run_mc_gain_sweep(args7{:})), ...
    'run_mc_gain_sweep:NoSuchCondition'), 'an unknown condition raises run_mc_gain_sweep:NoSuchCondition');
args7 = common;
args7{find(strcmp(common, 'protocol')) + 1} = struct('T_training', 5);
n_fail = n_fail + report(strcmp(error_id(@() run_mc_gain_sweep(args7{:})), ...
    'run_mc_gain_sweep:BadProtocolField'), 'an unknown protocol field raises run_mc_gain_sweep:BadProtocolField');

%% ---- 8. the second axis is the caller's: self-feedback ----
fprintf('\n=== 8. gain by self-feedback ===\n');
% The default axis first: a bias run's job files carry no trace of the
% option (the files of earlier runs must still match).
first = load(fullfile(job_files(1).folder, job_files(1).name), 'job');
n_fail = n_fail + report(isfield(first.job, 'bias') && ~isfield(first.job.signature, 'second_axis') ...
    && ~isempty(regexp(job_files(1).name, '^job_g[\dp]+_b[\dp]+_s3001\.mat$', 'once')) ...
    && strcmp(s.second_axis, 'S_c') && isequal(s.bias_levels, biases), ...
    'the bias run: job field bias, file job_g*_b*_s*, no second_axis in the signature');

ks = [0 1];
tmp8 = fullfile(tmp, 'selffb');
[out8, f8] = run_mc_gain_sweep(common{:}, 'gain_levels', gains, 'second_axis', 'self_feedback', ...
    'second_levels', ks, 'output_dir', tmp8);
L8 = load(f8);
T8 = L8.T;  s8 = L8.settings;
files8 = dir(fullfile(out8, 'jobs', 'job_*.mat'));
n_fail = n_fail + report(isequal(T8.Properties.VariableNames, strrep(want_cols, 'bias', 'self_feedback')) ...
    && height(T8) == 8 && ~any(T8.failed) && isequal(unique(T8.self_feedback)', ks) ...
    && isequal(unique(T8.gain)', gains) && numel(files8) == 4 ...
    && all(startsWith({files8.name}, 'job_g')) && all(contains({files8.name}, '_self_feedback')), ...
    '2 gains x 2 k x 1 seed x 2 conditions: the column and the job files are named by the property');
n_fail = n_fail + report(strcmp(s8.second_axis, 'self_feedback') && isequal(s8.self_feedback_levels, ks) ...
    && s8.preset_self_feedback == 0 && s8.preset_gain == 1 && ~isfield(s8, 'bias_levels') ...
    && ~isfield(s8, 'preset_bias'), ...
    'settings: second_axis, self_feedback_levels, preset_self_feedback 0 (the preset has none)');

% The swept value reached the model, and the conditions are paired: each
% condition's row equals a reservoir built by hand on the same seeds with
% self_feedback = 1 (the existing full-trajectory path), whose W is the
% SAME drawn matrix in both conditions with the diagonal lowered by R_W.
W_drawn = cell(1, numel(conds));
paired = true;
for c = 1:numel(conds)
    cond = conditions{strcmp(names, conds{c})};
    cond_args = struct2namevalue(rmfield(cond, 'name'));
    preset_args = struct2namevalue(d);
    esn = SRNN_ESN_reservoir('rng_seeds', [seed, seed + 1000], preset_args{:}, model_args{:}, ...
        'fs', fs, 'ode_solver', 'sra1', 'input_type', P.input_type, 'T_hold', P.T_hold, ...
        'T_wash', protocol.T_wash_sec * fs, 'T_train', protocol.T_train_sec * fs, ...
        'T_test', protocol.T_test_sec * fs, 'd_max', protocol.d_max_sec * fs, ...
        cond_args{:}, 'level_of_chaos', gains(1), 'self_feedback', ks(2), 'verbose', 'near-none');
    esn.build();
    W_drawn{c} = esn.W;
    shift = full(diag(esn.W) - diag(esn.cached_params.W));
    paired = paired && max(abs(shift - ks(2) * esn.R_gain1)) < 1e-12 && esn.R_gain1 > 0;
    [mc_ref, r2_ref] = esn.run_memory_capacity('readout_signal', P.readout_signal, ...
        'store_timeseries', false, 'verbose', false);
    row = T8(T8.gain == gains(1) & T8.self_feedback == ks(2) & strcmp(T8.condition, conds{c}), :);
    paired = paired && height(row) == 1 && isequal(row.R2_delay, r2_ref(:)') && row.mc_total == mc_ref;
end
paired = paired && isequal(W_drawn{1}, W_drawn{2});
n_fail = n_fail + report(paired, ...
    'k = 1: each condition is the hand-built reservoir with W_ii lowered by R_W, same drawn W in both');
at = @(k) T8.mc_total(T8.gain == gains(1) & T8.self_feedback == k);
n_fail = n_fail + report(~isequal(at(ks(1)), at(ks(2))), 'k = 0 and k = 1 give different memory capacities');

% The readers, the statistics and the figures take the axis from the result.
D8 = read_mc_gain_sweep(f8);
[D8j, info8] = read_mc_gain_sweep_jobs(tmp8);
key8 = {'gain', 'self_feedback', 'seed', 'condition'};
n_fail = n_fail + report(strcmp(D8.second_axis.property, 'self_feedback') && isequal(D8.second_levels, ks) ...
    && D8.preset_second == 0 && ~isfield(D8, 'biases') ...
    && info8.n_read == 4 && isequaln(sortrows(D8j.T, key8), sortrows(D8.T, key8)) ...
    && D8j.preset_second == 0 && strcmp(D8j.second_axis.column, 'self_feedback'), ...
    'read_mc_gain_sweep and read_mc_gain_sweep_jobs: the axis, its levels and the preset''s value');
Qk = mc_sweep_quantity('self_feedback', D8);
n_fail = n_fail + report(strcmp(Qk.label, 'self-feedback k (R_W)') && isequal(Qk.levels, ks) ...
    && strcmp(Qk.at(1), 'self-feedback k = 1 R_W') ...
    && strcmp(error_id(@() mc_sweep_quantity('bias', D8)), 'mc_sweep_quantity:NotInResult'), ...
    'mc_sweep_quantity labels it "self-feedback k (R_W)"; the bias is not a quantity of this result');
fig_dir = fullfile(tmp, 'selffb_figs');
drawn = make_block_C_figures(f8, fig_dir, 'sets', {'both', conds}, 'reference', conds{2}, ...
    'verbose', 'near-none');
S8 = drawn.both.S;
stat_md = fileread(fullfile(fig_dir, 'both', 'MC_gain_bias_statistics.md'));
n_fail = n_fail + report(all(cellfun(@(f) exist(f, 'file') == 2, drawn.both.files)) ...
    && any(contains(drawn.both.files, 'Fig_MC_gain_bias_maps')) ...
    && any(contains(drawn.both.files, 'Fig_MC_gain_curves')) ...
    && any(contains(drawn.both.files, 'Fig_MC_robustness')) ...
    && any(contains(drawn.both.files, 'Fig_MC_vs_lambda1')) ...
    && any(endsWith(drawn.both.files, '_by_self_feedback.csv')) ...
    && all(ismember({'peak_self_feedback', 'gain_oct_fitself_feedback', 'self_feedback_range_best', ...
        'fit_self_feedback'}, S8.Properties.VariableNames)) ...
    && ~any(contains(S8.Properties.VariableNames, 'bias')) && S8.fit_self_feedback(1) == 0 ...
    && contains(stat_md, 'gain and self-feedback') && contains(stat_md, '(R_W)') ...
    && ~contains(stat_md, 'bias') && ~contains(stat_md, 'S_c'), ...
    'make_block_C_figures draws the set; the statistics name self-feedback and R_W, not the bias');
maps = fig_mc_gain_bias_maps('data', D8, 'save', false, 'visible', false);
labels = arrayfun(@(a) strjoin(cellstr(a.YLabel.String), ' | '), findobj(maps.figs, 'Type', 'axes'), ...
    'UniformOutput', false);
sc = fig_mc_sweep_scatter('data', D8, 'x', 'self_feedback', 'y', 'mc_total', 'color_by', 'gain', ...
    'save', false, 'visible', false);
n_fail = n_fail + report(any(contains(labels, 'Self-feedback k (R_W)')) && ~any(contains(labels, 'Bias')) ...
    && height(sc.stats) == numel(conds), ...
    'the maps'' y axis reads "Self-feedback k (R_W)"; the scatter takes self_feedback as x');
close([maps.figs, sc.figs]);

% Options that do not fit the axis fail before any compute.
selffb = [common, {'gain_levels', gains, 'second_axis', 'self_feedback', 'output_dir', fullfile(tmp, 'never')}];
n_fail = n_fail + report( ...
    strcmp(error_id(@() run_mc_gain_sweep(selffb{:})), 'run_mc_gain_sweep:NoSecondLevels') ...
    && strcmp(error_id(@() run_mc_gain_sweep(selffb{:}, 'bias_levels', biases)), 'run_mc_gain_sweep:BiasLevelsOnOtherAxis') ...
    && strcmp(error_id(@() run_mc_gain_sweep(common{:}, 'second_levels', ks, 'bias_levels', biases)), 'run_mc_gain_sweep:TwoLevelOptions') ...
    && strcmp(error_id(@() run_mc_gain_sweep(common{:}, 'second_axis', 'no_such_property', 'second_levels', 1)), 'run_mc_gain_sweep:BadSecondAxis') ...
    && strcmp(error_id(@() run_mc_gain_sweep(common{:}, 'second_axis', 'level_of_chaos', 'second_levels', 1)), 'run_mc_gain_sweep:BadSecondAxis') ...
    && ~isfolder(fullfile(tmp, 'never')), ...
    'no levels, bias_levels on another axis, both level options, an unknown or reserved property: refused');

%% ---- Verdict ----
fprintf('\n');
if n_fail == 0
    fprintf('ALL CHECKS PASSED\n');
else
    error('test_mc_gain_sweep:Failed', '%d check(s) FAILED', n_fail);
end

%% ---- helpers ----
function id = error_id(fn)
id = '';
try
    fn();
catch ME
    id = ME.identifier;
end
end

function failed = report(ok, label)
if ok
    fprintf('[PASS] %s\n', label);
    failed = 0;
else
    fprintf('[FAIL] %s\n', label);
    failed = 1;
end
end

function [out_dir, seconds] = regression_rapid_run(out_dir)
% REGRESSION_RAPID_RUN The ultrarapid regression: the model, a sweep and a memory-capacity run, in about a minute.
%
%   out_dir = regression_rapid_run()          % a new folder under tempdir
%   [out_dir, seconds] = regression_rapid_run(out_dir)
%
% TWO TIERS OF REGRESSION (TR, 2026-11-10). The FULL tier
% (regression_fast_run) runs every analysis stage and every figure of the
% manuscript's bundle on the 100-neuron network; it takes about an hour, so
% it is for milestones. THIS is the ultrarapid tier, cheap enough to run
% after every edit to the model class. It does not go through
% run_all_paper_analyses or any stage function, because the stages size
% themselves by run mode and the smallest mode is the slow part. It calls
% the three classes the pipeline is built on directly:
%
%   1. SRNNCellTypePairs, the regression network (n = 100) in each of the
%      preset's three conditions, noise off (rk4) and noise on (sra1): the
%      full state trajectory, the unpacked time series (rates, synaptic
%      output, adaptation and depression variables), the top-K Lyapunov
%      results, and at three states the analytic Jacobian, its largest real
%      eigenvalue, and the matrix-free Jacobian product.
%   2. ParamSpaceAnalysis2, a 12-job serial sweep (two gains x two reps x
%      three conditions) set up from resolve_run_context as the sweep stages
%      set theirs up.
%   3. SRNN_ESN_reservoir, a short memory-capacity run in each condition on
%      the noise-free twin, with the synaptic readout.
%
% What it writes: OUT_DIR/rapid_model.mat, rapid_sweep.mat, rapid_mc.mat.
% regression_digest reduces them and regression_compare compares, exactly as
% for the full tier. Wall-clock fields are removed before saving.
%
% WHAT IT DOES NOT COVER: the analysis stage functions, the figures, and
% anything at the manuscript's network size. A change that passes here can
% still move a stage's output; the full tier is the check for that.
%
% DO NOT RETUNE the settings below. Like the regression preset, their job is
% to stay the same. If they must change, the reference changes with them and
% the commit says why.
%
% See also: regression_rapid_check, regression_digest, regression_fast_run
arguments
    out_dir (1,:) char = ''
end
t_all = tic;
if isempty(out_dir)
    out_dir = fullfile(tempdir, ['regression_rapid_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
end
if isfolder(out_dir) && numel(dir(out_dir)) > 2
    error('regression_rapid_run:NotEmpty', 'Output folder is not empty: %s', out_dir);
end
if ~isfolder(out_dir); mkdir(out_dir); end

P    = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
P_mc = 'regression_pairs_noStim_noise0_dualStd_3cond';
[~, ~, conditions] = srnn_param_preset(P);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));

%% 1. The model: three conditions x two integrators
run_args = {'rng_seeds', [5 6], 'T_range', [0 6], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', false, 'lya_dt', 0.05, ...
    'lya_T_interval', [2 6], 'lya_warmup', 2, 'verbose', 'near-none'};
integrators = struct('name', {'rk4_noise_off', 'sra1_noise_on'}, ...
    'args', {{'sigma_u_noise', 0, 'ode_solver', 'rk4'}, {'ode_solver', 'sra1'}});
model = struct();
for i = 1:numel(cond_names)
    for s = 1:numel(integrators)
        m = build_from_preset(P, cond_names{i}, run_args{:}, integrators(s).args{:});
        m.run();
        p = m.get_params();
        p.u_interpolant = m.u_interpolant;   % get_params leaves it to run()
        r = struct();
        r.W        = m.W;
        r.S0       = m.S0;
        r.S_out    = m.S_out;
        r.t_out    = m.t_out;
        r.plot_data = m.plot_data;
        r.lya_results = drop_timing(m.lya_results);
        r.lya_summary = drop_timing(m.lya_summary());
        % Jacobian, its leading eigenvalue and the matrix-free product, at
        % three states along the trajectory. Y is fixed, not random.
        N = size(m.S_out, 2);
        Y = cos((1:N)' * (1:4)) / sqrt(N);
        idx = round(size(m.S_out, 1) * [0.4 0.7 1.0]);
        r.state_index = idx;
        r.J = cell(1, numel(idx)); r.JY = cell(1, numel(idx)); r.dS = cell(1, numel(idx));
        r.max_real_eig = zeros(1, numel(idx));
        for k = 1:numel(idx)
            S = m.S_out(idx(k), :)';
            J = SRNNCellTypePairs.compute_Jacobian_fast(S, p);
            r.J{k} = J;
            r.JY{k} = SRNNCellTypePairs.jacobian_times(S, Y, p);
            r.max_real_eig(k) = max(real(eig(full(J))));
            r.dS{k} = SRNNCellTypePairs.dynamics_fast(m.t_out(idx(k)), S, p);
        end
        model.(cond_names{i}).(integrators(s).name) = r;
    end
end
save(fullfile(out_dir, 'rapid_model.mat'), 'model', '-v7');

%% 2. A sweep, set up as the sweep stages set theirs up
ctx = resolve_run_context('param_space', 'preset_name', P, 'run_mode', 'fast', 'verbose', 'near-none');
d = ctx.model_defaults;
d.n = 40; d.indegree = 8; d.F_tracks_network = true;
d.T_range = [0 4]; d.lya_T_interval = [2 4]; d.lya_warmup = 1;
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 6, 'note', 'regression_rapid', ...
    'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = d;
psa.output_dir     = fullfile(out_dir, 'sweep_scratch');
psa.set_conditions(ctx.conditions);
psa.add_grid_parameter('level_of_chaos', [0.5, 1.5]);
psa.add_grid_parameter('reps', 1:2);
psa.run();
sweep = struct();
for i = 1:numel(cond_names)
    rc = psa.results.(cond_names{i});
    sweep.(cond_names{i}) = cellfun(@drop_timing, rc, 'UniformOutput', false);
end
% Only the results are kept: the sweep's own files hold wall-clock times and
% the object, and are not what is being compared.
rmdir(psa.output_dir, 's');
save(fullfile(out_dir, 'rapid_sweep.mat'), 'sweep', '-v7');

%% 3. Memory capacity through the reservoir class
[preset_mc, ~, conditions_mc] = srnn_param_preset(P_mc);
fs = 200;
base_args = [struct2namevalue(preset_mc), ...
    {'fs', fs, 'ode_solver', 'sra1', 'rng_seeds', [3 4], ...
     'input_type', 'sample_hold', 'T_hold', 0.3, ...
     'T_wash', 5 * fs, 'T_train', 30 * fs, 'T_test', 15 * fs, 'd_max', 3 * fs}];
mc = struct();
for i = 1:numel(conditions_mc)
    c = conditions_mc{i};
    cargs = struct2namevalue(rmfield(c, 'name'));
    esn = SRNN_ESN_reservoir(base_args{:}, cargs{:});
    esn.build();
    [MC, R2, res] = esn.run_memory_capacity('readout_signal', 'synaptic', ...
        'verbose', false, 'store_timeseries', false);
    mc.(c.name) = struct('MC', MC, 'R2', R2, 'delay_s', res.delay_s, ...
        'hold_len', res.hold_len, 'W', esn.W, 'W_in', esn.W_in);
end
save(fullfile(out_dir, 'rapid_mc.mat'), 'mc', '-v7');

seconds = toc(t_all);
fprintf('regression_rapid_run: %.1f s -> %s\n', seconds, out_dir);
end

%% ------------------------------------------------------------------------
function s = drop_timing(s)
% Wall-clock fields differ between identical runs and are not results.
if ~isstruct(s); return; end
timing = {'seconds', 'run_duration', 'elapsed', 'wall_time'};
s = rmfield(s, intersect(fieldnames(s), timing));
end

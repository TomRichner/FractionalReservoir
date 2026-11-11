function [out_dir, seconds] = regression_rapid_run(out_dir)
% REGRESSION_RAPID_RUN The ultrarapid regression: the model, a sweep and a memory-capacity run, in about half a minute.
%
%   out_dir = regression_rapid_run()          % a new folder under tempdir
%   [out_dir, seconds] = regression_rapid_run(out_dir)
%
% THREE TIERS OF REGRESSION (TR, 2026-11-10):
%   ultrarapid  THIS                    ~30 s    the three classes, directly
%   stages      regression_stages_run   ~3 min   stage functions and figures
%   full        regression_fast_run     ~50 min  the whole bundle
%
% This is the tier cheap enough to run after every edit to the model class.
% It does not go through run_all_paper_analyses or any stage function,
% because the stages size themselves by run mode and the smallest mode is
% the slow part. It calls the three classes the pipeline is built on
% directly:
%
%   1. SRNNCellTypePairs, the regression network (n = 100) in each of the
%      preset's three conditions, noise off (rk4) and noise on (sra1): the
%      full state trajectory, the unpacked time series (rates, synaptic
%      output, adaptation and depression variables), the top-K Lyapunov
%      results, and at three states the analytic Jacobian, its largest real
%      eigenvalue, and the matrix-free Jacobian product. Also, on the
%      most-adapted condition, the transient gain in its three variants and
%      the finite-difference Jacobian; on a reduced network, the Benettin
%      and QR Lyapunov estimators (the main runs use top-K); and a
%      three-type network with every synapse option the paper's presets do
%      not use (facilitation, the geometric-mean combine, a route weight
%      scale, the zero floor), so those numbers are pinned as well.
%   2. ParamSpaceAnalysis2, a 12-job serial sweep (two gains x two reps x
%      three conditions) set up from resolve_run_context as the sweep stages
%      set theirs up.
%   3. SRNN_ESN_reservoir, a short memory-capacity run in each condition on
%      the noise-free twin, with the synaptic readout.
%
% What it writes: OUT_DIR/rapid_model.mat, rapid_model_extra.mat,
% rapid_sweep.mat, rapid_mc.mat.
% regression_digest reduces them and regression_compare compares, exactly as
% for the full tier. Wall-clock fields are removed before saving.
%
% WHAT IT DOES NOT COVER: the analysis stage functions, the figures, and
% anything at the manuscript's network size. A change that passes here can
% still move a stage's output; the stages tier and the full tier are the
% checks for that.
%
% DO NOT RETUNE the settings below. Like the regression preset, their job is
% to stay the same. If they must change, the reference changes with them and
% the commit says why.
%
% See also: regression_rapid_check, regression_stages_run, regression_digest,
%           regression_fast_run
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

        % 1b. On the most-adapted condition with noise on: the transient
        % gain in its three variants, and the finite-difference Jacobian.
        if i == numel(cond_names) && s == numel(integrators)
            extra = struct();
            i0 = find(m.t_out >= 3, 1);
            for variant = {'frozen_x', 'frozen_full', 'active'}
                [G, g_info] = SRNNCellTypePairs.transient_gain(m.S_out, m.t_out, i0, p, 0.5, ...
                    struct('variant', variant{1}));
                extra.transient_gain.(variant{1}) = struct('G', G, 'info', drop_timing(g_info));
            end
            S = m.S_out(idx(end), :)';
            Jfd = SRNNCellTypePairs.finite_difference_jacobian(S, p, 1e-6);
            extra.finite_difference.J = Jfd;
            extra.finite_difference.max_abs_vs_analytic = ...
                full(max(abs(Jfd - SRNNCellTypePairs.compute_Jacobian_fast(S, p)), [], 'all'));
        end
    end
end
save(fullfile(out_dir, 'rapid_model.mat'), 'model', '-v7');

%% 1c. The other two Lyapunov estimators, on a reduced network
% Benettin re-integrates the dynamics and QR integrates the variational
% equation with the Jacobian, each by its own route; the main runs above use
% top-K only. Small n because QR carries N^2 equations.
small = {'n', 20, 'indegree', 5, 'F_tracks_network', true, 'rng_seeds', [5 6], ...
    'T_range', [0 3], 'fs', 200, 'lya_T_interval', [1.5 3], 'lya_warmup', 1, 'verbose', 'near-none'};
estimators = struct();
for i = [1, numel(cond_names)]
    mb = build_from_preset(P, cond_names{i}, small{:}, 'ode_solver', 'sra1', 'lya_method', 'benettin');
    mb.run();
    mq = build_from_preset(P, cond_names{i}, small{:}, 'sigma_u_noise', 0, 'ode_solver', 'rk4', 'lya_method', 'qr');
    mq.run();
    estimators.(cond_names{i}) = struct('benettin', drop_timing(mb.lya_results), ...
        'qr', drop_timing(mq.lya_results));
end

%% 1d. Synapse options the paper's presets do not use
% Three cell types with facilitation, a geometric-mean combine, a route
% weight scale and the zero floor all live, so their numbers are pinned too.
sc = struct();
sc.E.E.std   = struct('tau_rec', [0.3 1], 'tau_rel', 0.2, 'combine', 'geomean');
sc.E.PV.stf  = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
sc.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
sc.E.SST.stf = struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
sc.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
sc.PV.PV.scale = 2;
sc.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
mo = SRNNCellTypePairs( ...
    'n', 30, 'indegree', 10, 'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], 'mu_tilde_relative', [3 -4 -2], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {[0.25 1], [], 2}, 'c', [0.3 0 0.2], 'synapse_config', sc, ...
    'std_zero_floor', true, 'rng_seeds', [5 6], 'T_range', [0 3], 'fs', 200, ...
    'ode_solver', 'rk4', 'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, ...
    'lya_T_interval', [1.5 3], 'lya_warmup', 1, 'store_full_state', true, 'verbose', 'near-none');
mo.build(); mo.run();
po = mo.get_params(); po.u_interpolant = mo.u_interpolant;
So = mo.S_out(end, :)';
No = numel(So);
options = struct('W', mo.W, 'S_out', mo.S_out, 'plot_data', mo.plot_data, ...
    'lya_results', drop_timing(mo.lya_results), ...
    'J', SRNNCellTypePairs.compute_Jacobian_fast(So, po), ...
    'JY', SRNNCellTypePairs.jacobian_times(So, cos((1:No)' * (1:3)) / sqrt(No), po), ...
    'dS', SRNNCellTypePairs.dynamics_fast(mo.t_out(end), So, po));
save(fullfile(out_dir, 'rapid_model_extra.mat'), 'extra', 'estimators', 'options', '-v7');

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

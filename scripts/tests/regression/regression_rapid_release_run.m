function [out_dir, seconds] = regression_rapid_release_run(out_dir)
% REGRESSION_RAPID_RELEASE_RUN The ultrarapid regression for release noise: four conditions with the gamma process on, both estimators, a lean memory-capacity run and a sweep, in under a minute.
%
%   out_dir = regression_rapid_release_run()          % a new folder under tempdir
%   [out_dir, seconds] = regression_rapid_release_run(out_dir)
%
% THE RELEASE NOISE'S OWN EXACT REFERENCE (FD-027, docs/Plan_release_noise.md,
% section 4). The other exact tiers (regression_rapid_run,
% regression_rapid_chain_run, regression_stages_run, regression_fast_run)
% hold no run with spontaneous release: they say "with release_rate and
% release_sigma at zero nothing moved" and nothing about the mechanism.
% This is their twin for runs with it on, with its own reference
% (scripts/tests/fixtures/regression_rapid_release_digest.mat), captured
% when the release noise was final (step N3). It is the bit-for-bit gate
% for any edit that should move no release-noise number: plumbing, the
% reservoir class, the sweep driver, display code.
%
% Like the other ultrarapid tiers it calls the three classes directly, no
% stage function. The networks are built here from the regression presets'
% defaults; the chain condition is the regression network's sfa3_std2 with
% its four routes switched to the chain, as in regression_rapid_chain_run.
%
%   1. SRNNCellTypePairs, integrator 'sra1'.
%      a. The four conditions of the paper on the regression network
%         (n = 100): no_adaptation (pass-through routes: the release is the
%         increment itself), sfa1_std1 (one pool), chain (the chained
%         pool) and sfa3_std2 (the product of two pools), with the Wiener
%         term off and release noise on in both cell types (r_s = 0.01,
%         sigma_gamma = 0.007). Per case: the full state trajectory, the
%         unpacked time series, the top-K Lyapunov results and summary,
%         the centered release increments, the realized release of every
%         route over every step (release_per_step: the mean part and what
%         each kick released), and at three states the analytic Jacobian
%         of the mean dynamics, its largest real eigenvalue, the
%         matrix-free product and the right-hand side.
%      b. chain_both: the chain condition with BOTH noise sources (Wiener
%         0.025 and the release noise), and a noise_seed_schedule that
%         switches both streams at 3 s.
%      c. chain_leak: the chain condition with the leak alone
%         (sigma_gamma = 0), the mean dynamics at r_s = 0.01.
%      d. mixed: four cell types (n = 40; regression_rapid_chain_run's
%         mixed network) with a different release in each presynaptic type
%         (r_s [0.01 0.02 0.005 0.01], sigma_gamma [0.007 0.01 0 0.005]: one
%         type has the leak alone) and the Wiener term on. Its routes are
%         a chain, a chain with a route scale, a two-pool product, a
%         geometric-mean product, single pools, an undepressed route and a
%         facilitating route without a pool: every form of the kick.
%   2. Benettin on the four conditions of 1a (the perturbed trajectory
%      gets the kicks of the same steps), and on chain_both.
%   3. ParamSpaceAnalysis2: an 8-job serial sweep at n = 40 with
%      release_sigma_all as a grid axis (0 and 0.014) x two reps x two
%      conditions (the chain and the single pool), r_s = 0.01 a model
%      default.
%   4. SRNN_ESN_reservoir: lean memory-capacity runs with the
%      'release_filtered' readout (the release through a tau_d filter,
%      sampled at the end of each hold) on the chain and the product
%      condition with release noise, and on the chain with both sources;
%      and the chain once more with the default 'synaptic' readout.
%
% What it writes: OUT_DIR/rapid_release_model.mat, rapid_release_extra.mat,
% rapid_release_sweep.mat, rapid_release_mc.mat. regression_digest reduces
% them and regression_compare compares, exactly as for the other tiers.
% Wall-clock fields and every field named params are removed before saving
% (regression_rapid_chain_run says why).
%
% DO NOT RETUNE the settings below. Their job is to stay the same; if they
% must change, the reference changes with them and the commit says why.
%
% See also: regression_rapid_release_check, regression_rapid_chain_run,
%           regression_tol_run, regression_digest, make_regression_reference
arguments
    out_dir (1,:) char = ''
end
t_all = tic;
if isempty(out_dir)
    out_dir = fullfile(tempdir, ['regression_rapid_release_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
end
if isfolder(out_dir) && numel(dir(out_dir)) > 2
    error('regression_rapid_release_run:NotEmpty', 'Output folder is not empty: %s', out_dir);
end
if ~isfolder(out_dir); mkdir(out_dir); end

P    = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
P_mc = 'regression_pairs_noStim_noise0_dualStd_3cond';
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
sc4 = struct();
sc4.E.E.std = chain; sc4.E.I.std = chain; sc4.I.E.std = chain; sc4.I.I.std = chain;
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));

R_S = 0.01;             % spontaneous release rate, both cell types
SIGMA_G = 0.007;        % release-noise amplitude, both cell types
release = {'release_rate', [R_S R_S], 'release_sigma', [SIGMA_G SIGMA_G]};

%% 1. The model class
run_args = {'rng_seeds', [5 6], 'T_range', [0 6], 'fs', 200, 'store_full_state', true, ...
    'ode_solver', 'sra1', 'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', false, 'lya_dt', 0.05, ...
    'lya_T_interval', [2 6], 'lya_warmup', 2, 'verbose', 'near-none'};
% name, preset, condition of the preset, extra arguments
cases = { ...
    'no_adaptation', P_mc, 'no_adaptation', release; ...
    'sfa1_std1',     P_mc, 'sfa1_std1',     release; ...
    'chain',         P_mc, 'sfa3_std2',     [{'synapse_config', sc4}, release]; ...
    'sfa3_std2',     P_mc, 'sfa3_std2',     release; ...
    'chain_both',    P,    'sfa3_std2',     [{'synapse_config', sc4, 'noise_seed_schedule', [3, 91]}, release]; ...
    'chain_leak',    P_mc, 'sfa3_std2',     {'synapse_config', sc4, 'release_rate', [R_S R_S]}; ...
    'mixed',         '',   '',              {}};
model = struct();
for q = 1:size(cases, 1)
    if strcmp(cases{q, 1}, 'mixed')
        m = mixed_model(run_args);
    else
        m = build_from_preset(cases{q, 2}, cases{q, 3}, run_args{:}, cases{q, 4}{:});
    end
    m.run();
    model.(cases{q, 1}) = model_results(m);
end
model = results_only(model);
save(fullfile(out_dir, 'rapid_release_model.mat'), 'model', '-v7');

%% 2. Benettin: the perturbed trajectory gets the same kicks
% The perturbation's direction comes from the global stream, which is
% seeded here per case.
benettin = struct();
for q = 1:5
    m = build_from_preset(cases{q, 2}, cases{q, 3}, run_args{:}, cases{q, 4}{:}, ...
        'lya_method', 'benettin', 'store_full_state', false);
    rng(17, 'twister');
    m.run();
    benettin.(cases{q, 1}) = m.lya_results;
end
benettin = results_only(benettin);
save(fullfile(out_dir, 'rapid_release_extra.mat'), 'benettin', '-v7');

%% 3. A sweep with the release noise as a grid axis
ctx = resolve_run_context('param_space', 'preset_name', P_mc, 'run_mode', 'fast', 'verbose', 'near-none');
c_one  = ctx.conditions{2};
c_full = ctx.conditions{3};
sweep_conditions = { ...
    struct('name', 'chain_all', 'tau_a', {c_full.tau_a}, 'synapse_config', sc4), ...
    struct('name', 'sfa1_std1', 'tau_a', {c_one.tau_a},  'synapse_config', c_one.synapse_config)};
d = ctx.model_defaults;
d.n = 40; d.indegree = 8; d.F_tracks_network = true;
d.T_range = [0 4]; d.lya_T_interval = [2 4]; d.lya_warmup = 1;
d.ode_solver = 'sra1';                   % the kicks need a stochastic integrator
d.release_rate = [R_S R_S];
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 4, 'note', 'regression_rapid_release', ...
    'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = d;
psa.output_dir     = fullfile(out_dir, 'sweep_scratch');
psa.set_conditions(sweep_conditions);
psa.add_grid_parameter('release_sigma_all', [0, 0.014]);
psa.add_grid_parameter('reps', 1:2);
psa.run();
sweep = struct();
for i = 1:numel(sweep_conditions)
    sweep.(sweep_conditions{i}.name) = psa.results.(sweep_conditions{i}.name);
end
rmdir(psa.output_dir, 's');
sweep = results_only(sweep);
save(fullfile(out_dir, 'rapid_release_sweep.mat'), 'sweep', '-v7');

%% 4. Memory capacity on the lean path, the release read through the tau_d filter
[preset_mc, ~, conditions_mc] = srnn_param_preset(P_mc);
fs = 200;
c_chain = conditions_mc{end};
c_chain.synapse_config = sc4;
mc_cases = { ...
    'chain',          conditions_mc{end}.name, c_chain,            {},                      'release_filtered'; ...
    'sfa3_std2',      conditions_mc{end}.name, conditions_mc{end}, {},                      'release_filtered'; ...
    'chain_both',     conditions_mc{end}.name, c_chain,            {'sigma_u_noise', 0.025}, 'release_filtered'; ...
    'chain_synaptic', conditions_mc{end}.name, c_chain,            {},                      'synaptic'};
base_args = struct2namevalue(preset_mc);
mc = struct();
for q = 1:size(mc_cases, 1)
    cargs = struct2namevalue(rmfield(mc_cases{q, 3}, 'name'));
    esn = SRNN_ESN_reservoir(base_args{:}, 'fs', fs, 'ode_solver', 'sra1', 'rng_seeds', [3 4], ...
        'input_type', 'sample_hold', 'T_hold', 0.3, 'verbose', 'near-none', ...
        'T_wash', 5 * fs, 'T_train', 30 * fs, 'T_test', 15 * fs, 'd_max', 3 * fs, ...
        'lya_method', 'none', 'lean_input', true, cargs{:}, release{:}, mc_cases{q, 4}{:});
    esn.build();
    [MC, R2, res] = esn.run_memory_capacity_lean('readout_signal', mc_cases{q, 5}, ...
        'segment_len', 700, 'verbose', false);
    mc.(mc_cases{q, 1}) = struct('MC', MC, 'R2', R2, 'results', res, 'W', esn.W, 'W_in', esn.W_in);
end
mc = results_only(mc);
save(fullfile(out_dir, 'rapid_release_mc.mat'), 'mc', '-v7');

seconds = toc(t_all);
fprintf('regression_rapid_release_run: %.1f s -> %s\n', seconds, out_dir);
end

%% ------------------------------------------------------------------------
function r = model_results(m)
% What one run of the model class leaves, and what the release noise did in it.
p = m.get_params();
p.u_interpolant = m.u_interpolant;   % get_params leaves it to run()
r = struct();
r.W        = m.W;
r.S0       = m.S0;
r.S_out    = m.S_out;
r.t_out    = m.t_out;
r.plot_data = m.plot_data;
r.lya_results = m.lya_results;
r.lya_summary = m.lya_summary();
% The centered increments the run consumed (drawn again from the seed: the
% run clears them), and the realized release of every route per step.
[rate, sigma] = m.release_noise_rows();
n_steps = numel(m.t_out) - 1;
dRc = release_increments(rate, sigma, m.fs, [1, n_steps], m.noise_base_seed(), ...
    m.noise_seed_schedule, m.t_out(1), n_steps, 'SRNNCellTypePairs') - rate .* (1 / m.fs);
r.dRc = dRc;
if ~any(sigma > 0); dRc = []; end
[r.release, ~, r.kicked] = SRNNCellTypePairs.release_per_step(m.S_out, dRc, p, 1 / m.fs);
% Jacobian of the mean dynamics, its leading eigenvalue, the matrix-free
% product and the right-hand side at three states. Y is fixed, not random.
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
% One kick written out: the state after it from a stored state and a fixed
% increment (a zero draw, a small and a large one in turn over the neurons).
S = m.S_out(idx(2), :)';
fixed = repmat([-1; 0.3; 40], ceil(m.n / 3), 1);
fixed = fixed(1:m.n) .* rate .* (1 / m.fs);
r.kick_state = SRNNCellTypePairs.release_kick(S, fixed, p);
end

function m = mixed_model(run_args)
% regression_rapid_chain_run's mixed network (four cell types; from the
% presynaptic type E a chain, a two-component product, a one-component and
% an undepressed route), written out again so the references do not depend
% on one another's file, with a different release in each presynaptic type
% and the Wiener term on.
chain = @(tau_rec, tau_rel, n2) struct('model', 'chain', 'tau_rec', tau_rec, 'tau_rel', tau_rel, 'n2', n2);
sc = struct();
sc.E.E.std     = chain([2 5], 0.25, 2);
sc.E.PV.std    = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
sc.E.SST.std   = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.PV.E.std    = chain([1 3], 0.3, 1.5);
sc.PV.E.scale  = 1.3;
sc.PV.PV.std   = struct('tau_rec', 0.6, 'tau_rel', 0.3);
sc.SST.E.std   = struct('tau_rec', [0.3 1], 'tau_rel', 0.2, 'combine', 'geomean');
sc.SST.PV.stf  = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
sc.VIP.SST.std = chain([0.5 8], 0.2, 0.5);
sc.VIP.VIP.std = chain([1.5 4], 0.4, 1e6);
m = SRNNCellTypePairs('n', 40, 'indegree', 10, 'n_cellTypes', 4, ...
    'cell_type_names', {'E', 'PV', 'SST', 'VIP'}, ...
    'f', [0.4 0.2 0.2 0.2], 'mu_tilde_relative', [3 -4 -2 -1.5], 'sigma_tilde_relative', [1 1 1 1], ...
    'tau_a', {[0.25 1], [], 2, 1}, 'c', [0.3 0 0.2 0.1], 'synapse_config', sc, ...
    'level_of_chaos', 1, 'input_config', struct('n_steps', 6, 'step_density', struct(), 'amp', 0.5, ...
        'no_stim_pattern', false(1, 6), 'intrinsic_drive', 0, 'positive_only', false), ...
    run_args{:}, 'sigma_u_noise', 0.025, ...
    'release_rate', [0.01 0.02 0.005 0.01], 'release_sigma', [0.007 0.01 0 0.005]);
m.build();
end

function v = results_only(v)
% Remove, at any depth of nested structs and cells, the wall-clock fields
% (they differ between identical runs) and every field named params (the
% model's internal packing; regression_rapid_chain_run's header).
drop = {'seconds', 'run_duration', 'elapsed', 'wall_time', 'params'};
if isstruct(v)
    v = rmfield(v, intersect(fieldnames(v), drop));
    names = fieldnames(v);
    for k = 1:numel(v)
        for j = 1:numel(names)
            x = v(k).(names{j});
            if isstruct(x) || iscell(x)
                v(k).(names{j}) = results_only(x);
            end
        end
    end
elseif iscell(v)
    for k = 1:numel(v)
        if isstruct(v{k}) || iscell(v{k})
            v{k} = results_only(v{k});
        end
    end
end
end

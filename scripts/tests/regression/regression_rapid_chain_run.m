function [out_dir, seconds] = regression_rapid_chain_run(out_dir)
% REGRESSION_RAPID_CHAIN_RUN The ultrarapid regression for chained depression: chain networks, a sweep and a memory-capacity run, in about half a minute.
%
%   out_dir = regression_rapid_chain_run()          % a new folder under tempdir
%   [out_dir, seconds] = regression_rapid_chain_run(out_dir)
%
% THE CHAIN'S OWN EXACT REFERENCE (Plan A, section 7.3, step 5). The three
% byte-exact tiers (regression_rapid_run, regression_stages_run,
% regression_fast_run) hold no route with chained vesicle-pool depression
% (std.model = 'chain', route_pool): they say "the networks the paper had
% are unchanged" and nothing about the chain. This is regression_rapid_run's
% twin for chain networks, with its own reference
% (scripts/tests/fixtures/regression_rapid_chain_digest.mat), captured when
% the chain was first merged. It is the bit-for-bit gate for any edit that
% should move no chain number: plumbing, the reservoir class, the sweep
% driver, and chained facilitation with facilitation off.
%
% Like the ultrarapid tier it calls the three classes directly, no stage
% function. The networks are built here from the regression presets'
% defaults plus a chain synapse_config; no preset is involved in WHICH
% routes are chains.
%
%   1. SRNNCellTypePairs.
%      a. chain: the regression network (n = 100) in its sfa3_std2
%         condition with all four routes switched to the chain (tau_rec
%         [2 5], tau_rel 0.25, n2 2), noise off (rk4) and on (sra1): the
%         full state trajectory, the unpacked time series, every route's
%         depression factor as display code gets it (std_factors), the top-K
%         Lyapunov results and summary, and at three states the analytic
%         Jacobian, its largest real eigenvalue, the matrix-free Jacobian
%         product and the right-hand side. With noise on, also the
%         transient gain in its three variants and the finite-difference
%         Jacobian.
%      b. mixed: four cell types (n = 40) in which one presynaptic type
%         has a chain, a two-component product, a one-component and an
%         undepressed route, beside a chain with a route scale, a
%         geometric-mean product, facilitation on a product route, and
%         chains with a small and a huge reserve (regression_tol_run's
%         mixed network at gain 1). Noise off and on; the same quantities.
%      c. estimators: Benettin (noise on) and QR (noise off) on a reduced
%         chain network (n = 12).
%      d. chaotic: the chain network at gain 2, noise off, where it is
%         clearly chaotic (lambda_1 about +1.2 1/s; its product twin about
%         +0.2 to +0.5), top-K and Benettin over [4, 12] s. The full-spectrum
%         QR is left out: on a chaotic network of this kind it took 14
%         minutes at n = 40.
%      e. WITH THE PER-NEURON SPREAD (FD-015; Plan A step 5 asked for the
%         chain's exact reference "with and without the spread"; added
%         2026-11-12, when the spread was built). Every case above except
%         the chaotic pair, again with the spread on, under keys of their
%         own: chain_spread (1a with std_spread [0.25 0.25], the
%         manuscript's setting; the ordering rule moves two E and two I
%         neurons) and mixed_spread (1b with std_spread [0.25 0.4 0.1 0.3]
%         and stf_spread on its one facilitated type), both integrators,
%         the same quantities plus the standard normals drawn, std_clipped
%         and the value each neuron has of every route parameter;
%         estimators_spread (1c); and in 2 and 3 a second sweep with the
%         spread as a model default (sweep_spread) and a second
%         memory-capacity run (mc.chain_spread). The cases without the
%         spread are computed first and exactly as before: adding these
%         changed no leaf the reference already held.
%   2. ParamSpaceAnalysis2: an 8-job serial sweep (two gains x two reps x
%      two conditions) at n = 40. The conditions are made here: chain_all
%      (the chain network) and chain_EE (one SFA timescale, a chain on
%      E->E beside one-component product routes). Then the same sweep
%      with the spread (1e).
%   3. SRNN_ESN_reservoir: a short memory-capacity run with the synaptic
%      readout on the noise-free chain network (four identical chains, so
%      the readout is defined). Then the same with the spread (1e).
%
% What it writes: OUT_DIR/rapid_chain_model.mat, rapid_chain_extra.mat,
% rapid_chain_sweep.mat, rapid_chain_mc.mat. regression_digest reduces
% them and regression_compare compares, exactly as for the other tiers.
%
% WHAT IS REMOVED BEFORE SAVING, at any depth: wall-clock fields, and ANY
% FIELD NAMED params. The existing ultrarapid digest holds the params
% struct that Benettin and QR results carry, the model's internal packing,
% which is why the route-interface refactor showed up there as structure
% (new keys, changed shapes) and not only as values. A reference should
% hold results, not the packing; the tolerance tier's parameter view is
% where the packed physics is checked.
%
% DO NOT RETUNE the settings below. Their job is to stay the same; if they
% must change, the reference changes with them and the commit says why.
%
% See also: regression_rapid_chain_check, regression_rapid_run,
%           regression_tol_run, regression_digest, make_regression_reference
arguments
    out_dir (1,:) char = ''
end
t_all = tic;
if isempty(out_dir)
    out_dir = fullfile(tempdir, ['regression_rapid_chain_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
end
if isfolder(out_dir) && numel(dir(out_dir)) > 2
    error('regression_rapid_chain_run:NotEmpty', 'Output folder is not empty: %s', out_dir);
end
if ~isfolder(out_dir); mkdir(out_dir); end

P    = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
P_mc = 'regression_pairs_noStim_noise0_dualStd_3cond';
full_cond = 'sfa3_std2';                 % three SFA timescales; its routes are replaced
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
sc4 = struct();
sc4.E.E.std = chain; sc4.E.I.std = chain; sc4.I.E.std = chain; sc4.I.I.std = chain;
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));

%% 1a, 1b. The chain network and the mixed network, two integrators each
run_args = {'rng_seeds', [5 6], 'T_range', [0 6], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', false, 'lya_dt', 0.05, ...
    'lya_T_interval', [2 6], 'lya_warmup', 2, 'verbose', 'near-none'};
integrators = struct('name', {'rk4_noise_off', 'sra1_noise_on'}, ...
    'args', {{'sigma_u_noise', 0, 'ode_solver', 'rk4'}, {'ode_solver', 'sra1'}});
model = struct();
extra = struct();
% The last two are the first two with the per-neuron spread on (section 1e
% of the header). They come after, and add keys; nothing of the first two
% depends on them.
nets = {'chain', 'mixed', 'chain_spread', 'mixed_spread'};
spread_chain = {'std_spread', [0.25 0.25]};
spread_mixed = {'std_spread', [0.25 0.4 0.1 0.3], 'stf_spread', [0 0 0.3 0]};
for q = 1:numel(nets)
    for s = 1:numel(integrators)
        switch q
            case 1
                m = build_from_preset(P, full_cond, run_args{:}, integrators(s).args{:}, 'synapse_config', sc4);
            case 2
                m = mixed_model(integrators(s).args{end}, run_args);
            case 3
                m = build_from_preset(P, full_cond, run_args{:}, integrators(s).args{:}, 'synapse_config', sc4, ...
                    spread_chain{:});
            case 4
                m = mixed_model(integrators(s).args{end}, [run_args, spread_mixed]);
        end
        m.run();
        p = m.get_params();
        p.u_interpolant = m.u_interpolant;   % get_params leaves it to run()
        r = struct();
        r.W        = m.W;
        r.S0       = m.S0;
        r.S_out    = m.S_out;
        r.t_out    = m.t_out;
        r.plot_data = m.plot_data;
        % Every route's depression factor as display code gets it
        % (std_factors -> std_route_factor): b1 on a chain route.
        r.std_factor = SRNNCellTypePairs.std_factors(m.plot_data.b, p);
        r.lya_results = m.lya_results;
        r.lya_summary = m.lya_summary();
        % Jacobian, its leading eigenvalue, the matrix-free product and the
        % right-hand side at three states. Y is fixed, not random.
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
        if q >= 3
            % What the spread drew, whom the ordering rule moved, and the
            % value each neuron has of every route parameter.
            r.std_draws = m.std_draws;
            r.stf_draws = m.stf_draws;
            r.std_clipped = m.std_clipped;
            r.realised = realised_parameters(m, p);
        end
        model.(nets{q}).(integrators(s).name) = r;

        % With noise on: the transient gain in its three variants, and the
        % finite-difference Jacobian.
        if s == numel(integrators)
            e = struct();
            i0 = find(m.t_out >= 3, 1);
            for variant = {'frozen_x', 'frozen_full', 'active'}
                [G, g_info] = SRNNCellTypePairs.transient_gain(m.S_out, m.t_out, i0, p, 0.5, ...
                    struct('variant', variant{1}));
                e.transient_gain.(variant{1}) = struct('G', G, 'info', g_info);
            end
            S = m.S_out(idx(end), :)';
            Jfd = SRNNCellTypePairs.finite_difference_jacobian(S, p, 1e-6);
            e.finite_difference.J = Jfd;
            e.finite_difference.max_abs_vs_analytic = ...
                full(max(abs(Jfd - SRNNCellTypePairs.compute_Jacobian_fast(S, p)), [], 'all'));
            extra.(nets{q}) = e;
        end
    end
end
model = results_only(model);
save(fullfile(out_dir, 'rapid_chain_model.mat'), 'model', '-v7');

%% 1c. The other two Lyapunov estimators, on a reduced chain network
% n = 12 because QR integrates N^2 variational equations with ode45 and is
% the slowest part of this function: 17 s at n = 20 (N = 160), which alone
% would be half the run.
small = {'n', 12, 'indegree', 4, 'F_tracks_network', true, 'rng_seeds', [5 6], ...
    'T_range', [0 3], 'fs', 200, 'lya_T_interval', [1.5 3], 'lya_warmup', 1, 'verbose', 'near-none', ...
    'synapse_config', sc4};
mb = build_from_preset(P, full_cond, small{:}, 'ode_solver', 'sra1', 'lya_method', 'benettin');
mb.run();
mq = build_from_preset(P, full_cond, small{:}, 'sigma_u_noise', 0, 'ode_solver', 'rk4', 'lya_method', 'qr');
mq.run();
estimators = struct('benettin', mb.lya_results, 'qr', mq.lya_results);
% ... and the same two with the spread on (1e).
mb = build_from_preset(P, full_cond, small{:}, 'ode_solver', 'sra1', 'lya_method', 'benettin', spread_chain{:});
mb.run();
mq = build_from_preset(P, full_cond, small{:}, 'sigma_u_noise', 0, 'ode_solver', 'rk4', 'lya_method', 'qr', ...
    spread_chain{:});
mq.run();
estimators_spread = struct('benettin', mb.lya_results, 'qr', mq.lya_results);

%% 1d. A chain network that is clearly chaotic: top-K and Benettin
% Gain 2 on the n = 100 chain network, noise off (measured 2026-11-12:
% lambda_1 = +1.23 by top-K and +1.23 by Benettin over [5, 20] s; the
% product network at the same gain +0.19 and +0.18).
chaos = {'rng_seeds', [5 6], 'T_range', [0 12], 'fs', 200, 'lya_T_interval', [4 12], 'lya_warmup', 3, ...
    'verbose', 'near-none', 'sigma_u_noise', 0, 'ode_solver', 'rk4', 'level_of_chaos', 2, ...
    'synapse_config', sc4};
mt = build_from_preset(P, full_cond, chaos{:}, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', false, 'lya_dt', 0.05);
mt.run();
mb = build_from_preset(P, full_cond, chaos{:}, 'lya_method', 'benettin');
mb.run();
chaotic = struct('topk', mt.lya_results, 'topk_summary', mt.lya_summary(), 'benettin', mb.lya_results);

extra = results_only(extra);
estimators = results_only(estimators);
estimators_spread = results_only(estimators_spread);
chaotic = results_only(chaotic);
save(fullfile(out_dir, 'rapid_chain_extra.mat'), 'extra', 'estimators', 'estimators_spread', 'chaotic', '-v7');

%% 2. A sweep with chain conditions, set up as the sweep stages set theirs up
ctx = resolve_run_context('param_space', 'preset_name', P, 'run_mode', 'fast', 'verbose', 'near-none');
c_one  = ctx.conditions{2};
c_full = ctx.conditions{3};
sc_EE = c_one.synapse_config;
sc_EE.E.E.std = chain;
sweep_conditions = { ...
    struct('name', 'chain_all', 'tau_a', {c_full.tau_a}, 'synapse_config', sc4), ...
    struct('name', 'chain_EE',  'tau_a', {c_one.tau_a},  'synapse_config', sc_EE)};
d = ctx.model_defaults;
d.n = 40; d.indegree = 8; d.F_tracks_network = true;
d.T_range = [0 4]; d.lya_T_interval = [2 4]; d.lya_warmup = 1;
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 4, 'note', 'regression_rapid_chain', ...
    'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = d;
psa.output_dir     = fullfile(out_dir, 'sweep_scratch');
psa.set_conditions(sweep_conditions);
psa.add_grid_parameter('level_of_chaos', [0.5, 1.5]);
psa.add_grid_parameter('reps', 1:2);
psa.run();
sweep = struct();
for i = 1:numel(sweep_conditions)
    sweep.(sweep_conditions{i}.name) = psa.results.(sweep_conditions{i}.name);
end
% Only the results are kept: the sweep's own files hold wall-clock times and
% the object, and are not what is being compared.
rmdir(psa.output_dir, 's');
sweep = results_only(sweep);

% The same sweep with the spread as a model default, as a preset carries it (1e).
d.std_spread = spread_chain{2};
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 4, 'note', 'regression_rapid_chain_spread', ...
    'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = d;
psa.output_dir     = fullfile(out_dir, 'sweep_scratch');
psa.set_conditions(sweep_conditions);
psa.add_grid_parameter('level_of_chaos', [0.5, 1.5]);
psa.add_grid_parameter('reps', 1:2);
psa.run();
sweep_spread = struct();
for i = 1:numel(sweep_conditions)
    sweep_spread.(sweep_conditions{i}.name) = psa.results.(sweep_conditions{i}.name);
end
rmdir(psa.output_dir, 's');
sweep_spread = results_only(sweep_spread);
save(fullfile(out_dir, 'rapid_chain_sweep.mat'), 'sweep', 'sweep_spread', '-v7');

%% 3. Memory capacity through the reservoir class, four identical chains
[preset_mc, ~, conditions_mc] = srnn_param_preset(P_mc);
fs = 200;
c = conditions_mc{end};
c.synapse_config = sc4;
base_args = struct2namevalue(preset_mc);
cargs = struct2namevalue(rmfield(c, 'name'));
esn = SRNN_ESN_reservoir(base_args{:}, 'fs', fs, 'ode_solver', 'sra1', 'rng_seeds', [3 4], ...
    'input_type', 'sample_hold', 'T_hold', 0.3, ...
    'T_wash', 5 * fs, 'T_train', 30 * fs, 'T_test', 15 * fs, 'd_max', 3 * fs, cargs{:});
esn.build();
[MC, R2, res] = esn.run_memory_capacity('readout_signal', 'synaptic', ...
    'verbose', false, 'store_timeseries', false);
mc = struct();
mc.chain = struct('MC', MC, 'R2', R2, 'delay_s', res.delay_s, ...
    'hold_len', res.hold_len, 'W', esn.W, 'W_in', esn.W_in);
% ... and with the spread on: four routes that stay identical neuron by
% neuron, so the synaptic readout is still defined (1e).
esn = SRNN_ESN_reservoir(base_args{:}, 'fs', fs, 'ode_solver', 'sra1', 'rng_seeds', [3 4], ...
    'input_type', 'sample_hold', 'T_hold', 0.3, ...
    'T_wash', 5 * fs, 'T_train', 30 * fs, 'T_test', 15 * fs, 'd_max', 3 * fs, cargs{:}, spread_chain{:});
esn.build();
[MC, R2, res] = esn.run_memory_capacity('readout_signal', 'synaptic', ...
    'verbose', false, 'store_timeseries', false);
mc.chain_spread = struct('MC', MC, 'R2', R2, 'delay_s', res.delay_s, ...
    'hold_len', res.hold_len, 'W', esn.W, 'W_in', esn.W_in, 'std_draws', esn.std_draws, ...
    'std_clipped', esn.std_clipped);
save(fullfile(out_dir, 'rapid_chain_mc.mat'), 'mc', '-v7');

seconds = toc(t_all);
fprintf('regression_rapid_chain_run: %.1f s -> %s\n', seconds, out_dir);
end

%% ------------------------------------------------------------------------
function m = mixed_model(solver, run_args)
% Four cell types; from the presynaptic type E a chain, a two-component
% product, a one-component and an undepressed route. The same network as
% regression_tol_run's mixed network at gain 1 (written out again here so
% the two references do not depend on one another's file).
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
if strcmp(solver, 'sra1'); noise = 0.025; else; noise = 0; end
m = SRNNCellTypePairs('n', 40, 'indegree', 10, 'n_cellTypes', 4, ...
    'cell_type_names', {'E', 'PV', 'SST', 'VIP'}, ...
    'f', [0.4 0.2 0.2 0.2], 'mu_tilde_relative', [3 -4 -2 -1.5], 'sigma_tilde_relative', [1 1 1 1], ...
    'tau_a', {[0.25 1], [], 2, 1}, 'c', [0.3 0 0.2 0.1], 'synapse_config', sc, ...
    'level_of_chaos', 1, 'input_config', struct('n_steps', 6, 'step_density', struct(), 'amp', 0.5, ...
        'no_stim_pattern', false(1, 6), 'intrinsic_drive', 0, 'positive_only', false), ...
    run_args{:}, 'sigma_u_noise', noise, 'ode_solver', solver);
m.build();
end

function realised = realised_parameters(m, p)
% Per route, every per-neuron array of the route model's parameters (the
% value each neuron has; `clipped`, whom the chain's ordering rule moved),
% under the model's own field names. Routes without state have none.
names = m.cell_type_names;
realised = struct();
for pre = 1:numel(names)
    for post = 1:numel(names)
        P = p.route{pre, post}.P;
        fields = fieldnames(P);
        route = struct('jittered', P.jittered);
        for k = 1:numel(fields)
            v = P.(fields{k});
            if (isnumeric(v) || islogical(v)) && size(v, 1) == P.npre && ~isempty(v) && ...
                    ~any(strcmp(fields{k}, {'npre', 'ns', 'n_b', 'n_g', 'jittered'}))
                route.(fields{k}) = v;
            end
        end
        realised.([names{pre} '_to_' names{post}]) = route;
    end
end
end

function v = results_only(v)
% Remove, at any depth of nested structs and cells, the wall-clock fields
% (they differ between identical runs) and every field named params (the
% model's internal packing; see the header).
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

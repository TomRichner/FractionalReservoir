function [R, seconds] = regression_tol_run(ref_data, suite)
% REGRESSION_TOL_RUN Build every case of the tolerance regression tier and return the raw arrays.
%
%   R = regression_tol_run()            % CAPTURE: choose the pointwise states now
%   R = regression_tol_run(ref_data)    % CHECK: evaluate at the states stored in ref_data
%   R = regression_tol_run([], 'chain') % the CHAIN suite (below), capture
%   R = regression_tol_run(ref_data, 'chain')
%   [R, seconds] = regression_tol_run(...)
%
% TWO SUITES, each with a reference of its own. SUITE 'base' (the default)
% is everything described from here to "THE CHAIN SUITE": the networks the
% model class had before Plan A touched it, whose reference is the anchor.
% SUITE 'chain' is the same two layers on networks with chained
% vesicle-pool depression (route_pool); it is described at the end of this
% header and its reference is regression_tol_chain.mat. The anchor holds no
% chain route and is never re-captured, which is why the chain's cases are
% a suite beside it and not more cases in it.
%
% THE TOLERANCE TIER (FD-020; docs/Plan_A_chained_vessicles.md, section 2.2)
% is the gate for a rewrite of the model's ARITHMETIC. The three byte-exact
% tiers (regression_rapid_run and friends) cannot gate that: a reordered
% floating-point expression changes the last bit, and a tolerance means
% nothing on a chaotic trajectory. This function builds the first two of
% the tier's three layers and returns them as plain arrays in one struct;
% regression_tol_compare compares two such structs element by element.
%
%   R.pointwise   single evaluations at STORED states: no trajectory is
%                 compared, so the regime (chaotic or not) does not matter.
%   R.stable      small networks that are STABLE with a clear margin, end to
%                 end: trajectories, Lyapunov estimates, transient gain, a
%                 sweep, memory capacity.
%
% (The third layer, the chaotic and bundle digests within a loose band, is
% regression_tol_band_check; it takes minutes and is run on demand.)
%
% LAYER 1, POINTWISE. Per case: the canonical parameter view (below), the
% states S (N x ns) and their times t, and at each state the right-hand
% side (dynamics_fast), the assembled Jacobian (compute_Jacobian_fast,
% sparse), the matrix-free product with a fixed N x 4 matrix
% (jacobian_times), and the rates and per-route synaptic outputs
% (unpack_and_compute_states). The cases:
%   pairs_noStim.<cond>.<integrator>   the regression network (n = 100),
%        chaotic, three conditions, noise off (rk4) and on (sra1)
%   pairs_steps.<cond>.<integrator>    the same under the step staircase,
%        which drives rates to both ends and depression deep
%   options_chaotic                    the three-cell-type network of
%        regression_rapid_run: facilitation, dual facilitation, the
%        geometric-mean combine, the zero floor, a route scale
%   stable_steps.<cond>.<integrator>   the stable networks of layer 2
%   options_stable                     the stable options network of layer 2
% In a CAPTURE the chaotic networks are run and the states picked from the
% trajectory: four evenly spaced, plus the time of the highest mean rate
% and the time of the deepest depression. In a CHECK nothing chaotic is
% run: the network is built and evaluated AT THE REFERENCE'S STATES, so a
% trajectory that diverged after a harmless reordering cannot fail this
% layer. The state vector's ordering is fixed by Plan A section 4.0, so
% states, right-hand sides and Jacobians compare element for element.
%
% THE CANONICAL PARAMETER VIEW. The params struct of get_params is NOT
% stored: its fields change in the refactor this tier gates (FD-021). What
% is stored is the physics, read through accessors that keep their
% signatures: the scaled weight matrix the dynamics read (params.W),
% tau_d, c_eff, tau_a_matrix, the setpoints, the nonlinearity's name and
% shape parameters, the zero-floor flag, and per route the struct returned
% by SRNNCellTypePairs.route_fields. A later model may add fields to that
% struct; the comparison reports leaves that exist only in the new run and
% does not fail on them.
%
% LAYER 2, STABLE END TO END. Presets regression_stable_steps_* (the
% regression network at n = 40 with the gain lowered to 0.2, under the step
% staircase; see their case blocks in srnn_param_preset for the margins):
%   stable.model.<cond>.<integrator>   S0, every 5th row of S_out and the
%        last, t_out, plot_data, the top-K Lyapunov results and lya_summary
%   stable.extra        on the most-adapted condition with noise on: the
%        transient gain in its three variants, the finite-difference Jacobian
%   stable.estimators   Benettin (noise on) and QR (noise off) on a reduced
%        network, first and last condition
%   stable.options      the options network at gain 0.2 under a staircase
%   stable.sweep        a 12-job serial ParamSpaceAnalysis2 sweep (two
%        gains x two reps x three conditions), set up as the stages do
%   stable.mc           SRNN_ESN_reservoir memory capacity per condition
% Wall-clock fields and the params struct that lya_results carries are
% removed before anything is returned.
%
% HOW STABLE "STABLE" IS (measured 2026-11-11 with the settings below;
% exponents in 1/s over [1.5, 5] s; noise off / noise on):
%                  top-K lambda_1   its largest     Benettin         QR, largest
%                                   LOCAL rate                       (noise off)
%   no_adaptation  -4.35 / -3.99    -2.98 / -2.78   -4.25 / -4.43    -4.22
%   sfa1_std1      -0.83 / -0.85    -0.69 / -0.68   -0.82 / -0.82    -0.56
%   sfa3_std2      -0.154 / -0.152  -0.099 / -0.095 -0.170 / -0.168  -0.070
%   options net    -0.77            -0.65           -0.76            -0.53
% No estimator finds a positive exponent anywhere (QR: none of 40, 200,
% 320 and 171), the leading direction's local growth rate never reaches
% zero, lya_K_auto never retries (K stays 10; 15 in the sweep), and no SFA
% ladder is clipped. The 12 sweep jobs have lambda_1 between -6.2 and
% -0.22. The adapting conditions sit as close to zero as they do because
% their leading exponent is the decay rate of the slowest adaptation
% variable, about -1/tau with tau up to 10 s: lowering the gain further
% does not move it, and the three estimators differ there because that
% slow band is nearly degenerate (docs/notes/Lyapunov_estimation_methods.md).
% Against round-off the margin is ten orders of magnitude. Meanwhile the
% staircase keeps the arithmetic busy: rates span 0 to 1 (a quarter of
% the no-adaptation samples at exactly 0, some at exactly 1), depression
% reaches b = 0.13, facilitation g = 1.65 in the options network.
% One consequence of stability: the transient gain never exceeds 1, so
% its peak is at t = 0 and only the decay curves G.worst(t), G.noise(t)
% carry information.
%
% EVERYTHING RUNS IN THE CLIENT. The sweep is serial (use_parallel false)
% and nothing here opens or uses a pool, so the result does not depend on
% whether one is up. It does depend, at round-off, on the client's thread
% count (see regression_tol_settings); R.info.threads records it.
%
% A CASE THAT THROWS. In a capture the error is rethrown: a reference with
% a hole in it is no reference. In a check it is recorded and the run goes
% on, so one case that no longer runs does not hide what the others have
% to say: the case's branch holds a single field ERROR with the message,
% the reference's leaves under it then count as missing (which fails), and
% R.info.errors lists the messages. (How this was found: a tau_rec off by
% 1e-6 on one route makes SRNN_ESN_reservoir refuse the synaptic readout,
% because the routes of one presynaptic type are no longer identical; the
% first version of this function died there and reported nothing.)
%
% IT USES ONLY THE PUBLIC SURFACE THAT THE REFACTOR KEEPS: build, run,
% get_params, the statics dynamics_fast, compute_Jacobian_fast,
% jacobian_times, finite_difference_jacobian, transient_gain, route_fields,
% unpack_and_compute_states, the result properties, ParamSpaceAnalysis2
% and SRNN_ESN_reservoir as regression_rapid_run uses them.
%
% DO NOT RETUNE the settings below. If they must change, the anchor
% changes with them, and the anchor describes the code BEFORE the model
% class was touched: it cannot be re-captured once that code is gone.
%
% THE CHAIN SUITE (Plan A, step 5; 2026-11-12). Chained depression is
% std.model = 'chain' on a route: a releasable pool b1 and a reserve b2 per
% presynaptic neuron, a 2 x 2 coupled linearisation, output factor b1
% (route_pool). The suite's networks are built here from the regression
% presets' defaults plus a chain synapse_config, as the options network is;
% no preset is added for them. Two kinds of network:
%   chain   the regression network in its sfa3_std2 condition with all four
%           routes switched to the chain (tau_rec [2 5], tau_rel 0.25,
%           n2 2): SFA on both types, per-neuron setpoints and ladders
%   mixed   FOUR cell types (E, PV, SST, VIP; n = 40) so that one
%           presynaptic type, E, has every kind of route side by side: a
%           chain (E->E), a two-component product (E->PV), one component
%           (E->SST) and no depression (E->VIP). Also a chain with a route
%           scale (PV->E), a geometric-mean product (SST->E), facilitation
%           on a product route (SST->PV), a chain with a small reserve
%           (VIP->SST, n2 0.5) and one with a huge reserve (VIP->VIP,
%           n2 1e6). The full-A and the diagonal branch of jacobian_times
%           both run inside every one of its Jacobian products.
% LAYER 1, POINTWISE (as above: view, S, t, dS, J, JY, rates, outputs):
%   chain_noStim.<integrator>   the chain network at n = 100, gain 1, no
%        input, as it comes (noise off: weakly chaotic, lambda_1 about
%        +0.1 where the product network is stable)
%   chain_steps.<integrator>    the same under the step staircase
%   mixed_gain1                 the mixed network at gain 1 under a staircase
%   chain_stable.<integrator>, mixed_stable.<integrator>   layer 2's networks
% The states are those of pick_states plus, for the chain routes: the
% first sample (every pool exactly 1), the time of the smallest b1, the
% time of the smallest b2, and TWO STATES THAT NO TRAJECTORY HERE REACHES,
% made from the mid-run state by overwriting the chain routes' rows: both
% pools deeply depleted (b1 in [0.01, 0.05], b2 in [0.02, 0.10]) and a
% nearly full reserve over a part-used pool (b2 in [1 - 1e-4, 1], b1 in
% [0.15, 0.85]). A single evaluation needs no trajectory to get there.
% Beyond the base suite's quantities a case holds std_factor, every
% route's depression factor at the states as display code gets it
% (std_factors, which calls std_route_factor: b1 on a chain route, not
% prod(b)), and its view holds params.std_chain and params.std_n2, which
% display code reads and route_fields does not. Layer 2's model runs hold
% std_factor along the run. No arithmetic reads those three, so without
% these leaves an error in them would pass.
% LAYER 2, STABLE END TO END (gain 0.2, staircase, as the base suite):
%   stable.model.chain.<integrator>, stable.model.mixed.<integrator>
%   stable.extra.chain, stable.extra.mixed   transient gain in its three
%        variants and the finite-difference Jacobian, noise on
%   stable.estimators.chain, .mixed   Benettin (noise on) and QR (noise
%        off) on reduced networks (n = 20)
%   stable.sweep        an 8-job serial ParamSpaceAnalysis2 sweep: two
%        gains x two reps x two conditions, chain_all (the chain network)
%        and chain_EE (one SFA timescale, a chain on E->E beside one-
%        component routes). The conditions are local to this file.
%   stable.mc.chain     memory capacity, synaptic readout, four identical
%        chain routes
% HOW STABLE: see the table in chain_suite below.
%
% See also: regression_tol_compare, regression_tol_check,
%           make_regression_tol_reference, regression_tol_settings,
%           regression_tol_band_check, regression_rapid_run
arguments
    ref_data = []
    suite (1,:) char {mustBeMember(suite, {'base', 'chain'})} = 'base'
end
if strcmp(suite, 'chain')
    [R, seconds] = chain_suite(ref_data);
    return;
end
t_all = tic;
capture = isempty(ref_data);
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));

P_chaos  = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
P_steps  = 'regression_pairs_steps5s_noise0p025_dualStd_3cond';
P_stable = struct('rk4_noise_off', 'regression_stable_steps_noise0_dualStd_3cond', ...
                  'sra1_noise_on', 'regression_stable_steps_noise0p025_dualStd_3cond');
[~, ~, conditions] = srnn_param_preset(P_chaos);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
integrators = struct('name', {'rk4_noise_off', 'sra1_noise_on'}, ...
    'args', {{'sigma_u_noise', 0, 'ode_solver', 'rk4'}, {'ode_solver', 'sra1'}});

R = struct('pointwise', struct(), 'stable', struct(), 'info', struct());
R.info.threads = maxNumCompThreads;
errors = {};

%% Layer 1a. Chaotic networks: states only, then single evaluations
chaos_args = {'rng_seeds', [5 6], 'T_range', [0 4], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'none', 'verbose', 'near-none'};
sets = struct('name', {'pairs_noStim', 'pairs_steps'}, 'preset', {P_chaos, P_steps});
for q = 1:numel(sets)
    for i = 1:numel(cond_names)
        for s = 1:numel(integrators)
            key = {'pointwise', sets(q).name, cond_names{i}, integrators(s).name};
            try
                m = build_from_preset(sets(q).preset, cond_names{i}, chaos_args{:}, integrators(s).args{:});
                if capture
                    m.run();
                    [S, t] = pick_states(m);
                else
                    [S, t] = stored_states(ref_data, key(2:end));
                end
                R = setfield(R, key{:}, pointwise_case(m, S, t));
            catch err
                [R, errors] = failed(R, errors, capture, err, key);
            end
        end
    end
end

key = {'pointwise', 'options_chaotic'};
try
    mo = options_model(1, false);
    if capture
        mo.run();
        [S, t] = pick_states(mo);
    else
        [S, t] = stored_states(ref_data, key(2:end));
    end
    R = setfield(R, key{:}, pointwise_case(mo, S, t));
catch err
    [R, errors] = failed(R, errors, capture, err, key);
end

%% Layer 2a. The stable networks: three conditions x two integrators
stable_args = {'rng_seeds', [5 6], 'T_range', [0 5], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', true, 'lya_dt', 0.05, ...
    'lya_T_interval', [1.5 5], 'lya_warmup', 1.5, 'verbose', 'near-none'};
for i = 1:numel(cond_names)
    for s = 1:numel(integrators)
        key_run = {'stable', 'model', cond_names{i}, integrators(s).name};
        key_pw  = {'pointwise', 'stable_steps', cond_names{i}, integrators(s).name};
        keys = {key_run, key_pw};
        with_extra = i == numel(cond_names) && s == numel(integrators);
        if with_extra; keys{end + 1} = {'stable', 'extra'}; end %#ok<AGROW>
        try
            % The preset states the noise; only the integrator is named here.
            m = build_from_preset(P_stable.(integrators(s).name), cond_names{i}, stable_args{:}, ...
                'ode_solver', integrators(s).args{end});
            m.run();
            R = setfield(R, key_run{:}, end_to_end(m));
            if capture
                [S, t] = pick_states(m);
            else
                [S, t] = stored_states(ref_data, key_pw(2:end));
            end
            R = setfield(R, key_pw{:}, pointwise_case(m, S, t));

            % 2b. On the most-adapted condition with noise on: the transient
            % gain in its three variants, and the finite-difference Jacobian.
            if with_extra
                p = m.get_params();
                p.u_interpolant = m.u_interpolant;
                extra = struct();
                i0 = find(m.t_out >= 2.5, 1);
                for variant = {'frozen_x', 'frozen_full', 'active'}
                    [G, g_info] = SRNNCellTypePairs.transient_gain(m.S_out, m.t_out, i0, p, 0.5, ...
                        struct('variant', variant{1}));
                    extra.transient_gain.(variant{1}) = struct('G', G, 'info', drop_volatile(g_info));
                end
                extra.finite_difference_J = SRNNCellTypePairs.finite_difference_jacobian( ...
                    m.S_out(end, :)', p, 1e-6);
                R.stable.extra = extra;
            end
        catch err
            [R, errors] = failed(R, errors, capture, err, keys{:});
        end
    end
end

%% Layer 2c. The other two Lyapunov estimators, on a reduced stable network
% QR carries N^2 equations, so n is halved again; first and last condition.
small = {'n', 20, 'indegree', 4, 'F_ref_n', 20, 'F_ref_indegree', 4, 'rng_seeds', [5 6], ...
    'T_range', [0 3], 'fs', 200, 'lya_T_interval', [1.5 3], 'lya_warmup', 1, 'verbose', 'near-none'};
for i = [1, numel(cond_names)]
    key = {'stable', 'estimators', cond_names{i}};
    try
        mb = build_from_preset(P_stable.sra1_noise_on, cond_names{i}, small{:}, ...
            'ode_solver', 'sra1', 'lya_method', 'benettin');
        mb.run();
        mq = build_from_preset(P_stable.rk4_noise_off, cond_names{i}, small{:}, ...
            'ode_solver', 'rk4', 'lya_method', 'qr');
        mq.run();
        R = setfield(R, key{:}, struct('benettin', drop_volatile(mb.lya_results), ...
            'qr', drop_volatile(mq.lya_results)));
    catch err
        [R, errors] = failed(R, errors, capture, err, key);
    end
end

%% Layer 2d. The options network, stable
key_run = {'stable', 'options'};
key_pw  = {'pointwise', 'options_stable'};
try
    mo = options_model(0.2, true);
    mo.run();
    R = setfield(R, key_run{:}, end_to_end(mo));
    if capture
        [S, t] = pick_states(mo);
    else
        [S, t] = stored_states(ref_data, key_pw(2:end));
    end
    R = setfield(R, key_pw{:}, pointwise_case(mo, S, t));
catch err
    [R, errors] = failed(R, errors, capture, err, key_run, key_pw);
end

%% Layer 2e. A sweep, set up as the sweep stages set theirs up
scratch = tempname;
try
    ctx = resolve_run_context('param_space', 'preset_name', P_stable.sra1_noise_on, ...
        'run_mode', 'fast', 'verbose', 'near-none');
    d = ctx.model_defaults;
    d.T_range = [0 4]; d.lya_T_interval = [2 4]; d.lya_warmup = 1;
    psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 6, 'note', 'regression_tol', ...
        'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
    psa.model_class    = ctx.model_class;
    psa.integer_params = ctx.integer_params;
    psa.model_defaults = d;
    psa.output_dir     = scratch;
    psa.set_conditions(ctx.conditions);
    psa.add_grid_parameter('level_of_chaos', [0.1, 0.2]);
    psa.add_grid_parameter('reps', 1:2);
    psa.run();
    sweep = struct();
    for i = 1:numel(cond_names)
        rc = psa.results.(cond_names{i});
        sweep.(cond_names{i}) = cellfun(@drop_volatile, rc, 'UniformOutput', false);
    end
    R.stable.sweep = sweep;
catch err
    [R, errors] = failed(R, errors, capture, err, {'stable', 'sweep'});
end
% Only the results are kept: the sweep's own files hold wall-clock times.
if isfolder(scratch); rmdir(scratch, 's'); end

%% Layer 2f. Memory capacity through the reservoir class
[preset_mc, ~, conditions_mc] = srnn_param_preset(P_stable.rk4_noise_off);
fs = 200;
base_args = [struct2namevalue(preset_mc), ...
    {'fs', fs, 'ode_solver', 'sra1', 'rng_seeds', [3 4], ...
     'input_type', 'sample_hold', 'T_hold', 0.3, ...
     'T_wash', 5 * fs, 'T_train', 30 * fs, 'T_test', 15 * fs, 'd_max', 3 * fs}];
for i = 1:numel(conditions_mc)
    c = conditions_mc{i};
    key = {'stable', 'mc', c.name};
    try
        cargs = struct2namevalue(rmfield(c, 'name'));
        esn = SRNN_ESN_reservoir(base_args{:}, cargs{:});
        esn.build();
        [MC, R2, res] = esn.run_memory_capacity('readout_signal', 'synaptic', ...
            'verbose', false, 'store_timeseries', false);
        R = setfield(R, key{:}, struct('MC', MC, 'R2', R2, 'delay_s', res.delay_s, ...
            'hold_len', res.hold_len, 'W', esn.W, 'W_in', esn.W_in));
    catch err
        [R, errors] = failed(R, errors, capture, err, key);
    end
end

% Only a check can have recorded errors (a capture rethrows), and a
% capture's output stays exactly what the anchor holds.
if ~capture
    R.info.errors = errors;
end
seconds = toc(t_all);
end

%% ------------------------------------------------------------------------
function [R, seconds] = chain_suite(ref_data)
% The chain suite; what it holds is in the file's header. Same structure,
% same helpers and same rules as the base suite above (a capture rethrows,
% a check records a case that threw and goes on).
%
% HOW STABLE the layer-2 networks are (measured 2026-11-12 with the
% settings below; exponents in 1/s; noise off / noise on):
%              top-K lambda_1     its largest       Benettin   QR, largest
%                                 LOCAL rate        (n = 20)   (n = 20)
%   chain      -0.167 / -0.166    -0.125 / -0.122   -0.218     -0.075
%   mixed      -0.326 / -0.325    -0.282 / -0.281   -0.430     -0.250
% No estimator finds a positive exponent (QR: none of 160 and 124), the
% leading direction's local rate never reaches zero, lya_K_auto never
% retries (K stays 10) and no SFA ladder is clipped. The 8 sweep jobs have
% lambda_1 between -0.42 and -0.21. As in the base suite the leading
% exponent is the decay rate of a slow variable (an SFA rung, or a
% reserve), so the gain does not move it, and the three estimators differ
% there for the reason given above. Under the staircase the pools are
% worked: in the stable chain network rates span 0 to 1, b1 falls to 0.10
% and b2 to 0.58; in the mixed one the small reserve of VIP->SST falls to
% 0.03 (its b1 to 0.07) while the huge one of VIP->VIP stays within 1e-6
% of 1. The pointwise layer's two made states go further.
%
% DO NOT RETUNE the settings below; if they must change, the chain
% reference changes with them and the commit says why.
t_all = tic;
capture = isempty(ref_data);
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));

P_chaos  = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
P_steps  = 'regression_pairs_steps5s_noise0p025_dualStd_3cond';
P_stable = struct('rk4_noise_off', 'regression_stable_steps_noise0_dualStd_3cond', ...
                  'sra1_noise_on', 'regression_stable_steps_noise0p025_dualStd_3cond');
full_cond = 'sfa3_std2';                 % three SFA timescales; its routes are replaced
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
sc4 = struct();
sc4.E.E.std = chain; sc4.E.I.std = chain; sc4.I.E.std = chain; sc4.I.I.std = chain;
integrators = struct('name', {'rk4_noise_off', 'sra1_noise_on'}, ...
    'args', {{'sigma_u_noise', 0, 'ode_solver', 'rk4'}, {'ode_solver', 'sra1'}});

R = struct('pointwise', struct(), 'stable', struct(), 'info', struct());
R.info.threads = maxNumCompThreads;
errors = {};

%% Layer 1a. The chain network at full gain: states only, then single evaluations
chaos_args = {'rng_seeds', [5 6], 'T_range', [0 4], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'none', 'verbose', 'near-none'};
sets = struct('name', {'chain_noStim', 'chain_steps'}, 'preset', {P_chaos, P_steps});
for q = 1:numel(sets)
    for s = 1:numel(integrators)
        key = {'pointwise', sets(q).name, integrators(s).name};
        try
            m = build_from_preset(sets(q).preset, full_cond, chaos_args{:}, integrators(s).args{:}, ...
                'synapse_config', sc4);
            if capture
                m.run();
                [S, t] = pick_chain_states(m);
            else
                [S, t] = stored_states(ref_data, key(2:end));
            end
            R = setfield(R, key{:}, chain_pointwise(m, S, t));
        catch err
            [R, errors] = failed(R, errors, capture, err, key);
        end
    end
end

key = {'pointwise', 'mixed_gain1'};
try
    mm = mixed_model(1, 'rk4', {});
    if capture
        mm.run();
        [S, t] = pick_chain_states(mm);
    else
        [S, t] = stored_states(ref_data, key(2:end));
    end
    R = setfield(R, key{:}, chain_pointwise(mm, S, t));
catch err
    [R, errors] = failed(R, errors, capture, err, key);
end

%% Layer 2a. The stable networks, chain and mixed, noise off and on
stable_args = {'rng_seeds', [5 6], 'T_range', [0 5], 'fs', 200, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', true, 'lya_dt', 0.05, ...
    'lya_T_interval', [1.5 5], 'lya_warmup', 1.5, 'verbose', 'near-none'};
nets = {'chain', 'mixed'};
for q = 1:numel(nets)
    for s = 1:numel(integrators)
        key_run = {'stable', 'model', nets{q}, integrators(s).name};
        key_pw  = {'pointwise', [nets{q} '_stable'], integrators(s).name};
        keys = {key_run, key_pw};
        with_extra = s == numel(integrators);
        if with_extra; keys{end + 1} = {'stable', 'extra', nets{q}}; end %#ok<AGROW>
        try
            solver = integrators(s).args{end};
            if q == 1
                % The preset states the noise; only the integrator is named here.
                m = build_from_preset(P_stable.(integrators(s).name), full_cond, stable_args{:}, ...
                    'ode_solver', solver, 'synapse_config', sc4);
            else
                m = mixed_model(0.2, solver, {});
            end
            m.run();
            e2e = end_to_end(m);
            % Every route's depression factor along the run, as display
            % code gets it (std_factors, which calls std_route_factor): b1
            % on a chain route, not prod(b).
            e2e.std_factor = SRNNCellTypePairs.std_factors(m.plot_data.b, m.get_params());
            R = setfield(R, key_run{:}, e2e);
            if capture
                [S, t] = pick_chain_states(m);
            else
                [S, t] = stored_states(ref_data, key_pw(2:end));
            end
            R = setfield(R, key_pw{:}, chain_pointwise(m, S, t));

            % 2b. With noise on: the transient gain in its three variants,
            % and the finite-difference Jacobian.
            if with_extra
                p = m.get_params();
                p.u_interpolant = m.u_interpolant;
                extra = struct();
                i0 = find(m.t_out >= 2.5, 1);
                for variant = {'frozen_x', 'frozen_full', 'active'}
                    [G, g_info] = SRNNCellTypePairs.transient_gain(m.S_out, m.t_out, i0, p, 0.5, ...
                        struct('variant', variant{1}));
                    extra.transient_gain.(variant{1}) = struct('G', G, 'info', drop_volatile(g_info));
                end
                extra.finite_difference_J = SRNNCellTypePairs.finite_difference_jacobian( ...
                    m.S_out(end, :)', p, 1e-6);
                R.stable.extra.(nets{q}) = extra;
            end
        catch err
            [R, errors] = failed(R, errors, capture, err, keys{:});
        end
    end
end

%% Layer 2c. The other two Lyapunov estimators, on reduced stable networks
small = {'n', 20, 'indegree', 4, 'F_ref_n', 20, 'F_ref_indegree', 4};
small_run = {'rng_seeds', [5 6], 'T_range', [0 3], 'fs', 200, 'lya_T_interval', [1.5 3], ...
    'lya_warmup', 1, 'verbose', 'near-none'};
for q = 1:numel(nets)
    key = {'stable', 'estimators', nets{q}};
    try
        if q == 1
            mb = build_from_preset(P_stable.sra1_noise_on, full_cond, small{:}, small_run{:}, ...
                'ode_solver', 'sra1', 'lya_method', 'benettin', 'synapse_config', sc4);
            mq = build_from_preset(P_stable.rk4_noise_off, full_cond, small{:}, small_run{:}, ...
                'ode_solver', 'rk4', 'lya_method', 'qr', 'synapse_config', sc4);
        else
            mb = mixed_model(0.2, 'sra1', [{'n', 20, 'indegree', 5}, small_run, {'lya_method', 'benettin'}]);
            mq = mixed_model(0.2, 'rk4',  [{'n', 20, 'indegree', 5}, small_run, {'lya_method', 'qr'}]);
        end
        mb.run();
        mq.run();
        R = setfield(R, key{:}, struct('benettin', drop_volatile(mb.lya_results), ...
            'qr', drop_volatile(mq.lya_results)));
    catch err
        [R, errors] = failed(R, errors, capture, err, key);
    end
end

%% Layer 2e. A sweep with chain conditions, set up as the sweep stages set theirs up
scratch = tempname;
try
    ctx = resolve_run_context('param_space', 'preset_name', P_stable.sra1_noise_on, ...
        'run_mode', 'fast', 'verbose', 'near-none');
    % Two conditions made here (no preset carries a chain yet): the chain
    % network, and one SFA timescale with a chain on E->E beside the
    % preset's one-component routes.
    c_one  = ctx.conditions{2};
    c_full = ctx.conditions{3};
    sc_EE = c_one.synapse_config;
    sc_EE.E.E.std = chain;
    sweep_conditions = { ...
        struct('name', 'chain_all', 'tau_a', {c_full.tau_a}, 'synapse_config', sc4), ...
        struct('name', 'chain_EE',  'tau_a', {c_one.tau_a},  'synapse_config', sc_EE)};
    d = ctx.model_defaults;
    d.T_range = [0 4]; d.lya_T_interval = [2 4]; d.lya_warmup = 1;
    psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 4, 'note', 'regression_tol_chain', ...
        'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
    psa.model_class    = ctx.model_class;
    psa.integer_params = ctx.integer_params;
    psa.model_defaults = d;
    psa.output_dir     = scratch;
    psa.set_conditions(sweep_conditions);
    psa.add_grid_parameter('level_of_chaos', [0.1, 0.2]);
    psa.add_grid_parameter('reps', 1:2);
    psa.run();
    sweep = struct();
    for i = 1:numel(sweep_conditions)
        rc = psa.results.(sweep_conditions{i}.name);
        sweep.(sweep_conditions{i}.name) = cellfun(@drop_volatile, rc, 'UniformOutput', false);
    end
    R.stable.sweep = sweep;
catch err
    [R, errors] = failed(R, errors, capture, err, {'stable', 'sweep'});
end
if isfolder(scratch); rmdir(scratch, 's'); end

%% Layer 2f. Memory capacity through the reservoir class, four identical chains
key = {'stable', 'mc', 'chain'};
try
    [preset_mc, ~, conditions_mc] = srnn_param_preset(P_stable.rk4_noise_off);
    fs = 200;
    c = conditions_mc{end};
    c.synapse_config = sc4;
    cargs = struct2namevalue(rmfield(c, 'name'));
    base_args = struct2namevalue(preset_mc);
    esn = SRNN_ESN_reservoir(base_args{:}, ...
        'fs', fs, 'ode_solver', 'sra1', 'rng_seeds', [3 4], ...
        'input_type', 'sample_hold', 'T_hold', 0.3, ...
        'T_wash', 5 * fs, 'T_train', 30 * fs, 'T_test', 15 * fs, 'd_max', 3 * fs, cargs{:});
    esn.build();
    [MC, R2, res] = esn.run_memory_capacity('readout_signal', 'synaptic', ...
        'verbose', false, 'store_timeseries', false);
    R = setfield(R, key{:}, struct('MC', MC, 'R2', R2, 'delay_s', res.delay_s, ...
        'hold_len', res.hold_len, 'W', esn.W, 'W_in', esn.W_in));
catch err
    [R, errors] = failed(R, errors, capture, err, key);
end

if ~capture
    R.info.errors = errors;
end
seconds = toc(t_all);
end

function m = mixed_model(gain, solver, extra_args)
% The mixed network of the chain suite: four cell types, and from the one
% presynaptic type E a chain, a two-component product, a one-component and
% an undepressed route (see the header). Under a step staircase always, so
% the routes' variables move at any gain; noise on with 'sra1', off with
% 'rk4'. EXTRA_ARGS are applied last.
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
args = {'n', 40, 'indegree', 10, 'n_cellTypes', 4, 'cell_type_names', {'E', 'PV', 'SST', 'VIP'}, ...
    'f', [0.4 0.2 0.2 0.2], 'mu_tilde_relative', [3 -4 -2 -1.5], 'sigma_tilde_relative', [1 1 1 1], ...
    'tau_a', {[0.25 1], [], 2, 1}, 'c', [0.3 0 0.2 0.1], 'synapse_config', sc, ...
    'level_of_chaos', gain, 'sigma_u_noise', noise, 'ode_solver', solver, ...
    'input_config', struct('n_steps', 6, 'step_density', struct(), 'amp', 0.5, ...
        'no_stim_pattern', false(1, 6), 'intrinsic_drive', 0, 'positive_only', false), ...
    'rng_seeds', [5 6], 'T_range', [0 5], 'fs', 200, ...
    'lya_method', 'topk', 'lya_K', 10, 'lya_K_auto', true, 'lya_dt', 0.05, ...
    'lya_T_interval', [1.5 5], 'lya_warmup', 1.5, 'store_full_state', true, 'verbose', 'near-none'};
m = SRNNCellTypePairs(args{:}, extra_args{:});
m.build();
end

function c = chain_pointwise(m, S, t)
% pointwise_case, and what only a network with chain routes has: the two
% params fields that say which routes are chains and their n2 (display
% code reads them; route_fields, and so the view, does not), and every
% route's depression factor at the states as std_route_factor gives it.
c = pointwise_case(m, S, t);
p = m.get_params();
c.view.std_chain = p.std_chain;
c.view.std_n2 = p.std_n2;
[~, ~, b] = SRNNCellTypePairs.unpack_and_compute_states(S', p);
c.std_factor = SRNNCellTypePairs.std_factors(b, p);
end

function [S, t] = pick_chain_states(m)
% States for the chain suite's pointwise layer, as columns: those of
% pick_states (four evenly spaced, the highest mean rate, the smallest b on
% any route), the first sample (every pool exactly 1), the times of the
% smallest b1 and of the smallest b2 on any chain route, and two states
% made from the mid-run state by overwriting the chain routes' rows:
%   depleted      b1 in [0.01, 0.05], b2 in [0.02, 0.10]
%   full reserve  b1 in [0.15, 0.85], b2 in [1 - 1e-4, 1]
% with fixed patterns over the neurons, not random. The made states carry
% the mid-run time.
[S, t] = pick_states(m);
p = m.get_params();
lay = p.state_layout;
nt = size(m.S_out, 1);
b1_min = inf(nt, 1);
b2_min = inf(nt, 1);
i_mid = round(nt / 2);
depleted = m.S_out(i_mid, :)';
reserve_full = depleted;
for pre = 1:p.n_cellTypes
    npre = p.n_per_type(pre);
    w1 = 0.5 + 0.5 * cos((1:npre)');
    w2 = 0.5 + 0.5 * sin(2 * (1:npre)');
    for post = 1:p.n_cellTypes
        if ~p.std_chain(pre, post); continue; end
        rows = lay.b{pre, post};             % component-major: b1 then b2
        rows1 = rows(1:npre);
        rows2 = rows(npre + 1:end);
        b1_min = min(b1_min, min(m.S_out(:, rows1), [], 2));
        b2_min = min(b2_min, min(m.S_out(:, rows2), [], 2));
        depleted(rows1) = 0.01 + 0.04 * w1;
        depleted(rows2) = 0.02 + 0.08 * w2;
        reserve_full(rows1) = 0.15 + 0.70 * w1;
        reserve_full(rows2) = 1 - 1e-4 * w2;
    end
end
[~, i1] = min(b1_min);
[~, i2] = min(b2_min);
idx = unique([1, i1, i2]);
S = [S, m.S_out(idx, :)', depleted, reserve_full];
t = [t, m.t_out(idx)', m.t_out(i_mid), m.t_out(i_mid)];
end

%% ------------------------------------------------------------------------
function [R, errors] = failed(R, errors, capture, err, varargin)
% A case that threw. In a CAPTURE that is fatal: a reference with a hole in
% it is no reference. In a CHECK the error is recorded and the run goes
% on, so that one case that no longer runs (a route check that now
% refuses, say) does not hide what every other case has to say. Each key
% path in VARARGIN that the case had not yet filled gets a one-field
% struct holding the message; the reference's leaves under it then count
% as missing, which fails the comparison, and R.info.errors lists the
% messages for regression_tol_check to print.
if capture
    rethrow(err);
end
text = sprintf('%s (%s)', err.message, err.identifier);
for k = 1:numel(varargin)
    path = varargin{k};
    node = R;
    present = true;
    for j = 1:numel(path)
        if isstruct(node) && isfield(node, path{j})
            node = node.(path{j});
        else
            present = false;
            break;
        end
    end
    if ~present
        R = setfield(R, path{:}, struct('ERROR', text));
        errors{end + 1} = sprintf('%s: %s', strjoin(path, '.'), text); %#ok<AGROW>
    end
end
end

%% ------------------------------------------------------------------------
function m = options_model(gain, staircase)
% The three-cell-type network of regression_rapid_run, section 1d: every
% synapse option the paper's presets leave off. GAIN = 1 without a
% staircase is that network exactly (chaotic); gain 0.2 with a staircase is
% its stable variant, driven so that rates and the route variables move.
sc = struct();
sc.E.E.std   = struct('tau_rec', [0.3 1], 'tau_rel', 0.2, 'combine', 'geomean');
sc.E.PV.stf  = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
sc.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
sc.E.SST.stf = struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
sc.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
sc.PV.PV.scale = 2;
sc.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
args = {'n', 30, 'indegree', 10, 'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], 'mu_tilde_relative', [3 -4 -2], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {[0.25 1], [], 2}, 'c', [0.3 0 0.2], 'synapse_config', sc, ...
    'std_zero_floor', true, 'level_of_chaos', gain, 'rng_seeds', [5 6], 'T_range', [0 3], 'fs', 200, ...
    'ode_solver', 'rk4', 'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, ...
    'lya_T_interval', [1.5 3], 'lya_warmup', 1, 'store_full_state', true, 'verbose', 'near-none'};
if staircase
    args = [args, {'input_config', struct('n_steps', 6, 'step_density', struct(), 'amp', 0.5, ...
        'no_stim_pattern', false(1, 6), 'intrinsic_drive', 0, 'positive_only', false)}];
end
m = SRNNCellTypePairs(args{:});
m.build();
end

function [S, t] = pick_states(m)
% States for the pointwise layer, as columns: four evenly spaced along the
% trajectory, the time of the highest mean rate, and the time of the
% deepest depression (the smallest b on any route), when the model has one.
nt = size(m.S_out, 1);
idx = round(nt * [0.25 0.5 0.75 1.0]);
[~, ~, b, r] = SRNNCellTypePairs.unpack_and_compute_states(m.S_out, m.get_params());
rates = struct2cell(r);
[~, i_rate] = max(mean(vertcat(rates{:}), 1));
idx(end + 1) = i_rate;
b_min = inf(1, nt);
pres = fieldnames(b);
for a = 1:numel(pres)
    posts = fieldnames(b.(pres{a}));
    for c = 1:numel(posts)
        v = b.(pres{a}).(posts{c});          % npre x n_b x nt
        if ~isempty(v)
            b_min = min(b_min, reshape(min(min(v, [], 1), [], 2), 1, nt));
        end
    end
end
if any(isfinite(b_min))
    [~, i_b] = min(b_min);
    idx(end + 1) = i_b;
end
idx = unique(idx);
S = m.S_out(idx, :)';
t = m.t_out(idx)';
end

function [S, t] = stored_states(ref_data, path)
node = ref_data.pointwise;
for k = 1:numel(path)
    node = node.(path{k});
end
S = node.S;
t = node.t;
end

function c = pointwise_case(m, S, t)
% Every single-evaluation quantity at the columns of S.
p = m.get_params();
p.u_interpolant = m.u_interpolant;       % get_params leaves it to run()
[N, ns] = size(S);
Y = cos((1:N)' * (1:4)) / sqrt(N);        % fixed, not random
c = struct();
c.view = parameter_view(m, p);
c.t = t;
c.S = S;
c.dS = zeros(N, ns);
c.J = cell(1, ns);
c.JY = zeros(N, 4, ns);
for k = 1:ns
    c.dS(:, k) = SRNNCellTypePairs.dynamics_fast(t(k), S(:, k), p);
    c.J{k} = SRNNCellTypePairs.compute_Jacobian_fast(S(:, k), p);
    c.JY(:, :, k) = SRNNCellTypePairs.jacobian_times(S(:, k), Y, p);
end
[~, ~, ~, r, theta] = SRNNCellTypePairs.unpack_and_compute_states(S', p);
c.rate = r;
c.synaptic_output = theta;
end

function v = parameter_view(m, p)
% The physics the model integrates, read through accessors that keep their
% signatures. NOT the params struct, whose fields change in the refactor.
v = struct();
v.n_per_type = m.n_per_type;
v.N = m.N_sys_eqs;
v.W_scaled = p.W;                        % what the dynamics read: route scale folded in
v.tau_d = m.tau_d;
v.c_eff = p.c_eff;
v.tau_a_matrix = p.tau_a_matrix;
v.activation = m.activation;
v.S_a = m.S_a;
v.S_c = m.S_c;
v.S_c_vec = m.S_c_vec;
v.std_zero_floor = m.std_zero_floor;
names = m.cell_type_names;
for pre = 1:numel(names)
    for post = 1:numel(names)
        v.route.([names{pre} '_to_' names{post}]) = SRNNCellTypePairs.route_fields(p, pre, post);
    end
end
end

function r = end_to_end(m)
% What a run produced. S_out is kept at every fifth sample plus the last:
% a difference anywhere along a trajectory is still present five samples
% later, and the fixture stays small.
nt = size(m.S_out, 1);
keep = unique([1:5:nt, nt]);
r = struct();
r.W = m.W;
r.S0 = m.S0;
r.t_kept = m.t_out(keep);
r.S_kept = m.S_out(keep, :);
r.plot_data = m.plot_data;
r.lya_results = drop_volatile(m.lya_results);
r.lya_summary = drop_volatile(m.lya_summary());
r.tau_a_clipped = m.tau_a_clipped;
end

function s = drop_volatile(s)
% Wall-clock fields differ between identical runs and are not results; the
% params struct some results carry is the model's internal packing, which
% the refactor changes (the canonical view stands in for it).
if ~isstruct(s); return; end
volatile = {'seconds', 'run_duration', 'elapsed', 'wall_time', 'params'};
s = rmfield(s, intersect(fieldnames(s), volatile));
end

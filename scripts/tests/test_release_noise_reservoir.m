%% test_release_noise_reservoir.m
% Release noise in the memory-capacity reservoir (FD-027, step N2):
% SRNN_ESN_reservoir.run_memory_capacity_lean with release noise on, the
% noise held a segment at a time, and the 'release_filtered' readout.
%
% Checks:
%   1. wiener_segments hands out wiener_increments' two tensors window by
%      window, bit for bit, with and without a seed schedule, and leaves
%      the global random stream alone.
%   2. filter_release on made-up input: a constant passes unchanged, a ramp
%      and an impulse follow the closed form, and cutting the samples into
%      pieces changes no bit.
%   3. With release noise the lean path does not depend on segment_len.
%      BIT FOR BIT (isequal), for every readout and with the Wiener term on
%      too: a segment is handed the increments of its own absolute steps,
%      which are the numbers a single pass consumes there, so every step,
%      and every filter update, is the same arithmetic on the same numbers.
%   4. The filtered readout against a direct computation written here from
%      the equations (the traced states of a single pass, the increments
%      drawn again from the seed, the closed-form release of a chain, a
%      two-pool product and an undepressed route, and the filter's closed
%      form), to round-off: release noise on, Wiener on, both, neither.
%   5. With both noises off it is a low-pass of the synaptic output
%      (MATLAB's filter, the same coefficients), and the two
%      memory-capacity paths give the same bits for it; with the Wiener
%      term alone too.
%   6. The realized release: released_by_kick is what release_kick put
%      into x (W dTheta / tau_d); a step's kick is never below the kick of
%      the smallest increment, -r_s dt; release_per_step is the trapezoid
%      of the mean output plus the kick and is never below -r_s dt times
%      half the largest change of the output factor over a step (the
%      bound derived in the comment there); the filter's output stays
%      above zero.
%   7. Route redundancy with release noise: two routes of one neuron keep
%      identical pools, synaptic output and kick release (isequal), and a
%      type whose routes differ is refused for the filtered readout.
%   8. Memory capacity end to end on a small network of the paper's
%      preset: the four conditions share their build (verify_shared_build),
%      all give a finite memory capacity with release noise, and the
%      Lyapunov run() of the same reservoir gives finite top-K exponents
%      and differs from the leak-only run (the kicks are in it).
%   9. The path that holds the whole trajectory (run_memory_capacity,
%      run_reservoir_esn) refuses release noise by name; an unknown
%      readout is refused before any integration.
%
% Prints [PASS]/[FAIL] per check and errors at the end if anything failed.
% About 30 s. Assumes setup_paths has run.

n_fail = 0;
fs = 200;
dt = 1 / fs;
r_s = 0.01;
sigma_g = 0.014;
net = { ...
    'n', 60, 'indegree', 30, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [3 -4; 3 -4], 'sigma_tilde_relative', [1 1; 1 1], ...
    'fs', fs, 'level_of_chaos', 1.2, 'S_c', 0.35, 'c', [0.15 0], 'rng_seeds', [42, 43], ...
    'tau_a', {[0.1 1.0 10], zeros(1, 0)}, 'ode_solver', 'sra1', 'lya_method', 'none'};
proto = {'input_type', 'sample_hold', 'T_hold', 0.15, ...
    'T_wash', 2*fs, 'T_train', 10*fs, 'T_test', 5*fs, 'd_max', 2*fs};
release_on = {'release_rate', [r_s r_s], 'release_sigma', [sigma_g sigma_g]};

% Two networks. CHAIN: a vesicle chain on all four routes. MIXED: a
% two-pool product on the E routes and no synapse dynamics on the I routes.
chain = struct('model', 'chain', 'tau_rec', [1 10], 'tau_rel', 0.25, 'n2', 3);
sc_chain = struct();
sc_chain.E.E.std = chain;  sc_chain.E.I.std = chain;
sc_chain.I.E.std = chain;  sc_chain.I.I.std = chain;
two = struct('tau_rec', [0.5 4], 'tau_rel', [0.25 1.5]);
sc_mixed = struct();
sc_mixed.E.E.std = two;  sc_mixed.E.I.std = two;
networks = struct('name', {'chain', 'mixed'}, 'synapse_config', {sc_chain, sc_mixed});

%% ---- 1. wiener_segments ----
fprintf('=== 1. wiener_segments ===\n');
n_rows = 7;  n_steps = 3000;
global_before = rng;
[xi1, xi2] = wiener_increments(n_rows, n_steps, 4321);
w = wiener_segments(n_rows, n_steps, 4321);
[p1, p2] = in_pieces(w, [1 599 600 3 1500 297]);
n_fail = n_fail + report(isequal(xi1, p1) && isequal(xi2, p2), ...
    'no schedule: six windows of unequal length are the two tensors, bit for bit');
schedule = [2.5 11; 7 4321; 9.005 99];
[xi1, xi2] = wiener_increments(n_rows, n_steps, 4321, schedule, 0, fs, 'test');
w = wiener_segments(n_rows, n_steps, 4321, schedule, 0, fs, 'test');
[p1, p2] = in_pieces(w, [1 999 801 1 1198]);
n_fail = n_fail + report(isequal(xi1, p1) && isequal(xi2, p2), ...
    'a three-row seed schedule: the spliced tensors, bit for bit');
n_fail = n_fail + report(isequal(global_before, rng) ...
    && strcmp(error_id(@() w.next(1)), 'wiener_segments:OutOfRange'), ...
    'the global stream is untouched; a window past the run''s last step is refused');

%% ---- 2. the filter on made-up input ----
fprintf('\n=== 2. filter_release ===\n');
tau = 0.1;
m = 400;
t = (0:m - 1) * dt;
flat = repmat([0.3; 0.7], 1, m);
Y = SRNN_ESN_reservoir.filter_release(flat, zeros(2, m), [], tau, dt);
n_fail = n_fail + report(max(abs(Y - flat), [], 'all') < 1e-14, ...
    sprintf('a constant input is passed unchanged (max difference %.1e)', max(abs(Y - flat), [], 'all')));
% A ramp theta = c + s t from rest on its first value: y = c + s (t - tau) + s tau exp(-t/tau).
slope = 2;
ramp = 0.3 + slope * t;
Y = SRNN_ESN_reservoir.filter_release(ramp, zeros(1, m), [], tau, dt);
want = 0.3 + slope * (t - tau) + slope * tau * exp(-t / tau);
err_ramp = max(abs(Y - want));
% An impulse of size q at the end of step 100 on zero input: q/tau exp(-(t - t_100)/tau).
kick = zeros(1, m);
kick(101) = 0.05;
Y = SRNN_ESN_reservoir.filter_release(zeros(1, m), kick, [], tau, dt);
want = [zeros(1, 100), (0.05 / tau) * exp(-(t(101:end) - t(101)) / tau)];
err_kick = max(abs(Y - want));
n_fail = n_fail + report(err_ramp < 1e-12 && err_kick < 1e-13, ...
    sprintf('a ramp and an impulse follow the closed form (max errors %.1e, %.1e): no step-size error', ...
    err_ramp, err_kick));
rng(5, 'twister');
theta = rand(3, m);  kicks = 0.01 * randn(3, m);
whole = SRNN_ESN_reservoir.filter_release(theta, kicks, [], tau, dt);
cut = [1 2 50 51 299 m];
pieces = zeros(3, 0);  state = [];
for k = 1:numel(cut)
    if k == 1; cols = 1:cut(1); else; cols = cut(k - 1) + 1:cut(k); end
    [part, state] = SRNN_ESN_reservoir.filter_release(theta(:, cols), kicks(:, cols), state, tau, dt);
    pieces = [pieces, part]; %#ok<AGROW>
end
n_fail = n_fail + report(isequal(whole, pieces), ...
    'cut into six pieces (one of a single sample): the same bits');

%% ---- 3. independent of segment_len, bit for bit ----
fprintf('\n=== 3. segment_len ===\n');
runs = struct('label', {'release noise', 'release and Wiener noise'}, ...
    'extra', {release_on, [release_on, {'sigma_u_noise', 0.02}]});
for v = 1:numel(networks)
    for k = 1:numel(runs)
        args = [net, proto, {'synapse_config', networks(v).synapse_config}, runs(k).extra];
        signals = {'release_filtered'};
        if v == 1 && k == 2
            signals = {'release_filtered', 'synaptic', 'rate'};
        end
        for sig = signals
            ref = lean_run(args, sig{1}, {'segment_len', 1e6});     % a single pass
            same = true;
            for seg = [23, 600]
                same = same && same_scores(ref, lean_run(args, sig{1}, {'segment_len', seg}));
            end
            n_fail = n_fail + report(same, sprintf( ...
                '%s, %s, %s readout: segments of 23 and 600 = a single pass, bit for bit (MC = %.4f)', ...
                networks(v).name, runs(k).label, sig{1}, ref.MC));
        end
    end
end

%% ---- 4 to 6. the readout and the realized release against the equations ----
fprintf('\n=== 4-6. the filtered readout, written out ===\n');
cases = struct('label', {'release noise', 'release and Wiener noise', 'Wiener noise', 'no noise'}, ...
    'extra', {release_on, [release_on, {'sigma_u_noise', 0.02}], {'sigma_u_noise', 0.02}, {}});
for v = 1:numel(networks)
    for k = 1:numel(cases)
        args = [net, proto, {'synapse_config', networks(v).synapse_config}, cases(k).extra];
        esn = SRNN_ESN_reservoir(args{:});
        esn.build();
        t_end = esn.t_ex(end);
        [~, ~, res] = esn.run_memory_capacity_lean('readout_signal', 'release_filtered', ...
            'segment_len', 173, 'trace_window', [0, t_end]);
        tr = res.trace;
        params = esn.cached_params;
        n_steps = numel(tr.t) - 1;
        has_release = any(esn.release_sigma > 0);

        % The increments, drawn again from the seed, centered.
        dRc = zeros(esn.n, n_steps);
        if has_release
            dRc = release_increments(repmat(r_s, esn.n, 1), sigma_g, fs, [1, n_steps], esn.rng_seeds(1) + 224737) - r_s * dt;
        end
        % What each step's kick released, from the pools after it.
        [kicked_E, F_E] = written_out_release(networks(v).name, 'E', tr, dRc(params.type_indices{1}, :));
        [kicked_I, F_I] = written_out_release(networks(v).name, 'I', tr, dRc(params.type_indices{2}, :));
        kicked = [kicked_E; kicked_I];
        theta = [tr.synaptic_output.E.E; tr.synaptic_output.I.E];

        % 4. The filter, by its closed form per step: with the input the
        % line theta0 + slope s, y(h) = yp(h) + (y(0) - yp(0)) exp(-h/tau),
        % yp(s) = theta0 + slope (s - tau); then the kick, / tau.
        tau = esn.tau_d;
        y = theta(:, 1);
        direct = zeros(size(theta));
        direct(:, 1) = y;
        for j = 1:n_steps
            slope = (theta(:, j + 1) - theta(:, j)) / dt;
            y = theta(:, j) + slope * (dt - tau) + (y - theta(:, j) + slope * tau) * exp(-dt / tau) ...
                + kicked(:, j) / tau;
            direct(:, j + 1) = y;
        end
        largest = max(abs(direct), [], 'all');
        err = max(abs(tr.release_filtered - direct), [], 'all') / largest;
        n_fail = n_fail + report(err < 1e-11, sprintf( ...
            '%s, %s: the readout is the written-out filter of the written-out release (relative difference %.1e)', ...
            networks(v).name, cases(k).label, err));

        if ~has_release
            % 5. A low-pass of the synaptic output: MATLAB's filter.
            a = exp(-dt / tau);
            c1 = 1 - (tau / dt) * (1 - a);
            c0 = (1 - a) - c1;
            low = zeros(size(theta));
            for i = 1:size(theta, 1)
                % Started at rest on its first value: filter the deviation.
                low(i, :) = theta(i, 1) + filter([c1 c0], [1 -a], theta(i, :) - theta(i, 1));
            end
            err = max(abs(tr.release_filtered - low), [], 'all') / largest;
            whole_path = SRNN_ESN_reservoir(args{:});
            whole_path.build();
            [~, ~, res_full] = whole_path.run_memory_capacity('readout_signal', 'release_filtered', 'verbose', false);
            n_fail = n_fail + report(err < 1e-11 && same_scores(res_full, res), sprintf( ...
                ['%s, %s: a low-pass of the synaptic output (relative difference %.1e); ' ...
                 'run_memory_capacity gives the same bits'], networks(v).name, cases(k).label, err));
            continue
        end

        % 6. The realized release.
        % The smallest increment is -r_s dt, and the kick is increasing in it.
        floor_E = written_out_release(networks(v).name, 'E', tr, repmat(-r_s * dt, size(kicked_E)));
        floor_I = written_out_release(networks(v).name, 'I', tr, repmat(-r_s * dt, size(kicked_I)));
        slack = 1e-15;
        ok_floor = all(kicked_E >= floor_E - slack, 'all') && all(kicked_I >= floor_I - slack, 'all');
        % A step's release, trapezoid of the mean output plus the kick.
        % theta >= F r_s at both ends and the kick >= -F_after r_s dt, so
        % the sum >= r_s dt (F_before - F_after) / 2: below zero only while
        % the factor recovers, and then by no more than r_s dt times half
        % its change over the step.
        release = (dt / 2) * (theta(:, 1:end - 1) + theta(:, 2:end)) + kicked;
        F = [F_E; F_I];
        bound = -r_s * dt * max(abs(diff(F, 1, 2)), [], 'all') / 2;
        mean_release = mean(release, 2) / dt;           % per unit time
        mean_theta = mean(theta, 2);
        n_fail = n_fail + report(ok_floor && min(release, [], 'all') >= bound - slack ...
            && min(tr.release_filtered, [], 'all') > 0 ...
            && max(abs(mean_release - mean_theta) ./ mean_theta) < 0.5, sprintf( ...
            ['%s, %s: no kick below the kick of -r_s dt; smallest release of a step %.2e ' ...
             '(bound %.2e; typical %.1e); the readout stays above zero (min %.3g)'], ...
            networks(v).name, cases(k).label, min(release, [], 'all'), bound, ...
            median(release, 'all'), min(tr.release_filtered, [], 'all')));
    end
end

% released_by_kick and release_per_step are those numbers, and the kick's.
fprintf('\n=== 6b. released_by_kick is what release_kick put into x ===\n');
for v = 1:numel(networks)
    args = [net, proto, {'synapse_config', networks(v).synapse_config}, release_on];
    esn = SRNN_ESN_reservoir(args{:});
    esn.build();
    params = esn.cached_params;
    rng(11, 'twister');
    S = esn.S0;
    for q = 1:2
        for p = 1:2
            idx = params.route{q, p}.idx;
            if p == 1
                pools = 0.2 + 0.7 * rand(numel(idx), 1);   % the same pools on a neuron's two routes
            end
            S(idx) = pools;
        end
    end
    step_dRc = release_increments(repmat(r_s, esn.n, 1), 0.2, fs, [1 1], 9) - r_s * dt;   % large draws
    step_dRc(1:5) = [-r_s * dt; 0; 1; 5; 1e-9];
    S_after = SRNNCellTypePairs.release_kick(S, step_dRc, params);
    kicked = SRNNCellTypePairs.released_by_kick(S_after.', step_dRc, params);
    x_rows = params.state_layout.x;
    want = zeros(esn.n, 1);
    for p = 1:2
        for q = 1:2
            want(params.type_indices{p}) = want(params.type_indices{p}) + ...
                params.W(params.type_indices{p}, params.type_indices{q}) * ...
                kicked.(params.cell_type_names{q}).(params.cell_type_names{p}) / params.tau_d;
        end
    end
    moved = S_after(x_rows) - S(x_rows);
    err = max(abs(moved - want)) / max(abs(moved));
    % Many steps at once = step by step; release_per_step composes them.
    S3 = [S.'; S_after.'; S_after.'];
    two_steps = [step_dRc, step_dRc / 2];
    [release, theta, kicked2] = SRNNCellTypePairs.release_per_step(S3, two_steps, params, dt);
    one = SRNNCellTypePairs.released_by_kick(S_after.', step_dRc / 2, params);
    none = SRNNCellTypePairs.released_by_kick(S3(2:end, :), [], params);
    ok = isequal(kicked2.E.I(:, 1), kicked.E.I) && isequal(kicked2.I.E(:, 2), one.I.E) ...
        && isequal(release.E.E, (dt / 2) * (theta.E.E(:, 1:2) + theta.E.E(:, 2:3)) + kicked2.E.E) ...
        && ~any(none.E.E, 'all') && isequal(size(none.I.I), [params.n_per_type(2), 2]);
    n_fail = n_fail + report(err < 1e-12 && ok, sprintf( ...
        ['%s: W dTheta / tau_d from the recovered release is the jump of x (relative difference %.1e); ' ...
         'two steps at once = one at a time; no increments, no release'], networks(v).name, err));
end

%% ---- 7. route redundancy ----
fprintf('\n=== 7. route redundancy with release noise ===\n');
args = [net, proto, {'synapse_config', sc_chain}, release_on, {'sigma_u_noise', 0.02}];
esn = SRNN_ESN_reservoir(args{:});
esn.build();
[~, ~, res] = esn.run_memory_capacity_lean('readout_signal', 'release_filtered', ...
    'trace_window', [0, esn.t_ex(end)]);
tr = res.trace;
n_steps = numel(tr.t) - 1;
dRc = release_increments(repmat(r_s, esn.n, 1), sigma_g, fs, [1, n_steps], esn.rng_seeds(1) + 224737) - r_s * dt;
n_fail = n_fail + report(isequal(tr.b.E.E, tr.b.E.I) && isequal(tr.b.I.E, tr.b.I.I) ...
    && isequal(tr.synaptic_output.E.E, tr.synaptic_output.E.I) ...
    && isequal(written_out_release('chain', 'E', tr, dRc(esn.cached_params.type_indices{1}, :)), ...
        chain_release(tr.b.E.I, dRc(esn.cached_params.type_indices{1}, :), 0.25)) ...
    && min(tr.b.E.E, [], 'all') >= 0, ...
    'the two routes of a neuron keep identical pools, synaptic output and kick release; no pool below zero');
sc_differ = sc_chain;
sc_differ.E.I.std.n2 = 4;
differ = SRNN_ESN_reservoir(net{:}, proto{:}, 'synapse_config', sc_differ, release_on{:});
differ.build();
n_fail = n_fail + report(strcmp(error_id(@() differ.run_memory_capacity_lean('readout_signal', 'release_filtered')), ...
    'SRNN_ESN_reservoir:AmbiguousSynapticOutput'), ...
    'a type whose routes differ has no one release per neuron: the filtered readout is refused');

%% ---- 8. memory capacity end to end, four conditions, and the Lyapunov run ----
fprintf('\n=== 8. the paper''s preset, four conditions ===\n');
preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
[d, ~, conditions] = srnn_param_preset(preset_name);
preset_args = struct2namevalue(d);
small = {'n', 50, 'indegree', 25, 'rng_seeds', [3001 4001], 'fs', fs, 'ode_solver', 'sra1', ...
    'input_type', 'sample_hold', 'T_hold', 0.3, 'verbose', 'near-none'};
mc_window = {'T_wash', 2*fs, 'T_train', 12*fs, 'T_test', 6*fs, 'd_max', 1.5*fs, 'lya_method', 'none', ...
    'lean_input', true};
lya_window = {'T_wash', 3*fs, 'T_train', 3*fs, 'T_test', 0, 'd_max', 1.5*fs, ...
    'lya_method', 'topk', 'lya_K', 3, 'lya_K_auto', false, 'lya_dt', 0.05, ...
    'lya_T_interval', [3, 6 - dt], 'lya_warmup', 3, ...
    'store_decimated_state', false, 'store_full_state', false};
n_cond = numel(conditions);
esn = cell(1, n_cond);
mc = nan(1, n_cond);
lle = nan(n_cond, 2);
for c = 1:n_cond
    cond_args = struct2namevalue(rmfield(conditions{c}, 'name'));
    esn{c} = SRNN_ESN_reservoir(preset_args{:}, small{:}, cond_args{:}, mc_window{:}, ...
        'release_rate_all', r_s, 'release_sigma_all', sigma_g);
    esn{c}.build();
    mc(c) = esn{c}.run_memory_capacity_lean('readout_signal', 'release_filtered');
    levels = [sigma_g, 0];          % with the kicks, and the leak alone
    for k = 1:2
        m = SRNN_ESN_reservoir(preset_args{:}, small{:}, cond_args{:}, lya_window{:}, ...
            'release_rate_all', r_s, 'release_sigma_all', levels(k));
        m.build();
        m.run();
        lle(c, k) = m.lya_results.LLE;
        left_over = ~isempty(m.noise_increments);
    end
end
verify_shared_build(esn, {'tau_a', 'synapse_config'}, {'W', 'W_in', 'u_scalar', 't_ex'});
n_fail = n_fail + report(all(isfinite(mc)) && all(mc > 0), sprintf( ...
    'four conditions share W, the input weights and the input; memory capacity with release noise: %s', ...
    strtrim(sprintf('%.3f ', mc))));
n_fail = n_fail + report(all(isfinite(lle), 'all') && all(lle(:, 1) ~= lle(:, 2)) && ~left_over, sprintf( ...
    'run() gives top-K exponents with release noise (lambda_1 %s; leak only %s), and clears the increments after', ...
    strtrim(sprintf('%+.3f ', lle(:, 1))), strtrim(sprintf('%+.3f ', lle(:, 2)))));

%% ---- 9. refusals ----
fprintf('\n=== 9. refusals ===\n');
args = [net, proto, {'synapse_config', sc_chain}, release_on];
esn = SRNN_ESN_reservoir(args{:});
esn.build();
leak = SRNN_ESN_reservoir(net{:}, proto{:}, 'synapse_config', sc_chain, 'release_rate', [r_s r_s]);
leak.build();
mc_leak = leak.run_memory_capacity('readout_signal', 'synaptic', 'verbose', false);
n_fail = n_fail + report(strcmp(error_id(@() esn.run_memory_capacity('verbose', false)), ...
        'SRNN_ESN_reservoir:ReleaseNoiseNeedsLeanPath') ...
    && strcmp(error_id(@() esn.run_reservoir_esn()), 'SRNN_ESN_reservoir:ReleaseNoiseNeedsLeanPath') ...
    && isfinite(mc_leak), ...
    'run_memory_capacity and run_reservoir_esn refuse release noise, naming the lean path; a leak alone runs there');
n_fail = n_fail + report(strcmp(error_id(@() esn.run_memory_capacity_lean('readout_signal', 'released')), ...
        'SRNN_ESN_reservoir:BadReadoutSignal') ...
    && strcmp(error_id(@() SRNN_ESN_reservoir.readout_from_states(struct(), struct(), ...
        esn.cached_params, 'release_filtered')), 'SRNN_ESN_reservoir:BadReadoutSignal'), ...
    'an unknown readout is refused before any integration; readout_from_states does not make the filtered one');

%% ---- Verdict ----
fprintf('\n');
if n_fail == 0
    fprintf('ALL CHECKS PASSED\n');
else
    error('test_release_noise_reservoir:Failed', '%d check(s) FAILED', n_fail);
end

%% ---- helpers ----
function [xi1, xi2] = in_pieces(w, counts)
xi1 = zeros(w.n_rows, 0);  xi2 = xi1;
for c = counts
    [a, b] = w.next(c);
    xi1 = [xi1, a];  xi2 = [xi2, b]; %#ok<AGROW>
end
end

function res = lean_run(args, readout_signal, opts)
esn = SRNN_ESN_reservoir(args{:});
esn.build();
[~, ~, res] = esn.run_memory_capacity_lean('verbose', false, 'readout_signal', readout_signal, opts{:});
end

function ok = same_scores(a, b)
ok = isequal(a.R2_d, b.R2_d) && a.MC == b.MC && isequal(a.weights, b.weights) ...
    && isequal(a.predictions, b.predictions) && isequal(a.t_pred, b.t_pred);
end

function [kicked, F] = written_out_release(network, type, tr, dRc)
% What the kick of each step released on the first route of a presynaptic
% type, from the traced pools AFTER the step and the step's centered
% increment, by the equations of the route (route_pool, route_product:
% "the release kick"); and the route's output factor F at every sample.
%   chain            tau_rel b1 expm1(dRc / tau_rel),            F = b1
%   two-pool product F tau_e expm1(dRc / tau_e), F = b1 b2, 1/tau_e = 1/tau_rel1 + 1/tau_rel2
%   no pool          dRc,                                         F = 1
if strcmp(network, 'chain')
    b = tr.b.(type).E;
    kicked = chain_release(b, dRc, 0.25);
    F = squeeze(b(:, 1, :));
elseif strcmp(type, 'E')
    b = tr.b.E.E;
    F = squeeze(b(:, 1, :) .* b(:, 2, :));
    tau_e = 1 / (1 / 0.25 + 1 / 1.5);
    kicked = F(:, 2:end) .* tau_e .* expm1(dRc ./ tau_e);
else
    kicked = dRc;
    F = ones(size(dRc, 1), size(dRc, 2) + 1);
end
end

function kicked = chain_release(b, dRc, tau_rel)
b1_after = squeeze(b(:, 1, 2:end));
kicked = tau_rel .* b1_after .* expm1(dRc ./ tau_rel);
end

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

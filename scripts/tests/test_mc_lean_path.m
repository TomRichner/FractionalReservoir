%% test_mc_lean_path.m
% SRNN_ESN_reservoir.run_memory_capacity_lean against run_memory_capacity
% on the same reservoir.
%
% The lean path integrates the memory-capacity task segment by segment and
% never holds the state trajectory; the existing path integrates it in one
% call. Both score with the same methods (mc_readout_plan,
% readout_from_states, score_readout).
%
% THE TOLERANCE IS ZERO, and why. A fixed-step integrator takes each step
% from the time grid it is handed and the state alone, so cutting the grid
% into segments that share their end points repeats the single call's
% arithmetic step for step: the states are the same bits. The readout
% matrix is then the same bits, and the scoring is the same code on the
% same matrix. So every comparison below is isequal, not a tolerance. A
% difference of one bit would mean a segment boundary changed a step.
%
% Checks:
%   1. sample_hold, 'rk4', rate and synaptic readouts, an adapting network
%      with depression: R2 curve, MC, weights, predictions, t_pred equal.
%   2. Segment lengths that do not divide the hold or the sequence (7, 61,
%      and one longer than the sequence) give the same bits.
%   3. 'sra1' with noise OFF, and with noise ON (the increments are indexed
%      by absolute time, so segments consume what the single call would).
%   4. 'white' input, scored at sample resolution (hold_len = 1).
%   5. rate_mean equals the time mean of the existing path's stored rates
%      over the training and test windows (to round-off: a sum taken in
%      segments is not the same floating-point sum as one taken at once).
%   6. trace_window returns the existing path's stored series over that
%      window, bit for bit, and nothing outside it.
%   7. The lean path leaves no trajectory on the object (S_out empty).
%   8. lean_input = true: build() makes no n x T input table, the lean
%      path gives the same bits, and the existing path refuses to run.
%
% Prints [PASS]/[FAIL] per check and errors at the end if anything failed.
% About 15 s. Assumes setup_paths has run.

n_fail = 0;
fs = 200;
T_hold = 0.15;
net = { ...
    'n', 60, 'indegree', 30, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [3 -4; 3 -4], 'sigma_tilde_relative', [1 1; 1 1], ...
    'fs', fs, 'level_of_chaos', 2, 'S_c', 0.35, ...
    'c', [0.15 0], 'rng_seeds', [42, 43]};
std_E = struct('tau_rec', 1.0, 'tau_rel', 0.25);
sc_on = struct();
sc_on.E.E.std = std_E;
sc_on.E.I.std = std_E;
adapt = {'tau_a', {[0.1 1.0 10], zeros(1, 0)}, 'synapse_config', sc_on};
proto = {'input_type', 'sample_hold', 'T_hold', T_hold, ...
    'T_wash', 4*fs, 'T_train', 20*fs, 'T_test', 10*fs, 'd_max', 3*fs};

%% ---- 1. rk4, both readouts ----
fprintf('=== 1. rk4, rate and synaptic readouts ===\n');
for sig = {'rate', 'synaptic'}
    [full, lean, esn_lean] = both_paths([net, adapt, proto, {'ode_solver', 'rk4'}], sig{1}, {});
    n_fail = n_fail + report(same_scores(full, lean), ...
        sprintf('%s readout: R2, MC, weights, predictions, t_pred bit-identical (MC = %.4f, max |dR2| = %.1e)', ...
        sig{1}, lean.MC, max(abs(full.R2_d - lean.R2_d))));
    n_fail = n_fail + report(isempty(esn_lean.S_out) && isempty(lean.lya_results) && lean.lean, ...
        sprintf('%s readout: the lean path leaves no trajectory on the object', sig{1}));
end

%% ---- 2. segment lengths ----
fprintf('\n=== 2. segment lengths ===\n');
for seg = [7, 61, 1e6]
    [full2, lean2] = both_paths([net, adapt, proto, {'ode_solver', 'rk4'}], 'synaptic', {'segment_len', seg});
    n_fail = n_fail + report(same_scores(full2, lean2), ...
        sprintf('segment_len = %g: bit-identical', seg));
end

%% ---- 3. sra1, noise off and on ----
fprintf('\n=== 3. sra1, noise off and on ===\n');
for sigma = [0, 0.02]
    [full3, lean3] = both_paths([net, adapt, proto, ...
        {'ode_solver', 'sra1', 'sigma_u_noise', sigma}], 'synaptic', {'segment_len', 100});
    n_fail = n_fail + report(same_scores(full3, lean3), ...
        sprintf('sra1, sigma_u_noise = %g: bit-identical (MC = %.4f)', sigma, lean3.MC));
end

%% ---- 4. white input at sample resolution ----
fprintf('\n=== 4. white input (hold_len = 1) ===\n');
[full4, lean4] = both_paths([net, adapt, {'ode_solver', 'rk4', 'input_type', 'white', ...
    'T_wash', 2*fs, 'T_train', 15*fs, 'T_test', 10*fs, 'd_max', 25}], 'rate', {'segment_len', 333});
n_fail = n_fail + report(lean4.hold_len == 1 && same_scores(full4, lean4), ...
    'white input: bit-identical at sample resolution');

%% ---- 5 and 6. rate_mean and the trace window ----
fprintf('\n=== 5-6. time-mean rate and trace window ===\n');
win = [25, 28];        % s, inside the test window (24 s to 34 s)
[full5, lean5] = both_paths([net, adapt, proto, {'ode_solver', 'rk4'}], 'synaptic', ...
    {'trace_window', win});
rows = lean5.rate_window(1):lean5.rate_window(2);
r_full = [full5.r.E; full5.r.I];
want = mean(r_full(:, rows), 2);
err = max(abs(lean5.rate_mean - want));
% The window is in whole holds (30 samples): 27 washout holds, then the
% 133 training and 66 test holds that fit the 226-hold sequence.
n_fail = n_fail + report(isequal(lean5.rate_window, [27*30 + 1, 226*30]) && err < 1e-13, ...
    sprintf('rate_mean is the time mean over samples %d..%d (max difference %.1e)', rows(1), rows(end), err));
tr = lean5.trace;
k = tr.rows;
ok = all(tr.t >= win(1) & tr.t <= win(2)) && numel(k) == round(diff(win) * fs) + 1 ...
    && isequal(tr.r.E, full5.r.E(:, k)) && isequal(tr.x.I, full5.x.I(:, k)) ...
    && isequal(tr.a.E, full5.a.E(:, :, k)) ...
    && isequal(tr.synaptic_output.E.I, full5.synaptic_output.E.I(:, k)) ...
    && isequal(tr.b.E.E, full5.b.E.E(:, :, k)) && isequal(tr.u_ex, full5.u_ex(:, k));
n_fail = n_fail + report(ok, sprintf('trace_window [%g %g] s: %d samples, the stored series bit for bit', ...
    win(1), win(2), numel(k)));

%% ---- 8. lean_input: the n x T input table is never built ----
fprintf('\n=== 8. lean_input ===\n');
args8 = [net, adapt, proto, {'ode_solver', 'sra1', 'lya_method', 'none'}];
full8 = SRNN_ESN_reservoir(args8{:});
full8.build();
[~, ~, res_full8] = full8.run_memory_capacity('verbose', false, 'readout_signal', 'synaptic');
lean8 = SRNN_ESN_reservoir(args8{:}, 'lean_input', true);
lean8.build();
[~, ~, res_lean8] = lean8.run_memory_capacity_lean('readout_signal', 'synaptic', 'trace_window', win);
n_fail = n_fail + report(isempty(lean8.u_ex) && isempty(lean8.u_interpolant), ...
    'lean_input = true: build() makes neither u_ex nor its interpolant');
n_fail = n_fail + report(same_scores(res_full8, res_lean8) ...
    && isequal(res_lean8.trace.u_ex, res_full8.u_ex(:, res_lean8.trace.rows)), ...
    'lean_input = true: scores and the traced input bit-identical to the existing path');
threw = '';
try
    lean8.run_memory_capacity('verbose', false);
catch ME
    threw = ME.identifier;
end
n_fail = n_fail + report(strcmp(threw, 'SRNN_ESN_reservoir:LeanInput'), ...
    sprintf('the existing path refuses a lean_input reservoir (%s)', threw));

%% ---- Verdict ----
fprintf('\n');
if n_fail == 0
    fprintf('ALL CHECKS PASSED\n');
else
    error('test_mc_lean_path:Failed', '%d check(s) FAILED', n_fail);
end

%% ---- helpers ----
function [full, lean, esn_lean] = both_paths(args, readout_signal, lean_opts)
% The same reservoir built twice from the same arguments, scored by the
% existing path (full trajectory kept) and by the lean path.
esn_full = SRNN_ESN_reservoir(args{:}, 'lya_method', 'none');
esn_full.build();
[~, ~, full] = esn_full.run_memory_capacity('verbose', false, 'readout_signal', readout_signal);
esn_lean = SRNN_ESN_reservoir(args{:}, 'lya_method', 'none');
esn_lean.build();
[~, ~, lean] = esn_lean.run_memory_capacity_lean('verbose', false, ...
    'readout_signal', readout_signal, lean_opts{:});
end

function ok = same_scores(full, lean)
ok = isequal(full.R2_d, lean.R2_d) && full.MC == lean.MC ...
    && isequal(full.weights, lean.weights) && isequal(full.predictions, lean.predictions) ...
    && isequal(full.t_pred, lean.t_pred) && isequal(full.u_test, lean.u_test) ...
    && isequal(full.delay_s, lean.delay_s) && isequal(full.d, lean.d);
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

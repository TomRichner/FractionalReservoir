%% test_sample_hold_mc.m
% A tiny memory-capacity run through SRNN_ESN_reservoir, and the
% sample-and-hold input (input_type = 'sample_hold' + hold-unit MC
% decimation).
%
% Verifies:
%   1. sample_hold: hold_len, R2 length (= floor(d_max/hold_len)), delay_s
%      axis, finite MC in a sensible range, and that verify_shared_build
%      passes (W, W_in and u_scalar identical across conditions).
%   2. The 'synaptic' readout runs on a network whose routes are redundant
%      per presynaptic type, and gives a different MC from the 'rate' one
%      (depression is engaged, so the two signals differ).
%   3. Reproducibility: the same reservoir built and run twice gives a
%      bit-identical R2 curve, with and without store_timeseries.
%   4. input_type = 'white' still runs at sample resolution
%      (hold_len = 1, numel(R2) = d_max).
%
% Prints [PASS]/[FAIL] per check and errors at the end if anything failed.
% Small networks on 'rk4', about 5 s. Assumes setup_paths has run.
%
% HISTORY. Written against the SRNNModel2-parented reservoir ('f', 0.6,
% 'n_a_E', 'n_b_E', 'c_E', 'tau_b_E_rec', ...). SRNN_ESN_reservoir was
% re-parented onto SRNNCellTypePairs on 2026-09-02 (commit 7966fc7) and this
% test was not ported, so from then until 2026-11-10 it died in the
% constructor ("Required constructor properties: n_cellTypes, ...") and
% checked nothing. Ported 2026-11-10: same network (n = 100, 60 % E,
% level_of_chaos 3, logistic at S_c = 0.35, the SRNNModel2 default weights
% mu = [3 -4] F, sigma = [1 1] F), with the adapting condition written as
% the Pairs equivalent of n_a_E = 3 / n_b_E = 1: a three-rung SFA ladder on
% E at total budget c = 0.15 (was c_E = 0.15/3 per rung) and one depression
% component on both E routes. Checks 2 and 3 are new.

n_fail = 0;

%% ---- Shared network ----
fs = 200;
T_hold = 0.15;                 % s  -> hold_len = 30 samples
hold_len_expect = round(T_hold * fs);
d_max = 3 * fs;                % 600 samples of horizon
d_max_eff_expect = floor(d_max / hold_len_expect);   % = 20

net = { ...
    'n', 100, 'indegree', 100, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [3 -4; 3 -4], 'sigma_tilde_relative', [1 1; 1 1], ...
    'fs', fs, 'level_of_chaos', 3, 'ode_solver', 'rk4', 'S_c', 0.35, ...
    'std_zero_floor', true, 'c', [0.15 0], 'rng_seeds', [42, 43]};

% Two conditions, the Pairs spelling of (n_a_E, n_b_E) = (0, 0) and (3, 1).
std_E = struct('tau_rec', 1.0, 'tau_rel', 0.25);
sc_on = struct();
sc_on.E.E.std = std_E;
sc_on.E.I.std = std_E;
cond_args = { ...
    {'tau_a', {zeros(1, 0), zeros(1, 0)}, 'synapse_config', struct()}, ...
    {'tau_a', {[0.1 1.0 10], zeros(1, 0)}, 'synapse_config', sc_on}};

%% ---- Test 1: sample_hold input + hold-unit MC ----
fprintf('=== Test 1: sample_hold ===\n');
base_sh = [net, {'input_type', 'sample_hold', 'T_hold', T_hold, ...
    'T_wash', 6*fs, 'T_train', 30*fs, 'T_test', 15*fs, 'd_max', d_max}];

esn = cell(1, numel(cond_args));
for i = 1:numel(cond_args)
    esn{i} = SRNN_ESN_reservoir(base_sh{:}, cond_args{i}{:});
    esn{i}.build();
end

verify_shared_build(esn, {'tau_a', 'synapse_config'}, {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});
fprintf('[PASS] verify_shared_build (W, W_in, u_scalar identical across conditions)\n');

MC = zeros(1, 2); R2 = cell(1, 2); res = cell(1, 2);
for i = 1:2
    [MC(i), R2{i}, res{i}] = esn{i}.run_memory_capacity('verbose', false);
end

n_fail = n_fail + report(res{1}.hold_len == hold_len_expect, ...
    sprintf('hold_len = %d (expect %d)', res{1}.hold_len, hold_len_expect));
n_fail = n_fail + report(numel(R2{1}) == d_max_eff_expect, ...
    sprintf('numel(R2) = %d (expect floor(%d/%d) = %d)', numel(R2{1}), d_max, hold_len_expect, d_max_eff_expect));
n_fail = n_fail + report(numel(res{1}.delay_s) == d_max_eff_expect, ...
    sprintf('numel(delay_s) = %d (expect %d)', numel(res{1}.delay_s), d_max_eff_expect));
n_fail = n_fail + report(abs(res{1}.delay_s(1) - hold_len_expect/fs) < 1e-9, ...
    sprintf('delay_s(1) = %.4f s (expect %.4f)', res{1}.delay_s(1), hold_len_expect/fs));
n_fail = n_fail + report(abs(res{1}.delay_s(end) - d_max_eff_expect*hold_len_expect/fs) < 1e-9, ...
    sprintf('delay_s(end) = %.4f s (expect %.4f)', res{1}.delay_s(end), d_max_eff_expect*hold_len_expect/fs));
n_fail = n_fail + report(numel(esn{1}.u_scalar) == (6+30+15)*fs, ...
    'u_scalar length = T_total (full staircase)');
for i = 1:2
    n_fail = n_fail + report(isfinite(MC(i)) && all(isfinite(R2{i})), ...
        sprintf('condition %d: MC = %.4f finite, R2 all finite', i, MC(i)));
    n_fail = n_fail + report(abs(MC(i) - sum(R2{i})) < 1e-12, ...
        sprintf('condition %d: MC is the sum of the R2 curve', i));
    % R2 is a squared correlation per delay, so MC lies in (0, number of
    % delays] for a reservoir that remembers anything at all.
    n_fail = n_fail + report(all(R2{i} >= 0 & R2{i} <= 1 + 1e-12) && MC(i) > 0, ...
        sprintf('condition %d: every R2 in [0, 1] and MC > 0', i));
end
n_fail = n_fail + report(MC(1) ~= MC(2), ...
    sprintf('adaptation changes the capacity (%.4f vs %.4f)', MC(1), MC(2)));

%% ---- Test 2: synaptic readout ----
fprintf('\n=== Test 2: synaptic readout ===\n');
esn_syn = SRNN_ESN_reservoir(base_sh{:}, cond_args{2}{:});
esn_syn.build();
[MC_syn, R2_syn] = esn_syn.run_memory_capacity('verbose', false, ...
    'readout_signal', 'synaptic', 'store_timeseries', false);
n_fail = n_fail + report(isfinite(MC_syn) && all(isfinite(R2_syn)) && numel(R2_syn) == d_max_eff_expect, ...
    sprintf('synaptic readout: MC = %.4f finite, %d delays', MC_syn, numel(R2_syn)));
n_fail = n_fail + report(MC_syn ~= MC(2), ...
    sprintf('synaptic and rate readouts differ (%.4f vs %.4f)', MC_syn, MC(2)));

%% ---- Test 3: reproducibility ----
fprintf('\n=== Test 3: reproducibility ===\n');
esn_again = SRNN_ESN_reservoir(base_sh{:}, cond_args{2}{:});
esn_again.build();
% store_timeseries = false is the lean path the run_memory_capacity analysis
% stage uses; the method's header promises MC and R2 are identical either
% way, so this one comparison checks both that promise and that a rebuild
% from the same seeds reproduces.
[MC_again, R2_again] = esn_again.run_memory_capacity('verbose', false, 'store_timeseries', false);
n_fail = n_fail + report(isequal(R2_again, R2{2}) && MC_again == MC(2), ...
    sprintf('a second build + lean run is bit-identical (max |dR2| = %.1e)', max(abs(R2_again - R2{2}))));

%% ---- Test 4: white noise at sample resolution ----
fprintf('\n=== Test 4: white (sample resolution) ===\n');
d_max_w = 25;                  % samples
esn_w = SRNN_ESN_reservoir(net{:}, cond_args{2}{:}, 'input_type', 'white', ...
    'T_wash', 2*fs, 'T_train', 15*fs, 'T_test', 15*fs, 'd_max', d_max_w);
esn_w.build();
[MC_w, R2_w, res_w] = esn_w.run_memory_capacity('verbose', false);

n_fail = n_fail + report(res_w.hold_len == 1, ...
    sprintf('hold_len = %d (expect 1)', res_w.hold_len));
n_fail = n_fail + report(numel(R2_w) == d_max_w, ...
    sprintf('numel(R2) = %d (expect d_max = %d)', numel(R2_w), d_max_w));
n_fail = n_fail + report(abs(res_w.delay_s(1) - 1/fs) < 1e-9, ...
    sprintf('delay_s(1) = %.5f s (expect %.5f)', res_w.delay_s(1), 1/fs));
n_fail = n_fail + report(isfinite(MC_w) && all(isfinite(R2_w)), ...
    sprintf('MC = %.4f finite, R2 all finite', MC_w));

%% ---- Summary ----
fprintf('\n=== Summary ===\n');
if n_fail == 0
    fprintf('All checks PASSED.\n');
else
    error('test_sample_hold_mc:Failed', '%d check(s) FAILED.', n_fail);
end

%% ---- Local helper ----
function failed = report(tf, msg)
    if tf
        fprintf('[PASS] %s\n', msg);
        failed = 0;
    else
        fprintf('[FAIL] %s\n', msg);
        failed = 1;
    end
end

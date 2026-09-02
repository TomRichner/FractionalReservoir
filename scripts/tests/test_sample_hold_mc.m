%% test_sample_hold_mc.m
% Sanity checks for the sample-and-hold memory-capacity feature added to
% SRNN_ESN_reservoir (input_type='sample_hold' + hold-unit MC decimation).
%
% Verifies:
%   1. sample_hold: hold_len, R2 length (= floor(d_max/hold_len)), delay_s axis,
%      finite MC, and that verify_shared_build passes (u_scalar identical across
%      conditions).
%   2. Backward-compat: input_type='white' still runs at sample resolution
%      (hold_len = 1, numel(R2) = d_max).
%
% Prints [PASS]/[FAIL] per check and errors at the end if anything failed.
% Uses small/fast networks (ode_rk4) so it runs in a few seconds.

clear; clc; close all;

n_fail = 0;

%% ---- Test 1: sample_hold input + hold-unit MC ----
fprintf('=== Test 1: sample_hold ===\n');
fs = 200;
T_hold = 0.15;                 % s  -> hold_len = 30 samples
hold_len_expect = round(T_hold * fs);
d_max = 3 * fs;                % 600 samples of horizon
d_max_eff_expect = floor(d_max / hold_len_expect);   % = 20

% Adaptation params are set explicitly (identical across conditions) so that
% verify_shared_build only sees n_a_E/n_b_E differ -- they are ignored by the
% dynamics when n_a_E=0 / n_b_E=0. Mirrors example_memory_capacity.m.
base_sh = { ...
    'n', 100, 'f', 0.6, 'fs', fs, 'level_of_chaos', 3, ...
    'ode_solver', 'rk4', 'S_c', 0.35, 'std_zero_floor', true, ...
    'n_a_I', 0, 'n_b_I', 0, 'c_E', 0.15/3, ...
    'tau_a_E', [0.1, 1.0, 10], 'tau_b_E_rec', 1.0, 'tau_b_E_rel', 0.25, ...
    'input_type', 'sample_hold', 'T_hold', T_hold, ...
    'T_wash', 6*fs, 'T_train', 30*fs, 'T_test', 15*fs, 'd_max', d_max, ...
    'rng_seeds', [42, 43]};

% Two conditions to also exercise verify_shared_build.
cond_args = {{'n_a_E', 0, 'n_b_E', 0}, {'n_a_E', 3, 'n_b_E', 1}};
esn = cell(1, numel(cond_args));
for i = 1:numel(cond_args)
    esn{i} = SRNN_ESN_reservoir(base_sh{:}, cond_args{i}{:});
    esn{i}.build();
end

verify_shared_build(esn, {'n_a_E', 'n_b_E'}, {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});
fprintf('[PASS] verify_shared_build (u_scalar identical across conditions)\n');

[MC, R2, res] = esn{1}.run_memory_capacity('verbose', false);

n_fail = n_fail + report(res.hold_len == hold_len_expect, ...
    sprintf('hold_len = %d (expect %d)', res.hold_len, hold_len_expect));
n_fail = n_fail + report(numel(R2) == d_max_eff_expect, ...
    sprintf('numel(R2) = %d (expect floor(%d/%d) = %d)', numel(R2), d_max, hold_len_expect, d_max_eff_expect));
n_fail = n_fail + report(numel(res.delay_s) == d_max_eff_expect, ...
    sprintf('numel(delay_s) = %d (expect %d)', numel(res.delay_s), d_max_eff_expect));
n_fail = n_fail + report(abs(res.delay_s(1) - hold_len_expect/fs) < 1e-9, ...
    sprintf('delay_s(1) = %.4f s (expect %.4f)', res.delay_s(1), hold_len_expect/fs));
n_fail = n_fail + report(abs(res.delay_s(end) - d_max_eff_expect*hold_len_expect/fs) < 1e-9, ...
    sprintf('delay_s(end) = %.4f s (expect %.4f)', res.delay_s(end), d_max_eff_expect*hold_len_expect/fs));
n_fail = n_fail + report(isfinite(MC) && all(isfinite(R2)), ...
    sprintf('MC = %.4f finite, R2 all finite', MC));
n_fail = n_fail + report(numel(esn{1}.u_scalar) == (6+30+15)*fs, ...
    'u_scalar length = T_total (full staircase)');

%% ---- Test 2: white noise backward-compat ----
fprintf('\n=== Test 2: white (backward-compat) ===\n');
d_max_w = 25;                  % samples
esn_w = SRNN_ESN_reservoir('n', 100, 'f', 0.6, 'fs', fs, 'level_of_chaos', 3, ...
    'ode_solver', 'rk4', 'S_c', 0.35, 'std_zero_floor', true, ...
    'input_type', 'white', ...
    'T_wash', 2*fs, 'T_train', 15*fs, 'T_test', 15*fs, 'd_max', d_max_w, ...
    'rng_seeds', [42, 43], 'n_a_E', 3, 'n_b_E', 1);
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

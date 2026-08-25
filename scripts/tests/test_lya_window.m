% test_lya_window.m
% Verify that the Lyapunov routines honour lya_T_interval and lya_warmup.
%
% Before this was fixed, SRNNModel2's Benettin used lya_T_interval only to trim
% its last sample and gated accumulation on t >= 0, so the documented "skip the
% first 15 seconds" had no effect; SRNNCellTypePairs honoured the interval, so
% the two classes measured different windows on the same nominal config. Both
% QR paths ignored it as well. The claims checked here are:
%
%   1. Accumulation starts at max(0, lya_T_interval(1)), not at t = 0.
%   2. Iteration starts lya_warmup seconds earlier, so the perturbation (or the
%      QR basis) has aligned by the time accumulation begins.
%   3. A warmup that would reach before T_range(1) is CLAMPED to the available
%      data and warned about, not dropped.
%   4. Nothing is accumulated past lya_T_interval(2).
%   5. Both model classes now agree on which samples count.
%   6. The exponent is genuinely warmup-sensitive, which is why lya_warmup
%      exists and why the default is not smaller.
%
% See also: SRNNModel2, SRNNCellTypePairs, test_benettin_vs_qr

fprintf('=== Testing the Lyapunov accumulation window and warmup ===\n\n');
all_passed = true;

%% Benettin accumulates only inside lya_T_interval
m = run_srnn2('T_range', [0 6], 'lya_T_interval', [3 6], 'lya_warmup', 2);
r = m.lya_results;

all_passed = check('iteration starts lya_warmup before the window', ...
    abs(r.t_lya(1) - 1) < 1e-9) && all_passed;
all_passed = check('accumulation starts at lya_T_interval(1)', ...
    abs(r.t_lya(find(~isnan(r.finite_lya), 1)) - 3) < 1e-9) && all_passed;
all_passed = check('every warmup sample is excluded from the estimate', ...
    all(isnan(r.finite_lya(r.t_lya < 3 - 1e-9)))) && all_passed;
all_passed = check('every in-window sample is included', ...
    all(~isnan(r.finite_lya(r.t_lya >= 3 - 1e-9)))) && all_passed;
all_passed = check('the warmup spans exactly lya_warmup/lya_dt samples', ...
    sum(isnan(r.finite_lya)) == round(2 / r.lya_dt)) && all_passed;
all_passed = check('no segment runs past lya_T_interval(2)', ...
    r.t_lya(end) + r.lya_dt <= 6 + 1e-9) && all_passed;

%% ...and the estimate really is the mean over that window
% LLE is the running mean of log(d_k/d0)/tau over the accumulated samples, so
% it must equal the plain mean of local_lya there.
in_window = r.t_lya >= 3 - 1e-9;
all_passed = check('LLE equals the mean local exponent over the window', ...
    abs(r.LLE - mean(r.local_lya(in_window))) < 1e-9) && all_passed;

%% A window starting before t = 0 still never counts the settling pre-window
m_neg = run_srnn2('T_range', [-4 6], 'lya_T_interval', [-2 6], 'lya_warmup', 1);
r_neg = m_neg.lya_results;
all_passed = check('a negative lya_T_interval(1) clamps accumulation to t = 0', ...
    abs(r_neg.t_lya(find(~isnan(r_neg.finite_lya), 1))) < 1e-9) && all_passed;

%% The warmup is clamped, not dropped, when the run starts too late
[wid, m_clamp] = capture_warning(@() run_srnn2( ...
    'T_range', [0 6], 'lya_T_interval', [1 6], 'lya_warmup', 5));
all_passed = check('an over-long warmup warns', ...
    strcmp(wid, 'SRNNModel:LyapunovWarmupClamped')) && all_passed;
all_passed = check('an over-long warmup clamps to T_range(1) rather than to zero', ...
    abs(m_clamp.lya_results.t_lya(1) - 0) < 1e-9) && all_passed;
all_passed = check('the clamped run still accumulates from the window start', ...
    abs(m_clamp.lya_results.t_lya(find(~isnan(m_clamp.lya_results.finite_lya), 1)) - 1) < 1e-9) ...
    && all_passed;

%% QR honours the same window
q = SRNNModel2('n', 8, 'indegree', 4, 'n_a_E', 0, 'n_a_I', 0, ...
    'n_b_E', 0, 'n_b_I', 0, 'T_range', [0 4], 'fs', 200, ...
    'lya_method', 'qr', 'lya_T_interval', [2 4], 'lya_warmup', 2, ...
    'store_full_state', true);
q.build();
evalc('q.run();');
rq = q.lya_results;
all_passed = check('QR iteration starts lya_warmup before the window', ...
    abs(rq.t_lya(1) - 0) < 1e-9) && all_passed;
all_passed = check('QR accumulation starts at lya_T_interval(1)', ...
    abs(rq.t_lya(find(~isnan(rq.finite_LE_spectrum_t(:, 1)), 1)) - 2) < 1e-9) && all_passed;

%% Both model classes agree on which samples count
p = run_pairs('T_range', [0 6], 'lya_T_interval', [3 6], 'lya_warmup', 2);
rp = p.lya_results;
all_passed = check('SRNNCellTypePairs starts iterating where SRNNModel2 does', ...
    abs(rp.t_lya(1) - r.t_lya(1)) < 1e-9) && all_passed;
all_passed = check('SRNNCellTypePairs starts accumulating where SRNNModel2 does', ...
    abs(rp.t_lya(find(~isnan(rp.finite_lya), 1)) - ...
        r.t_lya(find(~isnan(r.finite_lya), 1))) < 1e-9) && all_passed;
all_passed = check('both classes accumulate the same number of samples', ...
    sum(~isnan(rp.finite_lya)) == sum(~isnan(r.finite_lya))) && all_passed;

%% The exponent is warmup-sensitive, which is the point of the property
% Same accumulation window, different alignment time: a short warmup leaves the
% perturbation off the leading direction and biases the exponent. This is why
% the default is 5 s and why test_benettin_vs_qr raises it to 10.
lle_short = run_srnn2('T_range', [0 12], 'lya_T_interval', [8 12], ...
    'lya_warmup', 0.25).lya_results.LLE;
lle_long = run_srnn2('T_range', [0 12], 'lya_T_interval', [8 12], ...
    'lya_warmup', 8).lya_results.LLE;
fprintf('  (LLE with 0.25 s of warmup: %+.4f; with 8 s: %+.4f)\n', lle_short, lle_long);
all_passed = check('a too-short warmup measurably biases the exponent', ...
    abs(lle_short - lle_long) > 1e-4) && all_passed;

%% A window too narrow to hold a segment says so instead of reporting LLE = 0
% Segments still exist (they are laid out from the warmup start), but none
% lands inside the window, so nothing is accumulated. Returning a bare 0 there
% would read as "edge of chaos" rather than "not measured".
[wid_empty, m_empty] = capture_warning(@() run_srnn2( ...
    'T_range', [0 6], 'lya_T_interval', [5.995 6]));
all_passed = check('an unaccumulable window warns rather than reporting 0 silently', ...
    strcmp(wid_empty, 'SRNNModel:LyapunovWindowEmpty')) && all_passed;
all_passed = check('...and still returns a finite LLE for downstream code', ...
    isfinite(m_empty.lya_results.LLE)) && all_passed;

%% The segment grid itself, exercised directly
% lyapunov_sample_grid is what both methods and both classes share, so pin its
% behaviour without paying for a simulation.
dt_g = 0.005; deci_g = 4; tau_g = dt_g * deci_g;      % lya_dt = 0.02
t_g = (0:dt_g:6)';

[idx_g, tl_g, acc_g] = SRNN2TestAccess.lya_grid(t_g, dt_g, deci_g, tau_g, [3 6], 2);
all_passed = check('grid: accumulation start is max(0, interval(1))', ...
    acc_g == 3) && all_passed;
all_passed = check('grid: indices and times refer to the same samples', ...
    max(abs(t_g(idx_g) - tl_g)) == 0) && all_passed;
all_passed = check('grid: every segment fits inside the trajectory', ...
    all(idx_g + deci_g <= numel(t_g))) && all_passed;
all_passed = check('grid: no segment extends past interval(2)', ...
    all(tl_g + tau_g <= 6 + 1e-12)) && all_passed;

% A trajectory shorter than a single segment leaves nothing to iterate over.
[threw, err] = capture(@() SRNN2TestAccess.lya_grid( ...
    (0:dt_g:2*dt_g)', dt_g, deci_g, tau_g, [0 2*dt_g], 0));
all_passed = check('grid: a trajectory shorter than one segment errors and says so', ...
    threw && contains(err.identifier, 'NoLyapunovIntervals')) && all_passed;

% The Pairs twin must agree sample-for-sample with the SRNNModel2 one.
[idx_p, tl_p, acc_p] = SRNNPairsTestAccess.lya_grid(t_g, dt_g, deci_g, tau_g, [3 6], 2);
all_passed = check('grid: both classes produce an identical grid', ...
    isequal(idx_g, idx_p) && isequal(tl_g, tl_p) && acc_g == acc_p) && all_passed;

%% lya_dt is a property, with per-method defaults when left empty
all_passed = check('an empty lya_dt gives the benettin default of 0.02', ...
    abs(r.lya_dt - 0.02) < 1e-12) && all_passed;
all_passed = check('an empty lya_dt gives the qr default of 0.1', ...
    abs(rq.lya_dt - 0.1) < 1e-12) && all_passed;

m_dt = run_srnn2('T_range', [0 6], 'lya_T_interval', [3 6], 'lya_warmup', 2, ...
    'lya_dt', 0.05);
all_passed = check('an explicit lya_dt is honoured', ...
    abs(m_dt.lya_results.lya_dt - 0.05) < 1e-12) && all_passed;
all_passed = check('...and changes the segment spacing accordingly', ...
    abs(diff(m_dt.lya_results.t_lya(1:2)) - 0.05) < 1e-12) && all_passed;
all_passed = check('...and still reports lya_fs as its reciprocal', ...
    abs(m_dt.lya_results.lya_fs - 20) < 1e-12) && all_passed;

% fs = 200 here, so dt = 0.005: 0.007 is not a multiple, 0.01 is only 2*dt.
[threw_mult, err_mult] = capture(@() run_srnn2('T_range', [0 6], 'lya_dt', 0.007));
all_passed = check('a lya_dt that is not a multiple of dt errors', ...
    threw_mult && contains(err_mult.identifier, 'InvalidLyapunovStep')) && all_passed;
[threw_small, err_small] = capture(@() run_srnn2('T_range', [0 6], 'lya_dt', 0.01));
all_passed = check('a lya_dt below 3*dt errors', ...
    threw_small && contains(err_small.identifier, 'InvalidLyapunovStep')) && all_passed;
[threw_neg, err_neg] = capture(@() run_srnn2('T_range', [0 6], 'lya_dt', -0.02));
all_passed = check('a non-positive lya_dt errors', ...
    threw_neg && contains(err_neg.identifier, 'InvalidLyapunovStep')) && all_passed;

% The Pairs class must resolve and validate lya_dt identically.
p_dt = run_pairs('T_range', [0 6], 'lya_T_interval', [3 6], 'lya_warmup', 2, ...
    'lya_dt', 0.05);
all_passed = check('SRNNCellTypePairs honours lya_dt the same way', ...
    abs(p_dt.lya_results.lya_dt - 0.05) < 1e-12 && ...
    abs(rp.lya_dt - 0.02) < 1e-12) && all_passed;
[threw_p, err_p] = capture(@() run_pairs('T_range', [0 6], 'lya_dt', 0.007));
all_passed = check('SRNNCellTypePairs rejects a bad lya_dt the same way', ...
    threw_p && contains(err_p.identifier, 'InvalidLyapunovStep')) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function m = run_srnn2(varargin)
% A small, fast SRNNModel2 run with Benettin. ode_rk4 keeps it quick; the
% Benettin path uses full time grids, so the fixed-step solver is fine.
defaults = {'n', 24, 'indegree', 12, 'fs', 200, 'ode_solver', 'rk4', ...
    'lya_method', 'benettin', 'store_full_state', true};
m = SRNNModel2(defaults{:}, varargin{:});
m.build();
evalc('m.run();');
end

function p = run_pairs(varargin)
% The SRNNCellTypePairs counterpart of run_srnn2, at the same operating point.
defaults = {'n', 24, 'indegree', 12, 'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], ...
    'tau_a', {srnn_sfa_timescales(3), []}, 'fs', 200, 'ode_solver', 'rk4', ...
    'lya_method', 'benettin', 'store_full_state', true};
p = SRNNCellTypePairs(defaults{:}, varargin{:});
p.build();
evalc('p.run();');
end

function [wid, out] = capture_warning(fcn)
% Run fcn and report the identifier of the last warning it raised.
old = warning('off', 'all');
cleanup = onCleanup(@() warning(old));
lastwarn('', '');
out = fcn();
[~, wid] = lastwarn;
end

function [threw, err] = capture(fcn)
threw = false; err = [];
try
    fcn();
catch err
    threw = true;
end
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end

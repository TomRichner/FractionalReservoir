% test_release_noise.m
% The release noise (FD-027, step N1), below the model class: the gamma
% increments, and the route models' part (the release drive in the drift,
% the kick, its inverse). test_release_noise_model checks the class.
%
%   1. release_increments: moments (tolerance from the increments' own
%      kurtosis, 3 + 6/shape), never negative or non-finite, the exact leak
%      at sigma = 0, any range of steps = the single pass, the schedule's
%      splice, the global stream untouched, one seed's steps not another
%      seed's, the same events at two sigmas.
%   2. Every route variant, with a leak: the release drive, the
%      finite-difference linearization, the steady state; facilitation and
%      the calcium trace do not see the leak; r_s = 0 is the rate itself.
%   3. The kick: it is the exact flow of the RELEASE TERMS OF THE MODEL'S
%      OWN rhs (integrated here by ode45, with the output factor beside
%      it); kick(kick(s, a), b) = kick(s, a + b); `released` recovers the
%      amount from the state after; pools never negative, a large draw
%      included; what does not release does not move; one releasing pool
%      releases tau_rel times what it lost; a route with no pool passes
%      the increment through.
%   4. The realized release at a frozen state against its closed form,
%      mean and variance, where the exact update and the linear law
%      F dR differ by a third.
%   5. What is refused.
%
% See also: release_increments, route_model, route_product, route_pool,
%           test_release_noise_model

fprintf('=== Testing the release noise: increments and route models ===\n\n');
all_passed = true;

%% 1. The increments
fprintf('--- 1. release_increments ---\n');
fs = 200; dt = 1 / fs;
% Two shapes per step: the paper's (about 0.01, most of the mass near
% zero) and a smooth one.
cases = struct('rate', {0.01, 0.2}, 'sigma', {0.007, 0.01}, 'rows', {400, 100}, 'steps', {5000, 2000});
for k = 1:numel(cases)
    c = cases(k);
    stream_before = rng;
    dR = release_increments(c.rate, repmat(c.sigma, c.rows, 1), fs, [1, c.steps], 11);
    stream_after = rng;
    shape = c.rate ^ 2 * dt / c.sigma ^ 2;
    N = numel(dR);
    % Standard errors: of the mean, sd / sqrt(N); of the variance,
    % var * sqrt((kurtosis - 1) / N), kurtosis = 3 + 6 / shape.
    se_mean = c.sigma * sqrt(dt) / sqrt(N);
    se_var = c.sigma ^ 2 * dt * sqrt((2 + 6 / shape) / N);
    off_mean = (mean(dR(:)) - c.rate * dt) / se_mean;
    off_var = (var(dR(:)) - c.sigma ^ 2 * dt) / se_var;
    fprintf('  (shape %.4g per step, N = %.1e: mean off by %+.2f SE, variance by %+.2f SE; %.2f%% exact zeros)\n', ...
        shape, N, off_mean, off_var, 100 * mean(dR(:) == 0));
    all_passed = check(sprintf('shape %.3g: the mean is rate dt (5 SE)', shape), abs(off_mean) < 5) && all_passed;
    all_passed = check(sprintf('shape %.3g: the variance is sigma^2 dt (5 SE)', shape), abs(off_var) < 5) && all_passed;
    all_passed = check(sprintf('shape %.3g: finite and never negative', shape), ...
        all(isfinite(dR(:))) && all(dR(:) >= 0)) && all_passed;
    all_passed = check(sprintf('shape %.3g: the global random stream is untouched', shape), ...
        isequal(stream_before, stream_after)) && all_passed;
end

rate = [0.01; 0.01; 0.02; 0];
sigma = [0.007; 0; 0.01; 0];
whole = release_increments(rate, sigma, fs, [1, 300], 5);
all_passed = check('sigma = 0 is the exact leak, rate dt, and rate 0 is zero', ...
    all(whole(2, :) == rate(2) * dt) && all(whole(4, :) == 0)) && all_passed;
pieces = [release_increments(rate, sigma, fs, [1, 37], 5), ...
    release_increments(rate, sigma, fs, [38, 211], 5), ...
    release_increments(rate, sigma, fs, [212, 300], 5)];
all_passed = check('three ranges of steps = the single pass, bit for bit', isequal(whole, pieces)) && all_passed;
all_passed = check('a range asked for alone = those columns of the whole', ...
    isequal(release_increments(rate, sigma, fs, [120, 140], 5), whole(:, 120:140))) && all_passed;
all_passed = check('an empty range gives no columns', ...
    isequal(size(release_increments(rate, sigma, fs, [5, 4], 5)), [4 0])) && all_passed;
alone = release_increments(rate(1), sigma(1), fs, [1, 300], 5);
all_passed = check('a neuron''s increments do not depend on the other rows', ...
    isequal(whole(1, :), alone)) && all_passed;

% Seeds: consecutive seeds are unrelated streams, not shifted copies.
next = release_increments(rate, sigma, fs, [1, 300], 6);
all_passed = check('another seed gives other increments', ~isequal(next(1, :), whole(1, :))) && all_passed;
all_passed = check('seed + 1 is not the same stream one step on', ...
    ~isequal(next(1, 1:299), whole(1, 2:300)) && ~isequal(next(1, 2:300), whole(1, 1:299))) && all_passed;

% The schedule: spliced at the absolute step, each seed's own stream.
t0 = -1; n_steps = 300;
schedule = [0, 6; 0.25, 5];                 % steps 201.. from seed 6, 251.. from seed 5 again
spliced = release_increments(rate, sigma, fs, [1, n_steps], 5, schedule, t0, n_steps, 'T');
all_passed = check('a schedule is bit-identical before its switch', ...
    isequal(spliced(:, 1:200), whole(:, 1:200))) && all_passed;
all_passed = check('...holds the scheduled seed''s own columns after it', ...
    isequal(spliced(:, 201:250), next(:, 201:250))) && all_passed;
all_passed = check('...and a row naming the base seed returns to the base path', ...
    isequal(spliced(:, 251:300), whole(:, 251:300))) && all_passed;
all_passed = check('a windowed call with a schedule = the window of the whole', ...
    isequal(release_increments(rate, sigma, fs, [190, 260], 5, schedule, t0, n_steps, 'T'), ...
    spliced(:, 190:260))) && all_passed;
[threw, err] = capture(@() release_increments(rate, sigma, fs, [1, n_steps], 5, [0.0012, 6], t0, n_steps, 'T'));
all_passed = check('a switch off the time grid errors', ...
    threw && strcmp(err.identifier, 'T:InvalidNoiseSchedule')) && all_passed;

% One uniform per (neuron, step) whatever sigma is: the same events.
small = release_increments(0.01, 0.003, fs, [1, 2000], 5);
large = release_increments(0.01, 0.012, fs, [1, 2000], 5);
[~, order_small] = sort(small, 'descend');
[~, order_large] = sort(large, 'descend');
all_passed = check('two sigmas put their twenty largest increments on the same steps', ...
    isequal(order_small(1:20), order_large(1:20))) && all_passed;

%% 2. The route models with a leak
fprintf('\n--- 2. The release drive in every route variant ---\n');
variants = route_variants();
npre = 7;
r_s = 0.06;
rng(20261114, 'twister');
for v = 1:numel(variants)
    V = variants(v);
    R1 = compile(V, npre, r_s);
    R0 = compile(V, npre, 0);
    model = R1.model; P = R1.P; ns = P.ns;
    s = random_state(P, npre);
    r = 0.1 + 0.8 * rand(npre, 1);

    [ds, F, drive] = model.rhs(s, r, P);
    [ds0, F0, drive0] = model.rhs(s, r, R0.P);
    [A, B, F_lin, Fs, drive_lin] = model.linearize(s, r, P);
    ok = isequal(drive, r + r_s) && isequal(drive_lin, drive) && ...
        isequal(model.release_drive(r, P), drive) && isequal(drive0, r) && ...
        isequal(model.release_drive([r, 2 * r], P), [r, 2 * r] + r_s);
    all_passed = check(sprintf('%-22s the release drive is r + r_s, and r itself at r_s = 0', V.name), ok) && all_passed;
    all_passed = check(sprintf('%-22s the output factor does not depend on the leak', V.name), ...
        isequal(F, F0) && isequal(F_lin, F)) && all_passed;

    % The leak is in the depletion of the pools and nowhere else: the
    % g-block's derivative is what it is without a leak, bit for bit, and
    % the b-block's is the derivative at the rate r + r_s where the model
    % has no facilitation (the pools then see nothing but the drive).
    g_cols = P.n_b + 1:ns;
    ok = isequal(ds(:, g_cols), ds0(:, g_cols));
    if P.n_g == 0 && ns > 0
        shifted = model.rhs(s, r + r_s, R0.P);
        ok = ok && max(abs(ds - shifted), [], 'all') < 1e-14;
        ok = ok && any(ds(:) ~= ds0(:));
    end
    all_passed = check(sprintf('%-22s the leak reaches the pools only', V.name), ok) && all_passed;

    % The linearization, by central differences of rhs.
    if ns > 0
        h = 1e-6;
        A_fd = zeros(npre, ns, ns);
        Fs_fd = zeros(npre, ns);
        for k = 1:ns
            dk = zeros(npre, ns); dk(:, k) = h;
            [up, F_up] = model.rhs(s + dk, r, P);
            [dn, F_dn] = model.rhs(s - dk, r, P);
            A_fd(:, :, k) = (up - dn) / (2 * h);
            Fs_fd(:, k) = (F_up - F_dn) / (2 * h);
        end
        B_fd = (model.rhs(s, r + h, P) - model.rhs(s, r - h, P)) / (2 * h);
        worst = max([max(abs(A - A_fd), [], 'all'), max(abs(B - B_fd), [], 'all'), ...
            max(abs(Fs - Fs_fd), [], 'all')]);
        all_passed = check(sprintf('%-22s linearize = finite differences of rhs (%.1e)', V.name, worst), ...
            worst < 1e-7) && all_passed;

        ss = model.steady_state(r, P);
        residual = max(abs(model.rhs(ss.s, r, P)), [], 'all');
        all_passed = check(sprintf('%-22s steady_state zeroes rhs with the leak (%.1e)', V.name, residual), ...
            residual < 1e-13) && all_passed;
    end
end

%% 3. The kick
fprintf('\n--- 3. The kick ---\n');
for v = 1:numel(variants)
    V = variants(v);
    R1 = compile(V, npre, 1);               % release terms: rhs at leak 1 minus rhs at leak 0
    R0 = compile(V, npre, 0);
    model = R0.model; P = R0.P; ns = P.ns;
    s = random_state(P, npre);
    r = rand(npre, 1);
    dRc = [0.3; 0.02; 1.5; -0.004; 0; 1e-9; 0.08];

    [s_new, dTheta] = model.kick(s, dRc, P);

    % (a) The exact flow of the model's own release terms, with the output
    % factor integrated beside it. Per neuron, over a = 0 .. dRc.
    worst = 0;
    for i = 1:npre
        if dRc(i) == 0
            worst = max([worst, abs(dTheta(i)), max(abs(s_new(i, :) - s(i, :)))]);
            continue;
        end
        flow = @(~, y) release_flow(y, i, s, r, model, R1.P, R0.P);
        [~, y] = ode45(flow, [0, dRc(i) / 2, dRc(i)], [s(i, :), 0]', ...
            odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
        worst = max([worst, abs(y(end, 1:ns) - s_new(i, :)), abs(y(end, end) - dTheta(i))]);
    end
    all_passed = check(sprintf('%-22s kick = the flow of rhs''s release terms (%.1e)', V.name, worst), ...
        worst < 1e-10) && all_passed;

    % (b) Two kicks = one kick of the sum.
    a = abs(dRc); b = flipud(abs(dRc)) + 0.01;
    [s_a, th_a] = model.kick(s, a, P);
    [s_ab, th_b] = model.kick(s_a, b, P);
    [s_sum, th_sum] = model.kick(s, a + b, P);
    worst = max([max(abs(s_ab - s_sum), [], 'all'), max(abs(th_a + th_b - th_sum))]);
    all_passed = check(sprintf('%-22s kick(kick(s, a), b) = kick(s, a + b) (%.1e)', V.name, worst), ...
        worst < 1e-13) && all_passed;

    % (c) The inverse.
    back = model.released(s_new, dRc, P);
    worst = max(abs(back - dTheta) ./ max(abs(dTheta), 1e-300));
    all_passed = check(sprintf('%-22s released(state after, dRc) = the kick''s amount (%.1e rel)', V.name, worst), ...
        worst < 1e-12) && all_passed;

    % (d) Signs and bounds, with a draw far above anything likely.
    huge = [50; 5; 500; 0.5; 1e3; 20; 80];
    [s_huge, th_huge] = model.kick(s, huge, P);
    pools = 1:P.n_b;
    ok = all(th_huge > 0) && all(isfinite(th_huge)) && all(s_huge(:, pools) >= 0, 'all') && ...
        all(isfinite(s_huge(:))) && all(dTheta(dRc > 0) > 0) && all(dTheta(dRc < 0) < 0) && ...
        all(s_new(:, pools) > 0, 'all');
    all_passed = check(sprintf('%-22s release positive for a positive increment; pools never negative', V.name), ok) && all_passed;

    % (e) What does not release does not move.
    still = setdiff(1:ns, V.releasing);
    all_passed = check(sprintf('%-22s only the releasing pools move', V.name), ...
        isequal(s_new(:, still), s(:, still)) && isequal(s_huge(:, still), s(:, still)) && ...
        (isempty(V.releasing) || all(s_new(dRc ~= 0, V.releasing) ~= s(dRc ~= 0, V.releasing), 'all'))) && all_passed;

    % (f) One releasing pool: what is released is tau_rel times what the
    % pool lost (times the output-only facilitation factor, where the
    % product route has one).
    if isscalar(V.releasing)
        if strcmp(P.model, 'pool')
            tau_rel = P.tau_rel; extra = 1;
        else
            tau_rel = P.tau_drive(:, 1); extra = prod(s(:, P.n_b + 1:end), 2);
        end
        lost = extra .* tau_rel .* (s(:, 1) - s_new(:, 1));
        sizable = abs(dRc) > 1e-6;          % the loss is a difference: not below round-off
        worst = max(abs(lost(sizable) - dTheta(sizable)) ./ abs(dTheta(sizable)));
        all_passed = check(sprintf('%-22s release = tau_rel x the pool''s loss (%.1e rel)', V.name, worst), ...
            worst < 1e-10) && all_passed;
    end
    % (g) No pool: the increment times the output factor, exactly.
    if isempty(V.releasing)
        all_passed = check(sprintf('%-22s no pool: dTheta = F dRc exactly', V.name), ...
            isequal(dTheta, model.output(s, P) .* dRc)) && all_passed;
    end
end
passthrough = compile(variants(1), npre, 0.01);
[s_pt, th_pt] = passthrough.model.kick(zeros(npre, 0), dRc, passthrough.P);
all_passed = check('the pass-through route (no synapse dynamics): dTheta = dRc, no state', ...
    isequal(th_pt, dRc) && isequal(size(s_pt), [npre 0])) && all_passed;

%% 4. The realized release at a frozen state, against the closed form
fprintf('\n--- 4. The realized release: moments ---\n');
% One pool, b1 held at b: a step with no rate releases
%   b g r_s dt  (the leak, in the drift)  +  dTheta(dR - r_s dt)  (the kick)
% with dTheta = -tau b expm1(-g (dR - r_s dt)/tau). For dR ~ Gamma(shape a,
% scale th), E exp(-c dR) = (1 + c th)^-a, so with q = g/tau
%   mean     = b [ g r_s dt + tau (1 - e^{q r_s dt} (1 + q th)^-a) ]
%   variance = (b tau e^{q r_s dt})^2 [ (1 + 2 q th)^-a - (1 + q th)^-2a ]
% The linear law F dR would give b g r_s dt and (b g sigma)^2 dt. Chosen so
% that the two differ by a third (g sigma^2 / (r_s tau_rel) = 1.5).
rate = 0.05; sig = 0.05; tau_rel = 0.05;
b = 0.8; g = 1.5;                           % the frozen pool and release probability
rows = 500; steps = 8000;
cfg = struct('std', struct('tau_rec', 2, 'tau_rel', tau_rel), ...
    'stf', struct('model', 'chain', 'tau_dec', 1, 'tau_fac', 0.5, 'G', 2));
route = SRNNCellTypePairs.compile_route(cfg, struct('npre', rows, 'release_rate', rate));
s = repmat([b, g], rows, 1);
dR = release_increments(rate, repmat(sig, rows, 1), fs, [1, steps], 3);
release = zeros(rows, steps);
for j = 1:steps
    [~, dTheta] = route.model.kick(s, dR(:, j) - rate * dt, route.P);
    release(:, j) = b * g * rate * dt + dTheta;
end
shape = rate ^ 2 * dt / sig ^ 2; th = sig ^ 2 / rate; q = g / tau_rel;
mean_exact = b * (g * rate * dt + tau_rel * (1 - exp(q * rate * dt) * (1 + q * th) ^ -shape));
var_exact = (b * tau_rel * exp(q * rate * dt)) ^ 2 * ((1 + 2 * q * th) ^ -shape - (1 + q * th) ^ (-2 * shape));
mean_linear = b * g * rate * dt;
var_linear = (b * g * sig) ^ 2 * dt;
N = numel(release);
kurt = mean((release(:) - mean(release(:))) .^ 4) / var(release(:)) ^ 2;
se_mean = sqrt(var_exact / N);
se_var = var_exact * sqrt((kurt - 1) / N);
off_mean = (mean(release(:)) - mean_exact) / se_mean;
off_var = (var(release(:)) - var_exact) / se_var;
fprintf('  (mean: realized %.4e, closed form %.4e, linear law %.4e; %+.2f SE from the closed form)\n', ...
    mean(release(:)), mean_exact, mean_linear, off_mean);
fprintf('  (variance: realized %.4e, closed form %.4e, linear law %.4e; %+.2f SE)\n', ...
    var(release(:)), var_exact, var_linear, off_var);
all_passed = check('the mean realized release is the closed form''s (5 SE)', abs(off_mean) < 5) && all_passed;
all_passed = check('...and not the linear law''s, which is 20 SE or more away', ...
    abs(mean(release(:)) - mean_linear) / se_mean > 20) && all_passed;
all_passed = check('the variance of the realized release is the closed form''s (5 SE)', abs(off_var) < 5) && all_passed;
% The mean rate as dt -> 0 is b tau (r_s^2/sigma^2) log(1 + g sigma^2/(r_s tau)):
% below the linear law by about g sigma^2 / (2 r_s tau) where that is small.
ratio = @(e) log1p(e) ./ e;
e_paper = 0.007 ^ 2 / (0.01 * 0.25);
all_passed = check('at the paper''s values the exact update lowers the mean release by under 1%', ...
    abs((1 - ratio(e_paper)) - e_paper / 2) < e_paper ^ 2 && 1 - ratio(e_paper) < 0.01) && all_passed;
fprintf('  (paper''s values, sigma 0.007, r_s 0.01, tau_rel 0.25: lower by %.3f%%; g sigma^2/(2 r_s tau_rel) = %.3f%%)\n', ...
    100 * (1 - ratio(e_paper)), 100 * e_paper / 2);

%% 5. What is refused
fprintf('\n--- 5. Refusals ---\n');
[threw, err] = capture(@() release_increments(0, 0.01, fs, [1, 10], 1));
all_passed = check('sigma > 0 at rate 0 errors', ...
    threw && strcmp(err.identifier, 'release_increments:InvalidParams')) && all_passed;
[threw, err] = capture(@() release_increments(0.01, -0.01, fs, [1, 10], 1));
all_passed = check('a negative sigma errors', ...
    threw && strcmp(err.identifier, 'release_increments:InvalidParams')) && all_passed;
[threw, err] = capture(@() release_increments(0.01, 0.01, fs, [0, 10], 1));
all_passed = check('a step before the first errors', ...
    threw && strcmp(err.identifier, 'release_increments:InvalidParams')) && all_passed;
floor_cfg = struct('std', struct('tau_rec', [2 5], 'tau_rel', 0.25), 'stf', []);
[threw, err] = capture(@() SRNNCellTypePairs.compile_route(floor_cfg, ...
    struct('npre', 3, 'zero_floor', true, 'release_rate', 0.01, 'err_id_prefix', 'T', 'label', 'E->E')));
all_passed = check('the zero floor with a release rate errors, naming both', ...
    threw && strcmp(err.identifier, 'T:ReleaseZeroFloor') && contains(err.message, 'std_zero_floor') && ...
    contains(err.message, 'E->E')) && all_passed;
[threw, ~] = capture(@() SRNNCellTypePairs.compile_route(floor_cfg, ...
    struct('npre', 3, 'zero_floor', true, 'release_rate', 0)));
all_passed = check('...and without one it compiles as before', ~threw) && all_passed;
[threw, ~] = capture(@() SRNNCellTypePairs.compile_route(struct('std', [], 'stf', []), ...
    struct('npre', 3, 'zero_floor', true, 'release_rate', 0.01)));
all_passed = check('...as does a route without depression, where the floor is not in force', ~threw) && all_passed;
for bad = {-0.01, [0.01 0.02], NaN}
    [threw, ~] = capture(@() SRNNCellTypePairs.compile_route(floor_cfg, ...
        struct('npre', 3, 'release_rate', bad{1})));
    all_passed = check(sprintf('release_rate = %s on a route errors', mat2str(bad{1})), threw) && all_passed;
end

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function variants = route_variants()
% Every kind of route: its configuration, the spread it is compiled with,
% and which columns of its state a kick moves (the releasing pools).
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
coupled = struct('model', 'chain', 'tau_dec', 1, 'tau_fac', 0.5, 'G', 2);
with_trace = coupled; with_trace.tau_h = 0.3;
output_only = struct('tau_dec', [1 3], 'tau_fac', 0.5, 'G', [2 1.5]);
variants = struct('name', {}, 'std', {}, 'stf', {}, 'spread', {}, 'releasing', {});
variants(end + 1) = variant('no synapse dynamics', [], [], 0, []);
variants(end + 1) = variant('one pool (product)', struct('tau_rec', 2, 'tau_rel', 0.25), [], 0.25, 1);
variants(end + 1) = variant('product of two', struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4]), [], 0.25, [1 2]);
variants(end + 1) = variant('geometric mean of two', ...
    struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4], 'combine', 'geomean'), [], 0, [1 2]);
variants(end + 1) = variant('product of three', struct('tau_rec', [2 5 9], 'tau_rel', [0.25 0.4 0.3]), [], 0, [1 2 3]);
variants(end + 1) = variant('one pool x STF', struct('tau_rec', 2, 'tau_rel', 0.25), output_only, 0, 1);
variants(end + 1) = variant('output-only STF', [], output_only, 0, []);
variants(end + 1) = variant('chain', chain, [], 0.25, 1);
variants(end + 1) = variant('chain, coupled g, h', chain, with_trace, 0.25, 1);
variants(end + 1) = variant('single pool, coupled g', struct('tau_rec', 2, 'tau_rel', 0.25), coupled, 0, 1);
variants(end + 1) = variant('coupled g, h, no pool', [], with_trace, 0, []);
end

function V = variant(name, std_cfg, stf_cfg, spread, releasing)
V = struct('name', name, 'std', std_cfg, 'stf', stf_cfg, 'spread', spread, 'releasing', releasing);
end

function route = compile(V, npre, release_rate)
% The route as SRNNCellTypePairs compiles it (compile_route picks the
% model), with the variant's per-neuron spread and the given leak. The
% draws are fixed, so two compiles of one variant differ in the leak only.
z = cos((1:npre)' .* (1:10));
route = SRNNCellTypePairs.compile_route(struct('std', V.std, 'stf', V.stf), ...
    struct('npre', npre, 'label', V.name, 'release_rate', release_rate, ...
    'std_sigma', V.spread, 'std_z', z, 'stf_sigma', V.spread, 'stf_z', z));
end

function s = random_state(P, npre)
% A state inside every variable's range, away from rest: pools in (0, 1),
% facilitation above 1, a calcium trace in (0, 1).
s = 0.2 + 0.7 * rand(npre, P.ns);
if strcmp(P.model, 'pool')
    if P.has_g, s(:, P.n_b + 1) = 1 + 0.8 * rand(npre, 1); end
else
    s(:, P.n_b + 1:end) = 1 + 0.5 * rand(npre, P.n_g);
end
end

function dy = release_flow(y, i, s, r, model, P_one, P_zero)
% d[s; Theta]/da along a kick of neuron i: the release terms of the
% model's rhs (rhs is affine in the leak, so they are rhs at leak 1 minus
% rhs at leak 0, whatever the rate) and the output factor.
ns = size(s, 2);
state = s;
state(i, :) = y(1:ns)';
[with_leak, F] = model.rhs(state, r, P_one);
without = model.rhs(state, r, P_zero);
dy = [with_leak(i, :)' - without(i, :)'; F(i)];
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

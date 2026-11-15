% test_release_noise_model.m
% The release noise in SRNNCellTypePairs (FD-027, step N1). The increments
% and the route models' part are test_release_noise's; this file checks
% what the class does with them, on the four conditions of the
% chained-pool network (no adaptation: pass-through routes; one pool; the
% chain; the product of two) at n = 40.
%
%   1. Off means identical: one small run of each condition against the
%      fixture made from the code before the release noise existed
%      (release_noise_off_run), bit for bit, defaults and explicit zeros.
%   2. The properties: broadcasting, the scalar aliases, validation, what
%      is refused (a deterministic integrator, sigma without a rate, the
%      zero floor), and that analysis_run_config calls a preset with
%      release noise stochastic.
%   3. The switches. A leak-only run (sigma 0) is the mean dynamics: no
%      draw, no kick, any integrator. A release-only run (sigma_u_noise 0)
%      has kicks. Both sources together.
%   4. The Jacobian of the mean dynamics with the leak on: analytic against
%      finite differences, and jacobian_times against the assembled matrix.
%   5. The kick in the class: only the releasing pools and x move; x moves
%      by W dTheta / tau_d with the matrix the drift uses (gain, route
%      scale, self-coupling); pools never negative over a run, a large
%      draw included; two routes of one presynaptic type stay identical.
%   6. Increments by absolute step: a run in segments = the single pass;
%      a schedule is bit-identical up to its switch and differs after;
%      Benettin re-integrates with the same kicks.
%   7. Strong order 1 in the step (coarse increments made by summing fine
%      ones).
%
% See also: test_release_noise, release_noise_off_run, SRNNCellTypePairs

fprintf('=== Testing the release noise in SRNNCellTypePairs ===\n\n');
all_passed = true;
names = condition_names();
release = {'release_rate', 0.01, 'release_sigma', 0.007};

%% 1. Off means identical
fprintf('--- 1. Off = the outputs before the release noise existed ---\n');
fixture = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'release_noise_off_reference.mat');
reference = load(fixture);
defaults = release_noise_off_run();
zeros_named = release_noise_off_run({'release_rate', [0 0], 'release_sigma', 0});
for k = 1:numel(names)
    all_passed = check(sprintf('%-13s the properties left alone: bit for bit', names{k}), ...
        isequal(defaults.(names{k}), reference.(names{k}))) && all_passed;
    all_passed = check(sprintf('%-13s both set to zero: bit for bit', names{k}), ...
        isequal(zeros_named.(names{k}), reference.(names{k}))) && all_passed;
end

%% 2. The properties
fprintf('\n--- 2. Properties, aliases and refusals ---\n');
m = net(3);
all_passed = check('the default is off: zeros, one per cell type', ...
    isequal(m.release_rate, [0 0]) && isequal(m.release_sigma, [0 0])) && all_passed;
m = net(3, 'release_rate', 0.01, 'release_sigma', [0.007 0], 'ode_solver', 'sra1');
all_passed = check('a scalar is broadcast to every type; a row is kept', ...
    isequal(m.release_rate, [0.01 0.01]) && isequal(m.release_sigma, [0.007 0])) && all_passed;
m = net(3, 'release_rate_all', 0.02, 'release_sigma_all', 0.005, 'ode_solver', 'sra1');
all_passed = check('the scalar aliases write every type and read back as scalars', ...
    isequal(m.release_rate, [0.02 0.02]) && isequal(m.release_sigma, [0.005 0.005]) && ...
    m.release_rate_all == 0.02 && m.release_sigma_all == 0.005) && all_passed;
m.release_sigma = [0.005 0.001];
[threw, err] = capture(@() m.release_sigma_all);
all_passed = check('reading an alias when the types differ errors', ...
    threw && strcmp(err.identifier, 'SRNNCellTypePairs:ReleaseMismatch')) && all_passed;
[threw, err] = capture(@() net(3, 'release_sigma_all', [0.1 0.2]));
all_passed = check('an alias takes one number', ...
    threw && strcmp(err.identifier, 'SRNNCellTypePairs:InvalidParams')) && all_passed;

for solver = {'rk4', 'ode45', 'ode15s'}
    [threw, err] = capture(@() net(3, release{:}, 'ode_solver', solver{1}));
    all_passed = check(sprintf('release noise on ''%s'' errors, naming release_sigma', solver{1}), ...
        threw && strcmp(err.identifier, 'SRNNCellTypePairs:InvalidParams') && ...
        contains(err.message, 'release_sigma') && contains(err.message, 'sra1')) && all_passed;
end
m = net(3, release{:}, 'ode_solver', 'sra1');
m.ode_solver = 'rk4';
[threw, ~] = capture(@() m.build());
all_passed = check('...also when the solver is changed after construction', threw) && all_passed;
[threw, ~] = capture(@() net(3, 'release_rate', 0.01, 'ode_solver', 'rk4'));
all_passed = check('a leak alone (sigma 0) runs on a deterministic integrator', ~threw) && all_passed;
[threw, err] = capture(@() net(3, 'release_sigma', 0.007, 'ode_solver', 'sra1'));
all_passed = check('sigma > 0 without a rate errors', ...
    threw && contains(err.message, 'release_rate > 0')) && all_passed;
[threw, ~] = capture(@() net(3, 'release_rate', [0.01 0], 'release_sigma', [0.007 0.007], 'ode_solver', 'sra1'));
all_passed = check('...per cell type', threw) && all_passed;
for bad = {-0.01, [0.01 0.01 0.01], NaN, 'a'}
    [threw, ~] = capture(@() net(3, 'release_rate', bad{1}));
    all_passed = check(sprintf('release_rate = %s errors', mat2str(bad{1})), threw) && all_passed;
end
[threw, err] = capture(@() net(4, 'release_rate', 0.01, 'std_zero_floor', true));
all_passed = check('std_zero_floor with a release rate errors', ...
    threw && strcmp(err.identifier, 'SRNNCellTypePairs:ReleaseZeroFloor')) && all_passed;
[threw, ~] = capture(@() net(1, 'release_rate', 0.01, 'std_zero_floor', true));
all_passed = check('...but not where no route has depression (the floor is not in force)', ~threw) && all_passed;
m = net(4, 'release_rate', 0.01); m.build();
p = m.get_params(); p.u_interpolant = m.u_interpolant; p.std_zero_floor = true;
[threw, err] = capture(@() SRNNCellTypePairs.dynamics_fast(0, m.S0, p));
all_passed = check('...nor can the floor be switched on over the params struct', ...
    threw && strcmp(err.identifier, 'SRNNCellTypePairs:ReleaseZeroFloor')) && all_passed;

% The route's comparable data carries both, so routes of types that
% release differently are not "identical".
m = net(3, 'release_rate', [0.01 0.02], 'release_sigma', [0.007 0], 'ode_solver', 'sra1'); m.build();
p = m.get_params();
sig_EE = SRNNCellTypePairs.route_fields(p, 1, 1);
sig_EI = SRNNCellTypePairs.route_fields(p, 1, 2);
sig_IE = SRNNCellTypePairs.route_fields(p, 2, 1);
all_passed = check('route_fields carries the presynaptic type''s release rate and sigma', ...
    sig_EE.release_rate == 0.01 && sig_EE.release_sigma == 0.007 && ...
    sig_IE.release_rate == 0.02 && sig_IE.release_sigma == 0) && all_passed;
all_passed = check('two routes of one presynaptic type compare equal; routes of two types do not', ...
    isequal(sig_EE, sig_EI) && ~isequal(sig_EE, sig_IE) && ...
    ~SRNNCellTypePairs.routes_identical(p)) && all_passed;

% analysis_run_config: noise is on when either source is.
preset = struct('sigma_u_noise', 0);
quiet = analysis_run_config('sensitivity', 'fast', preset);
preset.release_rate = [0.01 0.01];
leak = analysis_run_config('sensitivity', 'fast', preset);
preset.release_sigma = [0.007 0];
noisy = analysis_run_config('sensitivity', 'fast', preset);
alias = analysis_run_config('sensitivity', 'fast', struct('release_sigma_all', 0.007));
all_passed = check('analysis_run_config: a preset with release noise is stochastic; a leak alone is not', ...
    ~quiet.is_stochastic && ~leak.is_stochastic && noisy.is_stochastic && alias.is_stochastic && ...
    strcmp(noisy.model.ode_solver, 'sra1') && strcmp(leak.model.ode_solver, quiet.model.ode_solver)) && all_passed;

%% 3. The switches
fprintf('\n--- 3. Leak only, release only, both ---\n');
for k = 1:numel(names)
    off = run_net(net(k, 'ode_solver', 'sra1'));
    leak = net(k, 'release_rate', 0.01, 'ode_solver', 'sra1');
    leak.build(); leak.build_noise_for_test();
    no_increments = isempty(leak.noise_increments);
    leak = run_net(leak);
    % The mean dynamics, integrated here with the plain scheme.
    p = leak.get_params(); p.u_interpolant = leak.u_interpolant;
    [~, S_mean] = sde_fixed_step(@(t, S) SRNNCellTypePairs.dynamics_fast(t, S, p), ...
        leak.t_ex, leak.S0, [], [], 'sra1');
    all_passed = check(sprintf('%-13s leak only: no increments, and the run is the mean dynamics', names{k}), ...
        no_increments && isequal(leak.S_out, S_mean)) && all_passed;
    all_passed = check(sprintf('%-13s ...which the leak moves away from the run without it', names{k}), ...
        ~isequal(leak.S_out, off.S_out)) && all_passed;

    only = net(k, release{:}, 'ode_solver', 'sra1');
    only.build(); only.build_noise_for_test();
    nz = only.noise_increments;
    ok = ~isempty(nz) && isempty(nz.xi1) && nz.sigma == 0 && isa(nz.kick, 'function_handle') && ...
        nz.kick_steps == numel(only.t_ex) - 1 && isequal(size(nz.release.dRc), [only.n, nz.kick_steps]);
    all_passed = check(sprintf('%-13s release only (sigma_u_noise 0): the noise struct has the kick and no Wiener tensor', names{k}), ok) && all_passed;
    only = run_net(only);
    x = only.cached_params.state_layout.x;
    all_passed = check(sprintf('%-13s ...and the kicks reach x from the first step', names{k}), ...
        all(only.S_out(2, x) ~= leak.S_out(2, x)) && all(isfinite(only.S_out(:)))) && all_passed;
    all_passed = check(sprintf('%-13s ...and the tensor is cleared after the run', names{k}), ...
        isempty(only.noise_increments)) && all_passed;

    both = net(k, release{:}, 'sigma_u_noise', 0.025, 'ode_solver', 'sra1');
    both.build(); both.build_noise_for_test();
    wiener = net(k, 'release_rate', 0.01, 'sigma_u_noise', 0.025, 'ode_solver', 'sra1');
    wiener.build(); wiener.build_noise_for_test();
    ok = isequal(both.noise_increments.xi1, wiener.noise_increments.xi1) && ...
        isequal(both.noise_increments.xi2, wiener.noise_increments.xi2) && ...
        ~isfield(wiener.noise_increments, 'kick') && ...
        isequal(both.noise_increments.release.dRc, nz.release.dRc);
    all_passed = check(sprintf('%-13s both sources: the Wiener path and the release increments are each what they are alone', names{k}), ok) && all_passed;
end
% One type with release noise, the other with a leak: only the first is kicked.
mixed = net(3, 'release_rate', [0.01 0.02], 'release_sigma', [0.007 0], 'ode_solver', 'sra1');
mixed.build(); mixed.build_noise_for_test();
dRc = mixed.noise_increments.release.dRc;
I = mixed.cached_params.type_indices{2};
all_passed = check('a type with a leak and no sigma has zero centered increments', ...
    all(dRc(I, :) == 0, 'all') && any(dRc(mixed.cached_params.type_indices{1}, :) ~= 0, 'all')) && all_passed;
all_passed = check('the centered increments are never below -r_s dt', ...
    min(dRc, [], 'all') >= -0.01 / mixed.fs) && all_passed;

%% 4. The Jacobian of the mean dynamics, with the leak
fprintf('\n--- 4. Jacobian with the leak ---\n');
jacobian_nets = [num2cell(1:numel(names)), {'facilitated'}];
for k = 1:numel(jacobian_nets)
    % A smooth activation: finite differences across the piecewise
    % activation's kinks are not a derivative.
    if ischar(jacobian_nets{k})
        make = @(leak) facilitated_net('release_rate', leak, 'activation', 'logistic', 'T_range', [0 1]);
        label = 'pools with coupled g, h';
    else
        make = @(leak) net(jacobian_nets{k}, 'release_rate', leak, 'activation', 'logistic', 'T_range', [0 1]);
        label = names{jacobian_nets{k}};
    end
    m = run_net(make([0.05 0.03]));
    p = m.get_params(); p.u_interpolant = m.u_interpolant;
    S = m.S_out(end, :)';
    J = full(SRNNCellTypePairs.compute_Jacobian_fast(S, p));
    J_fd = SRNNCellTypePairs.finite_difference_jacobian(S, p, 1e-6);
    worst = max(abs(J - J_fd), [], 'all') / max(abs(J), [], 'all');
    all_passed = check(sprintf('%-23s analytic = finite differences (%.1e of the largest entry)', label, worst), ...
        worst < 1e-8) && all_passed;
    Y = cos((1:numel(S))' .* [1 2 3]);
    JY = J * Y;
    worst = max(abs(SRNNCellTypePairs.jacobian_times(S, Y, p) - JY), [], 'all') / max(abs(JY), [], 'all');
    all_passed = check(sprintf('%-23s jacobian_times = the assembled Jacobian times Y (%.1e)', label, worst), ...
        worst < 1e-13) && all_passed;
    without = make(0);
    without.build();
    p0 = without.get_params();
    % The leak is in it wherever a route has state (it multiplies dF/ds
    % and empties the pools). Without synapse dynamics it is a constant
    % input, W r_s, and the Jacobian is the one without it.
    J_without = full(SRNNCellTypePairs.compute_Jacobian_fast(S, p0));
    if any(p.n_b_pairs(:) > 0) || any(p.n_g_pairs(:) > 0)
        all_passed = check(sprintf('%-23s ...and the leak is in it (it differs from the Jacobian without)', label), ...
            ~isequal(J, J_without)) && all_passed;
    else
        all_passed = check(sprintf('%-23s ...which without synapse dynamics is the Jacobian without the leak', label), ...
            isequal(J, J_without)) && all_passed;
    end
end

%% 5. The kick in the class
fprintf('\n--- 5. The kick: what moves, and by how much ---\n');
for k = 1:numel(names)
    % Gain, self-coupling and (where there are routes) a route scale, so
    % that the matrix the kick uses is not the drawn one.
    m = net(k, release{:}, 'ode_solver', 'sra1', 'level_of_chaos', 1.3, 'self_coupling', [-0.2 0.1]);
    if k > 1
        m.synapse_config.I.E.scale = 2;
    end
    m = run_net(m);
    p = m.cached_params;
    L = p.state_layout;
    S = m.S_out(end, :)';
    dRc = 0.02 * abs(cos(1:m.n))' - 0.01 / m.fs;
    S_new = SRNNCellTypePairs.release_kick(S, dRc, p);

    % Which states may move: x, and each route's releasing pools (column
    % 1 of a pool route's b-block; every column of a product route's).
    may_move = false(L.N, 1);
    may_move(L.x) = true;
    expected_x = zeros(m.n, 1);
    for pre = 1:2
        for post = 1:2
            R = p.route{pre, post};
            npre = p.n_per_type(pre);
            pools = L.b{pre, post};
            if strcmp(R.P.model, 'pool'), pools = pools(1:npre); end
            may_move(pools) = true;
            % The amount released, recovered from the state after.
            s_after = reshape(S_new(R.idx), npre, R.P.ns);
            dTheta = R.model.released(s_after, dRc(p.type_indices{pre}), R.P);
            expected_x(p.type_indices{post}) = expected_x(p.type_indices{post}) + ...
                p.W(p.type_indices{post}, p.type_indices{pre}) * dTheta;
        end
    end
    moved = S_new ~= S;
    all_passed = check(sprintf('%-13s only x and the releasing pools move (SFA, reserve, the rest: untouched)', names{k}), ...
        ~any(moved & ~may_move) && all(moved(L.x))) && all_passed;
    worst = max(abs(S_new(L.x) - S(L.x) - expected_x / m.tau_d)) / max(abs(expected_x / m.tau_d));
    all_passed = check(sprintf('%-13s x moves by W dTheta / tau_d, dTheta as `released` recovers it (%.1e)', names{k}, worst), ...
        worst < 1e-12) && all_passed;
    % The matrix is the drift's: the drawn one with the gain (in obj.W),
    % the route scale and the self-coupling.
    W_sim = m.W;
    if k > 1
        W_sim(p.type_indices{1}, p.type_indices{2}) = 2 * W_sim(p.type_indices{1}, p.type_indices{2});
    end
    coupling = zeros(m.n, 1);
    coupling(p.type_indices{1}) = -0.2; coupling(p.type_indices{2}) = 0.1;
    W_sim = W_sim + diag(coupling) * m.R_gain1;
    all_passed = check(sprintf('%-13s ...with the matrix the drift uses: gain, route scale, self-coupling', names{k}), ...
        max(abs(full(p.W) - full(W_sim)), [], 'all') < 1e-14 && ~isequal(full(p.W), full(m.W))) && all_passed;
    if k == 1
        all_passed = check(sprintf('%-13s pass-through routes: x moves by W dRc / tau_d', names{k}), ...
            max(abs(S_new(L.x) - S(L.x) - p.W * dRc / m.tau_d)) < 1e-13) && all_passed;
    end

    % Over the run, and after a draw far above anything likely: no pool
    % below zero, nothing non-finite.
    pools_all = [L.b{:}];
    S_huge = SRNNCellTypePairs.release_kick(S, 100 * ones(m.n, 1), p);
    ok = all(isfinite(m.S_out(:))) && all(isfinite(S_huge));
    if ~isempty(pools_all)
        ok = ok && min(m.S_out(:, pools_all), [], 'all') > 0 && min(S_huge(pools_all)) >= 0;
    end
    all_passed = check(sprintf('%-13s pools never negative over the run or after a draw of 100', names{k}), ok) && all_passed;

    % Two routes of one presynaptic type: the same state, every step.
    same = true;
    for pre = 1:2
        same = same && isequal(m.S_out(:, [L.b{pre, 1}, L.g{pre, 1}]), m.S_out(:, [L.b{pre, 2}, L.g{pre, 2}]));
    end
    all_passed = check(sprintf('%-13s the two routes of each presynaptic type stay isequal under the noise', names{k}), same) && all_passed;
end

%% 6. Increments by absolute step
fprintf('\n--- 6. Segments, the schedule, Benettin ---\n');
for k = [1 3 4]
    m = net(k, release{:}, 'sigma_u_noise', 0.025, 'ode_solver', 'sra1');
    m.build(); m.build_noise_for_test();
    p = m.get_params(); p.u_interpolant = m.u_interpolant;
    rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, p);
    solver = resolve_solver('sra1', m.noise_increments);
    [~, whole] = solver(rhs, m.t_ex, m.S0, []);
    cuts = [1, 57, 58 + 3, numel(m.t_ex)];
    pieces = m.S0';
    for c = 1:numel(cuts) - 1
        [~, part] = solver(rhs, m.t_ex(cuts(c):cuts(c + 1)), pieces(end, :)', []);
        pieces = [pieces; part(2:end, :)]; %#ok<AGROW>
    end
    all_passed = check(sprintf('%-13s integrated in three segments = the single pass, bit for bit', names{k}), ...
        isequal(whole, pieces)) && all_passed;
    [threw, err] = capture(@() solver(rhs, m.t_ex(end) + (0:3)' / m.fs, m.S0, []));
    all_passed = check(sprintf('%-13s a span past the last kick errors', names{k}), ...
        threw && strcmp(err.identifier, 'sde_fixed_step:NoiseOutOfRange')) && all_passed;
end

% The schedule, release only and with both sources.
for wiener = [0 0.025]
    base = run_net(net(3, release{:}, 'sigma_u_noise', wiener, 'ode_solver', 'sra1', 'noise_seed', 21));
    switched = run_net(net(3, release{:}, 'sigma_u_noise', wiener, 'ode_solver', 'sra1', 'noise_seed', 21, ...
        'noise_seed_schedule', [1.0, 77]));
    same_seed = run_net(net(3, release{:}, 'sigma_u_noise', wiener, 'ode_solver', 'sra1', 'noise_seed', 21, ...
        'noise_seed_schedule', [1.0, 21]));
    upto = find(abs(base.t_out - 1.0) < 1e-9);
    label = sprintf('sigma_u = %g:', wiener);
    all_passed = check(sprintf('%-15s a schedule is bit-identical up to and including its switch time', label), ...
        isequal(base.S_out(1:upto, :), switched.S_out(1:upto, :))) && all_passed;
    all_passed = check(sprintf('%-15s ...and differs from the next step on', label), ...
        ~isequal(base.S_out(upto + 1, :), switched.S_out(upto + 1, :))) && all_passed;
    all_passed = check(sprintf('%-15s a schedule naming the base seed reproduces the run', label), ...
        isequal(base.S_out, same_seed.S_out)) && all_passed;
end

% Benettin: finite, and reproducible, so the perturbed trajectory got the
% same kicks. The exponent is not the leak-only one (the path differs).
lya = {'lya_method', 'benettin', 'T_range', [0 4], 'lya_T_interval', [2 4], 'lya_warmup', 1, 'ode_solver', 'sra1'};
b1 = run_net(net(3, release{:}, lya{:}));
b2 = run_net(net(3, release{:}, lya{:}));
b0 = run_net(net(3, 'release_rate', 0.01, lya{:}));
fprintf('  (Benettin, chain condition: leak only %+.4f, release noise %+.4f)\n', ...
    b0.lya_results.LLE, b1.lya_results.LLE);
all_passed = check('Benettin returns a finite exponent under release noise', ...
    isfinite(b1.lya_results.LLE)) && all_passed;
all_passed = check('...reproducibly: the perturbed run saw the same increments', ...
    b1.lya_results.LLE == b2.lya_results.LLE && b1.lya_results.LLE ~= b0.lya_results.LLE) && all_passed;
% The perturbed segment gets the kicks: re-integrating a segment of the
% stored run from its own stored state returns the stored end state.
m = net(3, release{:}, 'ode_solver', 'sra1'); m.build(); m.build_noise_for_test();
p = m.get_params(); p.u_interpolant = m.u_interpolant;
solver = resolve_solver('sra1', m.noise_increments);
rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, p);
[~, whole] = solver(rhs, m.t_ex, m.S0, []);
[~, again] = solver(rhs, m.t_ex(101:105), whole(101, :)', []);
all_passed = check('a segment re-integrated from the stored state ends on the stored state', ...
    isequal(again(end, :), whole(105, :))) && all_passed;

%% 7. Strong order in the step
fprintf('\n--- 7. Strong order ---\n');
% The same gamma paths on coarser steps: increments add (the sum of m fine
% gamma increments is the coarse one, exactly), so the difference between a
% coarse run and the fine one is discretization error alone. The kick is
% exact for the release terms; what is first order is its splitting from
% the drift.
% At the values of the paper's figures (r_s 0.01, sigma 0.007), steps
% 1/200 s (the sweeps') to 1/1600 s, the reference 8 times finer still.
% Measured 2026-11-14: order 1.0 with pass-through routes and 0.9 to 1.0
% with the chain (6 and 8 paths; the increments are lumpy, a few large
% events per path, so the fit moves with the paths), and an error at
% 1/200 s of 0.2 to 0.3 % of the state's norm after 0.5 s.
T = 0.5; fs_fine = 12800; refinements = [64 32 16 8];
n_paths = 6;
for k = [1 3]
    m = net(k, release{:}, 'ode_solver', 'sra1', 'T_range', [0 T]);
    m.build();
    p = m.get_params(); p.u_interpolant = m.u_interpolant;
    rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, p);
    rate = 0.01 * ones(m.n, 1);
    err = zeros(numel(refinements), n_paths);
    for path = 1:n_paths
        dR = release_increments(rate, 0.007 * ones(m.n, 1), fs_fine, [1, T * fs_fine], 100 + path);
        S_ref = final_state(rhs, m.S0, dR, rate, fs_fine, T, p);
        for ri = 1:numel(refinements)
            w = refinements(ri);
            coarse = squeeze(sum(reshape(dR, m.n, w, []), 2));
            err(ri, path) = norm(final_state(rhs, m.S0, coarse, rate, fs_fine / w, T, p) - S_ref);
        end
    end
    rms_err = sqrt(mean(err .^ 2, 2));
    h_used = refinements(:) / fs_fine;
    fit = polyfit(log(h_used), log(rms_err), 1);
    fprintf('  (%s: rms error %s at steps 1/%s s, of a state of norm %.1f; order %.2f)\n', names{k}, ...
        mat2str(rms_err', 2), mat2str(fs_fine ./ refinements), norm(S_ref), fit(1));
    all_passed = check(sprintf('%-13s the release kicks are strong order 1 in the step (0.75 to 1.3)', names{k}), ...
        fit(1) > 0.75 && fit(1) < 1.3) && all_passed;
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
function names = condition_names()
[~, ~, conditions] = srnn_param_preset(base_preset());
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
end

function name = base_preset()
% The A-block network, noise-free: the tests switch each source on.
name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
end

function args = net_args(k)
[d, ~, conditions] = srnn_param_preset(base_preset());
d.n = 40;
d.indegree = 8;
d.F_ref_n = 40;
d.F_ref_indegree = 8;
args = [struct2namevalue(d), struct2namevalue(rmfield(conditions{k}, 'name')), ...
    {'T_range', [0 2], 'fs', 200, 'ode_solver', 'rk4', 'lya_method', 'none', ...
    'store_full_state', true, 'rng_seeds', [3 4]}];
end

function m = net(k, varargin)
% Condition k of the chained-pool network at n = 40.
args = net_args(k);
m = SRNNPairsTestAccess(args{:}, varargin{:});
end

function m = facilitated_net(varargin)
% The same network with vesicle-pool routes that facilitate: E routes a
% chain with release-coupled g and a calcium trace, I routes a single pool
% with release-coupled g.
args = net_args(3);
coupled = struct('model', 'chain', 'tau_dec', 1, 'tau_fac', 0.5, 'G', 2);
with_trace = coupled; with_trace.tau_h = 0.3;
E = struct('std', struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2), 'stf', with_trace);
I = struct('std', struct('tau_rec', 2, 'tau_rel', 0.25), 'stf', coupled);
config = struct('E', struct('E', E, 'I', E), 'I', struct('E', I, 'I', I));
m = SRNNPairsTestAccess(args{:}, 'synapse_config', config, varargin{:});
end

function m = run_net(m)
m.build();
evalc('m.run();');
end

function S_end = final_state(rhs, S0, dR, rate, fs, T, params)
% One 'sra1' run over [0, T] at step 1/fs with the given increments (one
% column per step), release noise alone; the state at T.
dRc = dR - rate ./ fs;
noise = struct('sigma', 0, 't0', 0, 'fs', fs, 'kick_steps', size(dRc, 2), ...
    'kick', @(S, j) SRNNCellTypePairs.release_kick(S, dRc(:, j), params));
[~, S] = sde_fixed_step(rhs, (0:round(T * fs))' / fs, S0, [], noise, 'sra1');
S_end = S(end, :)';
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

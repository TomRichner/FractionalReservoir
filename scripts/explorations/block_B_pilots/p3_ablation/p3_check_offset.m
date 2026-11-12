function err = p3_check_offset()
% P3_CHECK_OFFSET Check that a constant input of -abar stands in exactly for a frozen adaptation term abar.
%
%   err = P3_CHECK_OFFSET()
%
% The sign and exactness check of the mean-matched ablation's SFA part, on
% the bursting network with depression removed from both models so that
% only the adaptation term is under test.
%   Model A  SFA intact, no input. State: x, and adaptation variables a
%            chosen at random; abar_i = c_eff * sum_k a_ik.
%   Model B  SFA removed, constant input u_i = -abar_i. State: y = x - abar.
% With the term frozen, A's equation is tau dx/dt = -x + W theta(phi(x -
% abar)), and y = x - abar obeys tau dy/dt = -y - abar + W theta(phi(y)),
% which is B's. So B's dy/dt must equal A's dx/dt. Returns the largest
% absolute difference (round-off), and errors above 1e-12. Flipping the
% sign of the input gives a difference of order 2 abar / tau.
%
% See also: p3_ablation_args

no_std = struct('E', struct('E', struct(), 'I', struct()), 'I', struct('E', struct(), 'I', struct()));
mA = bursting_condition_model('original', 'model_args', {'synapse_config', no_std});
pA = mA.cached_params;
L = pA.state_layout;
rs = RandStream('mt19937ar', 'Seed', 7);
S = zeros(L.N, 1);
S(L.x) = 0.5 + 0.4 * randn(rs, mA.n, 1);
abar = zeros(mA.n, 1);
for q = 1:pA.n_cellTypes
    if pA.n_a(q) == 0, continue; end
    a = 0.6 * rand(rs, pA.n_per_type(q), pA.n_a(q));
    S(L.a{q}) = a(:);
    abar(pA.type_indices{q}) = pA.c_eff(q) .* sum(a, 2);
end
assert(any(abar > 0.01), 'The test needs a nonzero adaptation term.');
dA = mA.dynamics(0, S);

err = zeros(1, 2);
signs = [-1, 1];
for k = 1:2
    ic = struct('n_steps', 1, 'step_density', struct(), 'amp', 0, 'no_stim_pattern', true, ...
        'intrinsic_drive', signs(k) * abar, 'positive_only', false);
    mB = bursting_condition_model('original', 'input_config', ic, 'model_args', ...
        {'synapse_config', no_std, 'tau_a', {zeros(1, 0), zeros(1, 0)}, 'T_range', [0 110]});
    LB = mB.cached_params.state_layout;
    assert(LB.N == mB.n, 'Model B must have no adaptation or depression state.');
    dB = mB.dynamics(1, S(L.x) - abar);
    err(k) = max(abs(dB(LB.x) - dA(L.x)));
end
assert(err(1) < 1e-12, 'A constant input of -abar does not reproduce the frozen adaptation term (%.3g).', err(1));
assert(err(2) > 1e-3, 'The wrong sign was not told apart (%.3g).', err(2));
end

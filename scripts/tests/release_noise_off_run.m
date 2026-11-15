function out = release_noise_off_run(extra_args)
%RELEASE_NOISE_OFF_RUN One small run of each condition with the release noise off.
%
%   out = RELEASE_NOISE_OFF_RUN()             names no release property
%   out = RELEASE_NOISE_OFF_RUN(extra_args)   extra constructor name-value
%                                             pairs, e.g. {'release_rate', 0}
%
% What test_release_noise compares against
% fixtures/release_noise_off_reference.mat, which holds this function's
% output FROM THE CODE BEFORE THE RELEASE NOISE EXISTED (commit e859784,
% FD-027 step N1). The function names no release property, so it runs on
% that code and on today's; the comparison is by isequal, bit for bit. It
% is the statement "off means identical" on the four conditions of the
% chained-pool network, the pass-through routes of no_adaptation included,
% beside what the exact regression tiers already say on their networks.
%
% The network is the A-block preset (..._noise0p025_chainStd_4cond_mu7chain)
% at n = 40, with its Wiener noise on, so the stochastic integrator's path
% is the one exercised. Per condition:
%   S        the state every 0.2 s of a 2-s 'sra1' run (rows), with the end
%   LLE      Benettin's exponent (it re-integrates through the same solver)
%   dS       dynamics_fast at the final state
%   JY       compute_Jacobian_fast(S) * Y, and JtY = jacobian_times(S, Y),
%            for a fixed Y
%   theta_EE, theta_IE   two routes' synaptic output from
%            unpack_and_compute_states, at the stored times
%
% TO REMAKE THE FIXTURE (deliberately, as for any reference): check out
% the model code of e859784, delete the .mat, run
%   ref = release_noise_off_run(); save(<fixture>, '-struct', 'ref')
% and say in the commit why the expected values moved.
%
% See also: test_release_noise

if nargin < 1, extra_args = {}; end
[d, ~, conditions] = srnn_param_preset( ...
    'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain');
d.n = 40;
d.indegree = 8;
d.F_ref_n = 40;
d.F_ref_indegree = 8;
out = struct();
for k = 1:numel(conditions)
    condition = conditions{k};
    args = [struct2namevalue(d), struct2namevalue(rmfield(condition, 'name'))];
    m = SRNNCellTypePairs(args{:}, 'T_range', [0 2], 'fs', 200, 'ode_solver', 'sra1', ...
        'lya_method', 'benettin', 'lya_T_interval', [1 2], 'lya_warmup', 0.5, ...
        'store_full_state', true, 'rng_seeds', [3 4], extra_args{:});
    m.build();
    m.run();
    p = m.cached_params;
    p.u_interpolant = m.u_interpolant;
    S_end = m.S_out(end, :)';
    Y = cos((1:numel(S_end))' .* [1 2 3]);
    kept = 1:40:size(m.S_out, 1);
    [~, ~, ~, ~, theta] = SRNNCellTypePairs.unpack_and_compute_states(m.S_out(kept, :), p);
    r = struct();
    r.S = m.S_out(kept, :);
    r.LLE = m.lya_results.LLE;
    r.dS = SRNNCellTypePairs.dynamics_fast(m.t_out(end), S_end, p);
    r.JY = SRNNCellTypePairs.compute_Jacobian_fast(S_end, p) * Y;
    r.JtY = SRNNCellTypePairs.jacobian_times(S_end, Y, p);
    r.theta_EE = theta.E.E;
    r.theta_IE = theta.I.E;
    out.(condition.name) = r;
end
end

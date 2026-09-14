% test_route_scale.m - the per-route weight scale in SRNNCellTypePairs.
%
% synapse_config.<pre>.<post>.scale multiplies W(post, pre) inside get_params
% (params.W only; obj.W is untouched). It exists to match the steady-state
% synaptic output of a two-timescale STD route to a one-timescale one at a
% reference rate (the dualStdScaled presets).
%
% Checks, on the all-blocks 3-type network of test_jacobian_times:
%   1. Default: params.route_scale is all ones and params.W == obj.W.
%   2. A scaled route (E->SST x 2.5, PV->PV x 0.4): obj.W unchanged;
%      params.W's blocks scaled; the trajectory is bit-identical to the
%      unscaled model whose obj.W block was multiplied by hand.
%   3. The analytic Jacobian matches central finite differences with the
%      scale on (1e-5 relative), and jacobian_times equals J * Y (1e-12).
%   4. Bad scales (0, negative, vector, NaN) error
%      SRNNCellTypePairs:InvalidSynapseConfig; an unknown route field still errors.
%   5. The ESN route-redundancy check refuses two routes of one presynaptic
%      type with different scales and accepts equal ones.
%   6. The scaled paper preset builds, its dual routes carry the scale, its
%      single-timescale condition does not, and the usage preset carries the
%      matched tau_rel.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: SRNNCellTypePairs.compile_synapse_config, srnn_param_preset,
%           fig_STD_steady_state, test_jacobian_times

fprintf('=== Testing the per-route weight scale ===\n\n');
all_passed = true;

%% 1. default
rng(0, 'twister');
m0 = make_model(struct());
m0.build();
p0 = m0.get_params();
all_passed = check('default route_scale is ones and params.W == obj.W', ...
    isequal(p0.route_scale, ones(3)) && isequal(p0.W, m0.W)) && all_passed;

%% 2. scaled routes
sc = struct(); sc.E.SST.scale = 2.5; sc.PV.PV.scale = 0.4;
rng(0, 'twister');
m1 = make_model(sc);
m1.build();
p1 = m1.get_params();
ti = m1.type_indices;
W_hand = m0.W;
W_hand(ti{3}, ti{1}) = 2.5 * W_hand(ti{3}, ti{1});
W_hand(ti{2}, ti{2}) = 0.4 * W_hand(ti{2}, ti{2});
all_passed = check('obj.W unchanged by the scale; params.W blocks scaled', ...
    isequal(m1.W, m0.W) && isequal(p1.W, W_hand) && p1.route_scale(1, 3) == 2.5 && p1.route_scale(2, 2) == 0.4) && all_passed;
m1.run();
% obj.W is read-only, so the hand-scaled comparison is at the params level: the
% RHS with the scaled route equals the RHS of the unscaled params whose W block
% was multiplied by hand, exactly, at random states.
p_hand = p0; p_hand.W = W_hand;
p_hand.u_interpolant = m0.u_interpolant; p1.u_interpolant = m1.u_interpolant;
same = true;
for k = 1:20
    S = m1.S_out(randi(size(m1.S_out, 1)), :)';
    same = same && isequal(SRNNCellTypePairs.dynamics_fast(0, S, p1), SRNNCellTypePairs.dynamics_fast(0, S, p_hand));
end
all_passed = check('RHS bit-identical to the hand-scaled W at 20 states', same) && all_passed;

%% 3. Jacobian and jacobian_times with the scale on
p1.u_interpolant = m1.u_interpolant;
S = m1.S_out(end, :)';
J = SRNNCellTypePairs.compute_Jacobian_fast(S, p1);
Jfd = SRNNCellTypePairs.finite_difference_jacobian(S, p1, 1e-6);
e_fd = norm(full(J) - Jfd, 'fro') / norm(Jfd, 'fro');
Y = randn(numel(S), 5);
e_jt = norm(SRNNCellTypePairs.jacobian_times(S, Y, p1) - J * Y, 'fro') / norm(J * Y, 'fro');
all_passed = check(sprintf('analytic Jacobian vs finite differences with scale (rel %.1e)', e_fd), e_fd < 1e-5) && all_passed;
all_passed = check(sprintf('jacobian_times == J*Y with scale (rel %.1e)', e_jt), e_jt < 1e-12) && all_passed;

%% 4. validation
bad = {0, -1, [1 2], NaN};
ok = true;
for k = 1:numel(bad)
    s = struct(); s.E.E.scale = bad{k};
    ok = ok && throws_id(@() build_bad(s), 'SRNNCellTypePairs:InvalidSynapseConfig');
end
s = struct(); s.E.E.gain = 2;
ok = ok && throws_id(@() build_bad(s), 'SRNNCellTypePairs:InvalidSynapseConfig');
all_passed = check('bad scales and an unknown route field error InvalidSynapseConfig', ok) && all_passed;

%% 5. ESN route redundancy
dual = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
sc_eq = struct(); sc_eq.E.E = struct('std', dual, 'scale', 3); sc_eq.E.I = struct('std', dual, 'scale', 3);
sc_ne = sc_eq; sc_ne.E.I.scale = 1;
pe = esn_params(sc_eq); pn = esn_params(sc_ne);
all_passed = check('assert_route_redundancy accepts equal scales, refuses different ones', ...
    ~throws_prefix(@() SRNN_ESN_reservoir.assert_route_redundancy(pe), 'SRNN_ESN_reservoir:') && ...
    throws_prefix(@() SRNN_ESN_reservoir.assert_route_redundancy(pn), 'SRNN_ESN_reservoir:')) && all_passed;

%% 6. the presets
[~, ~, cs] = srnn_param_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdScaled_3cond_mu8p25');
[~, ~, cu] = srnn_param_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdUsage_3cond_mu8p25');
[~, ~, c0] = srnn_param_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25');
ms = build_from_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdScaled_3cond_mu8p25', 'sfa3_std2', 'n', 20, 'indegree', 5, 'T_range', [0 0.1]);
ps = ms.get_params();
s_val = cs{3}.synapse_config.E.E.scale;
r_ref = 0.25; rho = 0.125;
all_passed = check(sprintf('scaled preset: scale %.4g on all four dual routes = 1 + r_ref/rho at r_ref 0.25', s_val), ...
    all(ps.route_scale(:) == s_val) && abs(s_val - (1 + r_ref / rho)) < 1e-12 && ...
    isequal(cs{3}.synapse_config.E.E.std, c0{3}.synapse_config.E.E.std)) && all_passed;
all_passed = check('scaled preset: no_adaptation and sfa1_std1 identical to the unscaled preset', ...
    isequaln(cs{1}, c0{1}) && isequaln(cs{2}, c0{2})) && all_passed;
tr = cu{3}.synapse_config.E.E.std.tau_rel; trec = cu{3}.synapse_config.E.E.std.tau_rec;
rho_u = r_ref / (sqrt(1 + r_ref / rho) - 1);
theta_single = r_ref / (1 + r_ref / rho);
theta_usage  = r_ref / prod(1 + r_ref ./ (tr ./ trec));
all_passed = check(sprintf('usage preset: tau_rel = rho_u * tau_rec (rho_u %.4f) and steady state matches at r_ref (%.4f vs %.4f)', rho_u, theta_usage, theta_single), ...
    max(abs(tr - rho_u * trec)) < 1e-3 && abs(theta_usage - theta_single) < 1e-3 && ...
    isequaln(cu{1}, c0{1}) && isequaln(cu{2}, c0{2})) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL route scale TESTS PASSED ===\n');
else
    fprintf('=== SOME route scale TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function model = make_model(scales)
% The all-blocks network of test_jacobian_times, plus optional route scales.
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', [0.3 1], 'tau_rel', 0.2);
synapse_config.E.PV.stf = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
synapse_config.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
synapse_config.E.SST.stf = struct( ...
    'tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
synapse_config.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
synapse_config.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
pres = fieldnames(scales);
for i = 1:numel(pres)
    posts = fieldnames(scales.(pres{i}));
    for j = 1:numel(posts)
        synapse_config.(pres{i}).(posts{j}).scale = scales.(pres{i}).(posts{j}).scale;
    end
end
model = SRNNCellTypePairs( ...
    'n', 12, 'indegree', 4, ...
    'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], ...
    'mu_tilde_relative', [0.1 -0.2 -0.08], ...
    'sigma_tilde_relative', [0.01 0.02 0.015], ...
    'tau_a', {[0.25 1], [], 2}, ...
    'c', [0.05 0 0.03], 'synapse_config', synapse_config, ...
    'T_range', [0 0.5], 'fs', 200, 'lya_method', 'none', ...
    'store_full_state', true);
end

function build_bad(sc)
m = SRNNCellTypePairs('n', 6, 'indegree', 2, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [0.1 -0.1], 'sigma_tilde_relative', [0.01 0.01], ...
    'synapse_config', sc, 'T_range', [0 0.1], 'lya_method', 'none');
m.build();
end

function p = esn_params(sc)
m = SRNNCellTypePairs('n', 6, 'indegree', 2, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [0.1 -0.1], 'sigma_tilde_relative', [0.01 0.01], ...
    'synapse_config', sc, 'T_range', [0 0.1], 'lya_method', 'none');
m.build();
p = m.get_params();
end

function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
end

function ok = throws_prefix(fn, prefix)
ok = false;
try
    fn();
catch ME
    ok = startsWith(ME.identifier, prefix);
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

% test_std_geomean.m - the opt-in geometric-mean STD combine in SRNNCellTypePairs.
%
% synapse_config.<pre>.<post>.std.combine = 'geomean' replaces a route's
% depression product prod_m b_m in theta = r * f(b) * prod(g) by
% (prod_m b_m)^(1/M). With every pair at one rho = tau_rel/tau_rec the
% geomean steady state is the one-timescale 1/(1 + r/rho) at every constant
% rate (the product squares it at M = 2). The default ('product', or the
% field omitted) must be unchanged.
%
% Checks:
%   1. Default: params.std_geomean all false; omitted combine and explicit
%      'product' give bit-identical trajectories on the all-blocks network
%      of test_jacobian_times; the pure-STD route's synaptic output is
%      r .* prod(b) exactly.
%   2. Geomean route: synaptic output is r .* prod(b).^(1/M) exactly, the
%      product routes elsewhere are untouched.
%   3. Constant rate (three unconnected neurons, W = 0, no SFA, constant
%      drive): the dual-timescale geomean route's steady-state synaptic output
%      equals the single-timescale route's and r/(1 + r/rho) at every rate
%      (rho = 0.125: tau_rec [2 4], tau_rel [0.25 0.5]); the product route
%      gives the square of the factor.
%   4. The analytic Jacobian matches central finite differences with geomean
%      on (with and without std_zero_floor), and jacobian_times == J * Y.
%   5. A bad combine errors SRNNCellTypePairs:InvalidSynapseConfig; routes
%      that differ only in combine are not identical (routes_identical and
%      the ESN route-redundancy check).
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: SRNNCellTypePairs.combine_std_pair, test_route_scale,
%           test_jacobian_times

fprintf('=== Testing the geometric-mean STD combine ===\n\n');
all_passed = true;

%% 1. default unchanged
rng(0, 'twister');
m0 = make_model('');
m0.build(); m0.run();
p0 = m0.get_params();
rng(0, 'twister');
mp = make_model('product');
mp.build(); mp.run();
all_passed = check('default std_geomean is all false (omitted and explicit product)', ...
    isequal(p0.std_geomean, false(3)) && isequal(mp.get_params().std_geomean, false(3))) && all_passed;
all_passed = check('omitted combine and explicit product: bit-identical trajectories', ...
    isequal(m0.S_out, mp.S_out)) && all_passed;
pd0 = m0.plot_data;
B = pd0.b.E.E;
theta_hand = pd0.r.E .* reshape(prod(B, 2), size(B, 1), []);
all_passed = check('default E->E synaptic output == r .* prod(b) exactly', ...
    isequal(pd0.synaptic_output.E.E, theta_hand)) && all_passed;

%% 2. geomean route
rng(0, 'twister');
mg = make_model('geomean');
mg.build(); mg.run();
pg = mg.get_params();
pdg = mg.plot_data;
B = pdg.b.E.E;
theta_geo = pdg.r.E .* reshape(prod(B, 2), size(B, 1), []) .^ (1 / size(B, 2));
e_geo = max(abs(pdg.synaptic_output.E.E(:) - theta_geo(:)));
all_passed = check(sprintf('geomean E->E synaptic output == r .* prod(b).^(1/M) (max err %.1e)', e_geo), ...
    pg.std_geomean(1, 1) && nnz(pg.std_geomean) == 1 && e_geo < 1e-14) && all_passed;
all_passed = check('geomean changes the trajectory', ~isequal(mg.S_out, m0.S_out)) && all_passed;

%% 3. steady state at constant rate
tau_rec = [2 4]; tau_rel = [0.25 0.5]; rho = 0.125;
drive = [0.2; 0.4; 0.7];                    % logistic S_c 0.4: r ~ 0.31, 0.5, 0.77
ss_single  = steady_output(struct('tau_rec', 2, 'tau_rel', 0.25), drive);
ss_product = steady_output(struct('tau_rec', tau_rec, 'tau_rel', tau_rel), drive);
ss_geo     = steady_output(struct('tau_rec', tau_rec, 'tau_rel', tau_rel, 'combine', 'geomean'), drive);
r = ss_geo.r;
theta_formula = r ./ (1 + r ./ rho);
e_single  = max(abs(ss_geo.theta - ss_single.theta) ./ ss_single.theta);
e_formula = max(abs(ss_geo.theta - theta_formula) ./ theta_formula);
e_square  = max(abs(ss_product.theta - r ./ (1 + r ./ rho) .^ 2) ./ ss_product.theta);
all_passed = check(sprintf('geomean steady state == single timescale at r = %s (rel %.1e)', ...
    mat2str(round(r', 3)), e_single), e_single < 1e-8 && max(abs(ss_single.r - r)) < 1e-12) && all_passed;
all_passed = check(sprintf('geomean steady state == r/(1 + r/rho) (rel %.1e)', e_formula), e_formula < 1e-8) && all_passed;
all_passed = check(sprintf('product steady state == r/(1 + r/rho)^2 (rel %.1e)', e_square), e_square < 1e-8) && all_passed;

%% 4. Jacobian with geomean on
for zf = [false true]
    rng(0, 'twister');
    mj = make_model('geomean');
    mj.std_zero_floor = zf;
    mj.build(); mj.run();
    pj = mj.get_params();
    pj.u_interpolant = mj.u_interpolant;
    e_fd = 0; e_jt = 0;
    for k = [round(size(mj.S_out, 1) / 2), size(mj.S_out, 1)]
        S = mj.S_out(k, :)';
        J = SRNNCellTypePairs.compute_Jacobian_fast(S, pj);
        Jfd = SRNNCellTypePairs.finite_difference_jacobian(S, pj, 1e-6);
        e_fd = max(e_fd, norm(full(J) - Jfd, 'fro') / norm(Jfd, 'fro'));
        Y = randn(numel(S), 5);
        e_jt = max(e_jt, norm(SRNNCellTypePairs.jacobian_times(S, Y, pj) - J * Y, 'fro') / norm(J * Y, 'fro'));
    end
    all_passed = check(sprintf('geomean (zero floor %d): analytic Jacobian vs finite differences (rel %.1e)', zf, e_fd), ...
        e_fd < 1e-5) && all_passed;
    all_passed = check(sprintf('geomean (zero floor %d): jacobian_times == J*Y (rel %.1e)', zf, e_jt), ...
        e_jt < 1e-12) && all_passed;
end

%% 5. validation and route identity
ok = true;
bad = {'mean', 3, {'geomean'}, ''''};
for k = 1:numel(bad)
    s = struct(); s.E.E.std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
    s.E.E.std.combine = bad{k};   % assigned, not via struct(): a cell would expand
    ok = ok && throws_id(@() build_two_type(s), 'SRNNCellTypePairs:InvalidSynapseConfig');
end
s = struct(); s.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25, 'combin', 'geomean');
ok = ok && throws_id(@() build_two_type(s), 'SRNNCellTypePairs:InvalidSynapseConfig');
all_passed = check('bad combine values and a misspelt field error InvalidSynapseConfig', ok) && all_passed;
dual = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
dual_geo = dual; dual_geo.combine = 'geomean';
s_eq = struct(); s_eq.E.E.std = dual_geo; s_eq.E.I.std = dual_geo;
s_ne = struct(); s_ne.E.E.std = dual_geo; s_ne.E.I.std = dual;
pe = build_two_type(s_eq).get_params(); pn = build_two_type(s_ne).get_params();
all_passed = check('ESN redundancy accepts equal combines and refuses different ones', ...
    ~throws_prefix(@() SRNN_ESN_reservoir.assert_route_redundancy(pe), 'SRNN_ESN_reservoir:') && ...
    throws_prefix(@() SRNN_ESN_reservoir.assert_route_redundancy(pn), 'SRNN_ESN_reservoir:')) && all_passed;
s_all = struct(); s_all.E.E.std = dual; s_all.E.I.std = dual; s_all.I.E.std = dual; s_all.I.I.std = dual;
s_mix = s_all; s_mix.I.I.std = dual_geo;
all_passed = check('routes_identical is false when only combine differs', ...
    SRNNCellTypePairs.routes_identical(build_two_type(s_all).get_params()) && ...
    ~SRNNCellTypePairs.routes_identical(build_two_type(s_mix).get_params())) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL STD geomean TESTS PASSED ===\n');
else
    fprintf('=== SOME STD geomean TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function model = make_model(combine)
% The all-blocks network of test_jacobian_times; combine ('' = omitted) is
% set on the dual-timescale E->E route only.
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', [0.3 1], 'tau_rel', 0.2);
if ~isempty(combine)
    synapse_config.E.E.std.combine = combine;
end
synapse_config.E.PV.stf = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
synapse_config.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
synapse_config.E.SST.stf = struct( ...
    'tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
synapse_config.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
synapse_config.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
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

function out = steady_output(std_config, drive)
% Unconnected neurons (W = 0, no SFA) held at a constant drive: the rate is
% constant after ~tau_d, and b relaxes to its steady state within a few s.
sc = struct(); sc.E.E.std = std_config;
n = numel(drive);
m = SRNNCellTypePairs('n', n, 'indegree', 1, 'n_cellTypes', 1, 'cell_type_names', {'E'}, ...
    'f', 1, 'mu_tilde_relative', 0, 'sigma_tilde_relative', 0, ...
    'tau_a', {zeros(1, 0)}, 'c', 0, 'x0_std', 0, 'synapse_config', sc, ...
    'T_range', [0 60], 'fs', 100, 'ode_solver', 'ode45', 'lya_method', 'none');
m.input_config.amp = 0;
m.input_config.intrinsic_drive = drive;
m.build(); m.run();
out.r = m.plot_data.r.E(:, end);
out.theta = m.plot_data.synaptic_output.E.E(:, end);
end

function m = build_two_type(sc)
m = SRNNCellTypePairs('n', 6, 'indegree', 2, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [0.1 -0.1], 'sigma_tilde_relative', [0.01 0.01], ...
    'synapse_config', sc, 'T_range', [0 0.1], 'lya_method', 'none');
m.build();
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

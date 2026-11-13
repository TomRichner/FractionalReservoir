% test_self_feedback.m - negative self-feedback on the diagonal of W (SRNNCellTypePairs.self_feedback).
%
% self_feedback is one dimensionless strength k per cell type. get_params
% subtracts k R_W from the diagonal of params.W, R_W being the theoretical
% bulk radius at level_of_chaos = 1 (R_gain1), after the gain and after any
% route scale; obj.W stays the drawn matrix.
%
% Checks, on a 60-neuron cut of the mu7chain memory-capacity network:
%   1. Default (and k = 0, scalar or row): params.W is obj.W to the bit, and
%      the right-hand side, the Jacobian and jacobian_times are bit-identical
%      to a model that never set the property.
%   2. k = [0.5 1.5] by type: obj.W unchanged; the diagonal of params.W is
%      the drawn diagonal minus k_q R_W per neuron of type q; off-diagonal
%      entries untouched. A scalar k is broadcast to every type.
%   3. Independent of the gain: R_gain1 does not move with level_of_chaos
%      (R does), and at gain 2 the diagonal shift is the same k R_W. With a
%      route scale on E->E the scaled block's diagonal is scale * drawn - k R_W.
%   4. For no adaptation and for chained MTS, k > 0: the analytic Jacobian
%      matches central finite differences (1e-5 relative, kink rows
%      excluded), and jacobian_times equals J * Y (1e-12).
%   5. No adaptation, a fixed state: the Jacobian equals
%      (-I + (W - k R_W I) Phi') / tau_d built by hand from obj.W, and its
%      eigenvalues are that matrix's.
%   6. Bad values (wrong length, NaN, Inf) error SRNNCellTypePairs:InvalidParams.
%   7. The property round-trips through ParamSpaceAnalysis2's
%      resolved_defaults and effective_param.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: SRNNCellTypePairs.get_params, test_route_scale, test_jacobian_times

fprintf('=== Testing self-feedback on the diagonal of W ===\n\n');
all_passed = true;
K = [0.5 1.5];                      % by type: E, I

%% 1. default and k = 0
m0 = make_model('sfa3_stdc2');
p0 = m0.get_params();
all_passed = check('default: self_feedback is zeros(1, C) and params.W is obj.W', ...
    isequal(m0.self_feedback, [0 0]) && isequal(p0.W, m0.W)) && all_passed;
same = true;
for value = {0, [0 0]}
    mz = make_model('sfa3_stdc2', 'self_feedback', value{1});
    pz = mz.get_params();
    S = initial_state(mz);
    Y = probe_columns(numel(S));
    same = same && isequal(pz.W, p0.W) ...
        && isequal(rhs(mz, S), rhs(m0, S)) ...
        && isequal(SRNNCellTypePairs.compute_Jacobian_fast(S, pz), SRNNCellTypePairs.compute_Jacobian_fast(S, p0)) ...
        && isequal(SRNNCellTypePairs.jacobian_times(S, Y, pz), SRNNCellTypePairs.jacobian_times(S, Y, p0));
end
all_passed = check('k = 0 (scalar or row): W, RHS, Jacobian and jacobian_times bit-identical', same) && all_passed;

%% 2. the diagonal
m1 = make_model('sfa3_stdc2', 'self_feedback', K);
p1 = m1.get_params();
R_W = m1.R_gain1;
k_i = zeros(m1.n, 1);
for q = 1:2; k_i(m1.type_indices{q}) = K(q); end
off = @(A) A - diag(diag(A));
d_err = max(abs(full(diag(p1.W)) - (full(diag(m1.W)) - k_i * R_W)));
all_passed = check(sprintf('diag(params.W) = drawn diagonal - k_q R_W (max error %.1e); obj.W and off-diagonal unchanged', d_err), ...
    d_err <= 4 * eps(R_W * max(K)) && isequal(m1.W, m0.W) && isequal(off(p1.W), off(m0.W)) && R_W > 0) && all_passed;
ms = make_model('sfa3_stdc2', 'self_feedback', 0.75);
ps = ms.get_params();
all_passed = check('a scalar k is broadcast to every type', ...
    isequal(ms.self_feedback, [0.75 0.75]) && ...
    max(abs(full(diag(ps.W)) - (full(diag(ms.W)) - 0.75 * R_W))) <= 4 * eps(R_W)) && all_passed;

%% 3. independent of the gain and of the route scale
m2 = make_model('sfa3_stdc2', 'self_feedback', K, 'level_of_chaos', 2);
p2 = m2.get_params();
gain_ok = m2.R_gain1 == R_W && abs(m2.R - 2 * R_W) <= 4 * eps(2 * R_W) ...
    && isequal(m2.W, 2 * m0.W) ...
    && max(abs(full(diag(p2.W)) - (2 * full(diag(m0.W)) - k_i * R_W))) <= 8 * eps(R_W * max(K));
all_passed = check('gain 2: R doubles, R_gain1 and the diagonal shift do not', gain_ok) && all_passed;
sc = m0.synapse_config; sc.E.E.scale = 3;
m3 = make_model('sfa3_stdc2', 'self_feedback', K, 'synapse_config', sc);
p3 = m3.get_params();
iE = m3.type_indices{1};
d_E = full(diag(p3.W)); d_E = d_E(iE);
d0_E = full(diag(m0.W)); d0_E = d0_E(iE);
all_passed = check('route scale 3 on E->E: its diagonal is 3 * drawn - k_E R_W (the scale does not multiply k)', ...
    max(abs(d_E - (3 * d0_E - K(1) * R_W))) <= 8 * eps(R_W)) && all_passed;

%% 4. Jacobians with k > 0
for name = {'no_adaptation', 'sfa3_stdc2'}
    m = make_model(name{1}, 'self_feedback', K);
    m.run();
    p = m.get_params();
    p.u_interpolant = m.u_interpolant;
    S = m.S_out(end, :)';
    J = SRNNCellTypePairs.compute_Jacobian_fast(S, p);
    Jfd = SRNNCellTypePairs.finite_difference_jacobian(S, p, 1e-6);
    e_fd = norm(full(J) - Jfd, 'fro') / norm(Jfd, 'fro');
    Y = probe_columns(numel(S));
    e_jt = norm(SRNNCellTypePairs.jacobian_times(S, Y, p) - J * Y, 'fro') / norm(J * Y, 'fro');
    all_passed = check(sprintf('%s: analytic Jacobian vs finite differences (rel %.1e)', name{1}, e_fd), e_fd < 1e-5) && all_passed;
    all_passed = check(sprintf('%s: jacobian_times == J * Y (rel %.1e)', name{1}, e_jt), e_jt < 1e-12) && all_passed;
end

%% 5. J_eff = (-I + (W - k R_W I) Phi') / tau_d, no adaptation
m = make_model('no_adaptation', 'self_feedback', K);
m.run();
p = m.get_params();
x = m.S_out(end, :)';                       % no adaptation: the state is x
assert(numel(x) == m.n, 'no adaptation: the state should be x alone');
dphi = m.activation_function_derivative(x);
W_sim = full(m.W) - diag(k_i * m.R_gain1);   % from the DRAWN matrix, by hand
J_pred = (-eye(m.n) + W_sim * diag(dphi)) / m.tau_d;
J = full(SRNNCellTypePairs.compute_Jacobian_fast(x, p));
e_J = norm(J - J_pred, 'fro') / norm(J_pred, 'fro');
lam = sort_complex(eig(J)); lam_pred = sort_complex(eig(J_pred));
e_lam = max(abs(lam - lam_pred)) / max(abs(lam_pred));
% And against k = 0: the shift is -k_i R_W phi'_i / tau_d on the diagonal.
m_free = make_model('no_adaptation');
J_free = full(SRNNCellTypePairs.compute_Jacobian_fast(x, m_free.get_params()));
e_shift = max(abs(diag(J - J_free) + k_i .* m.R_gain1 .* dphi / m.tau_d)) / max(abs(diag(J_free)));
all_passed = check(sprintf('J = (-I + (W - k R_W I) Phi'')/tau_d (rel %.1e), eigenvalues (rel %.1e), diagonal shift vs k = 0 (rel %.1e)', ...
    e_J, e_lam, e_shift), e_J < 1e-13 && e_lam < 1e-10 && e_shift < 1e-13 ...
    && isequal(off(J), off(J_free))) && all_passed;

%% 6. validation
bad = {[1 2 3], NaN, [0 Inf]};
ok = true;
for k = 1:numel(bad)
    ok = ok && throws_id(@() make_model('no_adaptation', 'self_feedback', bad{k}), 'SRNNCellTypePairs:InvalidParams');
end
all_passed = check('wrong length, NaN and Inf are refused', ok) && all_passed;

%% 7. resolved_defaults
tmp = tempname;
cleanup = onCleanup(@() remove_dir(tmp));
[preset, model_class, conditions] = srnn_param_preset(preset_name());
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 4, 'note', 'self_feedback_test', 'verbose', 'near-none');
psa.model_class = model_class;
psa.use_parallel = false;
psa.output_dir = tmp;
psa.model_defaults = preset;
small = small_args();
for k = 1:2:numel(small); psa.model_defaults.(small{k}) = small{k + 1}; end
psa.model_defaults.self_feedback = K;
psa.model_defaults.lya_method = 'none';
psa.set_conditions(conditions(1));
psa.add_grid_parameter('level_of_chaos', [0.5 1]);
psa.run();
round_trip = isequal(psa.resolved_defaults.self_feedback, K) && isequal(psa.effective_param([], 'self_feedback'), K);
L = load(fullfile(psa.output_dir, 'psa_object.mat'));
loaded = L.(char(fieldnames(L)));
round_trip = round_trip && isequal(loaded.resolved_defaults.self_feedback, K);
all_passed = check('self_feedback is frozen into resolved_defaults, read by effective_param and saved', round_trip) && all_passed;

%% Banner
fprintf('\n');
if all_passed
    fprintf('=== ALL SELF-FEEDBACK TESTS PASSED ===\n');
else
    fprintf('=== SOME SELF-FEEDBACK TESTS FAILED ===\n');
    error('test_self_feedback:Failed', 'One or more self-feedback checks failed.');
end

%% ------------------------------------------------------------------------
function name = preset_name()
name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
end

function args = small_args()
% A cut of the network small and short enough for a unit test.
args = {'n', 60, 'indegree', 20, 'T_range', [0 3], 'fs', 200, 'ode_solver', 'rk4', ...
    'rng_seeds', [11 12], 'verbose', 'near-none'};
end

function m = make_model(condition, varargin)
% The preset, cut small, in the named condition, with overrides last; built.
[preset, ~, conditions] = srnn_param_preset(preset_name());
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond = conditions{strcmp(names, condition)};
preset_args = struct2namevalue(preset);
cond_args = struct2namevalue(rmfield(cond, 'name'));
small = small_args();
m = SRNNCellTypePairs(preset_args{:}, cond_args{:}, small{:}, 'lya_method', 'none', ...
    'store_full_state', true, varargin{:});
m.build();
end

function S = initial_state(m)
% A reproducible state off the rest point, the same for models of one size.
rng(5, 'twister');
S = 0.3 * rand(m.N_sys_eqs, 1);
end

function Y = probe_columns(n)
rng(6, 'twister');
Y = randn(n, 4);
end

function d = rhs(m, S)
p = m.get_params();
p.u_interpolant = m.u_interpolant;
d = SRNNCellTypePairs.dynamics_fast(0, S, p);
end

function v = sort_complex(v)
[~, order] = sortrows([real(v), imag(v)]);
v = v(order);
end

function ok = check(label, cond)
ok = logical(cond);
if ok; fprintf('  PASS  %s\n', label); else; fprintf('  FAIL  %s\n', label); end
end

function tf = throws_id(f, id)
tf = false;
try
    f();
catch err
    tf = strcmp(err.identifier, id);
end
end

function remove_dir(d)
if isfolder(d); rmdir(d, 's'); end
end

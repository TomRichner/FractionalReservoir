function test_rng_generator_independent()
% TEST_RNG_GENERATOR_INDEPENDENT A seed means one network, whatever generator the session has.
%
% Written 2026-11-12 (Plan C, fix F1). Until then the model seeded most of
% its draws with a bare rng(seed), which keeps whatever generator the
% global stream already has. An ordinary client session has the Mersenne
% twister; a parallel pool worker has Threefry. So the same seeds gave one
% network in a client and ANOTHER on a pool worker: a sweep job could not
% be rebuilt in a client to the bit without first setting the generator by
% hand, and only the stage functions that call rng(0, 'twister') before
% building were safe. Every seeding call in src/ now names the generator,
% rng(seed, 'twister'), and this test holds that in place.
%
% (a) A small SRNNCellTypePairs model (per-neuron setpoints, a per-neuron
%     SFA ladder, per-neuron depression and facilitation parameters, a
%     step input, noise on under 'sra1') and a small SRNN_ESN_reservoir
%     (input weights, a sample-and-hold input, noise on) are built and run
%     twice: once with the session's generator as it is, and once after
%     rng(0, 'threefry'), which is what a pool worker starts from. W, the
%     setpoints, the SFA ladders, the synapse draws, the input, the
%     initial state, the noise path and the trajectory must be
%     bit-identical (isequaln). The session's generator is restored at the
%     end, also when the test errors (onCleanup).
% (b) If this session holds a pool (pool_ensure), the same two models are
%     built and run on a worker inside a parfor and must equal the
%     client's. Each worker is first put on Threefry, where a fresh worker
%     starts, and the generator it built from is printed. With no pool,
%     (b) prints SKIP: run it after pool_set(2).
%
% Prints PASS/FAIL per check and a final banner, and errors at the end if
% anything failed. About 10 s. Assumes setup_paths has run.
%
% See also: test_pairs_reproducibility (same generator, same seeds),
%           test_shared_build_spread

fprintf('=== Testing that a model does not depend on the session''s random generator ===\n\n');
all_passed = true;
state_at_entry = rng;
restore = onCleanup(@() rng(state_at_entry));

%% (a) Client: the session's generator, then Threefry
fprintf('--- (a) client: generator ''%s'', then ''threefry'' ---\n', state_at_entry.Type);
kinds = {'pairs', 'esn'};
ref = cell(1, numel(kinds));
for k = 1:numel(kinds)
    rng(state_at_entry);
    ref{k} = fingerprint(kinds{k});
    rng(0, 'threefry');
    other = fingerprint(kinds{k});
    during = rng;
    all_passed = check(sprintf('%s: the second build really ran under another generator (''%s'')', ...
        kinds{k}, other.generator_at_entry), ...
        strcmp(other.generator_at_entry, 'threefry') || strcmp(state_at_entry.Type, 'threefry')) && all_passed;
    [same, what] = same_fingerprint(ref{k}, other);
    all_passed = check(sprintf('%s: bit-identical under ''%s'' and ''threefry''%s', ...
        kinds{k}, state_at_entry.Type, what), same) && all_passed;
    all_passed = check(sprintf('%s: the run did something (finite, states moved, noise path of %d steps)', ...
        kinds{k}, size(ref{k}.noise_xi1, 2)), ...
        all(isfinite(ref{k}.S_out), 'all') && any(ref{k}.S_out(end, :) ~= ref{k}.S_out(1, :)) && ...
        ~isempty(ref{k}.noise_xi1) && any(ref{k}.noise_xi1(:) ~= 0)) && all_passed;
    fprintf('  (%s: global stream after build and run is ''%s'')\n', kinds{k}, during.Type);
end
rng(state_at_entry);

%% (b) Pool worker against the client
fprintf('\n--- (b) pool worker against the client ---\n');
pool = pool_ensure();
if isempty(pool)
    fprintf('  SKIP: this session holds no pool (run pool_set(2) first)\n');
else
    n_kinds = numel(kinds);
    from_worker = cell(1, n_kinds);
    parfor k = 1:n_kinds
        % A fresh worker starts on Threefry, but one that has built a
        % model before is left on the twister (the model's seeding calls
        % set it and build() does not restore the stream). Put the worker
        % where a fresh one is, so a second run of this test in the same
        % pool still crosses generators.
        rng(0, 'threefry');
        from_worker{k} = fingerprint(kinds{k});
    end
    for k = 1:n_kinds
        [same, what] = same_fingerprint(ref{k}, from_worker{k});
        all_passed = check(sprintf('%s: built on a worker (generator ''%s'') equals the client''s%s', ...
            kinds{k}, from_worker{k}.generator_at_entry, what), same) && all_passed;
    end
end

fprintf('\n');
if all_passed
    fprintf('=== ALL generator-independence TESTS PASSED ===\n');
else
    fprintf('=== SOME generator-independence TESTS FAILED ===\n');
    error('test_rng_generator_independent:Failed', 'A model depends on the session''s random generator.');
end
end

%% ------------------------------------------------------------------------
function fp = fingerprint(kind)
% Build and run one small model and return everything random about it.
s = rng;
fp = struct('generator_at_entry', s.Type);

dual = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
fac  = struct('tau_dec', 1.0, 'tau_fac', 0.5, 'G', 2.0);
sc = struct();
for pre = 'EI'
    for post = 'EI'
        sc.(pre).(post).std = dual;
        sc.(pre).(post).stf = fac;
    end
end
common = {'n', 40, 'indegree', 8, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [3 -4; 3 -4], 'sigma_tilde_relative', ones(2, 2), ...
    'level_of_chaos', 1.5, 'S_c', 0.35, 'sigma_S_c', [0.1 0.1], ...
    'tau_a', {[0.25 1.58 10], [0.25 1.58 10]}, 'tau_a_spread', [0.25 0.25], ...
    'synapse_config', sc, 'std_spread', [0.25 0.25], 'stf_spread', [0.25 0.25], ...
    'rng_seeds', [5 6], 'fs', 100, 'ode_solver', 'sra1', 'sigma_u_noise', 0.025};

switch kind
    case 'pairs'
        m = SRNNCellTypePairs(common{:}, 'T_range', [0 3], 'store_full_state', true, ...
            'lya_method', 'none');
        evalc('m.build(); m.build_noise();');
        fp.noise_xi1 = m.noise_increments.xi1;
        fp.noise_xi2 = m.noise_increments.xi2;
        evalc('m.run();');
    case 'esn'
        m = SRNN_ESN_reservoir(common{:}, 'T_wash', 100, 'T_train', 200, 'T_test', 100, ...   % in samples: 4 s
            'd_max', 4, 'input_type', 'sample_hold', 'T_hold', 0.3);
        evalc('m.build(); m.build_noise();');
        fp.noise_xi1 = m.noise_increments.xi1;
        fp.noise_xi2 = m.noise_increments.xi2;
        evalc('m.run_reservoir_esn();');
        fp.W_in = m.W_in;
        fp.u_scalar = m.u_scalar;
end
fp.W = m.W;
fp.S_c_vec = m.S_c_vec;
fp.tau_a_matrix = m.tau_a_matrix;
fp.std_draws = m.std_draws;
fp.stf_draws = m.stf_draws;
fp.u_ex = m.u_ex;
fp.S0 = m.S0;
fp.S_out = m.S_out;
end

function [same, what] = same_fingerprint(A, B)
% Compare everything but the generator the build started from. `what`
% names the first field that differs.
same = true; what = '';
fields = setdiff(fieldnames(A), {'generator_at_entry'}, 'stable');
for k = 1:numel(fields)
    if ~isfield(B, fields{k}) || ~isequaln(A.(fields{k}), B.(fields{k}))
        same = false;
        what = sprintf(' (first difference: %s)', fields{k});
        return
    end
end
if isempty(A.W) || isempty(A.S_c_vec) || isempty(A.tau_a_matrix) || isempty(A.std_draws) || ...
        isempty(A.stf_draws) || isempty(A.S0)
    same = false;
    what = ' (a compared quantity is empty: the test compares nothing there)';
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

% TEST_SHARED_BUILD_SPREAD verify_shared_build on the conditions of a preset with a depression spread.
%
% run_memory_capacity pairs its trials by building every condition on one
% seed pair and asserting, with verify_shared_build, that everything but
% tau_a and synapse_config is identical. On 2026-11-12 the first
% memory-capacity run of a preset with std_spread failed there: std_clipped
% (how many neurons the chain's ordering rule moved) differs between a chain
% condition and the others, and was compared. No test had gone through that
% check with a spread: the exact chain tier's memory-capacity case calls the
% reservoir directly.
%
% Checks, on the mu7chain memory-capacity preset built as reservoirs the way
% run_memory_capacity builds them (a fifth of the size, a short protocol):
%   1. the four conditions pass the check (std_clipped is derived from the
%      condition's routes and is skipped);
%   2. the case is not vacuous: std_clipped does differ between conditions,
%      and the draws are equal;
%   3. the DRAWS are still compared: a condition built with another
%      std_seed, so that its neurons have other depression parameters, is
%      refused. That is what keeps a trial paired.
% Asserts; prints one PASS line.

preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
[d, cls, cnd] = srnn_param_preset(preset);
assert(strcmp(cls, 'SRNNCellTypePairs') && isequal(d.std_spread, [0.25 0.25]), ...
    'the mu7chain memory-capacity preset must carry std_spread');
d.n = 100; d.indegree = 20; d.F_ref_n = 100; d.F_ref_indegree = 20;   % small, same density
fs = 200;
base = [struct2namevalue(d), {'rng_seeds', [11 12], 'verbose', 'near-none', 'fs', fs, ...
    'ode_solver', 'sra1', 'T_wash', 2 * fs, 'T_train', 6 * fs, 'T_test', 3 * fs, 'd_max', fs}];

n_cond = numel(cnd);
esn = cell(1, n_cond);
names = cell(1, n_cond);
for i = 1:n_cond
    names{i} = cnd{i}.name;
    cargs = struct2namevalue(rmfield(cnd{i}, 'name'));
    esn{i} = SRNN_ESN_reservoir(base{:}, cargs{:});
    esn{i}.build();
end

%% 1. The four conditions of one trial pass the fairness check
verify_shared_build(esn, {'tau_a', 'synapse_config'}, {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});

%% 2. Not vacuous: std_clipped differs, the draws do not
clipped = cellfun(@(e) sum(e.std_clipped), esn);
i_chain = find(strcmp(names, 'sfa3_stdc2'));
assert(isscalar(i_chain) && clipped(i_chain) > 0, ...
    'the chain condition must have neurons moved by the ordering rule here (clipped: %s); pick another seed', mat2str(clipped));
assert(all(clipped(setdiff(1:n_cond, i_chain)) == 0), 'only the chain condition has an ordering rule');
assert(~isempty(esn{1}.std_draws), 'the draws must exist in every condition, with or without routes');
for i = 2:n_cond
    assert(isequal(esn{i}.std_draws, esn{1}.std_draws), 'std_draws differ between %s and %s', names{1}, names{i});
end

%% 3. Other draws in one condition are refused
cargs = struct2namevalue(rmfield(cnd{i_chain}, 'name'));
odd = SRNN_ESN_reservoir(base{:}, cargs{:}, 'std_seed', 987654);
odd.build();
tampered = esn; tampered{i_chain} = odd;
threw = false;
try
    verify_shared_build(tampered, {'tau_a', 'synapse_config'}, {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});
catch err
    threw = strcmp(err.identifier, 'verify_shared_build:Mismatch') && ...
        (contains(err.message, 'std_draws') || contains(err.message, 'std_seed'));
    if ~threw, rethrow(err); end
end
assert(threw, 'a condition with other depression draws must fail the fairness check');

fprintf('test_shared_build_spread: PASS\n');

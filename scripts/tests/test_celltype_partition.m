% test_celltype_partition.m - The model's type partition must agree with W's.
%
% SRNNCellTypePairs indexes its whole state vector by n_per_type/type_indices,
% while W is built by RMTBlocks, which partitions (n, f) itself. If the two
% partitions disagree, the model applies type-q parameters to a neuron whose
% outgoing synapses W built as type-q' -- a silent Dale's-law violation with no
% error and no warning.
%
% That is not hypothetical: until 2026-09-02 the model used
% RMTCellTypes.allocate_counts (largest remainder) while RMTBlocks used
% cumulative rounding. The two agree for every C = 2 model, which is why no
% paper preset was ever affected and why this went unnoticed -- and disagree for
% ~22% of C >= 3 draws.
%
% The sizes below are chosen to HIT that divergence. Every C >= 3 site that
% existed in the repo when this was found (4 examples, 5 tests) used round
% fractions and escaped by luck; none of them would have caught it.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs, RMTBlocks/set_types

all_passed = true;

fprintf('\n=== test_celltype_partition ===\n');

%% ------------------------------------------------------------------------
%  1. The two partition algorithms, compared directly.
%     Cumulative rounding is what RMTBlocks uses and what the model must use.
%% ------------------------------------------------------------------------
fprintf('\n--- 1. Partition agreement between model and RMTBlocks ---\n');

partition_cases = { ...
    300, [0.5 0.5];                    % the paper's shape
    300, [0.7 0.1 0.1 0.1];
    100, [1/3 1/3 1/3];                % DIVERGES under largest remainder
     10, [1/3 1/3 1/3];                % DIVERGES under largest remainder
     50, [0.25 0.25 0.25 0.25];        % DIVERGES under largest remainder
    128, [0.261 0.548 0.191]; ...
};

for k = 1:size(partition_cases, 1)
    n = partition_cases{k, 1};
    f = partition_cases{k, 2};
    C = numel(f);
    names = arrayfun(@(q) sprintf('T%d', q), 1:C, 'UniformOutput', false);

    m = SRNNCellTypePairs('n', n, 'n_cellTypes', C, 'cell_type_names', names, ...
        'f', f, 'mu_tilde_relative', zeros(C, C), ...
        'sigma_tilde_relative', ones(C, C), 'indegree', min(10, n - 1));

    rmt = RMTBlocks(n);
    rmt.set_types(f, zeros(C, C), ones(C, C));

    ok = isequal(m.n_per_type, rmt.n_per_type);
    all_passed = report(ok, sprintf('n=%3d f=%-22s model %s == RMTBlocks %s', ...
        n, mat2str(f, 3), mat2str(m.n_per_type), mat2str(rmt.n_per_type))) && all_passed;

    % type_indices must be the contiguous ranges implied by those counts
    ok = isequal(m.type_indices, rmt.type_indices);
    all_passed = report(ok, sprintf('n=%3d f=%-22s type_indices match', ...
        n, mat2str(f, 3))) && all_passed;
end

%% ------------------------------------------------------------------------
%  2. The sharp test: does W actually agree, neuron by neuron?
%     Dale's law makes W's columns sign-constant within a presynaptic type,
%     so the column signs recover W's TRUE partition independently.
%% ------------------------------------------------------------------------
fprintf('\n--- 2. W column signs vs the model''s type_indices ---\n');

sign_cases = { ...
    100, [1/3 1/3 1/3];
     10, [1/3 1/3 1/3];
     50, [0.25 0.25 0.25 0.25]; ...
};

for k = 1:size(sign_cases, 1)
    n = sign_cases{k, 1};
    f = sign_cases{k, 2};
    C = numel(f);
    names = [{'E'}, arrayfun(@(q) sprintf('I%d', q), 2:C, 'UniformOutput', false)];

    % Type 1 excitatory, all others inhibitory -> column sign identifies type 1.
    mu = repmat([3, -4 * ones(1, C - 1)], C, 1);

    m = SRNNCellTypePairs('n', n, 'n_cellTypes', C, 'cell_type_names', names, ...
        'f', f, 'mu_tilde_relative', mu, 'sigma_tilde_relative', ones(C, C), ...
        'indegree', min(20, n - 1));
    evalc('m.build()');                       % silence the build chatter

    col_is_exc = false(1, n);
    for j = 1:n
        v = nonzeros(m.W(:, j));
        col_is_exc(j) = ~isempty(v) && sum(v) > 0;
    end

    model_says_exc = false(1, n);
    model_says_exc(m.type_indices{1}) = true;

    ok = isequal(col_is_exc, model_says_exc);
    if ~ok
        bad = find(col_is_exc ~= model_says_exc);
        detail = sprintf(' -- disagree at neuron(s) %s', mat2str(bad(1:min(5, end))));
    else
        detail = '';
    end
    all_passed = report(ok, sprintf('n=%3d f=%-22s W signs match type_indices%s', ...
        n, mat2str(f, 3), detail)) && all_passed;
end

%% ------------------------------------------------------------------------
%  3. The empty-type guard. RMTBlocks permits a zero count (legal
%     connectivity); the model must reject it (a type with no neurons but its
%     own tau_a and synapse routes is a configuration error).
%% ------------------------------------------------------------------------
fprintf('\n--- 3. Empty-type guard ---\n');

threw = false;
try
    bad = SRNNCellTypePairs('n', 4, 'n_cellTypes', 3, ...
        'cell_type_names', {'A', 'B', 'C'}, 'f', [0.9 0.05 0.05], ...
        'mu_tilde_relative', zeros(3, 3), 'sigma_tilde_relative', ones(3, 3), ...
        'indegree', 2);
    evalc('bad.build()');
catch ME
    threw = strcmp(ME.identifier, 'SRNNCellTypePairs:InvalidParams') && ...
        contains(ME.message, 'at least one neuron');
end
all_passed = report(threw, 'n=4, f=[0.9 0.05 0.05] rejected with the empty-type error') ...
    && all_passed;

% ...and a valid partition of the same shape is accepted.
threw = false;
try
    good = SRNNCellTypePairs('n', 60, 'n_cellTypes', 3, ...
        'cell_type_names', {'A', 'B', 'C'}, 'f', [0.9 0.05 0.05], ...
        'mu_tilde_relative', zeros(3, 3), 'sigma_tilde_relative', ones(3, 3), ...
        'indegree', 10);
    evalc('good.build()');
catch
    threw = true;
end
all_passed = report(~threw, 'n=60, f=[0.9 0.05 0.05] accepted (counts [54 3 3])') ...
    && all_passed;

%% ------------------------------------------------------------------------
fprintf('\n');
if all_passed
    fprintf('test_celltype_partition: ALL TESTS PASSED\n');
else
    fprintf(2, 'test_celltype_partition: FAILURES ABOVE\n');
end

function ok = report(ok, msg)
if ok
    fprintf('  PASS  %s\n', msg);
else
    fprintf(2, '  FAIL  %s\n', msg);
end
end

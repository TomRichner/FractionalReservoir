% test_RMTCellTypes.m - Verify arbitrary-cell-type RMT connectivity.

clear; clc;

%% Largest-remainder allocation and deterministic population columns
rng(11);
rmt = RMTCellTypes(11, 1, [0.5 0.3 0.2], [0.2 -0.4 0.1], [0 0 0]);
assert(isequal(rmt.n_per_type, [6 3 2]));
assert(issparse(rmt.W));
assert(all(rmt.W(:, rmt.type_indices{1}) == 0.2, 'all'));
assert(all(rmt.W(:, rmt.type_indices{2}) == -0.4, 'all'));
assert(all(rmt.W(:, rmt.type_indices{3}) == 0.1, 'all'));

%% Reproducibility and sparse storage
rng(29);
rmt1 = RMTCellTypes(40, 0.2, [0.6 0.25 0.15], ...
    [0.1 -0.2 -0.05], [0.03 0.04 0.02]);
rng(29);
rmt2 = RMTCellTypes(40, 0.2, [0.6 0.25 0.15], ...
    [0.1 -0.2 -0.05], [0.03 0.04 0.02]);
assert(isequal(rmt1.W, rmt2.W));
assert(issparse(rmt1.W));
assert(abs(nnz(rmt1.W) / 40^2 - 0.2) < 0.05);

%% Validation
assert_throws(@() RMTCellTypes(10, 0.2, [0.6 0.3], [1 2], [1 1]));
assert_throws(@() RMTCellTypes(10, 0, [0.5 0.5], [1 -1], [1 1]));
assert_throws(@() RMTCellTypes(2, 1, [0.98 0.01 0.01], [1 1 1], [1 1 1]));

fprintf('test_RMTCellTypes: ALL TESTS PASSED\n');

function assert_throws(f)
threw = false;
try
    f();
catch
    threw = true;
end
assert(threw, 'Expected operation to throw an error.');
end

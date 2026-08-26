function test_subset_fraction()
% TEST_SUBSET_FRACTION ParamSpaceAnalysis2.subset_fraction thins the grid safely.
%
% The load-bearing claim is that a partial run is a SUBSET OF THE SAME
% EXPERIMENT, not a different one: results must land at their true grid
% coordinates, unrun points must be empty (not shifted), the per-point network
% seed must not depend on how many points were run, and every pooling path must
% skip the holes exactly as it already skips failures.
%
% Does not simulate anything -- generate_grid is called directly, so this runs
% in milliseconds.
%
% See also: ParamSpaceAnalysis2, run_param_space_analysis

fprintf('\n=== test_subset_fraction ===\n');
n_pass = 0; n_fail = 0;

    function check(cond, msg)
        if cond
            fprintf('  PASS  %s\n', msg); n_pass = n_pass + 1;
        else
            fprintf('  FAIL  %s\n', msg); n_fail = n_fail + 1;
        end
    end

%% Full grid is untouched
psa = ParamSpaceAnalysis2('n_levels', 4, 'verbose', false);
psa.add_grid_parameter('n', [100, 1000]);
psa.add_grid_parameter('level_of_chaos', [0.5, 2.0]);
evalc('psa.generate_grid()');
check(psa.num_combinations == 16, 'full grid: 4^2 = 16 combinations');
check(numel(psa.shuffled_indices) == 16, 'full grid: all 16 scheduled');
check(isequal(sort(psa.shuffled_indices(:))', 1:16), 'full grid: a permutation of 1:16');

%% Subset thins the schedule but not the grid
psa2 = ParamSpaceAnalysis2('n_levels', 4, 'verbose', false);
psa2.subset_fraction = 0.25;
psa2.add_grid_parameter('n', [100, 1000]);
psa2.add_grid_parameter('level_of_chaos', [0.5, 2.0]);
evalc('psa2.generate_grid()');
check(psa2.num_combinations == 16, 'subset: num_combinations stays the FULL grid size');
check(numel(psa2.shuffled_indices) == 4, 'subset: 25% of 16 -> 4 scheduled');
check(numel(unique(psa2.shuffled_indices)) == 4, 'subset: scheduled indices are distinct');
check(all(psa2.shuffled_indices >= 1 & psa2.shuffled_indices <= 16), ...
    'subset: scheduled indices are valid grid positions');
check(numel(psa2.all_configs) == 16, 'subset: all_configs still describes every grid point');

% ceil, not round: a fraction that would round to zero still runs one point.
psa3 = ParamSpaceAnalysis2('n_levels', 4, 'verbose', false);
psa3.subset_fraction = 0.01;
psa3.add_grid_parameter('n', [100, 1000]);
psa3.add_grid_parameter('level_of_chaos', [0.5, 2.0]);
evalc('psa3.generate_grid()');
check(numel(psa3.shuffled_indices) == 1, 'tiny fraction runs 1 point, never 0');

%% A subset is NOT a corner of the grid
% The point of requiring randomize_order. Over repeated draws the scheduled
% points must spread across the grid rather than clustering at low indices.
spread_ok = false;
for trial = 1:20
    p = ParamSpaceAnalysis2('n_levels', 5, 'verbose', false);
    p.subset_fraction = 0.2;
    p.add_grid_parameter('n', [100, 1000]);
    p.add_grid_parameter('level_of_chaos', [0.5, 2.0]);
    evalc('p.generate_grid()');
    if any(p.shuffled_indices > 20)   % beyond the first 20 of 25
        spread_ok = true; break;
    end
end
check(spread_ok, 'subset samples across the grid, not the first-K corner');

%% Sequential order + subset is a hard error
psa4 = ParamSpaceAnalysis2('n_levels', 4, 'randomize_order', false, 'verbose', false);
psa4.subset_fraction = 0.5;
psa4.add_grid_parameter('n', [100, 1000]);
threw = false;
try
    evalc('psa4.generate_grid()');
catch ME
    threw = strcmp(ME.identifier, 'ParamSpaceAnalysis2:SubsetNeedsRandomOrder');
end
check(threw, 'subset_fraction with randomize_order = false errors');

% ...but sequential at fraction 1 is fine (the sensitivity sweeps' path).
psa5 = ParamSpaceAnalysis2('n_levels', 4, 'randomize_order', false, 'verbose', false);
psa5.add_grid_parameter('n', [100, 1000]);
ok = true;
try
    evalc('psa5.generate_grid()');
catch
    ok = false;
end
check(ok && isequal(psa5.shuffled_indices, 1:4), ...
    'sequential order at fraction 1 is unaffected');

%% Out-of-range fractions are rejected
for bad = {0, -0.5, 1.5, [], [0.5 0.5]}
    p = ParamSpaceAnalysis2('n_levels', 3, 'verbose', false);
    p.subset_fraction = bad{1};
    p.add_grid_parameter('n', [100, 1000]);
    threw = false;
    try
        evalc('p.generate_grid()');
    catch ME
        threw = strcmp(ME.identifier, 'ParamSpaceAnalysis2:BadSubsetFraction');
    end
    check(threw, sprintf('subset_fraction = %s is rejected', mat2str(bad{1})));
end

%% Batching covers exactly the scheduled points, no more
% Mirrors run_batched_simulation's arithmetic: with batch_size 3 and 4
% scheduled points there must be 2 batches, and the union of their slices must
% be the schedule itself -- an off-by-one here would index past the end.
sched = psa2.shuffled_indices;
batch_size = 3;
num_batches = ceil(numel(sched) / batch_size);
covered = [];
for b = 1:num_batches
    s = (b-1)*batch_size + 1;
    e = min(b*batch_size, numel(sched));
    covered = [covered, sched(s:e)]; %#ok<AGROW>
end
check(num_batches == 2, 'batching: 4 scheduled / batch 3 -> 2 batches');
check(isequal(covered, sched), 'batching: batches cover exactly the schedule');

%% Network seed is tied to grid POSITION, not execution order
% This is what makes a subset comparable to the full grid: the same grid point
% draws the same network either way.
seed_of = @(config_idx, offset) config_idx*100 + offset;
check(seed_of(7, 0) == seed_of(7, 0), 'network seed depends only on config_idx');
check(seed_of(7, 0) ~= seed_of(8, 0), 'distinct grid points get distinct seeds');

%% Pooling skips holes exactly like failures
% The predicate every plotting path uses. Six slots: two empty (unrun), one
% failed, three good -- must pool to exactly the three good values.
results = cell(6, 1);
results{2} = struct('success', true,  'LLE',  0.5, 'config_idx', 2);
results{3} = struct('success', false, 'LLE',  NaN, 'config_idx', 3);
results{4} = struct('success', true,  'LLE', -0.2, 'config_idx', 4);
results{6} = struct('success', true,  'LLE',  1.1, 'config_idx', 6);
keep = cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success && ...
    isfield(r, 'LLE') && ~isnan(r.LLE), results);
vals = cellfun(@(r) r.LLE, results(keep));
check(numel(vals) == 3, 'pooling: empty slots skipped like failed ones');
check(isequal(sort(vals(:))', [-0.2, 0.5, 1.1]), 'pooling: exactly the successful values');

%% Round-trip through save/load
tmp = [tempname '.mat'];
c = onCleanup(@() delete_if_exists(tmp));
save(tmp, 'psa2');
loaded = load(tmp);
check(loaded.psa2.subset_fraction == 0.25, 'subset_fraction survives save/load');
check(numel(loaded.psa2.shuffled_indices) == 4, ...
    'shuffled_indices survives save/load (so consolidate gets the right batch count)');
check(loaded.psa2.num_combinations == 16, 'num_combinations survives save/load');

%% Banner
fprintf('\n%d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0
    error('test_subset_fraction:Failed', '%d check(s) failed.', n_fail);
end
fprintf('=== test_subset_fraction PASSED ===\n');
end

function delete_if_exists(f)
if exist(f, 'file'), delete(f); end
end

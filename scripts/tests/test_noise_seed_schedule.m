% test_noise_seed_schedule.m - SRNNCellTypePairs.noise_seed_schedule, the noise
% seed that changes at a given time, and wiener_increments behind it.
%
% Checks, on a small three-type network with SFA and depression (sra1):
%   1. An empty schedule draws the tensors build_noise always drew
%      (rng(seed, 'twister'); randn; randn), written out here independently.
%   2. A schedule naming the base seed reproduces the unscheduled run exactly.
%   3. Two copies that differ only after the switch are bit-identical up to
%      and including the state at the switch time, and differ at the very
%      next sample and at the end of the run.
%   4. Two copies with the same schedule are bit-identical throughout.
%   5. The splice is at the same absolute time index: after the switch the
%      tensors hold the scheduled seed's own columns; with two rows, each
%      segment holds its seed's columns.
%   6. The global random stream is untouched by the draw.
%   7. A schedule at sigma_u_noise = 0 changes nothing.
%   8. Malformed schedules error :InvalidNoiseSchedule: wrong shape, times not
%      increasing, a non-integer seed (all at construction), a time off the
%      grid, at the first sample, or at or beyond the last (at run).
%
% A few seconds. Prints PASS/FAIL per check and a final banner. Assumes
% setup_paths.
%
% See also: wiener_increments, check_noise_seed_schedule, SRNNCellTypePairs,
%           test_SRNNCellTypePairs_noise

fprintf('=== Testing noise_seed_schedule ===\n\n');
all_passed = true;

base_seed = 31; t_switch = 1; fs = 200; T = 2;
i_switch = round(t_switch * fs) + 1;            % sample index of the switch time

%% 1. Empty schedule: the historical draw
m = pairs('noise_seed', base_seed);
m.build(); m.build_noise_for_test();
n_steps = T * fs;
state = rng; rng(base_seed, 'twister');
ref1 = randn(m.n, n_steps); ref2 = randn(m.n, n_steps);
rng(state);
all_passed = check('empty schedule draws the historical tensors', ...
    isequal(m.noise_increments.xi1, ref1) && isequal(m.noise_increments.xi2, ref2)) && all_passed;

%% 2. A schedule naming the base seed changes nothing
plain = run_traj(pairs('noise_seed', base_seed));
same  = run_traj(pairs('noise_seed', base_seed, 'noise_seed_schedule', [t_switch, base_seed]));
all_passed = check('a schedule naming the base seed reproduces the unscheduled run', ...
    isequal(plain, same)) && all_passed;

%% 3. Copies that differ only after the switch
A = run_traj(pairs('noise_seed', base_seed, 'noise_seed_schedule', [t_switch, 101]));
B = run_traj(pairs('noise_seed', base_seed, 'noise_seed_schedule', [t_switch, 102]));
all_passed = check('pre-switch states bit-identical, the switch sample included', ...
    isequal(A(1:i_switch, :), B(1:i_switch, :))) && all_passed;
all_passed = check('...and equal to the unscheduled run there', ...
    isequal(A(1:i_switch, :), plain(1:i_switch, :))) && all_passed;
all_passed = check('the first sample after the switch differs', ...
    ~isequal(A(i_switch + 1, :), B(i_switch + 1, :))) && all_passed;
all_passed = check('the final states differ', ~isequal(A(end, :), B(end, :))) && all_passed;

%% 4. The same schedule reproduces
A2 = run_traj(pairs('noise_seed', base_seed, 'noise_seed_schedule', [t_switch, 101]));
all_passed = check('the same schedule gives the same trajectory', isequal(A, A2)) && all_passed;

%% 5. Same absolute time index
s = pairs('noise_seed', base_seed, 'noise_seed_schedule', [t_switch, 101]);
s.build(); s.build_noise_for_test();
state = rng; rng(101, 'twister');
alt1 = randn(s.n, n_steps); alt2 = randn(s.n, n_steps);
rng(state);
c = i_switch;                                   % first column of the new seed
all_passed = check('columns before the switch are the base seed''s', ...
    isequal(s.noise_increments.xi1(:, 1:c-1), ref1(:, 1:c-1)) && ...
    isequal(s.noise_increments.xi2(:, 1:c-1), ref2(:, 1:c-1))) && all_passed;
all_passed = check('columns from the switch on are the scheduled seed''s own', ...
    isequal(s.noise_increments.xi1(:, c:end), alt1(:, c:end)) && ...
    isequal(s.noise_increments.xi2(:, c:end), alt2(:, c:end))) && all_passed;

s2 = pairs('noise_seed', base_seed, 'noise_seed_schedule', [0.5, 101; 1.5, base_seed]);
s2.build(); s2.build_noise_for_test();
c1 = round(0.5 * fs) + 1; c2 = round(1.5 * fs) + 1;
all_passed = check('two rows: base, scheduled, base again, each on its own columns', ...
    isequal(s2.noise_increments.xi1(:, 1:c1-1), ref1(:, 1:c1-1)) && ...
    isequal(s2.noise_increments.xi1(:, c1:c2-1), alt1(:, c1:c2-1)) && ...
    isequal(s2.noise_increments.xi1(:, c2:end), ref1(:, c2:end)) && ...
    isequal(s2.noise_increments.xi2(:, c1:c2-1), alt2(:, c1:c2-1))) && all_passed;

%% 6. The global stream is untouched
rng(5, 'twister'); before = rng;
s.build_noise_for_test();
after = rng;
all_passed = check('the draw restores the global random stream', isequal(before, after)) && all_passed;

%% 7. Ignored without noise
d0 = run_traj(pairs('sigma_u_noise', 0, 'ode_solver', 'rk4'));
d1 = run_traj(pairs('sigma_u_noise', 0, 'ode_solver', 'rk4', 'noise_seed_schedule', [t_switch, 101]));
all_passed = check('a schedule at sigma_u_noise = 0 changes nothing', isequal(d0, d1)) && all_passed;

%% 8. Validation
id = 'SRNNCellTypePairs:InvalidNoiseSchedule';
all_passed = check('a 1 x 3 schedule errors at construction', ...
    throws_id(@() pairs('noise_seed_schedule', [1 2 3]), id)) && all_passed;
all_passed = check('times not increasing error at construction', ...
    throws_id(@() pairs('noise_seed_schedule', [1 5; 1 6]), id)) && all_passed;
all_passed = check('a non-integer seed errors at construction', ...
    throws_id(@() pairs('noise_seed_schedule', [1 5.5]), id)) && all_passed;
all_passed = check('a time off the grid errors at run', ...
    throws_id(@() run_traj(pairs('noise_seed_schedule', [1.0012 5])), id)) && all_passed;
all_passed = check('a switch at the first sample errors at run', ...
    throws_id(@() run_traj(pairs('noise_seed_schedule', [0 5])), id)) && all_passed;
all_passed = check('a switch at the last sample errors at run', ...
    throws_id(@() run_traj(pairs('noise_seed_schedule', [T 5])), id)) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function p = pairs(varargin)
% Three cell types, SFA on two and depression on one route, so the splice is
% checked on a state vector with more than the dendritic block.
syn = struct('A', struct('B', struct('std', struct('tau_rec', 0.5, 'tau_rel', 0.2))));
defaults = {'n', 24, 'indegree', 12, 'n_cellTypes', 3, ...
    'cell_type_names', {'A', 'B', 'C'}, 'f', [0.4 0.35 0.25], ...
    'mu_tilde_relative', [3 -4 -2], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {log_ladder(0.25, 10, 2), log_ladder(0.25, 10, 1), []}, ...
    'synapse_config', syn, 'T_range', [0 2], 'fs', 200, ...
    'sigma_u_noise', 0.3, 'ode_solver', 'sra1', ...
    'lya_method', 'none', 'store_full_state', true};
p = SRNNPairsTestAccess(defaults{:}, varargin{:});
end

function S = run_traj(p)
p.build();
evalc('p.run();');
S = p.S_out;
end

function tf = throws_id(fcn, id)
tf = false;
try
    fcn();
catch err
    tf = strcmp(err.identifier, id);
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

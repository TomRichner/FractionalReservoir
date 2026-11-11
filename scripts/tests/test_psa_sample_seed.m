% test_psa_sample_seed.m
% ParamSpaceAnalysis2.sample_seed: the optional seed for the random execution
% order, and so for which grid points a subset_fraction < 1 run simulates
% (FD-016).
%
% What is asserted:
%   1. two objects with the same seed draw identical shuffled_indices, and
%      the draw is the documented one (a private mt19937ar stream);
%   2. a different seed draws a different sample;
%   3. unset, two objects draw different samples (today's behaviour);
%   4. the seeded draw leaves the global random stream as it found it, and
%      the unset draw reseeds it, as it always did;
%   5. a bad seed is rejected, and a seed with sequential order warns;
%   6. two seeded RUNS give identical results at identical grid points;
%   7. the seed survives save and load, a loaded object redraws the same
%      order, and the summary records the seed when set and has no such
%      field when unset;
%   8. same_config does not compare it;
%   9. the config field reaches the PSA through resolve_run_context and
%      run_param_space_analysis.
%
% Tiny networks, short runs, serial sweeps; the one run_param_space_analysis
% call uses whatever pool the session holds (serial without one).
%
% See also: ParamSpaceAnalysis2, run_param_space_analysis, resolve_run_context

fprintf('=== test_psa_sample_seed ===\n');

P = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
scratch = fullfile(tempdir, sprintf('psa_sample_seed_%s', char(java.util.UUID.randomUUID())));
mkdir(scratch);
figs_before = findall(groot, 'Type', 'figure');
rng_at_start = rng;

ctx0 = resolve_run_context('param_space', 'preset_name', P, 'run_mode', 'fast', 'verbose', 'near-none');
assert(isfield(ctx0, 'param_space_seed') && isempty(ctx0.param_space_seed), ...
    'resolve_run_context must carry an empty param_space_seed by default');

%% 1-3. The draw itself (no simulation)
a = make_psa(ctx0, scratch); a.sample_seed = 11; a.generate_grid();
b = make_psa(ctx0, scratch); b.sample_seed = 11; b.generate_grid();
c = make_psa(ctx0, scratch); c.sample_seed = 12; c.generate_grid();
n_full = a.num_combinations;
n_run  = numel(a.shuffled_indices);
assert(n_full == 16 && n_run == 4, 'expected a 4-of-16 sample, got %d of %d', n_run, n_full);
assert(isequal(a.shuffled_indices, b.shuffled_indices), 'same seed, different sample');
assert(~isequal(a.shuffled_indices, c.shuffled_indices), 'different seeds drew the same sample');
expected = randperm(RandStream('mt19937ar', 'Seed', 11), n_full);
assert(isequal(a.shuffled_indices, expected(1:n_run)), ...
    'the seeded draw is not randperm on a private mt19937ar stream');
a.generate_grid();
assert(isequal(a.shuffled_indices, b.shuffled_indices), 'a redraw with the same seed changed the order');

% Unset: drawn from the clock. Two draws of 4 ordered points out of 16 agree
% by chance once in 43680, and rng('shuffle') called twice within its clock
% resolution would also agree; neither is allowed for here, and the pause
% makes the second impossible in practice.
u1 = make_psa(ctx0, scratch); u1.generate_grid();
pause(0.05);
u2 = make_psa(ctx0, scratch); u2.generate_grid();
assert(isempty(u1.sample_seed), 'sample_seed must default to empty');
assert(~isequal(u1.shuffled_indices, u2.shuffled_indices), ...
    'two unseeded objects drew the same sample (1 in 43680 by chance; rerun once before believing it)');
fprintf('  draw: same seed same sample, other seed other sample, unset differs\n');

%% 4. The global stream
rng(5, 'twister'); before = rng;
g = make_psa(ctx0, scratch); g.sample_seed = 3; g.generate_grid();
after = rng;
assert(isequal(before, after), 'the seeded draw changed the global random stream');
rng(5, 'twister'); r_ref = rand(1, 3);
rng(5, 'twister'); g.generate_grid(); r_got = rand(1, 3);
assert(isequal(r_ref, r_got), 'the seeded draw consumed numbers from the global stream');
% The generator is pinned: the same order under another global generator.
rng(5, 'philox'); h = make_psa(ctx0, scratch); h.sample_seed = 3; h.generate_grid();
assert(isequal(h.shuffled_indices, g.shuffled_indices), 'the seeded draw depends on the global generator');
% Unset keeps today's side effect: rng('shuffle') reseeds the global stream.
rng(5, 'twister'); before = rng;
u1.generate_grid();
assert(~isequal(before, rng), 'the unset draw no longer reseeds the global stream: the default path changed');
fprintf('  global stream: untouched when seeded, reseeded when unset\n');

%% 5. Validation
bad = {-1, 1.5, NaN, Inf, [1 2], 2^32, 'abc', {1}};
for k = 1:numel(bad)
    v = make_psa(ctx0, scratch); v.sample_seed = bad{k};
    assert(threw(@() v.generate_grid(), 'ParamSpaceAnalysis2:BadSampleSeed'), ...
        'bad sample_seed #%d was accepted by generate_grid', k);
end
% run() rejects it before it creates an output directory
v = make_psa(ctx0, fullfile(scratch, 'never_created')); v.sample_seed = -1;
assert(threw(@() v.run(), 'ParamSpaceAnalysis2:BadSampleSeed'), 'run() accepted a bad sample_seed');
assert(~isfolder(fullfile(scratch, 'never_created')), 'a rejected seed left an output directory behind');
% sequential order: the seed can do nothing, and says so
w = make_psa(ctx0, scratch); w.randomize_order = false; w.subset_fraction = 1; w.sample_seed = 4;
lastwarn('', '');
ws = warning('off', 'ParamSpaceAnalysis2:SampleSeedUnused');
w.generate_grid();
warning(ws);
[~, warn_id] = lastwarn;
assert(strcmp(warn_id, 'ParamSpaceAnalysis2:SampleSeedUnused'), 'no warning for a seed with sequential order');
assert(isequal(w.shuffled_indices, 1:w.num_combinations), 'sequential order was not sequential');
fprintf('  validation: %d bad values rejected, sequential order warns\n', numel(bad));

%% 6. Two seeded runs: same points, same numbers
r1 = make_psa(ctx0, fullfile(scratch, 'run1')); r1.sample_seed = 11; r1.run();
r2 = make_psa(ctx0, fullfile(scratch, 'run2')); r2.sample_seed = 11; r2.run();
assert(isequal(r1.shuffled_indices, a.shuffled_indices) && isequal(r2.shuffled_indices, a.shuffled_indices), ...
    'run() did not use the seeded order');
cond_names = cellfun(@(s) s.name, r1.conditions, 'UniformOutput', false);
n_ok = 0;
for i = 1:numel(cond_names)
    x1 = cellfun(@drop_timing, r1.results.(cond_names{i}), 'UniformOutput', false);
    x2 = cellfun(@drop_timing, r2.results.(cond_names{i}), 'UniformOutput', false);
    filled = find(cellfun(@isstruct, x1(:)'));
    assert(isequal(sort(filled), sort(r1.shuffled_indices)), ...
        '%s: results are not at the sampled grid points', cond_names{i});
    assert(all(cellfun(@(s) s.success, x1(filled))), '%s: a job failed', cond_names{i});
    assert(isequaln(x1, x2), '%s: two runs with the same seed gave different results', cond_names{i});
    n_ok = n_ok + numel(filled);
end
fprintf('  runs: %d jobs identical between two runs with the same seed\n', n_ok);

%% 7. Save, load, summary
L = load(fullfile(r1.output_dir, 'psa_object.mat'));
assert(isequal(L.psa.sample_seed, 11), 'sample_seed did not survive save and load');
assert(isequal(L.psa.shuffled_indices, r1.shuffled_indices), 'shuffled_indices did not survive save and load');
L.psa.generate_grid();
assert(isequal(L.psa.shuffled_indices, r1.shuffled_indices), 'a loaded object redrew a different order');
S = load(fullfile(r1.output_dir, 'param_space_summary.mat'));
assert(isfield(S.summary_data, 'sample_seed') && isequal(S.summary_data.sample_seed, 11), ...
    'the summary does not record the seed');

% Unset: the run still works, the object loads with an empty seed, and the
% summary has NO sample_seed field (its field set is what it was before the
% seed existed, which the full regression tier depends on).
r3 = make_psa(ctx0, fullfile(scratch, 'run3')); r3.run();
L3 = load(fullfile(r3.output_dir, 'psa_object.mat'));
assert(isempty(L3.psa.sample_seed), 'an unseeded object loaded with a seed');
S3 = load(fullfile(r3.output_dir, 'param_space_summary.mat'));
assert(~isfield(S3.summary_data, 'sample_seed'), 'an unseeded summary gained a sample_seed field');
% An object saved before the property existed loads as unset.
old = saveobj(r3); old = rmfield(old, 'sample_seed');
assert(isempty(ParamSpaceAnalysis2.loadobj(old).sample_seed), 'a pre-seed object did not load as unset');
fprintf('  save/load: seed kept, same order redrawn, summary records it only when set\n');

%% 8. same_config ignores it
[tf, reason] = r1.same_config(r3);
assert(tf, 'same_config refused two runs differing only in sample_seed: %s', reason);

%% 9. From the config to the PSA, through run_param_space_analysis
seed = 20261111;
ctx = resolve_run_context('param_space', 'preset_name', P, 'run_mode', 'fast', ...
    'verbose', 'near-none', 'output_dir', fullfile(scratch, 'stage'), 'param_space_seed', seed);
assert(isequal(ctx.param_space_seed, seed), 'resolve_run_context dropped param_space_seed');
ctx.n_levels = 2;                 % 2^7 = 128 grid points
ctx.cfg.n_grid_points = 2;        % of which two are simulated
ctx.model_defaults = short_run(ctx.model_defaults);
mkdir(ctx.output_dir);
stage_dir = run_param_space_analysis(ctx);
Ls = load(fullfile(stage_dir, 'psa_object.mat'));
assert(isequal(Ls.psa.sample_seed, seed), 'the seed did not reach the PSA through run_param_space_analysis');
expected = randperm(RandStream('mt19937ar', 'Seed', seed), Ls.psa.num_combinations);
assert(Ls.psa.num_combinations == 128 && isequal(Ls.psa.shuffled_indices, expected(1:2)), ...
    'the stage did not sample the seeded points');
Ss = load(fullfile(stage_dir, 'param_space_summary.mat'));
assert(isequal(Ss.summary_data.sample_seed, seed), 'the stage''s summary does not record the seed');
fprintf('  plumbing: ctx.param_space_seed reached the PSA and the summary\n');

%% Clean up
close(setdiff(findall(groot, 'Type', 'figure'), figs_before));
rng(rng_at_start);
rmdir(scratch, 's');
fprintf('PASS: test_psa_sample_seed\n');

%% ------------------------------------------------------------------------
function psa = make_psa(ctx, out_dir)
% A 4 x 4 grid of which a quarter is sampled, on a 40-neuron network.
psa = ParamSpaceAnalysis2('n_levels', 4, 'batch_size', 4, 'note', 'sample_seed', ...
    'verbose', 'near-none', 'use_parallel', false);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
d = short_run(ctx.model_defaults);
d.n = 40; d.indegree = 8;
psa.model_defaults = d;
psa.output_dir     = out_dir;
psa.set_conditions(ctx.conditions);
psa.add_grid_parameter('level_of_chaos', [0.5, 1.5]);
psa.add_grid_parameter(ctx.f_param, [0.4, 0.6]);
psa.subset_fraction = 0.25;
end

function d = short_run(d)
d.T_range = [0 2]; d.lya_T_interval = [1 2]; d.lya_warmup = 0.5;
end

function tf = threw(fn, id)
tf = false;
try
    fn();
catch err
    tf = strcmp(err.identifier, id);
end
end

function s = drop_timing(s)
% Wall-clock fields differ between identical runs and are not results.
if ~isstruct(s); return; end
timing = {'seconds', 'run_duration', 'elapsed', 'wall_time'};
s = rmfield(s, intersect(fieldnames(s), timing));
end

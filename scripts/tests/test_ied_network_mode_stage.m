%% test_ied_network_mode_stage.m
% The ied_stim bundle's network_mode stage and its appendix figure (Plan 2,
% package E): full_jacobian_network_mode, ied_network_mode_stage_plan,
% run_ied_network_mode_stage, fig_ied_network_mode.
%
% Checks:
%   1. full_jacobian_network_mode on a block-diagonal Jacobian whose
%      spectrum is known: the slow single-variable eigenvalue to the right
%      is passed over, the network mode found, its coherence 1 for a
%      population-uniform vector.
%   2. THE PLAN: sizes by mode (1 / 3 / 12 networks; 20 / 60 / 120 s;
%      3 / 6 / 6 discharges), the bundle's first and first stimulation
%      blocks; refused: an unknown knob, a seed not in the bundle, a window
%      longer than the block.
%   3. THE STAGE at short size (the modular network is fixed at 350
%      neurons, so short blocks instead: 20 s analysed, the mode every 2 s
%      over 4 s, one discharge at two lags), re-simulating the bundle's
%      blocks: index and file written, both blocks' series, Im >= 0, the
%      coherence a fraction, the event arrays sized by discharges x lags;
%      and against a bundle block file written by ied_stim_blocks with the
%      same settings, the re-simulation reproduces its mean E potential
%      (to the single-precision rounding of the saved signal).
%   4. THE FIGURE draws (PNG and _table.md), and errors on a run directory
%      without the stage.
%
% Asserts; prints one PASS line. About two minutes without a pool. Assumes
% setup_paths has run and a pool allocation is set (pool_set(0) will do).

t_all = tic;

%% 1. The dense network-mode rule
n = 4; is_x = [true(n, 1); false(n, 1)];
U = [ones(2, 1); zeros(2, 1)] / sqrt(2); U(:, 2) = [zeros(2, 1); ones(2, 1)] / sqrt(2);
Jxx = -10 * eye(n) + 0.5 * ones(n);                      % network mode: -10 + 2 = -8, along ones(4,1)
J = blkdiag(Jxx, diag([-0.05 -1 -2 -3]));                % slow single-variable rates, to the right
f = full_jacobian_network_mode(J, is_x, 0.1, U);
assert(abs(f.net - (-8)) < 1e-10 && abs(f.abscissa - (-0.05)) < 1e-12 && abs(f.net_coh - 1) < 1e-10);
assert(f.n_net_unstable == 0);

%% 2. The plan
sizes = {'essential', 1, 20, 3; 'fast', 3, 60, 6; 'production', 12, 120, 6};
for k = 1:size(sizes, 1)
    c = ied_stim_config(sizes{k, 1});
    p = ied_network_mode_stage_plan(c);
    assert(numel(p.runs) == sizes{k, 2} && p.window_s == min(sizes{k, 3}, c.sim.block_len) && ...
        p.n_events == sizes{k, 4} && isequal(p.blocks, [1, c.sim.n_no + 1]), 'sizes of %s', sizes{k, 1});
end
c = ied_stim_config('essential');
bad = {struct('foo', 1), struct('seeds', 5), struct('window_s', 1e4)};
for k = 1:numel(bad)
    cc = c; cc.network_mode = bad{k};
    assert(throws(@() ied_network_mode_stage_plan(cc), 'ied_network_mode_stage_plan:BadConfig'), 'refusal %d', k);
end

%% 3. The stage at short size, against the bundle's own block
run_dir = tempname; mkdir(run_dir);
cleanup = onCleanup(@() rmdir(run_dir, 's'));
c = ied_stim_config('essential', 'run_dir', run_dir);
c.sim.seeds = 42; c.sim.block_len = 20;
c.network_mode = struct('window_s', 4, 'full_dt', 2, 'n_events', 1, 'lags', [-0.1 0]);
sim_args = namedargs2cell(c.sim);
ied_stim_blocks('data_dir', run_dir, sim_args{:});
detail = run_ied_network_mode_stage(c);
assert(contains(detail, 'network_mode'));
L = load(fullfile(run_dir, 'network_mode', 'ied_netmode_index.mat'));
assert(height(L.runs) == 1 && L.runs.ok);
assert(L.runs.repro_max_diff < 1e-6, 'the stage re-simulates the bundle''s block: %g', L.runs.repro_max_diff);
M = load(fullfile(run_dir, 'network_mode', L.runs.file{1}));
N = M.N;
for b = 1:2
    assert(numel(N.block(b).t) == 2 && all(imag(N.block(b).net) >= 0));
    assert(all(N.block(b).net_coh >= 0 & N.block(b).net_coh <= 1 + 1e-9));
    assert(all(abs(sum(N.block(b).share, 2) - 1) < 1e-9));
end
assert(isequal(size(N.event.ev.net), [numel(N.event.t_event), 2]) && numel(N.event.t_ctrl) == numel(N.event.t_event));

%% 4. The figure
o = fig_ied_network_mode('data_dir', run_dir, 'out_dir', fullfile(run_dir, 'fig'));
assert(all(cellfun(@isfile, o.files)));
assert(throws(@() fig_ied_network_mode('data_dir', tempname, 'out_dir', fullfile(run_dir, 'x')), ''), ...
    'the figure errors without the stage');
fprintf('PASS test_ied_network_mode_stage (%.0f s)\n', toc(t_all));

function tf = throws(f, id)
tf = false;
try
    f();
catch err
    tf = isempty(id) || strcmp(err.identifier, id);
end
end

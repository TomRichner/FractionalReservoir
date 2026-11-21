%% test_outlier_structures_stage.m
% The opt-in stage 'outlier_structures' (FD-043): outlier_structures_stage_plan,
% ostruct_structural_check, run_outlier_structures_stage, ostruct_load,
% ostruct_mode_measures and the four figures, on tiny networks.
%
% Checks:
%   1. The plan in EVERY run mode, from mu7chain_config: the fifteen
%      networks in their groups, the three conditions, 1 seed at essential
%      and fast and 5 above; seeds outermost; the seed convention (new
%      weights per seed row for the free networks; for the two
%      fixed-weight networks rng_seeds fixed and seed row s >= 2 adding the
%      token -noise<seeds(s,1)>); the checks (one per free network and seed,
%      one per fixed network); the knobs taken; bad settings refused.
%   2. The noise token: 'reference-noise3' on weight seed [27 28] has the
%      weights of 'reference' on the same seed and another noise path.
%   3. The structural check at 80 neurons: an added real mode
%      (lrA1p5W0p0), an added rotating pair (lrA1p4W0p6) and a block-mean
%      network (EEx2p0) pass; a wrong expectation fails with its reason; a
%      network with no outlier intended passes as such.
%   4. END TO END, tiny: 80 neurons, 2 s settling, 2 s window, 3
%      calibration states; five networks (the three above, the local
%      outlier reference_w27, and IIx1p25, whose block-mean outlier lies
%      inside the bulk at this size and so FAILS the check and is skipped)
%      by two conditions, one seed. The two indexes, the skipped runs, the
%      check against the runs' own structural eigenvalue, the measures'
%      fields and one recomputed by hand.
%   5. Noise realizations of a fixed-weight network: two seeds, same
%      weights and structural mode, different trajectories, one check.
%   6. Every figure draws on the run (PNG and _table.md written), and
%      errors on an empty run directory.
%
% Asserts; prints one PASS line. About 30 s. Assumes setup_paths has run
% and a pool allocation is set (pool_set(0) will do: the runs then go one
% after another in the client).
%
% NOT TESTED HERE: the stage through run_all_paper_analyses (FD-046 wires
% it), a run on a pool, anything at 500 neurons.
%
% See also: outlier_structures_stage_plan, run_outlier_structures_stage, ostruct_load

%% 1. The plan in every run mode
sizes = {'essential', 1; 'fast', 1; 'fastmedium', 5; 'medium', 5; 'medium2', 5; 'production', 5};
assert(isequal(sizes(:, 1)', run_mode_names()), 'the size table below must cover run_mode_names');
conds3 = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
keys15 = {'lrA1p5W0p0', 'lrA1p4W0p15', 'lrA1p4W0p35', 'lrA1p4W0p6', 'lrA1p4W1p0', 'lrA2p5W0p0', 'lrA2p5W1p0', ...
    'reference_w27', 'ind200_w25', 'reference', 'EEx1p5', 'EIx0p9', 'IIx1p25', 'EEx1p6-EIx1p33-IEx1p33', 'EEx1p6-EIx1p31-IEx1p31'};
for m = 1:size(sizes, 1)
    [md, n_seeds] = sizes{m, :};
    c = mu7chain_config(md);
    pl = outlier_structures_stage_plan(c);
    assert(strcmp(pl.folder, 'outlier_structures') && isequal({pl.networks.key}, keys15) && ...
        isequal(pl.conditions, conds3) && pl.n_seeds == n_seeds && pl.T_settle == 60 && pl.window_s == 20 && ...
        pl.dt_sample == 0.05 && pl.n_cal == 8 && pl.dt_mode == 0.1 && pl.dt_full == 20 && pl.n_override == 0, '%s: the plan', md);
    assert(numel(pl.runs) == 15 * 3 * n_seeds && all([pl.runs.full]) && isequal(pl.seeds(1, :), [7 8]), '%s: the runs', md);
    assert(numel(pl.checks) == 13 * n_seeds + 2, '%s: one check per free network and seed, one per fixed network', md);
    assert(isequal(unique({pl.networks.group}), {'blockmean', 'local', 'rotation'}) && ...
        nnz(strcmp({pl.networks.group}, 'rotation')) == 7 && nnz(strcmp({pl.networks.group}, 'local')) == 3, '%s: groups', md);
    assert(all(arrayfun(@(r) r.realization == 1, pl.runs(1:45))), '%s: seeds outermost', md);
end
r = pl.runs;                                               % production: 5 seeds
fixed = strcmp({r.key}, 'reference_w27');
assert(all(arrayfun(@(x) isequal(x.seed, [27 28]), r(fixed))), 'a fixed weight seed is every run''s rng_seeds');
assert(isequal(unique({r(fixed).network}, 'stable'), {'reference', 'reference-noise1', 'reference-noise2', ...
    'reference-noise3', 'reference-noise4'}), 'the noise realizations are named by seeds(s, 1)');
assert(all(isnan([r(fixed & [r.realization] == 1).noise_seed])) && all([r(fixed & [r.realization] == 3).noise_seed] == 2));
free = strcmp({r.key}, 'EEx1p5');
assert(isequal(unique(cat(1, r(free).seed), 'rows'), sortrows(pl.seeds)) && all(strcmp({r(free).network}, 'EEx1p5')), ...
    'a free network gets the seed rows as rng_seeds');
assert(numel(unique({r.file})) == numel(r));
assert(strcmp(r(1).file, 'ospec_lrA1p5W0p0_no_adaptation_s7-8.mat'));

c = mu7chain_config('fastmedium');
c.outlier_structures_n_seeds = 2;
c.outlier_structures_seeds = [11 12; 13 14; 15 16];
c.outlier_structures_networks = {'EEx1p5', [], []; 'ind200', [25 26], complex(1.16, 0.64)};
c.outlier_structures_conditions = 'sfa1_std1';
c.outlier_structures_window_s = 5;
pl = outlier_structures_stage_plan(c);
assert(numel(pl.runs) == 4 && isequal(pl.seeds, [11 12; 13 14]) && isequal(pl.conditions, {'sfa1_std1'}) && ...
    pl.dt_full == 5 && isequal({pl.runs.network}, {'EEx1p5', 'ind200', 'EEx1p5', 'ind200-noise13'}) && ...
    numel(pl.checks) == 3, 'the knobs are taken');

base = mu7chain_config('essential');
refused = { ...
    @(b) setfield(b, 'outlier_structures_n_seed', 2), ...                                   an unknown knob
    @(b) setfield(b, 'outlier_structures_networks', {'EEx1.5', [], []}), ...                not a network name
    @(b) setfield(b, 'outlier_structures_networks', {'reference-noise2', [27 28], 1.3}), ... a noise token
    @(b) setfield(b, 'outlier_structures_networks', {'reference', [27 28], []}), ...         a local outlier without expect
    @(b) setfield(b, 'outlier_structures_networks', {'EEx1p5', [1 2 3], []}), ...           a bad weight seed
    @(b) setfield(b, 'outlier_structures_networks', {'EEx1p5', [], []; 'EEx1p5', [], []}), ... twice
    @(b) setfield(b, 'outlier_structures_networks', {'EEx1p5', []}), ...                    not m x 3
    @(b) setfield(b, 'outlier_structures_conditions', {'sfa9'}), ...                        not a condition
    @(b) setfield(b, 'outlier_structures_dt_mode', 0), ...                                  the mode not followed
    @(b) setfield(b, 'outlier_structures_n_override', 5), ...                               too small
    @(b) setfield(setfield(b, 'outlier_structures_seeds', [1 2]), 'outlier_structures_n_seeds', 2), ... too few rows
    @(b) setfield(setfield(b, 'outlier_structures_seeds', [1 2; 1 3]), 'outlier_structures_n_seeds', 2), ... noise names collide
    @(b) setfield(b, 'preset_name', b.mc_preset), ...                                       noise-free
    @(b) setfield(b, 'preset_name', '')};                                                  % no preset
for k = 1:numel(refused)
    assert(throws(@() outlier_structures_stage_plan(refused{k}(base)), 'outlier_structures_stage_plan:BadConfig'), 'refusal %d', k);
end
bad = base; bad.run_mode = 'huge';
assert(throws(@() outlier_structures_stage_plan(bad), 'outlier_structures_stage_plan:BadConfig'));

%% 2. The noise token
p = srnn_param_preset(base.preset_name);
size_args = {'n', 40, 'indegree', 8, 'F_ref_n', 40, 'F_ref_indegree', 8};
B = struct('preset_name', base.preset_name, 'condition', 'no_adaptation', 'net_args', {{}}, 'size_args', {size_args}, ...
    'added_structure', [], 'rng_seeds', [27 28], 'T_range', [0 0.4], 'fs', 500, 'lya_args', {{'lya_method', 'none'}});
B1 = B; B1.net_args = ospec_network('reference', p.mu_tilde_relative).net_args;
B2 = B; B2.net_args = ospec_network('reference-noise3', p.mu_tilde_relative).net_args;
assert(isequal(B2.net_args, {'noise_seed', 3}) && strcmp(ospec_network('reference-noise3').label, 'reference') && ...
    strcmp(ospec_network('ind200-noise3').target, 'rightmost'));
m1 = ospec_build_model(B1); m2 = ospec_build_model(B2);
w1 = ospec_weight_stats(m1.get_params()); w2 = ospec_weight_stats(m2.get_params());
assert(isequal(w1.Wf, w2.Wf), 'the noise token leaves the weights alone');
S1 = m1.run_sampled(1:50:201, 0.4); S2 = m2.run_sampled(1:50:201, 0.4);
assert(isequal(S1(:, 1), S2(:, 1)) && ~isequal(S1(:, end), S2(:, end)), 'the same start, another noise path');
clear m1 m2

%% 3. The structural check at 80 neurons
cfg3 = tiny_cfg();
cfg3.outlier_structures_networks = {'lrA1p5W0p0', [], []; 'lrA1p4W0p6', [], []; 'EEx2p0', [], []; ...
    'reference', [27 28], complex(3, 0); 'reference', [], NaN};
pl3 = outlier_structures_stage_plan(cfg3);
ch = arrayfun(@(k) ostruct_structural_check(pl3, pl3.checks(k).key, pl3.checks(k).seed), 1:numel(pl3.checks));
assert(isequal([ch.pass], [true true true false true]), 'the structural check');
assert(ch(1).lam0_im_R == 0 && abs(ch(1).expect_re - 1.5) < 1e-12 && ch(2).lam0_im_R > 0.3 && ...
    abs(complex(ch(2).expect_re, ch(2).expect_im) - complex(1.4, 0.6)) < 1e-12, 'added real mode and rotating pair');
assert(strcmp(ch(3).target, 'blockmean') && ch(3).lam0_re_R > 2 && ch(3).rel_err < 0.15, 'the block-mean outlier');
assert(contains(ch(4).reason, 'from the expected') && isempty(ch(1).reason), 'a failed check says why');
assert(strcmp(ch(5).reason, 'no outlier intended') && isnan(ch(5).expect_re) && isnan(ch(5).expect_im));

%% 4. End to end, tiny
tmp = tempname;
cleanup = onCleanup(@() rmdir_if(tmp));
cfg = tiny_cfg();
cfg.outlier_structures_networks = {'lrA1p5W0p0', [], []; 'lrA1p4W0p6', [], []; 'EEx2p0', [], []; ...
    'reference', [27 28], 1.33; 'IIx1p25', [], []};
plan = outlier_structures_stage_plan(cfg);
t0 = tic;
evalc('detail = run_outlier_structures_stage(cfg, tmp, plan);');
t_run = toc(t0);
assert(contains(detail, '8 runs') && contains(detail, '2 run(s) of 1 draw(s) skipped'), detail);
folder = fullfile(tmp, 'outlier_structures');
assert(isfile(fullfile(folder, 'ostruct_index.mat')) && isfile(fullfile(folder, 'ospec_index.mat')));
ospec = ospec_read_index(folder);
assert(height(ospec) == 8 && all(ospec.ok) && all(ospec.has_dense) && ~any(contains(ospec.file, 'IIx1p25')), ...
    'the machinery ran the checked runs only');
S = ostruct_load(tmp, 'keep_runs', true);
assert(height(S.runs) == 10 && nnz(S.runs.ok) == 8 && isequal(S.runs.checked, ~strcmp(S.runs.key, 'IIx1p25')) && ...
    all(startsWith(S.runs.message(~S.runs.checked), 'skipped')), 'the stage''s index and the skipped runs');
assert(height(S.checks) == 5 && isequal(S.checks.pass', [true true true true false]) && ...
    ~isempty(S.checks.reason{5}), 'the checks are recorded');
assert(isequal(S.settings.conditions, {'no_adaptation', 'sfa3_stdc2'}) && numel(S.networks) == 5);
assert(height(S.T) == 8 && numel(S.Rs) == 8 && numel(S.MM) == 8 && isequal(S.idx, find(S.runs.ok)));
for q = 1:8
    R = S.Rs{q};
    k = S.idx(q);
    chk = S.checks(strcmp(S.checks.key, S.runs.key{k}), :);
    assert(abs(R.ref.lam0 / R.ref.R_struct - complex(chk.lam0_re_R, chk.lam0_im_R)) < 1e-10, ...
        '%s: the check found the eigenvalue the run follows', S.runs.key{k});
    assert(isequal(R.seed, [S.runs.seed1(k), S.runs.seed2(k)]) && numel(R.mode.t) == 21 && numel(R.cal.t) == 3);
end
need = {'key', 'group', 'base', 'realization', 'network', 'condition', 'lam0_re_R', 'lam0_im_R', 'lambda_1', ...
    'mode_rho', 'mode_rho_uniform', 'mode_reined', 'mode_rho_abs', 'mode_reined_abs', 'mode_ov', 'mode_frac_rightmost', ...
    'rate_struct', 'rate_frozen', 'rate_full', 'rot_struct', 'rot_frozen', 'rot_full', 'slow_shift', 'cont_ov', ...
    'gain_ratio', 'mean_over_rms', 'shift_obs_re', 'shift_obs_im', 'shift_pred_re', 'shift_pred_im', 'on_top', ...
    'on_rest', 'frac_on_E', 'frac_on_I', 'right_rho', 'right_reined', 'full_abscissa', 'n_out103_mean'};
missing = setdiff(need, S.T.Properties.VariableNames);
assert(isempty(missing), 'measures lack %s', strjoin(missing, ', '));
assert(all(isfinite(S.T.mode_rho_uniform)) && all(isfinite(S.T.rate_full)) && all(isfinite(S.T.rate_frozen)) && ...
    all(S.T.frac_on_E >= 0 & S.T.frac_on_E <= 1) && all(S.T.mode_ov > 0 & S.T.mode_ov <= 1 + 1e-12), 'measures in range');
% One measure by hand: the detachment against the theoretical radius.
q = find(strcmp(S.T.key, 'EEx2p0') & strcmp(S.T.condition, 'sfa3_stdc2'));
R = S.Rs{q};
[~, at] = ismember(R.mode.t, R.t);
assert(abs(S.T.mode_rho(q) - median(real(R.mode.lam) ./ R.R_eff(at))) < 1e-12 && ...
    abs(S.T.rot_struct(q) - imag(R.ref.lam0) * median(R.mode.d_rms) / R.tau_d) < 1e-12, 'mode_rho by hand');
assert(strcmp(S.T.group{q}, 'blockmean') && strcmp(S.T.group{strcmp(S.T.key, 'reference_w27') & ...
    strcmp(S.T.condition, 'no_adaptation')}, 'local'));
S2 = ostruct_load(folder, 'keys', {'EEx2p0'}, 'measures', false);
assert(isempty(S2.T) && isequal(S2.idx, find(S2.runs.ok & strcmp(S2.runs.key, 'EEx2p0'))));
assert(throws(@() ostruct_load(tempname), 'ostruct_load:NoStage'));
assert(throws(@() ostruct_load(tmp, 'keys', {'nope'}), 'ostruct_load:Key'));

%% 5. Noise realizations on fixed weights
tmp5 = tempname;
cleanup5 = onCleanup(@() rmdir_if(tmp5));
cfg5 = tiny_cfg();
cfg5.outlier_structures_networks = {'reference', [27 28], 1.33};
cfg5.outlier_structures_conditions = {'no_adaptation'};
cfg5.outlier_structures_n_seeds = 2;
cfg5.outlier_structures_lya_K = 0;
plan5 = outlier_structures_stage_plan(cfg5);
evalc('run_outlier_structures_stage(cfg5, tmp5, plan5);');
S5 = ostruct_load(tmp5, 'keep_runs', true);
assert(height(S5.checks) == 1 && height(S5.T) == 2 && isequal(S5.T.realization', [1 2]) && ...
    isequal(S5.T.network', {'reference', 'reference-noise1'}) && isequal(unique(S5.T.key), {'reference_w27'}));
Ra = S5.Rs{1}; Rb = S5.Rs{2};
assert(isequal(Ra.ref.lam0, Rb.ref.lam0) && isequal(Ra.seed, Rb.seed) && ~isequal(Ra.cal.S, Rb.cal.S), ...
    'the same weights and structural mode, another noise realization');

%% 6. The figures
out_dir = fullfile(tmp, 'figs');
mkdir(out_dir);
figs = {@fig_ostruct_summary, {}, 'Fig_Ostruct_Summary'; @fig_ostruct_rotation, {}, 'Fig_Ostruct_Rotation'; ...
    @fig_ostruct_h1, {}, 'Fig_Ostruct_H1'; @fig_ostruct_plane, {'group', 'rotation'}, 'Fig_Ostruct_Plane_rotation'; ...
    @fig_ostruct_plane, {'group', 'local'}, 'Fig_Ostruct_Plane_local'; ...
    @fig_ostruct_plane, {'group', 'blockmean'}, 'Fig_Ostruct_Plane_blockmean'};
empty_run = tempname; mkdir(empty_run);
cleanup6 = onCleanup(@() rmdir_if(empty_run));
for i = 1:size(figs, 1)
    o = figs{i, 1}('run_dir', tmp, 'out_dir', out_dir, 'visible', false, 'verbose', 'near-none', figs{i, 2}{:});
    close(o.figs);
    assert(isfile(fullfile(out_dir, [figs{i, 3} '.png'])) && isfile(fullfile(out_dir, [figs{i, 3} '_table.md'])), ...
        '%s: written', figs{i, 3});
    txt = fileread(fullfile(out_dir, [figs{i, 3} '_table.md']));
    assert(contains(txt, '## The structural check') && contains(txt, 'median'), '%s: the sidecar', figs{i, 3});
    assert(throws_prefix(@() figs{i, 1}('run_dir', empty_run, 'save', false, 'visible', false, figs{i, 2}{:}), ...
        'resolve_data_file:'), '%s on an empty run directory', figs{i, 3});
end
txt = fileread(fullfile(out_dir, 'Fig_Ostruct_Plane_blockmean_table.md'));
assert(contains(txt, 'IIx1p25') && contains(txt, 'skipped'), 'the skipped network is named in its group''s sidecar');

fprintf('PASS: test_outlier_structures_stage -- the plan in every mode and its refusals, the noise token, the structural check (real, pair, block mean, a failure), the stage end to end with a skipped network (%.0f s), noise realizations, six figures.\n', t_run);

%% ------------------------------------------------------------------------
function cfg = tiny_cfg()
% The bundle's essential config with the stage shrunk to seconds.
cfg = mu7chain_config('essential');
cfg.verbose = 'near-none';
cfg.outlier_structures_n_override = 80;
cfg.outlier_structures_conditions = {'no_adaptation', 'sfa3_stdc2'};
cfg.outlier_structures_T_settle = 2;
cfg.outlier_structures_window_s = 2;
cfg.outlier_structures_lya_window_s = 1.5;
cfg.outlier_structures_n_cal = 3;
cfg.outlier_structures_n_null = 4;
cfg.outlier_structures_n_shuffle = 3;
cfg.outlier_structures_chunk_s = 0.7;
end

function tf = throws(fn, id)
try
    fn();
    tf = false;
catch err
    tf = strcmp(err.identifier, id);
    if ~tf; fprintf(2, 'expected %s, got %s: %s\n', id, err.identifier, err.message); end
end
end

function tf = throws_prefix(fn, prefix)
try
    fn();
    tf = false;
catch err
    tf = startsWith(err.identifier, prefix);
    if ~tf; fprintf(2, 'expected %s..., got %s: %s\n', prefix, err.identifier, err.message); end
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end

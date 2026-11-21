%% test_outlier_spectra_stage.m
% The opt-in stage 'outlier_spectra' (FD-040): outlier_spectra_stage_plan,
% run_outlier_spectra_stage, ospec_run and the readers ospec_read_index,
% ospec_load_run, ospec_rebuild_model, on a tiny network.
%
% Checks:
%   1. The plan in EVERY run mode, from mu7chain_config: five networks,
%      the three conditions, 2 seeds and a 20 s window at essential and
%      fast, 10 seeds and 200 s at 50 ms for reference and E-to-E x 1.5
%      above; seeds outermost, the exploration's seed [7 8] first; the
%      knobs are taken; bad settings are refused (an unknown knob, an
%      unknown network, a full network that is not a network, a condition
%      the preset does not state, more seeds than seed rows, a noise-free
%      preset, a size below 20, a sampling step below the integration
%      step, no preset).
%   2. END TO END, tiny: 60 neurons, 2 s of settling, a 2 s window, 3
%      calibration states, 1 seed, two networks (reference; E-to-E x 1.5
%      with the full-Jacobian passes), two conditions (no adaptation,
%      chained multiple timescales). The index lists exactly the files
%      written; every field of the contract (the header of ospec_load_run)
%      is there with its shape; the saved operating points are consistent
%      (the rate is the centered sigmoid at u; gain_route times W is the
%      effective connectivity; the frozen-gain series is the calibration
%      spectrum at the calibration states).
%   3. ospec_rebuild_model reproduces the saved spectra from the saved
%      states: the frozen-gain spectrum of every calibration state of
%      every run, and the dense full-Jacobian spectrum of one.
%   4. OspecNet.run_sampled is run(): the same states bit for bit and the
%      same Lyapunov summary, on the model of one run, with pieces that do
%      not divide the run.
%   5. The kinetic edge by neuron blocks is the dense eig of the slow
%      block; ospec_excursions on a made-up series (start, end, duration,
%      cut off at either end, NaN not beyond); ospec_measures runs on a
%      run with and one without the full passes.
%   6. A RERUN INTO THE NON-EMPTY FOLDER starts over, as the other stages
%      do: a smaller plan into the same folder replaces the index (it
%      lists the new plan's runs only), overwrites the run's file, and the
%      run is THE SAME RUN (everything but the wall-clock times equal, the
%      null draws included).
%   7. A run that fails does not stop the others: the index is written
%      with ok false and the message, then the stage errors.
%
% Asserts; prints one PASS line. About 20 s. Assumes setup_paths has run
% and a pool allocation is set (pool_set(0) will do; the stage then runs
% its runs one after another in the client).
%
% NOT TESTED HERE: the stage through run_all_paper_analyses (it is wired
% into paper_opt_in_stages and the bundles by FD-046), a run on a pool, and
% anything at 500 neurons.
%
% See also: outlier_spectra_stage_plan, run_outlier_spectra_stage, ospec_load_run

%% 1. The plan in every run mode
sizes = {'essential', 2, 20, 0.01; 'fast', 2, 20, 0.01; 'fastmedium', 10, 200, 0.05; ...
    'medium', 10, 200, 0.05; 'medium2', 10, 200, 0.05; 'production', 10, 200, 0.05};
assert(isequal(sizes(:, 1)', run_mode_names()), 'the size table below must cover run_mode_names');
nets = {'reference', 'EEx1p1', 'EEx1p25', 'EEx1p5', 'EEx1p75'};
conds3 = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
for m = 1:size(sizes, 1)
    [md, n_seeds, full_window, full_dt] = sizes{m, :};
    c = mu7chain_config(md);
    pl = outlier_spectra_stage_plan(c);
    assert(strcmp(pl.folder, 'outlier_spectra') && strcmp(pl.preset_name, c.preset_name) && ...
        isequal(pl.networks, nets) && isequal(pl.full_networks, {'reference', 'EEx1p5'}) && ...
        isequal(pl.conditions, conds3) && pl.n_seeds == n_seeds && pl.T_settle == 60 && ...
        pl.window_s == 20 && pl.dt_sample == 0.01 && pl.full_window_s == full_window && ...
        pl.full_dt_sample == full_dt && pl.n_cal == 12 && pl.n_override == 0 && pl.fs == 500, '%s: the plan', md);
    assert(numel(pl.runs) == 5 * 3 * n_seeds && isequal(pl.seeds(1, :), [7 8]) && isequal(pl.seeds(2, :), [1 2]), ...
        '%s: the runs', md);
    r = pl.runs;
    is_full = ismember({r.network}, {'reference', 'EEx1p5'});
    assert(isequal([r.full], is_full) && all([r(is_full).window_s] == full_window) && ...
        all([r(~is_full).window_s] == 20) && all([r(is_full).dt_sample] == full_dt) && ...
        all([r(~is_full).dt_sample] == 0.01), '%s: the long window goes to the full networks only', md);
    first = r(1:15);                                   % seeds outermost: the first 15 runs are seed [7 8]
    assert(all(arrayfun(@(x) isequal(x.seed, [7 8]), first)) && numel(unique({r.file})) == numel(r));
end
assert(isequal(pl.seeds, [7 8; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 8 9; 9 10; 10 11]), 'the ten default seeds');
assert(strcmp(pl.runs(1).file, 'ospec_reference_no_adaptation_s7-8.mat'));

c = mu7chain_config('fastmedium');
c.outlier_spectra_n_seeds = 3;
c.outlier_spectra_seeds = [11 12; 13 14; 15 16; 17 18];
c.outlier_spectra_networks = {'reference', 'EEx1p5-ind50'};
c.outlier_spectra_full_networks = {};
c.outlier_spectra_conditions = 'sfa1_std1';
c.outlier_spectra_full_window_s = 40;
pl = outlier_spectra_stage_plan(c);
assert(numel(pl.runs) == 6 && isequal(pl.seeds, [11 12; 13 14; 15 16]) && isempty(pl.full_networks) && ...
    ~any([pl.runs.full]) && all([pl.runs.window_s] == 20) && isequal(pl.conditions, {'sfa1_std1'}), 'the knobs are taken');

refused = { ...
    @(b) setfield(b, 'outlier_spectra_n_seed', 2), ...                           an unknown knob (a typo)
    @(b) setfield(b, 'outlier_spectra_networks', {'reference', 'EEx1.5'}), ...   not a network name
    @(b) setfield(b, 'outlier_spectra_full_networks', {'EEx2p0'}), ...           full, but not a network
    @(b) setfield(b, 'outlier_spectra_conditions', {'sfa9'}), ...                not a condition of the preset
    @(b) setfield(setfield(b, 'outlier_spectra_seeds', [1 2]), 'outlier_spectra_n_seeds', 2), ...   more seeds than seed rows
    @(b) setfield(b, 'preset_name', b.mc_preset), ...                            noise-free
    @(b) setfield(b, 'outlier_spectra_n_override', 5), ...                       too small
    @(b) setfield(b, 'outlier_spectra_dt_sample', 0.001), ...                    below the integration step
    @(b) setfield(b, 'outlier_spectra_x_thr', 1.5), ...                          not a fraction
    @(b) setfield(b, 'preset_name', '')};                                       % no preset
for k = 1:numel(refused)
    bad = refused{k}(mu7chain_config('essential'));
    assert(throws(@() outlier_spectra_stage_plan(bad), 'outlier_spectra_stage_plan:BadConfig'), 'setting %d', k);
end
bad = mu7chain_config('essential'); bad.run_mode = 'huge';
assert(throws(@() outlier_spectra_stage_plan(bad), 'outlier_spectra_stage_plan:BadConfig'));

% The network names build what the exploration's did (its function is on the path).
a = ospec_network('EEx1p5-ind50'); b = outlier_network('EEx1p5-ind50');
assert(isequal(a, b), 'ospec_network with the default block means is outlier_network');

%% 2. End to end, tiny
tmp = tempname;
cleanup = onCleanup(@() rmdir_if(tmp));
cfg = tiny_cfg();
plan = outlier_spectra_stage_plan(cfg);
assert(numel(plan.runs) == 4);
t0 = tic;
detail = run_outlier_spectra_stage(cfg, tmp, plan);
t_run = toc(t0);
assert(contains(detail, '4 runs'));

[runs, settings, folder] = ospec_read_index(tmp);
assert(strcmp(folder, fullfile(tmp, 'outlier_spectra')) && height(runs) == 4 && all(runs.ok));
on_disk = dir(fullfile(folder, 'ospec_*_s*.mat'));
names_on_disk = {on_disk.name};
assert(isequal(sort(names_on_disk(:)), sort(runs.file)), 'the index lists exactly the run files written');
assert(isequal(runs.Properties.VariableNames(1:8), {'network', 'condition', 'seed1', 'seed2', 'file', ...
    'has_full_series', 'has_dense', 'window_s'}), 'the columns of the contract');
assert(isequal(runs.network', {'reference', 'reference', 'EEx1p5', 'EEx1p5'}) && ...
    isequal(runs.condition', {'no_adaptation', 'sfa3_stdc2', 'no_adaptation', 'sfa3_stdc2'}) && ...
    all(runs.seed1 == 7) && all(runs.seed2 == 8) && isequal(runs.has_full_series', [false false true true]) && ...
    all(runs.window_s == 2));
assert(settings.n_override == 60 && settings.n_workers == 0 && ~isfield(settings, 'runs'));
[runs2, ~, folder2] = ospec_read_index(folder);          % the stage's folder itself is accepted
assert(isequal(runs2, runs) && strcmp(folder2, folder));
assert(throws(@() ospec_read_index(tempname), 'ospec_read_index:NoIndex'));

n = 60; n_s = 201; n_cal = 3; C = 2;
Rs = cell(1, 4);
for k = 1:4
    R = ospec_load_run(tmp, k);
    Rs{k} = R;
    tag = sprintf('%s | %s', R.network, R.condition);
    N = R.N;
    assert(strcmp(R.network, runs.network{k}) && strcmp(R.condition, runs.condition{k}) && ...
        isequal(R.seed, [7 8]) && isequal(R.rng_seeds, [7 8]) && strcmp(R.preset_name, cfg.preset_name) && ...
        R.n == n && R.window_s == 2 && R.has_full_series == runs.has_full_series(k) && ...
        R.has_dense == runs.has_dense(k) && R.stage_version == 1 && iscell(R.net_args), '%s: what the run is', tag);
    assert((strcmp(R.condition, 'no_adaptation') && N == n) || (strcmp(R.condition, 'sfa3_stdc2') && N > 4 * n), ...
        '%s: N = %d', tag, N);
    % The series
    assert(isequal(size(R.t), [1 n_s]) && abs(R.t(1) - 2) < 1e-9 && abs(R.t(end) - 4) < 1e-9 && abs(R.dt - 0.01) < 1e-12);
    shapes = {'ev_top', [3 n_s]; 'rate_top', [1 n_s]; 'R_eff', [1 n_s]; 'ratio_theory', [1 n_s]; ...
        're_keep', [30 n_s]; 'abs_keep', [30 n_s]; 'n_unstable', [1 n_s]; 'res', [3 n_s]; 'coh', [3 n_s]; ...
        'K_eig', [C n_s]; 'K_eff', [C C n_s]; 'dbar', [C C n_s]; 'drms', [C C n_s]; 'rate', [C n_s]; ...
        'frac_on', [C n_s]; 'vec_t', [1 21]; 'V_top', [n 3 21]};
    check_shapes(R, shapes, tag);
    assert(isequal(R.rate_top, (real(R.ev_top(1, :)) - 1) ./ R.tau_d) && isa(R.V_top, 'single') && ...
        max(R.res(:)) < 1e-8 && all(ismember(R.vec_t, R.t)), '%s: the series', tag);
    assert(all(isfinite([R.kappa, R.kappa_abs, R.agree_re, R.agree_abs])) && R.kappa > 0.5 && R.kappa < 3, ...
        '%s: the null-calibrated factors', tag);
    assert(isfield(R.lya, 'LLE') && isfinite(R.lya.LLE), '%s: lambda_1', tag);
    % The calibration states
    shapes = {'t', [1 n_cal]; 'idx', [1 n_cal]; 'S', [n_cal N]; 'gain', [n n_cal]; 'gain_route', [C n n_cal]; ...
        'slope', [n n_cal]; 'u', [n n_cal]; 'rate', [n n_cal]; 'cell_type', [n 1]; 'S_c', [n 1]; ...
        'ev_M', [n n_cal]; 'ev_uniform', [n n_cal]; 'R_eff', [1 n_cal]};
    check_shapes(R.cal, shapes, [tag ' cal']);
    assert(isequal(R.t(R.cal.idx), R.cal.t) && isequal(R.cal.idx, [1 101 201]) && ...
        isequal(R.cal.R_eff, R.R_eff(R.cal.idx)), '%s: the calibration states are states of the series', tag);
    assert(isequal(sort(unique(R.cal.cell_type))', 1:2) && sum(R.cal.cell_type == 1) == R.n_per_type(1));
    assert(max(abs(R.cal.rate - SRNNCellTypePairs.piecewiseSigmoid(R.cal.u, R.S_a, 0)), [], 'all') < 1e-12, ...
        '%s: the rate is the centered sigmoid at u', tag);
    assert(all(R.cal.gain(:) >= 0 | isnan(R.cal.gain(:))) && all(R.cal.slope(:) >= 0));
    % The rightmost real part of the series at a calibration state is that of the whole spectrum there.
    assert(max(abs(max(real(R.cal.ev_M), [], 1) - real(R.ev_top(1, R.cal.idx)))) < 1e-9, '%s: series against ev_M', tag);
    assert(isequal(size(R.null), [1 n_cal]) && isequal(size(R.null(1).null_re), [4 1]) && ...
        isequal(size(R.null(1).shuf_re), [3 1]) && isequal(R.null(2).ev_M, R.cal.ev_M(:, 2)));
    % The followed mode
    assert(isequal(size(R.mode.t), [1 21]) && isequal(size(R.mode.lam), [1 21]) && isfield(R.ref, 'lam0'));
    % The full Jacobian
    if R.has_full_series
        F = R.full;
        kk = F.k;
        assert(isequal(size(F.t), [1 21]) && abs(F.dt - 0.1) < 1e-12 && all(ismember(F.t, R.t)) && F.x_thr == 0.1);
        check_shapes(F, {'lam', [kk 21]; 'xfrac', [kk 21]; 'pr', [kk 21]; 'coh', [kk 21]; 'n_converged', [1 21]; ...
            'abscissa', [1 21]; 'net', [1 21]; 'net_found', [1 21]; 'bound', [1 21]; 'kin_edge', [1 21]}, [tag ' full']);
        assert(isequal(isnan(F.net), ~F.net_found) && all(imag(F.net(F.net_found)) >= 0), '%s: net', tag);
        % An eigenvalue beyond zero is a network mode: the two series have the same excursions.
        E1 = ospec_excursions(F.t, F.abscissa); E2 = ospec_excursions(F.t, real(F.net));
        assert(isequal(E1.t_start, E2.t_start) && isequal(E1.t_end, E2.t_end), '%s: excursions beyond zero', tag);
        assert(all(F.bound(~F.net_found) < 0), '%s: an unresolved state is bounded below zero', tag);
        if strcmp(R.condition, 'no_adaptation')
            % No slow variables: every eigenvalue is a network mode, and the
            % full Jacobian is the dendritic block.
            assert(all(F.net_found) && all(F.kin_edge == -Inf) && ...
                max(abs(F.abscissa - interp1(R.t, R.rate_top, F.t, 'nearest'))) < 1e-8, '%s: full = frozen', tag);
        else
            assert(all(F.kin_edge < 0) && all(isfinite(F.kin_edge)), '%s: the kinetic edge', tag);
        end
        assert(isequal(size(R.dense), [1 n_cal]) && isequal(size(R.dense(1).lam), [N 1]) && ...
            isa(R.dense(1).lam, 'single') && isequal(size(R.dense(1).xfrac), [N 1]) && ...
            all(isfield(R.dense, {'net', 'net_rank', 'kin_edge', 'cont', 'cont_ov', 'mode_cont', 'n_pos', 'abscissa'})), ...
            '%s: the dense spectra', tag);
    else
        assert(isempty(R.full) && isempty(R.dense), '%s: no full passes', tag);
    end
end

%% 3. The rebuilt model reproduces the saved spectra from the saved states
worst = 0;
for k = 1:4
    R = Rs{k};
    [model, params, ws] = ospec_rebuild_model(R);
    assert(isa(model, 'SRNNCellTypePairs') && ~model.has_run && params.N_sys_eqs == R.N);
    for q = 1:n_cal
        S = R.cal.S(q, :)';
        [M, d] = ospec_effective_matrix(S, params, ws);
        worst = max(worst, max(min(abs(R.cal.ev_M(:, q) - eig(M).'), [], 2)));
        % gain_route times W is M (the saved gains need no recovery).
        G = R.cal.gain_route(:, :, q); G(isnan(G)) = 0;
        assert(max(abs(ws.Wf .* G(R.cal.cell_type, :) - M), [], 'all') < 1e-12 && ...
            max(abs(G - d), [], 'all') == 0, 'gain_route');
    end
end
assert(worst < 1e-10, 'the frozen-gain spectra from the saved states: largest difference %.3g', worst);
R = Rs{4};
[~, params, ws] = ospec_rebuild_model(R);
lam = eig(full(SRNNCellTypePairs.compute_Jacobian_fast(R.cal.S(2, :)', params)));
d_dense = max(min(abs(double(R.dense(2).lam) - lam.'), [], 2));
assert(d_dense < 1e-3 * max(abs(lam)), 'the dense spectrum (saved in single): largest difference %.3g', d_dense);
R_broken = R; R_broken.N = R.N + 1;
assert(throws(@() ospec_rebuild_model(R_broken), 'ospec_rebuild_model:Size'));

%% 4. run_sampled is run()
m1 = ospec_build_model(R.build);
m1.store_full_state = true;
m1.run();
L1 = m1.lya_summary();
keep = [1, 7:13:numel(m1.t_out), numel(m1.t_out)];
m2 = ospec_build_model(R.build);
[S_keep, L2] = m2.run_sampled(keep, 0.37);              % pieces that do not divide the run
assert(isequal(S_keep, m1.S_out(keep, :)'), 'run_sampled must return the rows of run() bit for bit');
assert(isequaln(L1, L2), 'and the same Lyapunov summary');
assert(isequal(R.cal.S, m1.S_out(round(R.cal.t * R.build.fs) + 1, :)), 'the stage''s saved states are the trajectory''s');
clear m1 m2 S_keep

%% 5. The kinetic edge, the excursions, the measures
J = SRNNCellTypePairs.compute_Jacobian_fast(R.cal.S(1, :)', params);
[edge, band, by_blocks] = ospec_kinetic_edge(J, ws);
es = sort(real(eig(full(J(~ws.is_x, ~ws.is_x)))));
assert(by_blocks && abs(edge - es(end)) < 1e-10 && edge < 0 && all(isfinite(band)), 'the kinetic edge by blocks');
assert(abs(edge - R.dense(1).kin_edge) < 1e-12);
[~, params0, ws0] = ospec_rebuild_model(Rs{3});
assert(ospec_kinetic_edge(SRNNCellTypePairs.compute_Jacobian_fast(Rs{3}.cal.S(1, :)', params0), ws0) == -Inf);

t = 10 + 0.1 * (0:9);
E = ospec_excursions(t, [1 1 -1 2 NaN 3 3 3 -1 4]);
assert(isequal(E.t_start, t([1 4 6 10])) && isequal(E.t_end, t([2 4 8 10])) && ...
    max(abs(E.duration - [0.2 0.1 0.3 0.1])) < 1e-12 && isequal(E.cut_start, [true false false false]) && ...
    isequal(E.cut_end, [false false false true]) && isequal(E.censored, [true false false true]) && ...
    abs(E.fraction - 0.7) < 1e-12, 'ospec_excursions');
E = ospec_excursions(t, -ones(1, 10));
assert(isempty(E.duration) && E.fraction == 0);
E = ospec_excursions(t, 0.5 + zeros(1, 10), 1);
assert(isempty(E.duration));

for k = [1 4]
    mm = ospec_measures(Rs{k});
    assert(isequal(size(mm.ratio_theory), [1 n_s]) && isfinite(mm.row.ratio_theory_med) && ...
        isequal(mm.exc.zero.duration, mm.dwell.zero) && mm.row.window_s == 2, 'ospec_measures, run %d', k);
    assert(isnan(mm.row.X_full) == ~Rs{k}.has_full_series);
end

%% 6. A rerun into the non-empty folder starts over
file4 = fullfile(folder, runs.file{4});
before = dir(file4);
pause(1.1);                                             % file times have a one-second resolution on some systems
cfg2 = tiny_cfg();
cfg2.outlier_spectra_networks = {'EEx1p5'};
cfg2.outlier_spectra_conditions = {'sfa3_stdc2'};
plan2 = outlier_spectra_stage_plan(cfg2);
run_outlier_spectra_stage(cfg2, tmp, plan2);
runs_again = ospec_read_index(tmp);
assert(height(runs_again) == 1 && strcmp(runs_again.file{1}, runs.file{4}) && runs_again.ok, ...
    'the index of a rerun lists the new plan''s runs only');
after = dir(file4);
assert(after.datenum > before.datenum, 'the run''s file is written again');
assert(isfile(fullfile(folder, runs.file{1})), 'files the new plan does not name are left alone');
R_again = ospec_load_run(tmp, struct('network', 'EEx1p5', 'condition', 'sfa3_stdc2', 'seed', [7 8]));
assert(isequaln(strip_times(R_again), strip_times(Rs{4})), 'a rerun is the same run, null draws included');
assert(throws(@() ospec_load_run(tmp, 2), 'ospec_load_run:Row'));

%% 7. A run that fails does not stop the others
plan3 = outlier_spectra_stage_plan(tiny_cfg());
plan3.runs = plan3.runs([1 3]);
plan3.runs(1).network = 'no-such-network';
plan3.runs(1).file = 'ospec_broken.mat';
tmp3 = tempname;
cleanup3 = onCleanup(@() rmdir_if(tmp3));
evalc('failed = throws(@() run_outlier_spectra_stage(tiny_cfg(), tmp3, plan3), ''run_outlier_spectra_stage:RunsFailed'');');
assert(failed, 'a stage with a failed run errors');
runs3 = ospec_read_index(tmp3);
assert(isequal(runs3.ok', [false true]) && contains(runs3.message{1}, 'ospec_network:Name') && isempty(runs3.message{2}) && ...
    isfile(fullfile(tmp3, 'outlier_spectra', runs3.file{2})), 'the index records the failure and the finished run');
assert(throws(@() ospec_load_run(tmp3, 1), 'ospec_load_run:NoFile'));

fprintf('PASS: test_outlier_spectra_stage -- the plan in every mode, refusals, the contract, the rebuilt model (largest spectral difference %.2g), run_sampled = run(), a rerun, a failed run (stage %.0f s).\n', ...
    worst, t_run);

%% ------------------------------------------------------------------------
function cfg = tiny_cfg()
% The bundle's essential config with the stage shrunk to seconds.
cfg = mu7chain_config('essential');
cfg.verbose = 'near-none';
cfg.outlier_spectra_n_override = 60;
cfg.outlier_spectra_networks = {'reference', 'EEx1p5'};
cfg.outlier_spectra_full_networks = {'EEx1p5'};
cfg.outlier_spectra_conditions = {'no_adaptation', 'sfa3_stdc2'};
cfg.outlier_spectra_n_seeds = 1;
cfg.outlier_spectra_T_settle = 2;
cfg.outlier_spectra_window_s = 2;
cfg.outlier_spectra_full_window_s = 2;
cfg.outlier_spectra_lya_window_s = 1.5;
cfg.outlier_spectra_n_cal = 3;
cfg.outlier_spectra_n_null = 4;
cfg.outlier_spectra_n_shuffle = 3;
cfg.outlier_spectra_chunk_s = 0.7;
end

function check_shapes(s, shapes, tag)
for k = 1:size(shapes, 1)
    assert(isfield(s, shapes{k, 1}), '%s: no field %s', tag, shapes{k, 1});
    assert(isequal(size(s.(shapes{k, 1})), shapes{k, 2}), '%s: %s is %s, expected %s', tag, shapes{k, 1}, ...
        mat2str(size(s.(shapes{k, 1}))), mat2str(shapes{k, 2}));
end
end

function R = strip_times(R)
% A run without its wall-clock fields.
R = rmfield(R, 'seconds');
R.opts = rmfield(R.opts, 'out_dir');
if ~isempty(R.full); R.full = rmfield(R.full, 'seconds'); end
if ~isempty(R.dense); R.dense = rmfield(R.dense, 'seconds'); end
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

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end

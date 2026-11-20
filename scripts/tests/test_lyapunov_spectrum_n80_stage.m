%% test_lyapunov_spectrum_n80_stage.m
% The ranked Lyapunov spectrum of the 80-neuron network (FD-032): the opt-in
% stage 'lyapunov_spectrum_n80' (lyapunov_spectrum_n80_stage_plan,
% run_lyapunov_spectrum_n80_stage), the new window arguments of
% run_lyapunov_spectrum, and the figure fig_lyapunov_ranked_spectrum.
%
% Checks:
%   1. The plan in EVERY run mode, from mu7chain_config: the 80-neuron
%      preset, K 240 with no retry doubling (K_max = K), T 140 s, window
%      [80 140] s, warm-up 20 s, 500 Hz, noise on only, the bundle's
%      conditions, one seed at essential and fast and three above, the
%      folder. Overrides are taken; an unknown field, a bad window, a
%      missing or noise-free preset are refused, and through
%      run_all_paper_analyses before the run directory exists.
%   2. run_lyapunov_spectrum's DEFAULT path is what it was: a call naming
%      none of the new arguments records T, window, warm-up, fs and K_max
%      as the run mode's table gives them (T, [T/2 T], T/4, 400 Hz, 2K),
%      and returns exactly what a call naming those values explicitly
%      returns (a 30-neuron network, 'fast', one seed).
%   3. END TO END, tiny: a bundle naming only the stage, the 80-neuron
%      network at T = 4 s, window [2 4] s, warm-up 1 s, K = 100, two seeds,
%      two conditions; run_all_paper_analyses writes
%      <run_dir>/lyapunov_spectrum_n80/, the manifest records the stage, K is
%      capped at the state dimension (no adaptation: N = 80, all 80 taken;
%      one timescale: N > 100, 100 taken), the window and warm-up reached the
%      model (the Lyapunov segments start at t1 - warm-up and end at t2), the
%      integrator is sra1 at 500 Hz with the preset's noise; then
%      make_all_paper_figures draws the registered entry, which writes its
%      PNG and table, and the figure errors on a run without the stage.
%
% Asserts; prints one PASS line. About half a minute. Assumes setup_paths
% has run and a pool allocation is set (pool_set(0) will do): both entry
% points refuse without one.
%
% See also: lyapunov_spectrum_n80_stage_plan, run_lyapunov_spectrum,
%           fig_lyapunov_ranked_spectrum, test_lyapunov_spectrum_stage

root = fileparts(which('setup_paths'));
n80 = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain_n80';

%% 1. The plan in every run mode
seeds = {'essential', 1; 'fast', 1; 'fastmedium', 3; 'medium', 3; 'medium2', 3; 'production', 3};
assert(isequal(seeds(:, 1)', run_mode_names()), 'the seed table below must cover run_mode_names');
for m = 1:size(seeds, 1)
    [md, n_seeds] = seeds{m, :};
    c = mu7chain_config(md);
    assert(ismember('lyapunov_spectrum_n80', paper_stage_plan(c).stages), '%s does not run the stage', md);
    pl = lyapunov_spectrum_n80_stage_plan(c);
    want_conds = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
    if ~strcmp(md, 'essential'); want_conds = want_conds(1:3); end   % FD-030
    assert(strcmp(pl.folder, 'lyapunov_spectrum_n80') && strcmp(pl.preset_name, n80));
    assert(pl.n_seeds == n_seeds && pl.K == 240 && pl.T == 140 && isequal(pl.lya_T_interval, [80 140]) && ...
        pl.lya_warmup == 20 && pl.fs == 500 && isequal(pl.conditions, want_conds), '%s: the plan', md);
    assert(strcmp(arg(pl.args, 'preset_name'), n80) && strcmp(arg(pl.args, 'run_mode'), md) && ...
        arg(pl.args, 'n_seeds') == n_seeds && arg(pl.args, 'K') == 240 && arg(pl.args, 'K_max') == 240 && ...
        arg(pl.args, 'T') == 140 && isequal(arg(pl.args, 'lya_T_interval'), [80 140]) && ...
        arg(pl.args, 'lya_warmup') == 20 && arg(pl.args, 'fs') == 500 && ...
        arg(pl.args, 'noise_off_too') == false && isequal(arg(pl.args, 'conditions'), want_conds), ...
        '%s: the arguments', md);
end
% The A + C bundle runs it too, after block A's stages and before block C's.
ac = paper_stage_plan(mu7chain_AC_config('fast'));
k = find(strcmp(ac.stages, 'lyapunov_spectrum_n80'));
assert(isscalar(k) && strcmp(ac.stages{k + 1}, 'mc_gain_sweep'));

c = mu7chain_config('essential');
c.lyapunov_spectrum_n80 = struct('n_seeds', 2, 'K', 30, 'T', 6, 'lya_T_interval', [3 6], 'lya_warmup', 2);
pl = lyapunov_spectrum_n80_stage_plan(c);
assert(pl.n_seeds == 2 && pl.K == 30 && arg(pl.args, 'K_max') == 30 && pl.T == 6 && ...
    isequal(pl.lya_T_interval, [3 6]) && pl.lya_warmup == 2);

probe = fullfile(root, 'data', '__lyapunov_n80_preflight_probe__');
assert(~isfolder(probe));
refused = { ...
    @(b) setfield(b, 'lyapunov_spectrum_n80', struct('seeds', 2)), ...                     an unknown field
    @(b) setfield(b, 'lyapunov_spectrum_n80', struct('lya_T_interval', [80 150])), ...    window past T
    @(b) setfield(b, 'lyapunov_spectrum_n80', struct('lya_warmup', 90)), ...               warm-up before 0
    @(b) setfield(b, 'lyapunov_spectrum_n80_preset', ''), ...                              no preset
    @(b) setfield(b, 'lyapunov_spectrum_n80_preset', b.mc_preset)};                       % noise-free
for k = 1:numel(refused)
    bad = mu7chain_config('essential');
    bad.run_dir = probe;
    bad = refused{k}(bad);
    assert(throws(@() lyapunov_spectrum_n80_stage_plan(bad), 'lyapunov_spectrum_n80_stage_plan:BadConfig'), ...
        'setting %d', k);
    assert(throws(@() run_all_paper_analyses(bad), 'lyapunov_spectrum_n80_stage_plan:BadConfig'), ...
        'setting %d was not refused by the entry point', k);
end
assert(~isfolder(probe), 'a refused run must not create its run directory');

%% 2. The default path of run_lyapunov_spectrum
tmp_d = tempname; tmp_e = tempname;
cleanup_de = onCleanup(@() cellfun(@rmdir_if, {tmp_d, tmp_e}));
base = {'preset_name', mu7chain_config('fast').preset_name, 'run_mode', 'fast', 'verbose', 'near-none', ...
    'n_override', 30, 'n_seeds', 1, 'K', 5, 'noise_off_too', false, 'conditions', {'sfa1_std1'}};
f_d = run_lyapunov_spectrum(base{:}, 'out_dir', tmp_d);
f_e = run_lyapunov_spectrum(base{:}, 'out_dir', tmp_e, 'T', 20, 'lya_T_interval', [10 20], ...
    'lya_warmup', 5, 'fs', 400, 'K_max', 10);
Dd = load(f_d); De = load(f_e);
S = Dd.settings;
assert(S.T == 20 && isequal(S.lya_T_interval, [10 20]) && S.lya_warmup == 5 && S.fs == 400 && ...
    S.K == 5 && S.K_max == 10 && S.lya_dt == 0.05, 'the defaults are the run mode''s table');
assert(isequal(rmfield(Dd.settings, 'minutes'), rmfield(De.settings, 'minutes')));
assert(isequaln(strip_times(Dd.results), strip_times(De.results)), ...
    'naming the defaults explicitly must give the same results');
assert(Dd.results(1).noise_on.K_used >= 5, 'K at or below N is not touched by the cap');

%% 3. End to end, tiny
tmp = tempname; tmpf = tempname; empty_run = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {tmp, tmpf, empty_run}));
conds = {'no_adaptation', 'sfa1_std1'};
cfg = mu7chain_config('essential');
cfg.stages = {'lyapunov_spectrum_n80'};
cfg.sweeps = {};
cfg.eig_heatmap_examples = [];
cfg.run_conditions = conds;
cfg.conditions = conds;
cfg.run_dir  = tmp;
cfg.fig_root = tmpf;
cfg.verbose  = 'near-none';
win = [2 4]; warm = 1;
cfg.lyapunov_spectrum_n80 = struct('n_seeds', 2, 'K', 100, 'T', 4, 'lya_T_interval', win, 'lya_warmup', warm);
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
cfg.figures = cfg.figures(strcmp(names, 'fig_lyapunov_ranked_spectrum'));
assert(isscalar(cfg.figures), 'the essential registry holds the entry');

t0 = tic;
run_dir = run_all_paper_analyses(cfg);
t_run = toc(t0);
assert(strcmp(run_dir, tmp));
M = load(fullfile(tmp, 'run_manifest.mat')); M = M.run_manifest;
assert(isequal(M.stages_run, {'lyapunov_spectrum_n80'}));
mat = fullfile(tmp, 'lyapunov_spectrum_n80', 'lyapunov_spectrum_data.mat');
assert(isfile(mat), 'the stage writes its own folder');
D = load(mat);
S = D.settings;
assert(strcmp(S.preset_name, n80) && S.n == 80 && S.T == 4 && isequal(S.lya_T_interval, win) && ...
    S.lya_warmup == warm && S.fs == 500 && S.K == 100 && S.K_max == 100 && S.n_seeds == 2 && ...
    isequal(S.variants, {'noise_on'}));
assert(isequal({D.results.name}, conds));
N_seen = zeros(1, 2);
for i = 1:2
    runs = D.results(i).noise_on;
    assert(numel(runs) == 2 && isequal(runs(1).seeds, [1 2]) && isequal(runs(2).seeds, [2 3]));
    for s = 1:2
        r = runs(s);
        assert(r.n == 80 && strcmp(r.ode_solver, 'sra1') && r.noise_on && r.sigma_u_noise > 0);
        assert(r.K_used == min(100, r.N) && numel(r.LE_spectrum) == r.K_used, ...
            '%s: K_used %d with N %d', conds{i}, r.K_used, r.N);
        % The segments start at t1 - warm-up and the last starts before t2.
        assert(abs(r.t_lya(1) - (win(1) - warm)) < 0.05 + 1e-9 && r.t_lya(end) < win(2) && ...
            r.t_lya(end) >= win(2) - 0.1, '%s: Lyapunov segments %g to %g s', conds{i}, r.t_lya(1), r.t_lya(end));
    end
    N_seen(i) = runs(1).N;
end
assert(N_seen(1) == 80 && N_seen(2) > 100, 'no adaptation has N = 80 (capped, all taken); one timescale N > K');

t0 = tic;
R = make_all_paper_figures(cfg);
t_fig = toc(t0);
assert(isscalar(R) && R.ok, 'the figure entry failed: %s', R.err);
out_dir = fullfile(tmpf, 'fig_lyapunov_ranked_spectrum');
assert(isfile(fullfile(out_dir, 'Fig_Lyapunov_Ranked_Spectrum.png')) && ...
    isfile(fullfile(out_dir, 'Fig_Lyapunov_Ranked_Spectrum_table.md')), 'the PNG and the table');
tbl = fileread(fullfile(out_dir, 'Fig_Lyapunov_Ranked_Spectrum_table.md'));
assert(contains(tbl, 'bin width') && contains(tbl, 'Freedman-Diaconis') && contains(tbl, '| 80 | 80 |'), ...
    'the table states N, K used and the bin rule');
mkdir(empty_run);
assert(throws_prefix(@() fig_lyapunov_ranked_spectrum('run_dir', empty_run, 'save', false, 'visible', false), ...
    'resolve_data_file:'), 'a run without the stage is an error');

fprintf('PASS: test_lyapunov_spectrum_n80_stage -- the plan in every mode, refusals, the default path, end to end (run %.0f s, figure %.0f s).\n', ...
    t_run, t_fig);

%% ------------------------------------------------------------------------
function v = arg(args, name)
% The value a name-value list gives name (the last one, as a call takes it).
k = find(strcmp(args(1:2:end), name), 1, 'last');
assert(~isempty(k), 'no argument %s', name);
v = args{2 * k};
end

function R = strip_times(R)
% The results without the wall-clock fields.
for i = 1:numel(R)
    for v = {'noise_on', 'noise_off'}
        if ~isfield(R, v{1}) || isempty(R(i).(v{1})); continue; end
        runs = R(i).(v{1});
        R(i).(v{1}) = rmfield(runs, intersect(fieldnames(runs), {'seconds', 'lya_seconds'}));
    end
end
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
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end

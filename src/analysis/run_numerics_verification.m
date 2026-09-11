function mat_file = run_numerics_verification(cfg)
% RUN_NUMERICS_VERIFICATION Measure the precision of the paper's numerics.
%
%   mat_file = RUN_NUMERICS_VERIFICATION()
%   mat_file = RUN_NUMERICS_VERIFICATION('preset_name', p, 'run_mode', 'fast', ...
%                                        'out_dir', d)
%
% The COMPUTE half of three supplemental figures (fig_numerics_verification,
% variants 'solver', 'lya_method' and 'ensemble'). Four sub-experiments, run for
% every adaptation regime the preset states, on several network seeds each:
%
%   A. NOISE-FREE RESHOOT, full-size network. A reference trajectory is
%      integrated with ode45 at RelTol = AbsTol = ref_tol (1e-10). Then, for
%      each step size on fs_ladder, SRA1 is repeatedly RESET to the reference
%      state and integrated for one short segment; the distance from the
%      reference at the segment's end is the discretisation error. This is
%      Benettin's construction with the perturbed trajectory replaced by "the
%      same trajectory under a different integrator", and the reset is what
%      keeps chaos out of the measurement: the reference need only be a true
%      trajectory LOCALLY, which 1e-10 gives even where no trajectory is
%      globally shadowable. Two segment lengths are recorded -- two steps (the
%      shortest span the fixed-step integrators accept; the local truncation
%      error) and lya_dt = 0.02 s (the error over one Benettin segment, the
%      number to compare with the perturbation scale Benettin tracks).
%      At sigma = 0 SRA1's drift is Ralston's second-order RK, so the 0.02 s
%      error should fall 4x per halving of dt. The preset's piecewise
%      activation has a derivative jump at its breakpoints, which ode45 steps
%      around and a fixed-step scheme cannot; crossings contribute first-order
%      error, so a measured slope between 1 and 2 is a finding, not a defect.
%
%   B. NOISY RESHOOT, same network, the same construction with the reference
%      being SRA1 on a grid 8x finer than the top of the ladder, ON THE SAME
%      BROWNIAN PATH. Coarser runs consume the fine path rebuilt by
%      coarsen_noise (increments summed, the I(1,0) area carried with its base
%      correction), which sde_fixed_step's absolute-time indexing then pairs
%      with the right steps. The difference is the strong error; SRA1 is
%      strong order 1.5 for additive noise, so expect ~2.8x per halving. This
%      is the only meaningful precision check WITH noise, since ode45 cannot
%      run the SDE at all. 8x separation from the reference is the minimum
%      test_sde_integrators found necessary: closer, and the reference's own
%      error correlates with the tested run's and flatters the slope.
%
%   L. LLE AGREEMENT, full-size network, noise-free: Benettin with ode45 for
%      both fiducial and perturbed runs vs Benettin with SRA1 for both.
%      Because Benettin uses ONE integrator for both trajectories, most of
%      the integration error cancels in their difference, so the exponents
%      should agree more closely than the raw trajectory error suggests. A
%      short window of x traces from both free runs is kept for an overlay;
%      where the LLE is positive the two must diverge after a few Lyapunov
%      times whatever the precision, which is why A exists.
%
%   C. BENETTIN vs QR on a REDUCED network built from the same preset physics
%      (n_small neurons, F_tracks_network = true so the spectral radius
%      matches the full network). QR integrates an N x N variational system
%      per segment and is out of the question at the preset's ~4000 states;
%      it had only ever been checked on SRNNModel2 (test_benettin_vs_qr). Both
%      methods run on the same fiducial trajectory (same seeds); the largest
%      QR exponent should match the Benettin LLE.
%
% TRIALS, AND WHICH SUB-EXPERIMENTS GET HOW MANY. Each trial is a new network:
% rng_seeds = [k, k+1], which draws W, the stimulus, the initial state and the
% per-neuron setpoints. Two counts:
%
%   n_trials_reshoot  for A and B. The reshoot error is a local quantity with
%                     hundreds of restarts per trial; it varied by ~30%
%                     between seeds and needs only a handful.
%   n_trials_lle      for L and C. The 10 s finite-time LLE of the intermittent
%                     single-timescale regime scatters by +-0.2 from one
%                     trajectory to the next with EITHER integrator (2026-09-10
%                     reports), so a single seed cannot separate integrator
%                     bias from scatter. Paired per-trial values over 25-30
%                     seeds can: bias is a consistent sign, scatter is not.
%
% PARALLEL. The three blocks per condition are parfor loops over trials (L
% over trial x integrator). Every job builds its own model inside the worker
% from the preset and the seed; nothing is shared but the settings struct.
% n_workers defaults to min(12, cores), leaving headroom on a 14-core box;
% per-worker memory peaks at about 1.5 GB (an ode45 run holding 20 s of 4000
% states, or a reshoot holding the reference rows plus a 12800 Hz noise
% tensor), so 12 workers fit in 64 GB. Worker output is interleaved.
%
% Every model is an SRNNNumericsProbe -- SRNNCellTypePairs with the noise
% tensor injectable and a public segment integrator, nothing else. The
% reference trajectories are integrated in one-second chunks and only the rows
% on the top-of-ladder grid inside the reshoot window are kept, so the
% full-rate state (32 KB a row at 4000 states) never sits in memory whole.
%
% THE RESHOOT WINDOW STAYS INSIDE THE MIDDLE THIRD OF THE RUN. The stimulus is
% the preset's three-step pattern (off / on / off); its per-neuron amplitudes
% are drawn once from rng_seeds, but the step EDGES land at T/3 and 2T/3
% rounded to the nearest SAMPLE, and the linear interpolant ramps over one
% sample -- so runs at different fs see slightly different input around each
% edge. The first version of this stage let the window straddle the 2T/3 edge,
% and a single edge-crossing segment dominated the RMS: the noise-free error
% fell 3x from 400 to 800 Hz and then 23x to 1600 Hz, the one rate sharing the
% reference's grid. Between the edges the input is constant on every grid, so
% the window is [T/3 + 0.3, 2T/3 - 0.01] and every run in a sub-experiment
% shares T_range.
%
% Output: <out_dir>/numerics_verification_data.mat with results(i).reshoot
% (1 x n_trials_reshoot), results(i).lle and results(i).qr (1 x n_trials_lle)
% and results(i).summary (per-trial scalars). A standalone run lands in
% <root>/data/numerics_verification; the pipeline passes
% <run_dir>/numerics_verification (run_all_paper_analyses).
%
% See also: fig_numerics_verification, SRNNNumericsProbe, coarsen_noise,
%           sde_fixed_step, test_benettin_vs_qr, run_eig_heatmap

arguments
    cfg.preset_name      (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25'
    cfg.run_mode         (1,:) char    = 'production'
    cfg.out_dir          (1,:) char    = ''
    cfg.n_trials_reshoot (1,1) double  = 0     % 0 -> per run_mode
    cfg.n_trials_lle     (1,1) double  = 0     % 0 -> per run_mode
    cfg.n_workers        (1,1) double  = 0     % 0 -> min(12, cores)
    cfg.license_wait_s   (1,1) double  = 4 * 3600   % how long to wait for a PCT licence seat before going serial
    cfg.license_poll_s   (1,1) double  = 20         % seconds between licence checks
end

setup_paths();

%% Cost table -- BEFORE any compute, so a bad mode fails in milliseconds
% test_run_modes asserts the :badMode identifier and that the message names
% every mode in run_mode_names().
switch cfg.run_mode
    case 'fast'
        T_free = 6;   T_lle = 10;  T_noisy = 3;    n_small = 30;  T_small = [-8, 4];   max_restarts = 400;  n_reshoot = 1;  n_lle = 2;
    case {'medium', 'medium2'}
        T_free = 12;  T_lle = 20;  T_noisy = 4.5;  n_small = 40;  T_small = [-10, 10]; max_restarts = 800;  n_reshoot = 5;  n_lle = 25;
    case 'production'
        T_free = 20;  T_lle = 40;  T_noisy = 4.5;  n_small = 60;  T_small = [-10, 20]; max_restarts = 1500; n_reshoot = 5;  n_lle = 30;
    otherwise
        error('run_numerics_verification:badMode', ...
            'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.n_trials_reshoot > 0; n_reshoot = cfg.n_trials_reshoot; end
if cfg.n_trials_lle     > 0; n_lle     = cfg.n_trials_lle;     end

P = struct();                        % everything a worker needs, broadcast once
P.preset_name  = cfg.preset_name;
P.fs_ladder    = [400, 800, 1600];   % the paper runs at 400; the rest test convergence
P.fs_top       = max(P.fs_ladder);
P.fs_noisy_ref = 8 * P.fs_top;       % 12800: the noisy reference grid
P.ref_tol      = 1e-10;              % ode45 RelTol = AbsTol for the noise-free reference
P.seg_long     = 0.02;               % s; = Benettin's default lya_dt
P.seg_short    = 2;                  % steps; the shortest span sde_fixed_step accepts
P.fs_lle       = 400;                % the paper's rate, for L
P.T_free       = T_free;
P.T_noisy      = T_noisy;
P.T_lle        = T_lle;
P.n_small      = n_small;
P.T_small      = T_small;
P.lya_small_warmup = -T_small(1);    % iterate from T_small(1), accumulate from 0
P.max_restarts = max_restarts;
P.win_free  = [T_free  / 3 + 0.3, 2 * T_free  / 3 - 0.01];   % inside the stim-on third, see header
P.win_noisy = [T_noisy / 3 + 0.3, 2 * T_noisy / 3 - 0.01];

if isempty(cfg.out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'numerics_verification');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir); mkdir(out_dir); end

[preset, model_class, conditions] = srnn_param_preset(cfg.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_numerics_verification:WrongClass', ...
        ['Preset ''%s'' is written for %s; this stage probes SRNNCellTypePairs ' ...
         '(via SRNNNumericsProbe) only.'], cfg.preset_name, model_class);
end
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
titles     = cellfun(@(n) pretty(n), cond_names, 'UniformOutput', false);
n_cond     = numel(cond_names);
P.sigma_preset = 0;
if isfield(preset, 'sigma_u_noise'); P.sigma_preset = preset.sigma_u_noise; end

fprintf('[numerics_verification] preset=%s run_mode=%s\n', cfg.preset_name, cfg.run_mode);
fprintf('  A: noise-free reshoot, ode45 @ %g vs sra1 @ fs %s, T = %g s      x %d trials\n', ...
    P.ref_tol, mat2str(P.fs_ladder), T_free, n_reshoot);
fprintf('  B: noisy reshoot, sra1 @ %d vs sra1 @ fs %s on one path, T = %g s  x %d trials\n', ...
    P.fs_noisy_ref, mat2str(P.fs_ladder), T_noisy, n_reshoot);
fprintf('  L: Benettin LLE, ode45 vs sra1 @ %d Hz, T = %g s                     x %d trials\n', ...
    P.fs_lle, T_lle, n_lle);
fprintf('  C: Benettin vs QR, n = %d, T = %s                              x %d trials\n', ...
    n_small, mat2str(T_small), n_lle);

pool = ensure_pool(cfg.n_workers, cfg.license_wait_s, cfg.license_poll_s);
fprintf('  parallel pool: %d workers\n', pool.NumWorkers);

t_stage = tic;
res = struct('name', cond_names, 'title', titles, 'reshoot', [], 'lle', [], 'qr', [], 'summary', []);

for i = 1:n_cond
    cname = cond_names{i};
    fprintf('\n=== %d/%d %s ===\n', i, n_cond, titles{i});

    %% A + B. reshoot, one trial per worker -------------------------------
    t0 = tic;
    fprintf('  [A/B] reshoot on %d seed(s) ...\n', n_reshoot);
    R = cell(1, n_reshoot);
    parfor k = 1:n_reshoot
        R{k} = reshoot_trial(P, cname, [k, k + 1]);
    end
    res(i).reshoot = [R{:}];
    fprintf('  [A/B] done in %.0f s\n', toc(t0));

    %% L. LLE agreement, trial x integrator, one job per worker ---------------
    t0 = tic;
    fprintf('  [L] Benettin ode45 vs sra1 on %d seed(s) ...\n', n_lle);
    solvers = {'ode45', 'sra1'};
    J = cell(1, 2 * n_lle);
    parfor j = 1:2 * n_lle
        k = ceil(j / 2);
        s = solvers{2 - mod(j, 2)};        % j odd -> ode45, even -> sra1
        J{j} = lle_trial(P, cname, [k, k + 1], s);
    end
    L = cell(1, n_lle);
    for k = 1:n_lle
        L{k} = struct('seeds', [k, k + 1], 'ode45', J{2 * k - 1}, 'sra1', J{2 * k});
    end
    res(i).lle = [L{:}];
    fprintf('  [L] done in %.0f s\n', toc(t0));

    %% C. Benettin vs QR on the reduced network, one trial per worker --------
    t0 = tic;
    fprintf('  [C] Benettin vs QR on n = %d, %d seed(s) ...\n', n_small, n_lle);
    Q = cell(1, n_lle);
    parfor k = 1:n_lle
        Q{k} = qr_trial(P, cname, [k, k + 1]);
    end
    res(i).qr = [Q{:}];
    fprintf('  [C] done in %.0f s\n', toc(t0));

    res(i).summary = summarise(res(i).reshoot, res(i).lle, res(i).qr, P.fs_lle);
    print_summary(res(i).summary, titles{i});
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'model_class', model_class, 'n', preset.n, 'sigma_u_noise', P.sigma_preset, ...
    'fs_ladder', P.fs_ladder, 'fs_noisy_ref', P.fs_noisy_ref, 'ref_tol', P.ref_tol, ...
    'seg_long', P.seg_long, 'seg_short_steps', P.seg_short, ...
    'T_free', T_free, 'win_free', P.win_free, ...
    'T_noisy', T_noisy, 'win_noisy', P.win_noisy, ...
    'T_lle', T_lle, 'fs_lle', P.fs_lle, ...
    'n_small', n_small, 'T_small', T_small, 'lya_small_warmup', P.lya_small_warmup, ...
    'max_restarts', max_restarts, 'block_names', {{'x', 'a', 'b'}}, ...
    'n_trials_reshoot', n_reshoot, 'n_trials_lle', n_lle, 'n_workers', pool.NumWorkers, ...
    'minutes', toc(t_stage) / 60);

condition_titles = titles;  % saved name
results = res;
mat_file = fullfile(out_dir, 'numerics_verification_data.mat');
save(mat_file, 'results', 'cond_names', 'condition_titles', 'settings', '-v7.3');
fprintf('\nSaved: %s  (%.1f min)\n', mat_file, settings.minutes);
end

%% ------------------------------------------------------------------------
function pool = ensure_pool(n_workers, wait_s, poll_s)
% A pool with n_workers (default min(12, cores)); reuse one of the right size.
%
% The Parallel Computing Toolbox is a 15-seat network licence here, and on
% 2026-09-11 every seat was taken (none by us). A pool needs ONE seat, for the
% client, not one per worker. So: ask the licence server how many seats are
% free (scripts/tools/pct_licenses.ps1, whose exit code is the free count),
% and while it says none, wait poll_s and ask again, for up to wait_s. When a
% seat is free -- or the script cannot answer -- try parpool; a failed attempt
% (someone else got the seat first) just goes back to polling. Past wait_s,
% warn and return a stand-in with NumWorkers = 1: the parfor loops then run
% serially, which is slow but correct.
if n_workers <= 0
    n_workers = min(12, feature('numcores'));
end
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers == n_workers
    return
end
if ~isempty(pool); delete(pool); end

% The licence polling below shells out to lmutil at a path on ONE machine and
% ties up a workstation for hours; it is only meant to run there (TR,
% 2026-09-11). Anywhere else, stop rather than guess.
host = getenv('COMPUTERNAME');
if ~strcmpi(host, 'R5456622')
    error('run_numerics_verification:WrongHost', ...
        ['This stage''s parallel-pool and licence handling is written for ' ...
         'R5456622; this is %s. Run it there, or edit ensure_pool.'], host);
end

t_wait = tic;
n_polls = 0;
while true
    free = pct_free_seats();
    if isnan(free) || free > 0
        try
            pool = parpool(parallel.defaultProfile, n_workers);
            return
        catch ME
            last_err = ME.message;
        end
    else
        last_err = sprintf('licence server reports 0 free seats');
    end
    if toc(t_wait) > wait_s
        warning('run_numerics_verification:NoPool', ...
            ['No parallel pool after %.0f min (%s). Running the trial loops ' ...
             'SERIALLY; expect roughly n_workers times the wall time.'], ...
            toc(t_wait) / 60, strtok(last_err, newline));
        pool = struct('NumWorkers', 1);
        return
    end
    n_polls = n_polls + 1;
    if n_polls == 1 || mod(n_polls, 15) == 0
        fprintf('  [pool] waiting for a PCT licence seat (%s); polling every %d s, %.0f min so far\n', ...
            strtok(last_err, newline), poll_s, toc(t_wait) / 60);
    end
    pause(poll_s);
end
end

function free = pct_free_seats()
% Free Parallel Computing Toolbox seats on the network licence, or NaN if the
% question cannot be answered (not Windows, script or lmutil missing).
free = NaN;
if ~ispc; return; end
script = fullfile(fileparts(which('setup_paths')), 'scripts', 'tools', 'pct_licenses.ps1');
if ~isfile(script); return; end
[status, ~] = system(sprintf( ...
    'powershell -NoProfile -ExecutionPolicy Bypass -File "%s" -Quiet', script));
if status >= 0 && status < 200      % the script exits with the free count; 254/255 are its errors
    free = status;
end
end

function T = reshoot_trial(P, cname, seeds)
% Sub-experiments A and B on one network. Returns seeds / free / noisy.
T = struct('seeds', seeds, 'free', [], 'noisy', []);

% A. noise-free: ode45 reference at the top of the ladder
ref = build_probe(P.preset_name, cname, 'rng_seeds', seeds, ...
    'sigma_u_noise', 0, 'ode_solver', 'ode45', 'fs', P.fs_top, ...
    'T_range', [0, P.T_free], 'lya_method', 'none');
ref.ode_opts = odeset('RelTol', P.ref_tol, 'AbsTol', P.ref_tol, 'MaxStep', 1 / P.fs_top);
[t_ref, S_ref] = integrate_reference(ref, 1, P.win_free);
S_free_end = S_ref(end, :)';           % settled state; seeds the noisy reference
blocks = state_blocks(ref.cached_params.state_layout);

free = cell(1, numel(P.fs_ladder));
for j = 1:numel(P.fs_ladder)
    fs = P.fs_ladder(j);
    m  = P.fs_top / fs;
    probe = build_probe(P.preset_name, cname, 'rng_seeds', seeds, ...
        'sigma_u_noise', 0, 'ode_solver', 'sra1', 'fs', fs, ...
        'T_range', [0, P.T_free], 'lya_method', 'none');
    probe.arm_noise();                 % no-op at sigma = 0; keeps the call shape uniform
    free{j} = reshoot(probe, t_ref(1:m:end), S_ref(1:m:end, :), ...
        blocks, P.seg_short, P.seg_long, P.max_restarts);
end
T.free = [free{:}];
fprintf('    [A] seeds %s: sra1 @ %d Hz |err| over %g s = %.3e, slope %.2f\n', ...
    mat2str(seeds), P.fs_lle, P.seg_long, ...
    T.free([T.free.fs] == P.fs_lle).err_long_total_rms, ...
    fit_slope([T.free.fs], [T.free.err_long_total_rms]));

% B. noisy: sra1 reference on the seeded path, coarser runs on the same path
if P.sigma_preset > 0
    ref = build_probe(P.preset_name, cname, 'rng_seeds', seeds, ...
        'ode_solver', 'sra1', 'fs', P.fs_noisy_ref, ...
        'T_range', [0, P.T_noisy], 'lya_method', 'none');
    ref.arm_noise();                                   % seeded draw at 12800 Hz
    nz_fine = ref.consumed_noise;
    [t_ref, S_ref] = integrate_reference(ref, P.fs_noisy_ref / P.fs_top, P.win_noisy, S_free_end);
    ref.disarm_noise();

    noisy = cell(1, numel(P.fs_ladder));
    for j = 1:numel(P.fs_ladder)
        fs = P.fs_ladder(j);
        m  = P.fs_noisy_ref / fs;
        [xi1_c, xi2_c] = coarsen_noise(nz_fine.xi1, nz_fine.xi2, 1 / P.fs_noisy_ref, m);
        probe = build_probe(P.preset_name, cname, 'rng_seeds', seeds, ...
            'ode_solver', 'sra1', 'fs', fs, ...
            'T_range', [0, P.T_noisy], 'lya_method', 'none');
        probe.set_noise(struct('xi1', xi1_c, 'xi2', xi2_c, 't0', 0, 'fs', fs, ...
            'sigma', 0, 'idx', []));                   % sigma/idx filled by the probe
        probe.arm_noise();
        step = P.fs_top / fs;
        noisy{j} = reshoot(probe, t_ref(1:step:end), S_ref(1:step:end, :), ...
            blocks, P.seg_short, P.seg_long, P.max_restarts);
        probe.disarm_noise();
    end
    T.noisy = [noisy{:}];
    fprintf('    [B] seeds %s: sra1 @ %d Hz |err| over %g s = %.3e, slope %.2f\n', ...
        mat2str(seeds), P.fs_lle, P.seg_long, ...
        T.noisy([T.noisy.fs] == P.fs_lle).err_long_total_rms, ...
        fit_slope([T.noisy.fs], [T.noisy.err_long_total_rms]));
end
end

function out = lle_trial(P, cname, seeds, solver)
% Sub-experiment L for one network and one integrator.
model = build_probe(P.preset_name, cname, 'rng_seeds', seeds, ...
    'sigma_u_noise', 0, 'ode_solver', solver, 'fs', P.fs_lle, ...
    'T_range', [0, P.T_lle], 'lya_method', 'benettin', ...
    'lya_T_interval', [P.T_lle / 2, P.T_lle], 'plot_deci', 4);
if strcmp(solver, 'ode45')
    model.ode_opts = odeset('RelTol', P.ref_tol, 'AbsTol', P.ref_tol, 'MaxStep', 1 / P.fs_lle);
end
t0 = tic;
evalc('model.run();');
pd = model.plot_data;
E  = model.cell_type_names{1};
out = struct('LLE', model.lya_results.LLE, ...
    't_lya', model.lya_results.t_lya, 'local_lya', model.lya_results.local_lya, ...
    't', pd.t, 'x_examples', pd.x.(E)(1:min(3, end), :), ...
    'mean_rate', mean(pd.r.(E)(:)), 'seconds', toc(t0));
fprintf('    [L] seeds %s %-5s: LLE = %+.4f (%.0f s)\n', mat2str(seeds), solver, out.LLE, out.seconds);
end

function Q = qr_trial(P, cname, seeds)
% Sub-experiment C for one reduced network: Benettin then QR on the same
% fiducial trajectory.
Q = struct('seeds', seeds, 'n', P.n_small, 'T_range', P.T_small, 'benettin', [], 'qr', []);
for method = {'benettin', 'qr'}
    mth = method{1};
    model = build_probe(P.preset_name, cname, 'rng_seeds', seeds, ...
        'n', P.n_small, 'indegree', max(2, round(0.2 * P.n_small)), ...
        'F_tracks_network', true, ...
        'sigma_u_noise', 0, 'ode_solver', 'ode45', 'fs', P.fs_lle, ...
        'T_range', P.T_small, 'lya_method', mth, ...
        'lya_T_interval', [0, P.T_small(2)], 'lya_warmup', P.lya_small_warmup);
    t0 = tic;
    evalc('model.run();');
    r = model.lya_results;
    if strcmp(mth, 'benettin')
        Q.benettin = struct('LLE', r.LLE, 't_lya', r.t_lya, 'local_lya', r.local_lya, 'seconds', toc(t0));
    else
        Q.qr = struct('LE_spectrum', r.LE_spectrum, 't_lya', r.t_lya, ...
            'local_LE_spectrum_t', r.local_LE_spectrum_t, ...
            'N_sys_eqs', model.N_sys_eqs, 'seconds', toc(t0));
    end
end
fprintf('    [C] seeds %s: Benettin %+.4f, QR %+.4f of %d (%.0f s)\n', mat2str(seeds), ...
    Q.benettin.LLE, Q.qr.LE_spectrum(1), numel(Q.qr.LE_spectrum), Q.qr.seconds);
end

function probe = build_probe(preset_name, condition_name, varargin)
% build_from_preset's twin that instantiates SRNNNumericsProbe. Same
% precedence: preset < condition < overrides, ode_solver among the overrides.
[preset, ~, conditions] = srnn_param_preset(preset_name);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond  = rmfield(conditions{strcmp(names, condition_name)}, 'name');
args  = [struct2namevalue(preset), struct2namevalue(cond), varargin];
probe = SRNNNumericsProbe(args{:});
evalc('probe.build();');
end

function [t_keep, S_keep] = integrate_reference(probe, keep_every, window, S0)
% Integrate the probe's whole time vector in one-second chunks with its own
% integrator, keeping only every keep_every-th row inside window. The
% first kept row is aligned to the probe's grid origin, so subsampling the
% result by any divisor of keep_every's cofactors stays on-grid.
t  = probe.t_ex;
nt = numel(t);
if nargin < 4 || isempty(S0); S0 = probe.S0; end
keep = find(mod(0:nt - 1, keep_every) == 0 & t(:)' >= window(1) - 1e-9 & t(:)' <= window(2) + 1e-9);
S_keep = zeros(numel(keep), numel(S0));
t_keep = t(keep);
chunk = round(probe.fs);          % one second of steps
S  = S0(:);
i0 = 1;
p  = 1;
while i0 < nt
    i1 = min(i0 + chunk, nt);
    [~, Sc] = probe.integrate_segment(t(i0:i1), S);
    last = i1 - 1; if i1 == nt; last = nt; end
    sel = keep(keep >= i0 & keep <= last);
    S_keep(p:p + numel(sel) - 1, :) = Sc(sel - i0 + 1, :);
    p  = p + numel(sel);
    S  = Sc(end, :)';
    i0 = i1;
end
end

function blocks = state_blocks(layout)
% Index sets for the three state families; an absent family is empty.
blocks = {layout.x, [layout.a{:}], [layout.b{:}]};
end

function R = reshoot(probe, t_ref, S_ref, blocks, seg_short, seg_long, max_restarts)
% Reset-and-integrate against a reference on the probe's own grid.
%   t_ref / S_ref  reference rows, spaced 1/probe.fs, on probe.t_ex
%   seg_short      segment length in STEPS (>= 2)
%   seg_long       segment length in SECONDS (an integer number of steps)
% Returns per-restart per-block RMS errors for both segment lengths, the
% reference's per-block RMS over the window (for a relative error), and the
% totals.
fs = probe.fs;
L_long = round(seg_long * fs);
if abs(L_long - seg_long * fs) > 1e-9
    error('run_numerics_verification:SegmentOffGrid', ...
        'seg_long = %g s is not an integer number of steps at fs = %g.', seg_long, fs);
end
n_ref = size(S_ref, 1);
t_ex  = probe.t_ex;

R = struct('fs', fs, 'seg_short_steps', seg_short, 'seg_long', seg_long);
R.block_rms = cellfun(@(b) rms_or_nan(S_ref(:, b)), blocks);

for which = {'short', 'long'}
    if strcmp(which{1}, 'short'); L = seg_short; else; L = L_long; end
    starts = 1:L:(n_ref - L);
    if numel(starts) > max_restarts
        starts = starts(round(linspace(1, numel(starts), max_restarts)));
    end
    n_s = numel(starts);
    err_blocks = nan(n_s, numel(blocks));
    err_total  = nan(n_s, 1);
    for k = 1:n_s
        s  = starts(k);
        i0 = round(t_ref(s) * fs) + 1;                 % row of t_ref(s) in t_ex
        if abs(t_ex(i0) - t_ref(s)) > 1e-9
            error('run_numerics_verification:ReferenceOffGrid', ...
                'Reference time %g is not on the probe''s grid at fs = %g.', t_ref(s), fs);
        end
        [~, Sseg] = probe.integrate_segment(t_ex(i0:i0 + L), S_ref(s, :)');
        d = Sseg(end, :) - S_ref(s + L, :);
        err_total(k)     = sqrt(mean(d.^2));
        err_blocks(k, :) = cellfun(@(b) rms_or_nan(d(b)), blocks);
    end
    R.(['t_' which{1}])          = t_ref(starts);
    R.(['err_' which{1}])        = err_blocks;
    R.(['err_' which{1} '_total']) = err_total;
    R.(['err_' which{1} '_rms'])   = sqrt(mean(err_blocks.^2, 1, 'omitnan'));
    R.(['err_' which{1} '_total_rms']) = sqrt(mean(err_total.^2));
end
end

function v = rms_or_nan(x)
if isempty(x); v = NaN; else; v = sqrt(mean(x(:).^2)); end
end

function s = pretty(name)
st = manuscript_style();
if st.condition_title.isKey(name)
    s = st.condition_title(name);
else
    s = strrep(name, '_', ' ');
end
end


function S = summarise(reshoot, lle, qr, fs_lle)
% Per-trial scalars for the ensemble figure and the report. LLE and QR rows
% are indexed by the n_trials_lle seeds; reshoot rows by the n_trials_reshoot
% seeds (the first n_trials_reshoot of the same sequence).
S = struct('n_trials_lle', numel(lle), 'n_trials_reshoot', numel(reshoot), ...
    'seeds_lle', vertcat(lle.seeds), 'seeds_reshoot', vertcat(reshoot.seeds));
S.lle_ode45   = arrayfun(@(t) t.ode45.LLE, lle);
S.lle_sra1    = arrayfun(@(t) t.sra1.LLE, lle);
S.qr_benettin = arrayfun(@(t) t.benettin.LLE, qr);
S.qr_lambda1  = arrayfun(@(t) t.qr.LE_spectrum(1), qr);
j = find([reshoot(1).free.fs] == fs_lle, 1); if isempty(j); j = 1; end
S.fs_paper        = reshoot(1).free(j).fs;
S.err_free_paper  = arrayfun(@(t) t.free(j).err_long_total_rms, reshoot);
S.slope_free      = arrayfun(@(t) fit_slope([t.free.fs], [t.free.err_long_total_rms]), reshoot);
if ~isempty(reshoot(1).noisy)
    S.err_noisy_paper = arrayfun(@(t) t.noisy(j).err_long_total_rms, reshoot);
    S.slope_noisy     = arrayfun(@(t) fit_slope([t.noisy.fs], [t.noisy.err_long_total_rms]), reshoot);
else
    S.err_noisy_paper = nan(1, numel(reshoot));
    S.slope_noisy     = nan(1, numel(reshoot));
end
end

function p = fit_slope(fs, err)
c = polyfit(log(1 ./ fs), log(err), 1);
p = c(1);
end

function print_summary(S, title)
d  = S.lle_sra1 - S.lle_ode45;
dq = S.qr_lambda1 - S.qr_benettin;
fprintf('\n  summary for %s:\n', title);
fprintf('    LLE over %d seeds: ode45 %+.3f +- %.3f | sra1 %+.3f +- %.3f | paired sra1-ode45 %+.3f +- %.3f\n', ...
    S.n_trials_lle, mean(S.lle_ode45), std(S.lle_ode45), mean(S.lle_sra1), std(S.lle_sra1), mean(d), std(d));
fprintf('    reduced net over %d seeds: Benettin %+.3f +- %.3f | QR %+.3f +- %.3f | paired QR-Benettin %+.4f +- %.4f\n', ...
    S.n_trials_lle, mean(S.qr_benettin), std(S.qr_benettin), mean(S.qr_lambda1), std(S.qr_lambda1), mean(dq), std(dq));
fprintf('    reshoot @ %d Hz over %d seeds: free %.2e +- %.1e (slope %.2f +- %.2f)\n', S.fs_paper, ...
    S.n_trials_reshoot, mean(S.err_free_paper), std(S.err_free_paper), mean(S.slope_free), std(S.slope_free));
if all(isfinite(S.err_noisy_paper))
    fprintf('                                  noisy %.2e +- %.1e (slope %.2f +- %.2f)\n', ...
        mean(S.err_noisy_paper), std(S.err_noisy_paper), mean(S.slope_noisy), std(S.slope_noisy));
end
end

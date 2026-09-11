function mat_file = run_numerics_verification(cfg)
% RUN_NUMERICS_VERIFICATION Measure the precision of the paper's numerics.
%
%   mat_file = RUN_NUMERICS_VERIFICATION()
%   mat_file = RUN_NUMERICS_VERIFICATION('preset_name', p, 'run_mode', 'fast', ...
%                                        'out_dir', d)
%
% The COMPUTE half of two supplemental figures (fig_numerics_verification,
% variants 'solver' and 'lya_method'). Three sub-experiments, run for every
% adaptation regime the preset states, all from the same network seed:
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
% TRIALS. Every sub-experiment is repeated n_trials times (per run mode, or the
% n_trials argument), each trial a new network: rng_seeds = [k, k+1], which
% draws W, the stimulus, the initial state and the per-neuron setpoints. Within
% a trial all sub-experiments share that network. The finite-time LLE of an
% intermittent regime scatters by +-0.2 over 10 s from one trajectory to the
% next with EITHER integrator (see the 2026-09-10 reports), so a single-trial
% integrator comparison cannot distinguish bias from scatter; the paired
% per-trial values in results(i).summary can.
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
% Output: <out_dir>/numerics_verification_data.mat. A standalone run lands in
% <root>/data/numerics_verification; the pipeline passes
% <run_dir>/numerics_verification (run_all_paper_analyses).
%
% See also: fig_numerics_verification, SRNNNumericsProbe, coarsen_noise,
%           sde_fixed_step, test_benettin_vs_qr, run_eig_heatmap

arguments
    cfg.preset_name (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25'
    cfg.run_mode    (1,:) char    = 'production'
    cfg.out_dir     (1,:) char    = ''
    cfg.n_trials    (1,1) double  = 0     % 0 -> per run_mode; each trial is a new network seed
end

setup_paths();

%% Cost table -- BEFORE any compute, so a bad mode fails in milliseconds
% test_run_modes asserts the :badMode identifier and that the message names
% every mode in run_mode_names().
switch cfg.run_mode
    case 'fast'
        T_free = 6;   T_lle = 10;  T_noisy = 3;    n_small = 30;  T_small = [-8, 4];   max_restarts = 400;  n_trials = 1;
    case {'medium', 'medium2'}
        T_free = 12;  T_lle = 20;  T_noisy = 4.5;  n_small = 40;  T_small = [-10, 10]; max_restarts = 800;  n_trials = 3;
    case 'production'
        T_free = 20;  T_lle = 40;  T_noisy = 4.5;  n_small = 60;  T_small = [-10, 20]; max_restarts = 1500; n_trials = 5;
    otherwise
        error('run_numerics_verification:badMode', ...
            'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.n_trials > 0; n_trials = cfg.n_trials; end
fs_ladder    = [400, 800, 1600];     % the paper runs at 400; the rest test convergence
fs_top       = max(fs_ladder);
fs_noisy_ref = 8 * fs_top;           % 12800: the noisy reference grid
ref_tol      = 1e-10;                % ode45 RelTol = AbsTol for the noise-free reference
seg_long     = 0.02;                 % s; = Benettin's default lya_dt
seg_short    = 2;                    % steps; the shortest span sde_fixed_step accepts
fs_lle       = 400;                  % the paper's rate, for L
lya_small_warmup = -T_small(1);      % iterate from T_small(1), accumulate from 0
win_free  = [T_free  / 3 + 0.3, 2 * T_free  / 3 - 0.01];   % inside the stim-on third, see header
win_noisy = [T_noisy / 3 + 0.3, 2 * T_noisy / 3 - 0.01];

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
sigma_preset = 0;
if isfield(preset, 'sigma_u_noise'); sigma_preset = preset.sigma_u_noise; end

fprintf('[numerics_verification] preset=%s run_mode=%s\n', cfg.preset_name, cfg.run_mode);
fprintf('  A: noise-free reshoot, ode45 @ %g vs sra1 @ fs %s, T = %g s\n', ...
    ref_tol, mat2str(fs_ladder), T_free);
fprintf('  B: noisy reshoot, sra1 @ %d vs sra1 @ fs %s on one path, T = %g s\n', ...
    fs_noisy_ref, mat2str(fs_ladder), T_noisy);
fprintf('  L: Benettin LLE, ode45 vs sra1 @ %d Hz, T = %g s\n', fs_lle, T_lle);
fprintf('  C: Benettin vs QR, n = %d, T = %s\n', n_small, mat2str(T_small));

t_stage = tic;
res = struct('name', cond_names, 'title', titles, 'trials', [], 'summary', []);

for i = 1:n_cond
    cname = cond_names{i};
    fprintf('\n=== %d/%d %s ===\n', i, n_cond, titles{i});

    trials = cell(1, n_trials);
    for k = 1:n_trials
        seeds = [k, k + 1];                 % one network per trial, shared by all sub-experiments
        trial = struct('seeds', seeds, 'free', [], 'noisy', [], 'lle', [], 'qr', []);
        fprintf('  -- trial %d/%d, rng_seeds %s --\n', k, n_trials, mat2str(seeds));

        %% A. noise-free reshoot ------------------------------------------------
        fprintf('  [A] ode45 reference at %d Hz, tol %g ...', fs_top, ref_tol);
        t0 = tic;
        ref = build_probe(cfg.preset_name, cname, 'rng_seeds', seeds, ...
            'sigma_u_noise', 0, 'ode_solver', 'ode45', 'fs', fs_top, ...
            'T_range', [0, T_free], 'lya_method', 'none');
        ref.ode_opts = odeset('RelTol', ref_tol, 'AbsTol', ref_tol, 'MaxStep', 1 / fs_top);
        [t_ref, S_ref] = integrate_reference(ref, 1, win_free);
        S_free_end = S_ref(end, :)';           % settled state; seeds the noisy reference
        layout = ref.cached_params.state_layout;
        blocks = state_blocks(layout);
        fprintf(' %.0f s\n', toc(t0));

        free = cell(1, numel(fs_ladder));
        for j = 1:numel(fs_ladder)
            fs = fs_ladder(j);
            m  = fs_top / fs;
            probe = build_probe(cfg.preset_name, cname, 'rng_seeds', seeds, ...
                'sigma_u_noise', 0, 'ode_solver', 'sra1', 'fs', fs, ...
                'T_range', [0, T_free], 'lya_method', 'none');
            probe.arm_noise();                 % no-op at sigma = 0; keeps the call shape uniform
            free{j} = reshoot(probe, t_ref(1:m:end), S_ref(1:m:end, :), ...
                blocks, seg_short, seg_long, max_restarts);
            fprintf('  [A] sra1 @ %4d Hz: |err| over %g s = %.3e (x %.2e, a %.2e, b %.2e)\n', ...
                fs, seg_long, free{j}.err_long_total_rms, free{j}.err_long_rms);
        end
        trial.free = [free{:}];

        %% L. LLE agreement, free-running ---------------------------------------
        lle = struct();
        for solver = {'ode45', 'sra1'}
            s = solver{1};
            fprintf('  [L] Benettin with %s ...', s);
            t0 = tic;
            model = build_probe(cfg.preset_name, cname, 'rng_seeds', seeds, ...
                'sigma_u_noise', 0, 'ode_solver', s, 'fs', fs_lle, ...
                'T_range', [0, T_lle], 'lya_method', 'benettin', ...
                'lya_T_interval', [T_lle / 2, T_lle], 'plot_deci', 4);
            if strcmp(s, 'ode45')
                model.ode_opts = odeset('RelTol', ref_tol, 'AbsTol', ref_tol, 'MaxStep', 1 / fs_lle);
            end
            evalc('model.run();');
            pd = model.plot_data;
            E  = model.cell_type_names{1};
            lle.(s) = struct('LLE', model.lya_results.LLE, ...
                't_lya', model.lya_results.t_lya, 'local_lya', model.lya_results.local_lya, ...
                't', pd.t, 'x_examples', pd.x.(E)(1:min(3, end), :), ...
                'mean_rate', mean(pd.r.(E)(:)));
            fprintf(' LLE = %+.4f (%.0f s)\n', lle.(s).LLE, toc(t0));
        end
        trial.lle = lle;

        %% B. noisy reshoot -----------------------------------------------------
        if sigma_preset > 0
            fprintf('  [B] sra1 reference at %d Hz on the seeded path ...', fs_noisy_ref);
            t0 = tic;
            ref = build_probe(cfg.preset_name, cname, 'rng_seeds', seeds, ...
                'ode_solver', 'sra1', 'fs', fs_noisy_ref, ...
                'T_range', [0, T_noisy], 'lya_method', 'none');
            ref.arm_noise();                                   % seeded draw at 12800 Hz
            nz_fine = ref.consumed_noise;
            [t_ref, S_ref] = integrate_reference(ref, fs_noisy_ref / fs_top, win_noisy, S_free_end);
            ref.disarm_noise();
            fprintf(' %.0f s\n', toc(t0));

            noisy = cell(1, numel(fs_ladder));
            for j = 1:numel(fs_ladder)
                fs = fs_ladder(j);
                m  = fs_noisy_ref / fs;
                [xi1_c, xi2_c] = coarsen_noise(nz_fine.xi1, nz_fine.xi2, 1 / fs_noisy_ref, m);
                probe = build_probe(cfg.preset_name, cname, 'rng_seeds', seeds, ...
                    'ode_solver', 'sra1', 'fs', fs, ...
                    'T_range', [0, T_noisy], 'lya_method', 'none');
                probe.set_noise(struct('xi1', xi1_c, 'xi2', xi2_c, 't0', 0, 'fs', fs, ...
                    'sigma', 0, 'idx', []));                   % sigma/idx filled by the probe
                probe.arm_noise();
                step = fs_top / fs;
                noisy{j} = reshoot(probe, t_ref(1:step:end), S_ref(1:step:end, :), ...
                    blocks, seg_short, seg_long, max_restarts);
                probe.disarm_noise();
                fprintf('  [B] sra1 @ %4d Hz: |err| over %g s = %.3e (x %.2e, a %.2e, b %.2e)\n', ...
                    fs, seg_long, noisy{j}.err_long_total_rms, noisy{j}.err_long_rms);
            end
            clear nz_fine xi1_c xi2_c
            trial.noisy = [noisy{:}];
        else
            fprintf('  [B] skipped: preset is deterministic (sigma_u_noise = 0)\n');
        end

        %% C. Benettin vs QR on the reduced network ------------------------------
        qrc = struct('n', n_small, 'T_range', T_small);
        for method = {'benettin', 'qr'}
            mth = method{1};
            fprintf('  [C] %s on n = %d ...', mth, n_small);
            t0 = tic;
            model = build_probe(cfg.preset_name, cname, 'rng_seeds', seeds, ...
                'n', n_small, 'indegree', max(2, round(0.2 * n_small)), ...
                'F_tracks_network', true, ...
                'sigma_u_noise', 0, 'ode_solver', 'ode45', 'fs', fs_lle, ...
                'T_range', T_small, 'lya_method', mth, ...
                'lya_T_interval', [0, T_small(2)], 'lya_warmup', lya_small_warmup);
            evalc('model.run();');
            r = model.lya_results;
            if strcmp(mth, 'benettin')
                qrc.benettin = struct('LLE', r.LLE, 't_lya', r.t_lya, 'local_lya', r.local_lya);
                fprintf(' LLE = %+.4f (%.0f s)\n', r.LLE, toc(t0));
            else
                qrc.qr = struct('LE_spectrum', r.LE_spectrum, 't_lya', r.t_lya, ...
                    'local_LE_spectrum_t', r.local_LE_spectrum_t, ...
                    'N_sys_eqs', model.N_sys_eqs);
                fprintf(' largest = %+.4f of %d (%.0f s)\n', r.LE_spectrum(1), ...
                    numel(r.LE_spectrum), toc(t0));
            end
        end
        trial.qr = qrc;
        trials{k} = trial;
    end
    res(i).trials  = [trials{:}];
    res(i).summary = summarise(res(i).trials, fs_lle);
    print_summary(res(i).summary, titles{i});
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'model_class', model_class, 'n', preset.n, 'sigma_u_noise', sigma_preset, ...
    'fs_ladder', fs_ladder, 'fs_noisy_ref', fs_noisy_ref, 'ref_tol', ref_tol, ...
    'seg_long', seg_long, 'seg_short_steps', seg_short, ...
    'T_free', T_free, 'win_free', win_free, ...
    'T_noisy', T_noisy, 'win_noisy', win_noisy, ...
    'T_lle', T_lle, 'fs_lle', fs_lle, ...
    'n_small', n_small, 'T_small', T_small, 'lya_small_warmup', lya_small_warmup, ...
    'max_restarts', max_restarts, 'block_names', {{'x', 'a', 'b'}}, 'n_trials', n_trials, ...
    'minutes', toc(t_stage) / 60);

condition_titles = titles;  % saved name
results = res;
mat_file = fullfile(out_dir, 'numerics_verification_data.mat');
save(mat_file, 'results', 'cond_names', 'condition_titles', 'settings', '-v7.3');
fprintf('\nSaved: %s  (%.1f min)\n', mat_file, settings.minutes);
end

%% ------------------------------------------------------------------------
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

function S = summarise(trials, fs_lle)
% Per-trial scalars, one row per trial, for the ensemble figure and the report.
n = numel(trials);
S = struct('n_trials', n, 'seeds', vertcat(trials.seeds));
S.lle_ode45   = arrayfun(@(t) t.lle.ode45.LLE, trials);
S.lle_sra1    = arrayfun(@(t) t.lle.sra1.LLE, trials);
S.qr_benettin = arrayfun(@(t) t.qr.benettin.LLE, trials);
S.qr_lambda1  = arrayfun(@(t) t.qr.qr.LE_spectrum(1), trials);
j = find([trials(1).free.fs] == fs_lle, 1); if isempty(j); j = 1; end
S.fs_paper        = trials(1).free(j).fs;
S.err_free_paper  = arrayfun(@(t) t.free(j).err_long_total_rms, trials);
S.slope_free      = arrayfun(@(t) fit_slope([t.free.fs], [t.free.err_long_total_rms]), trials);
if ~isempty(trials(1).noisy)
    S.err_noisy_paper = arrayfun(@(t) t.noisy(j).err_long_total_rms, trials);
    S.slope_noisy     = arrayfun(@(t) fit_slope([t.noisy.fs], [t.noisy.err_long_total_rms]), trials);
else
    S.err_noisy_paper = nan(1, n);
    S.slope_noisy     = nan(1, n);
end
end

function p = fit_slope(fs, err)
c = polyfit(log(1 ./ fs), log(err), 1);
p = c(1);
end

function print_summary(S, title)
fprintf('\n  summary for %s over %d trial(s):\n', title, S.n_trials);
fprintf('    LLE ode45 : %s\n', mat2str(S.lle_ode45, 4));
fprintf('    LLE sra1  : %s\n', mat2str(S.lle_sra1, 4));
fprintf('    QR net    : Benettin %s | QR %s\n', mat2str(S.qr_benettin, 4), mat2str(S.qr_lambda1, 4));
fprintf('    reshoot @ %d Hz: free %s (slopes %s)\n', S.fs_paper, ...
    mat2str(S.err_free_paper, 3), mat2str(S.slope_free, 3));
if all(isfinite(S.err_noisy_paper))
    fprintf('                     noisy %s (slopes %s)\n', mat2str(S.err_noisy_paper, 3), mat2str(S.slope_noisy, 3));
end
end

function mat_file = run_transient_gain(cfg)
% RUN_TRANSIENT_GAIN Transient amplification, adaptation frozen vs active.
%
%   mat_file = RUN_TRANSIENT_GAIN('preset_name', p, 'run_mode', m, 'out_dir', d)
%
% The compute half of fig_transient_gain and fig_transient_gain_excursions.
% For every adaptation regime the preset states, on several network seeds,
% the x-in / x-out transient gain G(t) = ||P_x Phi(t) P_x'||_2 of the tangent
% flow at sampled states of the noisy trajectory (SRNNCellTypePairs.transient_gain),
% three ways per state:
%   frozen_x     Phi = expm(J_xx t)  -- the dendritic block fixed at the state:
%                the conventional rate-network picture, adaptation FROZEN
%   frozen_full  Phi = expm(J t)     -- the full Jacobian fixed at the state:
%                adaptation's feedback, no drift
%   active       Phi from J(S(t)) along the stored trajectory
% and, from the same n x n block, the worst-case gain, the noise-average
% gain, the gain along the E/I difference mode (balanced amplification), the
% E/I sum mode and the leading Lyapunov direction at that state, plus the
% optimal input direction at the peak with its E fraction, participation
% ratio and alignment with the leading direction.
%
% TWO SAMPLING RULES per trajectory, both inside [T/2 + 5, T - horizon]
% (5 s of history for the leading direction's warm-up, room for the horizon):
%   'regular'   n_regular states evenly spaced in time
%   'onset' / 'quiet'   up to n_excursion states each at the ONSETS of
%               local-Lyapunov-rate excursions (>= 4 positive 0.05-s segments
%               after >= 4 negative) and at the MIDPOINTS of quiet stretches
%               (>= 20 negative segments), from the top-K run's leading local
%               rate (SRNNCellTypePairs.excursion_samples). A regime with none
%               of a class reports 0 of it; nothing errors.
% The onset/quiet contrast asks whether the network's own intermittent
% divergence is non-normal amplification triggered by its fluctuations.
%
% Run modes (the network keeps the preset's n):
%   fast         T = 20 s, 1 seed,  6 regular, 4 onset + 4 quiet, horizon 1 s
%   medium(2)    T = 40 s, 2 seeds, 12 regular, 8 + 8,             horizon 2 s
%   production   T = 60 s, 3 seeds, 20 regular, 10 + 10,           horizon 3 s
% Each trial (condition x seed) runs top-K (K 15, cap 30, lya_dt 0.05) over
% [T/2, T] with a 5 s warm-up, with the preset's noise on. One parfor over the
% trials; each worker builds its own network and walks its samples serially,
% so no trajectory crosses a process boundary. Cost at n = 500, fs 400: the
% two full-state variants are ~2 x horizon x 400 jacobian_times calls on 500
% columns (~25 s per second of horizon each); frozen_x is a dense 500 x 500
% product, seconds.
%
% Builds happen INSIDE the parfor, so the generator is pinned to 'twister'
% first (workers default to 'threefry').
%
% See also: fig_transient_gain, fig_transient_gain_excursions,
%           SRNNCellTypePairs.transient_gain, run_lyapunov_spectrum (template)

arguments
    cfg.preset_name (1,:) char   = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25'
    cfg.run_mode    (1,:) char   = 'production'
    cfg.out_dir     (1,:) char   = ''
    cfg.verbose                    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (or a logical); see verbose_level
    cfg.n_seeds     (1,1) double = 0     % 0 -> per run_mode
    cfg.n_regular   (1,1) double = 0     % 0 -> per run_mode
    cfg.n_excursion (1,1) double = 0     % 0 -> per run_mode (onsets AND quiets, each)
    cfg.horizon_s   (1,1) double = 0     % 0 -> per run_mode
    cfg.n_override  (1,1) double = 0     % 0 -> the preset's n (tests use a small one)
    cfg.n_workers   (1,1) double = 0     % 0 -> min(12, cores)
    cfg.report_dt   (1,1) double = 0.02
end

setup_paths();

%% Cost table -- BEFORE any compute, so a bad mode fails in milliseconds
switch cfg.run_mode
    case 'fast'                 , T = 20; n_seeds = 1; n_regular = 6;  n_exc = 4;  horizon = 1;
    case {'medium', 'medium2'}  , T = 40; n_seeds = 2; n_regular = 12; n_exc = 8;  horizon = 2;
    case 'production'           , T = 60; n_seeds = 3; n_regular = 20; n_exc = 10; horizon = 3;
    otherwise
        error('run_transient_gain:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.n_seeds > 0;     n_seeds   = cfg.n_seeds;     end
if cfg.n_regular > 0;   n_regular = cfg.n_regular;   end
if cfg.n_excursion > 0; n_exc     = cfg.n_excursion; end
if cfg.horizon_s > 0;   horizon   = cfg.horizon_s;   end

P = struct();
P.preset_name = cfg.preset_name;
P.verbose     = cfg.verbose;
P.T           = T;
P.T_range     = [0, T];
P.lya_T_interval = [T / 2, T];
P.lya_warmup  = 5;
P.lya_dt      = 0.05;
P.K           = 15;
P.K_max       = 30;
P.fs          = 400;
P.n_override  = cfg.n_override;
P.horizon     = horizon;
P.report_dt   = cfg.report_dt;
P.n_regular   = n_regular;
P.n_exc       = n_exc;
P.min_pos     = 4;      % onset: >= 4 positive segments (0.2 s) after >= 4 negative
P.min_quiet   = 20;     % quiet: >= 20 negative segments (1 s)
P.dir_warmup  = 5;      % s of K = 1 propagation for the leading direction
P.variants    = {'frozen_x', 'frozen_full', 'active'};
P.t_lo        = T / 2 + P.dir_warmup;
P.t_hi        = T - horizon;

if isempty(cfg.out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'transient_gain');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir); mkdir(out_dir); end

[preset, model_class, conditions] = srnn_param_preset(cfg.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_transient_gain:WrongClass', ...
        'Preset ''%s'' is written for %s; this stage runs SRNNCellTypePairs only.', ...
        cfg.preset_name, model_class);
end
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
titles     = cellfun(@(n) pretty(n), cond_names, 'UniformOutput', false);
n_cond     = numel(cond_names);

vprintf(cfg.verbose, 'minimal', '[transient_gain] preset=%s run_mode=%s: %d conditions x %d seeds, T = %g s, %d regular + up to %d onset + %d quiet samples, horizon %g s, {%s}\n', ...
    cfg.preset_name, cfg.run_mode, n_cond, n_seeds, T, n_regular, n_exc, n_exc, horizon, strjoin(P.variants, ', '));

pool = ensure_pool(cfg.n_workers);
vprintf(cfg.verbose, 'verbose', '  parallel pool: %d workers\n', pool.NumWorkers);

%% One parfor over every (condition, seed) trial
[ci, si] = ndgrid(1:n_cond, 1:n_seeds);
trial_cond = ci(:); trial_seed = si(:);
n_trials = numel(trial_cond);
t_stage = tic;
R = cell(n_trials, 1);
parfor j = 1:n_trials
    R{j} = gain_trial(P, cond_names{trial_cond(j)}, [trial_seed(j), trial_seed(j) + 1]);
end
% The raw trials FIRST: 2026-09-13 a medium run lost 51 min of compute to a
% concatenation error in the assembly below. Whatever happens next, this is
% on disk.
if ~isfolder(out_dir); mkdir(out_dir); end
trials_file = fullfile(out_dir, 'transient_gain_trials.mat');
save(trials_file, 'R', 'trial_cond', 'trial_seed', 'cond_names', 'P', '-v7.3');

%% Assemble per condition
res = struct('name', cond_names, 'title', titles);
for i = 1:n_cond
    Ri = [R{trial_cond == i}];
    res(i).trials  = rmfield(Ri, 'samples');
    res(i).samples = vertcat(Ri.samples);   % columns of different lengths per seed: NOT [Ri.samples]
    res(i).n_onset_found = sum([Ri.n_onset_found]);
    res(i).n_quiet_found = sum([Ri.n_quiet_found]);
    res(i).N = Ri(1).N;
    res(i).n = Ri(1).n;
    vprintf(cfg.verbose, 'verbose', '  %-14s %d samples (%d regular, %d onset, %d quiet; found %d onsets, %d quiets) | %s\n', ...
        cond_names{i}, numel(res(i).samples), nnz(strcmp({res(i).samples.kind}, 'regular')), ...
        nnz(strcmp({res(i).samples.kind}, 'onset')), nnz(strcmp({res(i).samples.kind}, 'quiet')), ...
        res(i).n_onset_found, res(i).n_quiet_found, gmax_txt(res(i).samples, P.variants));
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'model_class', model_class, 'n', preset.n, 'n_override', cfg.n_override, ...
    'T', T, 'lya_T_interval', P.lya_T_interval, 'lya_warmup', P.lya_warmup, 'lya_dt', P.lya_dt, ...
    'K', P.K, 'K_max', P.K_max, 'fs', P.fs, 'n_seeds', n_seeds, 'n_regular', n_regular, ...
    'n_excursion', n_exc, 'horizon_s', horizon, 'report_dt', P.report_dt, ...
    'min_pos_segments', P.min_pos, 'min_quiet_segments', P.min_quiet, 'dir_warmup_s', P.dir_warmup, ...
    'sample_window', [P.t_lo, P.t_hi], 'variants', {P.variants}, ...
    'directions', {{'ei_diff', 'ei_sum', 'lyap'}}, ...
    'n_workers', pool.NumWorkers, 'minutes', toc(t_stage) / 60);
condition_titles = titles;   % saved name
results = res;               % saved name
mat_file = fullfile(out_dir, 'transient_gain_data.mat');
save(mat_file, 'results', 'cond_names', 'condition_titles', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '[transient_gain] %.1f min -> %s\n', settings.minutes, mat_file);
end

%% ------------------------------------------------------------------------
function r = gain_trial(P, cname, seeds)
% One network, one seed: build, run top-K, sample, propagate every variant.
rng(0, 'twister');                   % workers default to threefry; a seed must mean one network
args = {'rng_seeds', seeds, 'fs', P.fs, 'T_range', P.T_range, ...
    'lya_method', 'topk', 'lya_K', P.K, 'lya_K_auto', true, 'lya_K_max', P.K_max, ...
    'lya_dt', P.lya_dt, 'lya_T_interval', P.lya_T_interval, 'lya_warmup', P.lya_warmup, ...
    'store_full_state', true};
if P.n_override > 0
    args = [args, {'n', P.n_override, 'indegree', max(2, round(0.2 * P.n_override)), ...
        'F_tracks_network', true}];
end
t0 = tic;
m = build_from_preset(P.preset_name, cname, 'verbose', P.verbose, args{:});
m.run();
run_seconds = toc(t0);
S = m.lya_summary();
params = m.get_params();
S_out = m.S_out; t_out = m.t_out;
n = m.n; nE = params.n_per_type(1);
lr = S.local_rate_lead(:); tl = S.t_lya_lead(:);

% The sample times: regular, then onsets and quiets inside the window.
t_reg = linspace(P.t_lo, P.t_hi, P.n_regular);
E = SRNNCellTypePairs.excursion_samples(lr, tl, P.min_pos, P.min_quiet);
t_on = E.onset_t(E.onset_t >= P.t_lo & E.onset_t <= P.t_hi);
t_qu = E.quiet_t(E.quiet_t >= P.t_lo & E.quiet_t <= P.t_hi);
n_onset_found = numel(t_on); n_quiet_found = numel(t_qu);
t_on = pick_evenly(t_on, P.n_exc);
t_qu = pick_evenly(t_qu, P.n_exc);
t_all = [t_reg(:); t_on(:); t_qu(:)];
kinds = [repmat({'regular'}, numel(t_reg), 1); repmat({'onset'}, numel(t_on), 1); repmat({'quiet'}, numel(t_qu), 1)];

dirs0 = struct('ei_diff', [ones(nE, 1) / sqrt(nE); -ones(n - nE, 1) / sqrt(n - nE)], ...
    'ei_sum', ones(n, 1) / sqrt(n));
n_var = numel(P.variants);
samples = repmat(empty_sample(), numel(t_all), 1);
for s = 1:numel(t_all)
    ts = tic;
    t_s = t_all(s);
    i0 = find(t_out >= t_s - 1e-9, 1);
    k_l = find(tl <= t_s + 1e-9, 1, 'last');
    if isempty(k_l); lr_s = NaN; else; lr_s = lr(k_l); end
    [v_lyap, v_full] = SRNNCellTypePairs.leading_direction_at(S_out, t_out, i0, params, P.dir_warmup);
    B = SRNNCellTypePairs.leading_vector_blocks(v_full, params.state_layout);
    dirs = dirs0; dirs.lyap = v_lyap;
    smp = empty_sample();
    smp.t_sample = t_s; smp.kind = kinds{s}; smp.seeds = seeds; smp.local_rate_at_sample = lr_s;
    smp.lyap_block_fracs = [B.x, B.sfa, B.std, B.stf];
    % The frozen operating point's own rates, beside lambda_1 of the trajectory:
    % alpha_xx (the eig stage's quantity), omega_xx, and alpha_full, the
    % spectral abscissa of the FULL J at this state -- "the trajectory is more
    % stable than any of its states" needs this number. eigs on the sparse
    % nonsymmetric N x N can fail to converge; then NaN, never an error.
    J = SRNNCellTypePairs.compute_Jacobian_fast(S_out(i0, :)', params);
    Jxx = full(J(params.state_layout.x, params.state_layout.x));
    smp.alpha_xx = max(real(eig(Jxx)));
    smp.omega_xx = max(eig((Jxx + Jxx') / 2));
    try
        ws = warning('off', 'MATLAB:eigs:NotAllEigsConverged');
        ev = eigs(J, 6, 'largestreal', 'MaxIterations', 600, 'Display', false);
        warning(ws);
        smp.alpha_full = max(real(ev));       % max omits the NaN of unconverged ones
    catch
        smp.alpha_full = NaN;
    end
    for v = 1:n_var
        opts = struct('variant', P.variants{v}, 'report_dt', P.report_dt, 'directions', dirs);
        if strcmp(P.variants{v}, 'frozen_x'); opts.Jxx = Jxx; end
        [G, info] = SRNNCellTypePairs.transient_gain(S_out, t_out, i0, params, P.horizon, opts);
        if v == 1
            nt = numel(G.t);
            smp.t = G.t(:)';
            for f = {'G_worst', 'G_noise', 'G_ei_diff', 'G_ei_sum', 'G_lyap'}
                smp.(f{1}) = nan(n_var, nt);
            end
        end
        smp.G_worst(v, :)   = G.worst(:)';
        smp.G_noise(v, :)   = G.noise(:)';
        smp.G_ei_diff(v, :) = G.dir.ei_diff(:)';
        smp.G_ei_sum(v, :)  = G.dir.ei_sum(:)';
        smp.G_lyap(v, :)    = G.dir.lyap(:)';
        smp.G_max(v) = info.G_max;
        smp.t_peak(v) = info.t_peak;
        smp.align_opt_lyap(v) = info.cos_opt.lyap;
        smp.cos_opt_ei_diff(v) = info.cos_opt.ei_diff;
        smp.cos_opt_ei_sum(v) = info.cos_opt.ei_sum;
        smp.frac_E_opt(v) = info.frac_E;
        smp.participation_opt(v) = info.participation;
        switch P.variants{v}
            case 'active';   smp.v_opt_active = info.v_opt;
            case 'frozen_x'; smp.v_opt_frozen_x = info.v_opt;
        end
    end
    smp.seconds = toc(ts);
    samples(s) = smp;
end

r = struct('seeds', seeds, 'N', m.N_sys_eqs, 'n', n, 'LLE', S.LLE, 'K_used', S.K_used, ...
    'frac_local_positive', S.frac_local_positive, 'mean_positive_excursion_s', S.mean_positive_excursion_s, ...
    'n_onset_found', n_onset_found, 'n_quiet_found', n_quiet_found, ...
    'n_onset_used', numel(t_on), 'n_quiet_used', numel(t_qu), ...
    'sigma_u_noise', m.sigma_u_noise, 'ode_solver', m.ode_solver, ...
    'run_seconds', run_seconds, 'seconds', toc(t0), 'samples', samples);
vprintf(P.verbose, 'verbose', '  %-14s seed %d: N %d, lambda_1 %+.3f, %d samples (%d onsets, %d quiets found) in %.0f s\n', ...
    cname, seeds(1), r.N, r.LLE, numel(samples), n_onset_found, n_quiet_found, r.seconds);
end

function s = empty_sample()
s = struct('t_sample', NaN, 'kind', '', 'seeds', [NaN NaN], 'local_rate_at_sample', NaN, ...
    't', [], 'G_worst', [], 'G_noise', [], 'G_ei_diff', [], 'G_ei_sum', [], 'G_lyap', [], ...
    'G_max', nan(1, 3), 't_peak', nan(1, 3), 'align_opt_lyap', nan(1, 3), ...
    'cos_opt_ei_diff', nan(1, 3), 'cos_opt_ei_sum', nan(1, 3), ...
    'frac_E_opt', nan(1, 3), 'participation_opt', nan(1, 3), ...
    'v_opt_active', [], 'v_opt_frozen_x', [], 'lyap_block_fracs', nan(1, 4), 'alpha_xx', NaN, 'omega_xx', NaN, 'alpha_full', NaN, 'seconds', NaN);
end

function t = pick_evenly(t, k)
% Up to k entries spread evenly over the candidates (all of them if fewer).
t = t(:);
if numel(t) > k
    t = t(unique(round(linspace(1, numel(t), k))));
end
end

function s = gmax_txt(samples, variants)
reg = samples(strcmp({samples.kind}, 'regular'));
parts = cell(1, numel(variants));
for v = 1:numel(variants)
    g = arrayfun(@(x) x.G_max(v), reg);
    parts{v} = sprintf('%s %.2f', variants{v}, median(g));
end
s = ['median G_max (regular): ' strjoin(parts, ', ')];
end

function pool = ensure_pool(n_workers)
% Use the pool that is up; else try to start one; else run serially (warn).
pool = gcp('nocreate');
if ~isempty(pool); return; end
if n_workers <= 0; n_workers = min(12, feature('numcores')); end
try
    pool = parpool(parallel.defaultProfile, n_workers);
catch ME
    warning('run_transient_gain:NoPool', ...
        ['Could not start a parallel pool (%s). Running the trials SERIALLY. ' ...
         'Secure a pool first with wait_for_parpool.'], strtok(ME.message, newline));
    pool = struct('NumWorkers', 1);
end
end

function s = pretty(name)
st = manuscript_style();
if st.condition_title.isKey(name); s = st.condition_title(name); else; s = strrep(name, '_', ' '); end
end

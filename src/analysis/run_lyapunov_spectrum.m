function mat_file = run_lyapunov_spectrum(cfg)
% RUN_LYAPUNOV_SPECTRUM The top-K Lyapunov spectrum of the paper's network.
%
%   mat_file = RUN_LYAPUNOV_SPECTRUM('preset_name', p, 'run_mode', m, 'out_dir', d)
%
% The compute half of fig_lyapunov_spectrum. For every adaptation regime the
% preset states, on several network seeds, the K largest Lyapunov exponents
% by the discrete QR method (lya_method 'topk', matrix-free Jacobian product,
% retry doubling K while the Kaplan-Yorke dimension is unresolved), TWICE per
% seed: with the preset's noise (the paper's physics) and with noise off (the
% deterministic attractor, for the Engelken-style comparison -- fluctuating
% input lowers both the entropy rate and the attractor dimension). Per run it
% stores the spectrum, n_positive, h_KS, D_KY with its resolved flag, the
% conditioning and convergence diagnostics, the leading vector's block
% fractions and the first ten finite-time exponents' convergence curves.
%
% This is what the sweeps' per-job scalars (K = 15, retry to 30) cannot give:
% the SHAPE of the spectrum, a resolved D_KY in the strongly chaotic regime,
% and a longer accumulation window for the mid-spectrum exponents.
%
% Run modes (the network keeps the preset's n: shrinking it would change the
% thing being measured, cf. run_eig_heatmap):
%   fast         T = 20 s, 1 seed,  K = 50   (smoke: minutes)
%   medium(2)    T = 40 s, 3 seeds, K = 200
%   production   T = 100 s, 5 seeds, K = 300
% Accumulation over the second half of T after a warm-up of T/4;
% re-orthonormalisation every 0.05 s; lya_K_max = 2K.
%
% Builds happen INSIDE the parfor, so the generator is pinned to 'twister'
% first (workers default to 'threefry'; UserNotes.md) -- a seed then means
% the same network as on the client.
%
% See also: fig_lyapunov_spectrum, lyapunov_topk, SRNNCellTypePairs.lya_summary,
%           run_numerics_verification (the stage this one is modelled on)

arguments
    cfg.preset_name (1,:) char   = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25'
    cfg.run_mode    (1,:) char   = 'production'
    cfg.out_dir     (1,:) char   = ''
    cfg.verbose                    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (or a logical); see verbose_level
    cfg.n_seeds     (1,1) double = 0     % 0 -> per run_mode
    cfg.K           (1,1) double = 0     % 0 -> per run_mode
    cfg.n_override  (1,1) double = 0     % 0 -> the preset's n (tests use a small one)
    cfg.noise_off_too (1,1) logical = true
    cfg.n_workers   (1,1) double = 0     % 0 -> min(12, cores)
end

setup_paths();

%% Cost table -- BEFORE any compute, so a bad mode fails in milliseconds
switch cfg.run_mode
    case 'fast'                 , T = 20;  n_seeds = 1; K = 50;
    case {'medium', 'medium2'}  , T = 40;  n_seeds = 3; K = 200;
    case 'production'           , T = 100; n_seeds = 5; K = 300;
    otherwise
        error('run_lyapunov_spectrum:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.n_seeds > 0; n_seeds = cfg.n_seeds; end
if cfg.K > 0;       K = cfg.K; end

P = struct();
P.preset_name = cfg.preset_name;
P.verbose     = cfg.verbose;
P.T_range     = [0, T];
P.lya_T_interval = [T / 2, T];
P.lya_warmup  = T / 4;
P.lya_dt      = 0.05;
P.K           = K;
P.K_max       = 2 * K;
P.fs          = 400;
P.n_override  = cfg.n_override;
P.n_finite    = 10;                  % finite-time curves kept for the first 10 exponents

if isempty(cfg.out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'lyapunov_spectrum');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir); mkdir(out_dir); end

[preset, model_class, conditions] = srnn_param_preset(cfg.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_lyapunov_spectrum:WrongClass', ...
        'Preset ''%s'' is written for %s; this stage runs SRNNCellTypePairs only.', ...
        cfg.preset_name, model_class);
end
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
titles     = cellfun(@(n) pretty(n), cond_names, 'UniformOutput', false);
n_cond     = numel(cond_names);
sigma_preset = 0;
if isfield(preset, 'sigma_u_noise'); sigma_preset = preset.sigma_u_noise; end
variants = {'noise_on'};
if cfg.noise_off_too && sigma_preset > 0; variants{end + 1} = 'noise_off'; end
if sigma_preset == 0; variants = {'noise_off'}; end

vprintf(cfg.verbose, 'minimal', '[lyapunov_spectrum] preset=%s run_mode=%s: %d conditions x %d seeds x {%s}, T = %g s, K = %d (cap %d)\n', ...
    cfg.preset_name, cfg.run_mode, n_cond, n_seeds, strjoin(variants, ', '), T, K, P.K_max);

pool = ensure_pool(cfg.n_workers);
vprintf(cfg.verbose, 'verbose', '  parallel pool: %d workers\n', pool.NumWorkers);

t_stage = tic;
res = struct('name', cond_names, 'title', titles);
for i = 1:n_cond
    cname = cond_names{i};
    for v = 1:numel(variants)
        variant = variants{v};
        t0 = tic;
        R = cell(1, n_seeds);
        parfor k = 1:n_seeds
            R{k} = spectrum_trial(P, cname, [k, k + 1], strcmp(variant, 'noise_on'));
        end
        res(i).(variant) = [R{:}];
        r1 = res(i).(variant)(1);
        vprintf(cfg.verbose, 'verbose', '  %-14s %-9s %d seeds in %5.0f s | N = %d, K_used %d, lambda_1 %+.4f, n_pos %d, h_KS %.2f bit/s, D_KY %s\n', ...
            cname, variant, n_seeds, toc(t0), r1.N, r1.K_used, r1.LLE, r1.n_positive, r1.h_KS_bits, dky_txt(r1));
    end
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'model_class', model_class, 'n', preset.n, 'n_override', cfg.n_override, ...
    'T', T, 'lya_T_interval', P.lya_T_interval, 'lya_warmup', P.lya_warmup, ...
    'lya_dt', P.lya_dt, 'K', K, 'K_max', P.K_max, 'fs', P.fs, 'n_seeds', n_seeds, ...
    'sigma_u_noise', sigma_preset, 'variants', {variants}, ...
    'n_workers', pool.NumWorkers, 'minutes', toc(t_stage) / 60);
condition_titles = titles;   % saved name
results = res;               % saved name
mat_file = fullfile(out_dir, 'lyapunov_spectrum_data.mat');
save(mat_file, 'results', 'cond_names', 'condition_titles', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '[lyapunov_spectrum] %.1f min -> %s\n', settings.minutes, mat_file);
end

%% ------------------------------------------------------------------------
function r = spectrum_trial(P, cname, seeds, noise_on)
% One network, one variant: build, run, keep the spectrum and its quality.
rng(0, 'twister');                   % workers default to threefry; a seed must mean one network
args = {'rng_seeds', seeds, 'fs', P.fs, 'T_range', P.T_range, ...
    'lya_method', 'topk', 'lya_K', P.K, 'lya_K_auto', true, 'lya_K_max', P.K_max, ...
    'lya_dt', P.lya_dt, 'lya_T_interval', P.lya_T_interval, 'lya_warmup', P.lya_warmup, ...
    'store_full_state', true};
if ~noise_on
    args = [args, {'sigma_u_noise', 0, 'ode_solver', 'rk4'}];
end
if P.n_override > 0
    args = [args, {'n', P.n_override, 'indegree', max(2, round(0.2 * P.n_override)), ...
        'F_tracks_network', true}];
end
t0 = tic;
m = build_from_preset(P.preset_name, cname, 'verbose', P.verbose, args{:});
m.run();
S = m.lya_summary();
lr = m.lya_results;
r = S;
r = rmfield(r, {'local_rate_lead', 't_lya_lead'});
r.LE_spectrum = lr.LE_spectrum(:)';
r.finite_LE_spectrum_t = lr.finite_LE_spectrum_t(:, 1:min(P.n_finite, lr.K));
r.t_lya = lr.t_lya(:);
r.local_rate_lead = S.local_rate_lead(:)';
r.N = m.N_sys_eqs;
r.n = m.n;
r.seeds = seeds;
r.noise_on = noise_on;
r.sigma_u_noise = m.sigma_u_noise;
r.ode_solver = m.ode_solver;
r.seconds = toc(t0);
r.lya_seconds = lr.seconds;
end

function s = dky_txt(r)
if r.D_KY_resolved; s = sprintf('%.1f', r.D_KY); else; s = sprintf('> %d (unresolved)', r.K_used); end
end

function pool = ensure_pool(n_workers)
% Use the pool that is up; else try to start one; else run serially (warn).
pool = gcp('nocreate');
if ~isempty(pool); return; end
if n_workers <= 0; n_workers = min(12, feature('numcores')); end
try
    pool = parpool(parallel.defaultProfile, n_workers);
catch ME
    warning('run_lyapunov_spectrum:NoPool', ...
        ['Could not start a parallel pool (%s). Running the seed loops SERIALLY. ' ...
         'Secure a pool first with wait_for_parpool.'], strtok(ME.message, newline));
    pool = struct('NumWorkers', 1);
end
end

function s = pretty(name)
st = manuscript_style();
if st.condition_title.isKey(name); s = st.condition_title(name); else; s = strrep(name, '_', ' '); end
end

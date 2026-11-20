% EXPLORE_SMALL_NETWORK The 80-neuron mu7chain network under three adaptation conditions, 15 s each, top-20 Lyapunov exponents.
%
% WHY. At 500 neurons the leading few dozen Lyapunov exponents of a stable
% adapting network all sit in a band of slow single-neuron adaptation
% variables, so a top-K estimate never gets past it. An 80-neuron version
% of the same network is small enough for later work to take the full
% spectrum and the complex (Ritz) spectrum of the Jacobian projected on the
% leading Lyapunov directions. This script is the first step only: build
% the 80-neuron network (small_network_model, the one place it is
% defined), run three conditions for 15 s each, print diagnostics and plot
% each with the class's own plot method. Nothing here is tuned to produce
% an effect.
%
% THE NETWORK: the mu7chain bundle's network with only its size changed,
% n = 80, indegree = 16 (alpha = 0.2), F_ref_n = 80, F_ref_indegree = 16,
% so F = 0.18634 and the theoretical bulk radius R = 4.799 (as at 500
% neurons). See small_network_model for the rest.
%
% THE CONDITIONS: no_adaptation, sfa1_std1 (one SFA timescale, one
% depression pool), sfa3_stdc2 (three SFA timescales, chained vesicle pool).
%
% THE RUN: 15 s at 500 Hz, sra1, network seed [7 8]; top-20 Lyapunov
% exponents (topk, lya_dt 0.05 s) accumulated over 5 to 15 s after a 5 s
% warm-up. The conditions run in a parfor over a pool of 3 opened with the
% project's pool module (pool_set / pool_release); when pool_set fails the
% error is printed and the parfor runs serially.
%
% WHAT IT PRINTS: the network (n, alpha, F, R; spectral radius, spectral
% abscissa and largest singular value of the realized W; the norm of its
% block-mean part); per condition the state dimension, the 20 exponents at
% the end of the window, how many are positive, the Kaplan-Yorke dimension,
% the lambda_1 drift (end of window minus three quarters of it, a
% convergence check), the mean rate of each cell type over 5 to 15 s, the
% share of neuron-time at the floor of the nonlinearity (r = 0), at its
% ceiling (r = 1) and between, the number of neurons silent over the whole
% window, and wall time.
%
% THE RATES are r, the output of the nonlinearity BEFORE depression, read
% from the model's stored decimated output (plot_data) at fs / plot_deci =
% 100 Hz. The piecewise nonlinearity returns exactly 0 below its lower knee
% and exactly 1 above its upper knee, so the floor and ceiling tests are
% exact comparisons.
%
% FIGURES: one per condition, the class's m.plot() (input, x, r, synaptic
% output, SFA, depression factor, and the 20 local exponents in the
% Lyapunov panel), the window named with the condition. Saved as PNG to
% figs/explorations/small_network/ (gitignored) and left open.
%
% PATH: run setup_paths once in the session first, with the current folder
% at the repository root. setup_paths adds scripts/ with genpath, so this
% folder is on the path provided setup_paths ran after the folder existed;
% otherwise run setup_paths again (or addpath this folder).
%
% Takes well under a minute on a laptop.
%
% WHAT TO EXPECT with the settings below, on the code as of 2026-11-19
% (two runs from a cleared workspace gave these to the bit):
%   W:         spectral radius 5.3691, spectral abscissa 4.3948, largest
%              singular value 22.196, block-mean part 21.09, R = 4.7987
%   lambda_1:  no_adaptation +0.7939 (1 positive of 20, NOT settled: it
%              moved by -0.34 over the last quarter of the window),
%              sfa1_std1 -0.4635, sfa3_stdc2 -0.1185 (0 positive in both)
%   top 20:    sfa1_std1 all within -0.46 to -0.58; sfa3_stdc2 19 within
%              -0.12 to -0.15: a band of slow adaptation variables, so
%              top-20 does not get past it at 80 neurons either
%   rates:     no_adaptation is a slow collective oscillation (about 1.3 s
%              a cycle) with half the neuron-time at the floor and a third
%              at the ceiling, unlike the 500-neuron network; with
%              adaptation almost all neuron-time is between the two
% A different lambda_1 means the model, the preset or this script changed.
%
% See also: small_network_model, build_from_preset, pool_set, pool_release

%% Settings
conditions   = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
rng_seeds    = [7 8];          % network seed pair
T_range      = [0 15];         % s
fs           = 500;            % Hz
lya_K        = 20;             % top-K Lyapunov exponents
lya_dt       = 0.05;           % s, re-orthonormalization interval
lya_warmup   = 5;              % s
lya_T_interval = [5 15];       % s, accumulation window
stats_window = [5 15];         % s, window for rates and floor / ceiling shares
plot_deci    = 5;              % stored output at fs / plot_deci = 100 Hz
n_pool       = 3;              % pool workers (one per condition)
save_figs    = true;
out_dir      = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', 'small_network');

%% Pool (the project's pool module; never parpool directly)
max_workers = 0;               % 0: the parfor runs serially
try
    pool_set(n_pool, 'label', 'small_network');
    p = gcp('nocreate');
    if isempty(p)
        fprintf('[pool] pool_set returned but no pool is up; running serially\n');
    else
        fprintf('[pool] pool up with %d workers\n', p.NumWorkers);
        max_workers = p.NumWorkers;
    end
catch err
    fprintf('[pool] pool_set FAILED: %s\n  (%s)\n  running the conditions serially\n', ...
        err.message, err.identifier);
end

%% Run the conditions
n_cond = numel(conditions);
models = cell(1, n_cond);
wall_s = zeros(1, n_cond);
parfor (k = 1:n_cond, max_workers)
    t_start = tic;
    mk = small_network_model(conditions{k}, 'seed', rng_seeds, ...
        'T_range', T_range, 'fs', fs, 'lya_method', 'topk', ...
        'lya_K', lya_K, 'lya_dt', lya_dt, ...
        'lya_T_interval', lya_T_interval, 'lya_warmup', lya_warmup, ...
        'model_args', {'plot_deci', plot_deci});
    mk.run();
    wall_s(k) = toc(t_start);
    models{k} = mk;
end

%% Release the pool and confirm it is gone
pool_release();
if isempty(gcp('nocreate'))
    fprintf('[pool] released; no pool is up\n');
else
    fprintf('[pool] WARNING: a pool is still up after pool_release\n');
end

%% Diagnostics: the network
m1 = models{1};
W = m1.W;
same_W = all(cellfun(@(mk) isequal(mk.W, W), models));
ev = eig(full(W));
F = m1.default_val;
% The block-mean part of W: each (post type, pre type) block replaced by
% its mean (zeros of the mask included), the realized counterpart of the
% rank-one E-minus-I mean.
W_mean = zeros(size(W));
types = m1.cell_type_names;
for q = 1:numel(types)
    for s = 1:numel(types)
        iq = m1.type_indices{q};
        is = m1.type_indices{s};
        W_mean(iq, is) = mean(full(W(iq, is)), 'all');
    end
end
fprintf('\n=== The network (seed [%d %d]) ===\n', rng_seeds);
fprintf('  n = %d, indegree = %d, alpha = %.3f, F = %.5f, theoretical R = %.4f\n', ...
    m1.n, m1.indegree, m1.alpha, F, m1.R);
fprintf('  realized W: spectral radius %.4f, spectral abscissa %.4f, largest singular value %.4f\n', ...
    max(abs(ev)), max(real(ev)), norm(full(W)));
fprintf('  block-mean part of W: norm %.3f (expected about N alpha F 7 = %.3f)\n', ...
    norm(W_mean), m1.n * m1.alpha * F * 7);
fprintf('  connections: %d (expected n^2 alpha = %.0f)\n', nnz(W), m1.n^2 * m1.alpha);
fprintf('  W identical across the three conditions: %d\n', same_W);

%% Diagnostics: per condition
rows = struct([]);
for k = 1:n_cond
    mk = models{k};
    lr = mk.lya_results;
    S = mk.lya_summary();
    pd = mk.plot_data;
    in_win = pd.t >= stats_window(1) & pd.t <= stats_window(2);
    r_all = [];
    rate_txt = '';
    for q = 1:numel(types)
        rq = pd.r.(types{q})(:, in_win);
        rate_txt = [rate_txt, sprintf('%s %.4f  ', types{q}, mean(rq, 'all'))]; %#ok<AGROW>
        r_all = [r_all; rq]; %#ok<AGROW>
    end
    share_floor   = mean(r_all == 0, 'all');
    share_ceiling = mean(r_all == 1, 'all');
    share_between = 1 - share_floor - share_ceiling;
    n_silent      = sum(all(r_all == 0, 2));
    n_saturated   = sum(all(r_all == 1, 2));

    fprintf('\n=== %s ===\n', conditions{k});
    fprintf('  state dimension %d, wall time %.1f s\n', mk.N_sys_eqs, wall_s(k));
    fprintf('  lambda_1 = %+.6f, positive exponents %d of %d, D_KY = %.3f (resolved %d), lambda_1 drift %+.4f\n', ...
        lr.LLE, lr.n_positive, lr.K, lr.D_KY, lr.D_KY_resolved, S.lambda_1_drift);
    fprintf('  exponents (1/s):');
    fprintf(' %+.4f', lr.LE_spectrum(:));
    fprintf('\n');
    fprintf('  mean rate %d-%d s: %s\n', stats_window, rate_txt);
    fprintf('  neuron-time at floor (r = 0) %.3f, at ceiling (r = 1) %.3f, between %.3f\n', ...
        share_floor, share_ceiling, share_between);
    fprintf('  neurons silent over the whole window: %d of %d; saturated over it: %d\n', ...
        n_silent, mk.n, n_saturated);

    rows(k).condition = conditions{k};
    rows(k).lambda_1 = lr.LLE;
    rows(k).spectrum = lr.LE_spectrum(:)';
    rows(k).n_positive = lr.n_positive;
    rows(k).D_KY = lr.D_KY;
    rows(k).share_floor = share_floor;
    rows(k).share_ceiling = share_ceiling;
    rows(k).n_silent = n_silent;
end

%% Figures: the class's own plot, one per condition, left open
if save_figs && ~isfolder(out_dir)
    mkdir(out_dir);
end
for k = 1:n_cond
    fig = models{k}.plot();
    set(fig, 'Name', sprintf('Small network (n = 80): %s', conditions{k}), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 + 40 * k, 60, 900, 1000]);
    tl = findobj(fig, 'Type', 'tiledlayout');
    title(tl, sprintf('n = 80, %s, seed [%d %d], \\lambda_1 = %+.3f', ...
        strrep(conditions{k}, '_', '\_'), rng_seeds, models{k}.lya_results.LLE));
    if save_figs
        png = fullfile(out_dir, sprintf('small_network_%s.png', conditions{k}));
        exportgraphics(fig, png, 'Resolution', 150);
        fprintf('saved %s\n', png);
    end
end

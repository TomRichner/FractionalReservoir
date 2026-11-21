function [E, file] = explore_eigenvalue_density(opts)
%EXPLORE_EIGENVALUE_DENSITY Does adaptation flatten the eigenvalue density of the frozen-gain connectivity? Saved runs only.
%
%   E = explore_eigenvalue_density()
%   [E, file] = explore_eigenvalue_density('runs_dir', d, 'out_dir', d2)
%
% THE QUESTION (Tom, 2026-11-20). Harris et al. (2023) show that the bulk
% of a sparse E-I matrix is denser at the center than at the rim, which
% leaves fewer eigenvalues near the stability boundary. Does adaptation
% change that, against no adaptation, and one timescale against several?
%
% THE OBJECT. M = W diag(d), the frozen-gain effective connectivity: each
% neuron's outgoing weights times its gain d_j (slope of the nonlinearity
% times the available vesicle fraction). Its bulk is that of a random
% matrix whose column j has variance proportional to d_j^2, so the DENSITY
% inside the bulk, not only its radius, follows the distribution of the
% gains (edens_radial_theory). In the paper's network W alone has equal
% E and I variances, so its own density is uniform (Harris's g = 0): all
% the non-uniformity of M comes from the gains.
%
% WHAT THIS DOES, for every saved run of the outliers exploration
% (data/explorations/outliers/run_*.mat; nothing is simulated):
%   1. recovers the per-neuron gains at the run's 12 calibration states
%      (edens_run_gains) and at every gain_step-th sampled state;
%   2. takes the stored spectrum of M at those 12 states, the stored null
%      draws (the same gains on a shuffled, mean-free W: the random-matrix
%      ensemble for these gains, sampled) and the stored spectrum under a
%      uniform gain, and forms the cumulative radial distribution and the
%      radial density, in units of the predicted bulk radius R_eff;
%   3. evaluates the theory with the measured gains
%      (edens_radial_theory_blocks) at the same states;
%   4. takes the measures of edens_measures on all of them, and the band
%      counts on the full Jacobian's spectrum (the dense pass);
%   5. does the band counts for the switching run of the paper's Figure 3
%      (full-Jacobian spectra at 40 states per condition), when its file
%      is there.
%
% Output E: struct with fields
%   runs      struct array, one per run (see the assembly below). The
%             checks of the gain recovery are kept with each run:
%             check_ev, check_zeros and check_moments_cal at the 12
%             calibration states (everything in the figures and the table
%             comes from those), check_moments_t and frac_t_failed for the
%             states along the window (used only for p0_t, cv_nz_t, pr_t)
%   switching struct for the switching run ([] when the file is absent)
%   harris    the density of W alone (edens_harris_density)
%   settings  edens_measures('settings'), the grids
% saved as eigenvalue_density.mat in out_dir (default
% data/explorations/eigenvalue_density).
%
% See also: fig_eigenvalue_density, edens_check_theory

arguments
    opts.runs_dir (1,:) char = ''
    opts.out_dir  (1,:) char = ''
    opts.switching_file (1,:) char = ''
    opts.gain_step (1,1) double = 20       % every 20th sampled state: 0.2 s apart
    opts.verbose  (1,1) logical = true
end
root = fileparts(which('setup_paths'));
if isempty(opts.runs_dir); opts.runs_dir = fullfile(root, 'data', 'explorations', 'outliers'); end
if isempty(opts.out_dir);  opts.out_dir  = fullfile(root, 'data', 'explorations', 'eigenvalue_density'); end
if isempty(opts.switching_file)
    opts.switching_file = fullfile(root, 'data', 'fd034_mu7chain_essential_switchK100_t0', ...
        'local_lyapunov', 'local_lyapunov_spectra.mat');
end
if ~isfolder(opts.out_dir); mkdir(opts.out_dir); end

S = edens_measures('settings');
u = 0:0.01:1.4;                          % radius in units of R_eff, for the cumulative curves
u2_edges = 0:0.04:1.44;                  % bins of (r / R_eff)^2: equal areas, so a uniform disk is flat
u_theory = sqrt(linspace(0, 1, 101));    % where the theory is solved (uniform in area)

files = dir(fullfile(opts.runs_dir, 'run_*.mat'));
runs = cell(1, numel(files));
for f = 1:numel(files)
    L = load(fullfile(files(f).folder, files(f).name), 'R');
    R = L.R;
    G = edens_run_gains(R);
    n = R.n; n_cal = numel(R.cal);
    C = numel(G.ti);

    cdf_M = zeros(n_cal, numel(u)); cdf_null = cdf_M; cdf_unigain = cdf_M; cdf_theory = cdf_M;
    dens_M = zeros(n_cal, numel(u2_edges) - 1); dens_null = dens_M; dens_theory = dens_M;
    R_eff = zeros(1, n_cal); p0 = R_eff; cv_nz = R_eff; pr = R_eff; mean_gain = R_eff; rms_gain = R_eff;
    names = {'shell', 'r50', 'n_band', 'n_unstable'};
    mM = zeros(4, n_cal); mNull = mM; mTheory = mM; mUniform = mM; mUnigain = mM;
    for k = 1:n_cal
        d = G.d(:, k);
        Tk = edens_radial_theory_blocks(G.s2, G.type_of, d, []);
        Tk = edens_radial_theory_blocks(G.s2, G.type_of, d, u_theory * Tk.R);
        Re = Tk.R;
        R_eff(k) = Re;
        g = gain_stats(d);
        p0(k) = g.p0; cv_nz(k) = g.cv_nz; pr(k) = g.pr; mean_gain(k) = g.mean; rms_gain(k) = g.rms;

        ev_M = R.cal(k).ev_M; ev_null = R.cal(k).ev_null; ev_u = R.cal(k).ev_uniform;
        cdf_M(k, :) = radial_cdf(ev_M, Re, u);
        cdf_unigain(k, :) = radial_cdf(ev_u, Re, u);
        cdf_theory(k, :) = interp1([u_theory, 1.5], [Tk.n_less, 1], u);
        dens_M(k, :) = radial_density(ev_M, Re, u2_edges);
        th = interp1([u_theory.^2, 2], [Tk.n_less, 1], u2_edges);
        th(1) = Tk.p0;                                    % leave the atom out of the first bin
        dens_theory(k, :) = diff(th) ./ diff(u2_edges);
        nd = size(ev_null, 2);
        m_null = zeros(4, nd);
        for q = 1:nd
            cdf_null(k, :) = cdf_null(k, :) + radial_cdf(ev_null(:, q), Re, u) / nd;
            dens_null(k, :) = dens_null(k, :) + radial_density(ev_null(:, q), Re, u2_edges) / nd;
            m_null(:, q) = as_column(edens_measures('spectrum', ev_null(:, q), Re), names);
        end
        mM(:, k) = as_column(edens_measures('spectrum', ev_M, Re), names);
        mNull(:, k) = mean(m_null, 2);
        mTheory(:, k) = as_column(edens_measures('radial', Tk.r, Tk.n_less, n), names);
        mUniform(:, k) = as_column(edens_measures('uniform', Re, n), names);
        mUnigain(:, k) = as_column(edens_measures('spectrum', ev_u, Re), names);
    end

    % The full Jacobian at the same states (the dense pass): all its
    % eigenvalues, and those carried mostly by the dendritic variables.
    full = struct('n_band', nan(1, n_cal), 'n_unstable', nan(1, n_cal), ...
        'n_band_x', nan(1, n_cal), 'n_unstable_x', nan(1, n_cal));
    if ~isempty(R.dense)
        for k = 1:n_cal
            lam = double(R.dense(k).lam); xf = double(R.dense(k).xfrac);
            a = edens_measures('rates', lam);
            b = edens_measures('rates', lam(xf >= 0.5));
            full.n_band(k) = a.n_band; full.n_unstable(k) = a.n_unstable;
            full.n_band_x(k) = b.n_band; full.n_unstable_x(k) = b.n_unstable;
        end
    end

    % The gains along the whole window, for how their distribution moves.
    idx_t = 1:opts.gain_step:numel(R.t);
    Gt = edens_run_gains(R, idx_t);
    gt = arrayfun(@(k) gain_stats(Gt.d(:, k)), 1:numel(idx_t));

    by_type = zeros(C, 3);                    % p0, mean of the nonzero gains, cv of the nonzero gains
    for b = 1:C
        gb = gain_stats(reshape(G.d(G.ti{b}, :), [], 1));
        by_type(b, :) = [gb.p0, gb.mean_nz, gb.cv_nz];
    end

    k_ex = ceil(n_cal / 2);
    runs{f} = struct('file', files(f).name, 'network', R.network, 'condition', R.condition, ...
        'seed', R.seed(1), 'n', n, 'N', R.N, 'R_struct', G.R_struct, 's2', G.s2, 'mu', G.mu, ...
        'lambda_1', lya_first(R), ...
        'cdf_M', mean(cdf_M, 1), 'cdf_null', mean(cdf_null, 1), 'cdf_unigain', mean(cdf_unigain, 1), ...
        'cdf_theory', mean(cdf_theory, 1), ...
        'dens_M', mean(dens_M, 1), 'dens_null', mean(dens_null, 1), 'dens_theory', mean(dens_theory, 1), ...
        'gains_cal', G.d, 'R_eff', R_eff, 'p0', p0, 'cv_nz', cv_nz, 'pr', pr, ...
        'mean_gain', mean_gain, 'rms_gain', rms_gain, 'gain_by_type', by_type, ...
        'measure_names', {names}, 'm_M', mM, 'm_null', mNull, 'm_theory', mTheory, ...
        'm_uniform', mUniform, 'm_unigain', mUnigain, 'full', full, ...
        'p0_t', [gt.p0], 'cv_nz_t', [gt.cv_nz], 'pr_t', [gt.pr], 'rms_t', [gt.rms], 't', Gt.t, ...
        'ev_example', R.cal(k_ex).ev_M, 'R_eff_example', R_eff(k_ex), ...
        'check_ev', max(G.check_ev), 'check_moments_cal', max(G.check_moments), ...
        'check_moments_t', max(Gt.check_moments), 'frac_t_failed', mean(Gt.check_moments > 1e-6), ...
        'check_zeros', max(abs(G.check_zeros)), 'n_unseen', max([G.n_unseen, Gt.n_unseen]));
    if opts.verbose
        fprintf('%-38s p0 %.3f  cv %.2f  PR %.3f | shell/0.36: M %.2f null %.2f theory %.2f | band M %.1f | unstable M %.1f, full %.1f | check %.1e\n', ...
            files(f).name, mean(p0), mean(cv_nz), mean(pr), mean(mM(1, :)) / 0.36, mean(mNull(1, :)) / 0.36, ...
            mean(mTheory(1, :)) / 0.36, mean(mM(3, :)), mean(mM(4, :)), mean(full.n_unstable), runs{f}.check_ev);
    end
end
runs = [runs{:}];

% W alone: Harris et al. (2023) with the preset's statistics. mu_tilde and
% sigma_tilde are the class's normalized ones (per sqrt(n) of the paper's),
% so they go in times sqrt(n) and the radius comes out in units of W.
harris = [];
if ~isempty(runs)
    L = load(fullfile(files(1).folder, files(1).name), 'R');
    model = build_from_preset(L.R.preset, L.R.condition, 'rng_seeds', L.R.seed, ...
        'T_range', [0 1], 'verbose', 'near-none');
    p = model.get_params();
    sc = sqrt(p.n) * p.level_of_chaos;
    f_e = p.n_per_type(1) / p.n;
    harris = edens_harris_density(p.alpha, f_e, sc * p.mu_tilde(1, 1), sc * p.mu_tilde(1, 2), ...
        sc * p.sigma_tilde(1, 1), sc * p.sigma_tilde(1, 2));
    harris.R_class = model.R;
end

% The switching run of Figure 3: the full Jacobian's eigenvalues at 40
% states per condition (the per-neuron gains were not saved, only their rms).
switching = [];
if isfile(opts.switching_file)
    D = load(opts.switching_file, 'cond_names', 'evals_by_cond', 'disk_by_cond', 'state_times');
    n_c = numel(D.cond_names);
    switching = struct('cond_names', {D.cond_names}, 'n_band', zeros(1, n_c), 'n_unstable', zeros(1, n_c), ...
        'N', zeros(1, n_c), 'rms_gain', zeros(1, n_c), 'uniform_band', zeros(1, n_c), 'uniform_unstable', zeros(1, n_c));
    n_states = numel(D.state_times);
    for c = 1:n_c
        a = edens_measures('rates', D.evals_by_cond{c});
        dk = D.disk_by_cond{c};
        switching.n_band(c) = a.n_band / n_states;
        switching.n_unstable(c) = a.n_unstable / n_states;
        switching.N(c) = numel(D.evals_by_cond{c}) / n_states;
        switching.rms_gain(c) = mean(dk.rms_gain);
        mu_ = arrayfun(@(g) edens_measures('uniform', dk.R * g, dk.n), dk.rms_gain);
        switching.uniform_band(c) = mean([mu_.n_band]);
        switching.uniform_unstable(c) = mean([mu_.n_unstable]);
    end
end

E = struct('runs', runs, 'switching', switching, 'harris', harris, ...
    'settings', struct('measures', S, 'u', u, 'u2_edges', u2_edges, 'gain_step', opts.gain_step));
file = fullfile(opts.out_dir, 'eigenvalue_density.mat');
save(file, 'E');
if opts.verbose; fprintf('saved %s\n', file); end
end

function g = gain_stats(d)
% Descriptors of one gain vector. pr is the participation ratio of the
% squared gains, mean(d^2)^2 / mean(d^4): the density at the rim of the
% bulk relative to a uniform disk (Ahmadian et al. 2015, Eq. 2.18) when
% the weights have one variance.
d = d(:); nz = d(d ~= 0);
g = struct('p0', mean(d == 0), 'mean', mean(d), 'rms', sqrt(mean(d.^2)), ...
    'mean_nz', mean(nz), 'cv_nz', std(nz, 1) / mean(nz), 'pr', mean(d.^2)^2 / mean(d.^4));
end

function c = radial_cdf(ev, R, u)
a = sort(abs(ev(:))) / R;
c = sum(a <= u, 1) / numel(a);
end

function h = radial_density(ev, R, u2_edges)
% Density relative to a uniform disk of radius R, in equal-area bins, the
% exact zeros (the zero-gain columns) left out.
a2 = (abs(ev(:)) / R).^2;
a2 = a2(a2 > 1e-16);
h = histcounts(a2, u2_edges) ./ numel(ev) ./ diff(u2_edges);
end

function v = as_column(m, names)
v = cellfun(@(f) m.(f), names(:));
end

function l = lya_first(R)
l = NaN;
if isfield(R, 'lya') && isfield(R.lya, 'LLE') && ~isempty(R.lya.LLE)
    l = R.lya.LLE;
end
end

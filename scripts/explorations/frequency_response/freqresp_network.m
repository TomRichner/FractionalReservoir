function [file, A] = freqresp_network(network, condition, seed, opts)
%FREQRESP_NETWORK The network frequency response of one run: the spectrum of M(i w) over a grid of frequencies, at its calibration states.
%
%   [file, A] = freqresp_network(network, condition, seed, Name, Value, ...)
%
% At each of the run's 12 calibration states (freqresp_states) and each
% frequency f of the grid, with z = 2 pi i f:
%   the gains d(a, j; z) of every neuron with its slow variables free
%   (freqresp_gain, closed form), with adaptation alone free and with
%   depression alone free; from them the predicted bulk radius and the
%   eigenvalues of the mean part (freqresp_radius); and, with both free,
%   all eigenvalues of M(z) = W .* d.
% The condition without adaptation has no slow variables: its M does not
% depend on z and is computed once per state.
%
% THE EXACTNESS CHECK (opts.exact_states, default states 1 and 7; skipped
% without slow variables). All eigenvalues z_k of the dense full Jacobian
% are computed. For the n_lead rightmost, the one nearest the continuation
% of the outlier that outlier_run saved, and n_sample drawn at random:
%   sigma  smallest singular value of M(z_k) - (1 + tau_d z_k) I, over the
%          largest singular value of M(z_k)  (zero when 1 + tau_d z_k is an
%          eigenvalue of M(z_k))
%   eigd   the distance from 1 + tau_d z_k to the nearest eigenvalue of
%          M(z_k), over |1 + tau_d z_k|
%   dslow  the distance from z_k to the nearest eigenvalue of the slow
%          block J_ss (the elimination divides by z - J_ss: an eigenvalue
%          of J that is ALSO one of J_ss, a slow mode of a neuron whose
%          slope is zero, is not a root of the reduced problem and does
%          not need to be)
% and at the same states M(z) from the Jacobian's blocks (freqresp_schur)
% is compared with the closed form at exact_f frequencies: the relative
% difference of the gains and the off-diagonal remainder.
%
% THE BULK IN THE COMPLEX PLANE (opts.plane, default true). The predicted
% radius R(z) on a grid of complex z = sigma + i w. The bulk of the full
% Jacobian's network eigenvalues is where R(z) >= |1 + tau_d z|; its edge
% is the curve R(z) = |1 + tau_d z|, the frequency-dependent version of
% "radius 1".
%
% Output file: data/explorations/frequency_response/network_<network>_<condition>_s<a>-<b>.mat,
% variable A:
%   f (1 x nf, Hz), tau_d, t_cal, n, network, condition, seed
%   R        struct of n_states x nf: .both, .sfa, .std  (predicted radius)
%   R_inf    n_states x 1: the frozen radius (z -> infinity)
%   R_zero   struct of n_states x 1 at z = 0: .both, .sfa, .std
%   K_eig    C x nf x n_states complex, both free; K_inf, K_zero  C x n_states
%   ev       n x nf x n_states single complex: eigenvalues of M(2 pi i f)
%            (nf = 1 without adaptation)
%   ev_inf   n x n_states single complex: eigenvalues of the frozen M
%   pop      struct: slope, rate (n x n_states), d_inf (C x n x n_states),
%            route_diff (largest |d(1, j) - d(2, j)| over the largest |d|:
%            zero when a neuron has one gain for all its routes)
%   exact    struct array over exact_states (fields above, with z, kind)
%   plane    struct: sigma, w (rad/s), R (n_states x n_sigma x n_w single)
%
% See also: freqresp_gain, freqresp_radius, explore_frequency_response

arguments
    network   (1,:) char
    condition (1,:) char
    seed      (1,2) double = [7 8]
    opts.f    (1,:) double = logspace(-3, 2, 61)
    opts.exact_states (1,:) double = [1 7]
    opts.exact_f (1,:) double = [0.001 0.03 0.5 3 30]
    opts.n_lead   (1,1) double = 20
    opts.n_sample (1,1) double = 200
    opts.plane (1,1) logical = true
    opts.sigma (1,:) double = -12:0.1:8
    opts.w     (1,:) double = 0:0.25:30
    opts.out_dir (1,:) char = ''
    opts.verbose (1,1) logical = true
end

t0 = tic;
root = fileparts(which('setup_paths'));
out_dir = opts.out_dir;
if isempty(out_dir)
    out_dir = fullfile(root, 'data', 'explorations', 'frequency_response');
end
file = fullfile(out_dir, sprintf('network_%s_%s_s%d-%d.mat', network, condition, seed(1), seed(2)));
if isfile(file)
    if nargout >= 2; L = load(file); A = L.A; end
    return;
end

[~, Q] = freqresp_states(network, condition, seed, out_dir);
[~, params, ws] = freqresp_model(network, condition, seed);
n = params.n; C = params.n_cellTypes; tau_d = params.tau_d;
has_slow = params.state_layout.N > n;
n_states = size(Q.S_cal, 2);
f = opts.f; nf = numel(f);
z = 2i * pi * f;
only_sfa = struct('sfa', true, 'std', false);
only_std = struct('sfa', false, 'std', true);

R = struct('both', zeros(n_states, nf), 'sfa', zeros(n_states, nf), 'std', zeros(n_states, nf));
R_inf = zeros(n_states, 1);
R_zero = struct('both', zeros(n_states, 1), 'sfa', zeros(n_states, 1), 'std', zeros(n_states, 1));
K_eig = complex(zeros(C, nf, n_states));
K_inf = complex(zeros(C, n_states)); K_zero = K_inf;
nf_ev = nf; if ~has_slow; nf_ev = 1; end
ev = complex(zeros(n, nf_ev, n_states, 'single'));
ev_inf = complex(zeros(n, n_states, 'single'));
pop = struct('slope', zeros(n, n_states), 'rate', zeros(n, n_states), ...
    'd_inf', zeros(C, n, n_states), 'route_diff', 0);

for k = 1:n_states
    S = Q.S_cal(:, k);
    [d, d_inf, parts] = freqresp_gain(S, params, [z, 0]);
    [Rk, Kk] = freqresp_radius(d, ws);
    R.both(k, :) = Rk(1:nf);  R_zero.both(k) = Rk(end);
    K_eig(:, :, k) = Kk(:, 1:nf);  K_zero(:, k) = Kk(:, end);
    [R_inf(k), K_inf(:, k)] = freqresp_radius(complex(d_inf), ws);
    Rk = freqresp_radius(freqresp_gain(S, params, [z, 0], only_sfa), ws);
    R.sfa(k, :) = Rk(1:nf);  R_zero.sfa(k) = Rk(end);
    Rk = freqresp_radius(freqresp_gain(S, params, [z, 0], only_std), ws);
    R.std(k, :) = Rk(1:nf);  R_zero.std(k) = Rk(end);

    pop.slope(:, k) = parts.slope; pop.rate(:, k) = parts.rate; pop.d_inf(:, :, k) = d_inf;
    if C > 1
        pop.route_diff = max(pop.route_diff, ...
            max(abs(d(1, :, :) - d(2, :, :)), [], 'all') / max(abs(d), [], 'all'));
    end

    ev_inf(:, k) = single(eig(freqresp_matrix(d_inf, ws)));
    for q = 1:nf_ev
        ev(:, q, k) = single(eig(freqresp_matrix(d(:, :, q), ws)));
    end
    say(opts, '  [%s %s] state %d of %d: spectra done (%.0f s)\n', network, condition, k, n_states, toc(t0));
end

%% Exactness
exact = [];
if has_slow
    cont_ref = [];
    ref_file = fullfile(root, 'data', 'explorations', 'outliers', outlier_run_name(network, condition, seed, ''));
    if isfile(ref_file)
        L = load(ref_file, 'R');
        cont_ref = [L.R.dense.cont];
    end
    exact = cell(1, numel(opts.exact_states));
    for e = 1:numel(opts.exact_states)
        k = opts.exact_states(e);
        exact{e} = exactness(Q.S_cal(:, k), params, ws, opts, cont_ref, k);
        exact{e}.state = k;
        say(opts, '  [%s %s] exactness at state %d: sigma median %.2g, max %.2g (off the slow block: max %.2g); closed form against Schur %.2g, off-diagonal %.2g (%.0f s)\n', ...
            network, condition, k, median(exact{e}.sigma, 'omitnan'), max(exact{e}.sigma), ...
            max(exact{e}.sigma(exact{e}.dslow > 1e-6)), max(exact{e}.closed_diff), max(exact{e}.offdiag), toc(t0));
    end
    exact = [exact{:}];
end

%% The predicted radius over the complex plane
plane = [];
if opts.plane
    [SG, WW] = ndgrid(opts.sigma, opts.w);
    zz = complex(SG(:), WW(:)).';
    Rp = zeros(n_states, numel(zz), 'single');
    chunk = 4000;
    for k = 1:n_states
        if ~has_slow
            Rp(k, :) = single(R_inf(k));
            continue;
        end
        for c0 = 1:chunk:numel(zz)
            sel = c0:min(c0 + chunk - 1, numel(zz));
            Rp(k, sel) = single(freqresp_radius(freqresp_gain(Q.S_cal(:, k), params, zz(sel)), ws));
        end
    end
    plane = struct('sigma', opts.sigma, 'w', opts.w, ...
        'R', reshape(Rp, n_states, numel(opts.sigma), numel(opts.w)));
    say(opts, '  [%s %s] plane done (%.0f s)\n', network, condition, toc(t0));
end

A = struct('network', network, 'condition', condition, 'seed', seed, 'n', n, ...
    'tau_d', tau_d, 't_cal', Q.t_cal, 'f', f, 'R', R, 'R_inf', R_inf, 'R_zero', R_zero, ...
    'K_eig', K_eig, 'K_inf', K_inf, 'K_zero', K_zero, 'ev', ev, 'ev_inf', ev_inf, ...
    'pop', pop, 'exact', exact, 'plane', plane, 'seconds', toc(t0));
save(file, 'A');
end

function E = exactness(S, params, ws, opts, cont_ref, k)
tau_d = params.tau_d;
n = params.n;
J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
zk = eig(full(J));
[~, order] = sort(real(zk), 'descend');
zk = zk(order);
slow = find(~ws.is_x);
z_slow = eig(full(J(slow, slow)));

n_lead = min(opts.n_lead, numel(zk));
pick = (1:n_lead)';
kind = ones(n_lead, 1);                        % 1 leading, 2 continuation, 3 random
if ~isempty(cont_ref)
    [~, ic] = min(abs(zk - cont_ref(k)));
    pick(end + 1) = ic; kind(end + 1) = 2;
end
stream = RandStream('twister', 'Seed', 1);
rest = setdiff((1:numel(zk))', pick);
extra = rest(randperm(stream, numel(rest), min(opts.n_sample, numel(rest))));
pick = [pick; extra]; kind = [kind; 3 * ones(numel(extra), 1)];

m = numel(pick);
sigma = zeros(m, 1); eigd = zeros(m, 1); dslow = zeros(m, 1);
I = eye(n);
for q = 1:m
    zq = zk(pick(q));
    Mz = freqresp_matrix(freqresp_gain(S, params, zq), ws);
    target = 1 + tau_d * zq;
    sv = svd(Mz - target * I);
    sigma(q) = sv(end) / norm(Mz);
    eigd(q) = min(abs(eig(Mz) - target)) / abs(target);
    dslow(q) = min(abs(zq - z_slow));
end

zc = 2i * pi * opts.exact_f;
closed_diff = zeros(1, numel(zc)); offdiag = zeros(1, numel(zc));
ok = ws.col_sq > 0;
for q = 1:numel(zc)
    [~, dn, offdiag(q)] = freqresp_schur(J, params, ws, zc(q));
    dc = freqresp_gain(S, params, zc(q));
    closed_diff(q) = max(abs(dc(ok) - dn(ok))) / max(abs(dn(ok)));
end

E = struct('z', zk(pick), 'kind', kind, 'sigma', sigma, 'eigd', eigd, 'dslow', dslow, ...
    'n_eig', numel(zk), 'exact_f', opts.exact_f, 'closed_diff', closed_diff, 'offdiag', offdiag);
end

function say(opts, varargin)
if opts.verbose; fprintf(varargin{:}); end
end

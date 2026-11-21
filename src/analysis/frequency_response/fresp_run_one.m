function A = fresp_run_one(st, opts)
%FRESP_RUN_ONE The frequency response of one run (network, condition, seed) at its calibration states.
%
%   A = fresp_run_one(st, opts)
%
%   st    one run's calibration states and model parameters (fresp_load_states)
%   opts  the stage's settings (frequency_response_stage_plan, plan.opts)
%
% THE OBJECT. Linearize about a state and let every variable go as
% exp(z t). A neuron's adaptation and depression variables are driven by
% its own rate only, so they are eliminated neuron by neuron: the mode sees
% M(z) = W .* d(z), each presynaptic neuron's column scaled by its complex
% gain (fresp_gain), and z is an eigenvalue of the full Jacobian exactly
% when 1 + tau_d z is an eigenvalue of M(z). The frozen-gain matrix is
% z -> infinity; z = 0 is the fully adapted network.
%
% AT EVERY CALIBRATION STATE, on the frequency grid opts.f (z = 2 pi i f):
%   the predicted bulk radius of M(z) and the eigenvalues of its mean part
%   (fresp_radius), with adaptation and depression free ('both'), with
%   adaptation alone and with depression alone; the same at z = 0 and for
%   the frozen gains; and the largest real part among the eigenvalues of
%   M(z) (one dense n x n eig per frequency). A condition without slow
%   variables has an M that does not depend on z: computed once per state.
%   With opts.plane_all_states, also the predicted abscissa of the bulk
%   (fresp_plane_edge) from the radius on the grid opts.sigma + i opts.w.
%
% AT THE CHECK STATES (opts.exact_states; empty = the first and the middle
% calibration state, 1 and 7 of 12; the shown state is always one of them):
%   every eigenvalue of the dense full Jacobian, the exactness and locality
%   checks (fresp_exactness), the predicted radius over the complex grid,
%   and, where the frozen mean part has an eigenvalue outside the frozen
%   bulk, the root of the mean part alone (fresp_outlier_root), started at
%   the leftover mode.
%
% AT THE SHOWN STATE (opts.shown_state; empty = the last check state), for
% a condition with slow variables: single-neuron responses on the grid
% opts.f_single, G(f) = d(1, j; 2 pi i f) / d_inf(1, j) (the gain on the
% route to the excitatory type over the frozen gain), closed form and from
% the Jacobian's blocks for that neuron alone (fresp_neuron_schur), for
% three neurons chosen among those with nonzero slope:
%   'E median'  the excitatory neuron of median rate
%   'I median'  the inhibitory neuron of median rate
%   'E driven'  the excitatory neuron at the 90th percentile of rate
%
% A fields
%   network, condition, seed, n, N, tau_d, has_slow, t_cal (1 x n_states),
%   f (1 x nf), cont (1 x n_states: the leftover mode of the outlier_spectra
%   stage, 1/s, Im >= 0), exact_states, shown_state
%   R        struct of n_states x nf: .both, .sfa, .std (predicted radius)
%   R_inf    n_states x 1: the frozen radius;  R_zero  struct of
%            n_states x 1 at z = 0
%   K_eig    C x nf x n_states complex, both free; K_inf, K_zero  C x n_states
%   max_re   n_states x nf: rightmost real part of eig M(2 pi i f);
%            max_re_inf  n_states x 1 for the frozen M
%   edge     1 x n_states: predicted abscissa of the bulk (1/s; NaN where
%            not computed)
%   pop      struct of 1 x n_states: frac_on (share of neurons with nonzero
%            slope), slope_mean, d_inf_rms, d_inf_mean; route_diff (largest
%            |d(1, j) - d(2, j)| over the largest |d|: zero when a neuron
%            has one gain for all its routes)
%   structural   struct: K_struct_eig, R_struct (realized W), K_theory_eig,
%            R_theory (the ensemble), has_outlier (fresp_weight_stats)
%   exact    struct array over the check states: state, t, and the fields
%            of fresp_exactness without lam (empty checks without slow
%            variables)
%   plane    struct array over the check states: state, t, lam (N x 1
%            single complex, the full Jacobian's eigenvalues), abscissa,
%            R (n_sigma x n_w single), R_inf, edge, cont, root (NaN where
%            not defined), K_inf1
%   single   struct array (empty without slow variables): label, type,
%            neuron (NaN when the type has no neuron with nonzero slope),
%            slope, rate, F, G.both / .sfa / .std (1 x nf_single complex),
%            G0, max_diff (largest |closed - numeric| over |numeric|),
%            summary (fresp_summary of each)
%   seconds  struct: curves, plane, exact, single, all
%
% See also: frequency_response_stage_plan, run_frequency_response_stage,
%           fresp_load_states, fresp_gain, fresp_exactness

t_all = tic;
params = st.params;
ws = fresp_weight_stats(params);
n = params.n; C = params.n_cellTypes; tau_d = params.tau_d;
N = params.state_layout.N;
has_slow = N > n;
n_states = size(st.S, 2);
f = opts.f; nf = numel(f);
z = 2i * pi * f;
only_sfa = struct('sfa', true, 'std', false);
only_std = struct('sfa', false, 'std', true);

exact_states = opts.exact_states;
if isempty(exact_states); exact_states = unique([1, ceil((n_states + 1) / 2)]); end
shown_state = opts.shown_state;
if isempty(shown_state); shown_state = exact_states(end); end
exact_states = unique([exact_states(:)', shown_state]);
if any(exact_states > n_states)
    error('fresp_run_one:State', ...
        'The run %s / %s / seed [%d %d] has %d calibration states; state(s) %s asked for.', ...
        st.network, st.condition, st.seed(1), st.seed(2), n_states, mat2str(exact_states(exact_states > n_states)));
end

%% The curves against frequency, at every calibration state
t0 = tic;
R = struct('both', zeros(n_states, nf), 'sfa', zeros(n_states, nf), 'std', zeros(n_states, nf));
R_inf = zeros(n_states, 1);
R_zero = struct('both', zeros(n_states, 1), 'sfa', zeros(n_states, 1), 'std', zeros(n_states, 1));
K_eig = complex(zeros(C, nf, n_states));
K_inf = complex(zeros(C, n_states)); K_zero = K_inf;
max_re = zeros(n_states, nf); max_re_inf = zeros(n_states, 1);
pop = struct('frac_on', zeros(1, n_states), 'slope_mean', zeros(1, n_states), ...
    'd_inf_rms', zeros(1, n_states), 'd_inf_mean', zeros(1, n_states), 'route_diff', 0);
for k = 1:n_states
    S = st.S(:, k);
    [d, d_inf, parts] = fresp_gain(S, params, [z, 0]);
    [Rk, Kk] = fresp_radius(d, ws);
    R.both(k, :) = Rk(1:nf);  R_zero.both(k) = Rk(end);
    K_eig(:, :, k) = Kk(:, 1:nf);  K_zero(:, k) = Kk(:, end);
    [R_inf(k), K_inf(:, k)] = fresp_radius(complex(d_inf), ws);
    Rk = fresp_radius(fresp_gain(S, params, [z, 0], only_sfa), ws);
    R.sfa(k, :) = Rk(1:nf);  R_zero.sfa(k) = Rk(end);
    Rk = fresp_radius(fresp_gain(S, params, [z, 0], only_std), ws);
    R.std(k, :) = Rk(1:nf);  R_zero.std(k) = Rk(end);

    pop.frac_on(k) = mean(parts.slope > 0);
    pop.slope_mean(k) = mean(parts.slope);
    pop.d_inf_rms(k) = sqrt(mean(d_inf(:).^2));
    pop.d_inf_mean(k) = mean(d_inf(:));
    if C > 1
        pop.route_diff = max(pop.route_diff, ...
            max(abs(d(1, :, :) - d(2, :, :)), [], 'all') / max(abs(d), [], 'all'));
    end

    max_re_inf(k) = max(real(eig(fresp_matrix(d_inf, ws))));
    if has_slow
        for q = 1:nf
            max_re(k, q) = max(real(eig(fresp_matrix(d(:, :, q), ws))));
        end
    else
        max_re(k, :) = max_re_inf(k);
    end
end
sec.curves = toc(t0);

%% The predicted radius over the complex plane
t0 = tic;
edge = nan(1, n_states);
plane_states = exact_states;
if opts.plane_all_states; plane_states = 1:n_states; end
R_plane = cell(1, n_states);
for k = plane_states
    Rp = plane_radius(st.S(:, k), params, ws, opts, has_slow, R_inf(k));
    edge(k) = fresp_plane_edge(Rp, opts.sigma, opts.w, tau_d);
    if ismember(k, exact_states); R_plane{k} = Rp; end
end
sec.plane = toc(t0);

%% The check states: the dense full Jacobian, exactness, locality, the mean-part root
t0 = tic;
exact = cell(1, numel(exact_states)); plane = exact;
for e = 1:numel(exact_states)
    k = exact_states(e);
    S = st.S(:, k);
    E = fresp_exactness(S, params, ws, struct('n_lead', opts.n_lead, 'n_sample', opts.n_sample, ...
        'exact_f', opts.exact_f, 'cont', st.cont(k)));
    root = NaN;
    if has_slow && real(K_inf(1, k)) > R_inf(k) && isfinite(st.cont(k))
        root = fresp_outlier_root(S, params, ws, st.cont(k));
    end
    plane{e} = struct('state', k, 't', st.t(k), 'lam', single(E.lam), 'abscissa', E.abscissa, ...
        'R', R_plane{k}, 'R_inf', R_inf(k), 'edge', edge(k), 'cont', st.cont(k), ...
        'root', root, 'K_inf1', K_inf(1, k));
    E = rmfield(E, 'lam');
    E.state = k; E.t = st.t(k);
    exact{e} = E;
end
exact = [exact{:}]; plane = [plane{:}];
sec.exact = toc(t0);

%% Single-neuron responses at the shown state
t0 = tic;
sn = [];
if has_slow
    sn = single_neurons(st.S(:, shown_state), params, ws, opts.f_single);
end
sec.single = toc(t0);
sec.all = toc(t_all);

A = struct('network', st.network, 'condition', st.condition, 'seed', st.seed, 'n', n, 'N', N, ...
    'tau_d', tau_d, 'has_slow', has_slow, 't_cal', st.t, 'f', f, 'cont', st.cont, ...
    'exact_states', exact_states, 'shown_state', shown_state, ...
    'R', R, 'R_inf', R_inf, 'R_zero', R_zero, 'K_eig', K_eig, 'K_inf', K_inf, 'K_zero', K_zero, ...
    'max_re', max_re, 'max_re_inf', max_re_inf, 'edge', edge, 'pop', pop, ...
    'structural', struct('K_struct_eig', ws.K_struct_eig, 'R_struct', ws.R_struct, ...
        'K_theory_eig', ws.K_theory_eig, 'R_theory', ws.R_theory, ...
        'has_outlier', ws.has_structural_outlier), ...
    'exact', exact, 'plane', plane, 'single', sn, 'seconds', sec);
end

%% ------------------------------------------------------------------------
function Rp = plane_radius(S, params, ws, opts, has_slow, R_inf)
% The predicted radius on the grid sigma + i w (n_sigma x n_w, single).
n_sigma = numel(opts.sigma); n_w = numel(opts.w);
if ~has_slow
    Rp = repmat(single(R_inf), n_sigma, n_w);
    return;
end
[SG, WW] = ndgrid(opts.sigma, opts.w);
zz = complex(SG(:), WW(:)).';
Rp = zeros(1, numel(zz), 'single');
chunk = 2000;                               % C x n x chunk complex gains at a time
for c0 = 1:chunk:numel(zz)
    sel = c0:min(c0 + chunk - 1, numel(zz));
    Rp(sel) = single(fresp_radius(fresp_gain(S, params, zz(sel)), ws));
end
Rp = reshape(Rp, n_sigma, n_w);
end

function own = single_neurons(S, params, ws, f)
% The three neurons, each among its type's neurons with nonzero slope.
picks = {'E median', 1, 0.5; 'I median', 2, 0.5; 'E driven', 1, 0.9};
modes = struct('both', struct('sfa', true, 'std', true), ...
    'sfa', struct('sfa', true, 'std', false), 'std', struct('sfa', false, 'std', true));
names = fieldnames(modes);
z = [2i * pi * f, 0];
[~, d_inf, parts] = fresp_gain(S, params, 0);
J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
d_mode = struct();
for m = 1:numel(names)
    d_mode.(names{m}) = fresp_gain(S, params, z, modes.(names{m}));
end

own = cell(1, size(picks, 1));
nan_row = complex(nan(1, numel(f)));
for p = 1:size(picks, 1)
    q = picks{p, 2};
    o = struct('label', picks{p, 1}, 'type', q, 'neuron', NaN, 'slope', NaN, 'rate', NaN, 'F', NaN, ...
        'G', struct('both', nan_row, 'sfa', nan_row, 'std', nan_row), ...
        'G0', struct('both', NaN, 'sfa', NaN, 'std', NaN), 'max_diff', NaN, 'summary', struct());
    idx = [];
    if q <= params.n_cellTypes; idx = params.type_indices{q}; end
    % A neuron with nonzero slope and a synapse onto the excitatory type.
    active = idx(parts.slope(idx) > 0 & d_inf(1, idx)' ~= 0);
    if ~isempty(active)
        [~, order] = sort(parts.rate(active));
        j = active(order(max(1, round(picks{p, 3} * numel(active)))));
        o.neuron = j; o.slope = parts.slope(j); o.rate = parts.rate(j); o.F = parts.F(1, j);
        o.max_diff = 0;
        for m = 1:numel(names)
            G = reshape(d_mode.(names{m})(1, j, :), 1, []) ./ d_inf(1, j);
            Gn = fresp_neuron_schur(J, params, ws, j, z, modes.(names{m})) ./ d_inf(1, j);
            o.G.(names{m}) = G(1:end - 1);
            o.G0.(names{m}) = real(G(end));
            o.max_diff = max(o.max_diff, max(abs(G - Gn) ./ abs(Gn)));
            o.summary.(names{m}) = fresp_summary(f, G(1:end - 1), real(G(end)));
        end
    end
    own{p} = o;
end
own = [own{:}];
end

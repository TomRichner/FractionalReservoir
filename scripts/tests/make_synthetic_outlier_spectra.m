function index_file = make_synthetic_outlier_spectra(run_dir, opts)
% MAKE_SYNTHETIC_OUTLIER_SPECTRA A tiny outlier_spectra output (FD-040's data contract), for the tests of the stages that read it.
%
%   index_file = MAKE_SYNTHETIC_OUTLIER_SPECTRA(run_dir)
%   index_file = MAKE_SYNTHETIC_OUTLIER_SPECTRA(run_dir, 'n', 80, 'T', 3, ...)
%
% A test fixture, not a stage: it writes <run_dir>/outlier_spectra/ with an
% index ospec_index.mat (table runs: network, condition, seed1, seed2,
% file, has_full_series, window_s) and one run_*.mat per network,
% condition and seed holding a struct R in the layout of
% docs/Plan_spectral_supplement.md, section 5, as far as the
% eigenvalue_density stage reads it (FD-042). Written while FD-040's stage
% was being built in parallel, so it stands in for that stage's output and
% is to be checked against it at the merge.
%
% THE NETWORK is the mu7chain preset at n neurons (indegree 0.2 n, F_ref at
% the same size: the size rule of small_network_model), built through
% build_from_preset; 'EEx1p5' sets mu_EE_relative to 1.5 times the
% preset's 7. Each run integrates T seconds with the preset's noise (sra1,
% fs Hz) and takes n_cal states evenly spaced over [T/3, T]. At each state,
% as scripts/explorations/outliers/outlier_effective_matrix.m and
% outlier_null.m do:
%   M = I + tau_d J_xx, the frozen-gain effective connectivity, read off
%   the class's own Jacobian, and the per-neuron gain d_j recovered from it
%   column by column (least squares against W, exact); the fixture checks
%   that M = W diag(d) holds (one gain per neuron) and errors otherwise;
%   ev_M = eig(M); ev_null, n_null draws of the same gains on W with its
%   block means removed and its entries permuted within each block;
%   ev_uniform = eig(W rms(d)); R_eff from the block variances;
%   the dense spectrum of the full Jacobian with each eigenvector's
%   fraction of norm^2 in the dendritic rows (xfrac).
% Also R.cal.u (argument of the nonlinearity relative to the neuron's
% setpoint), rate, cell_type, S, t, and R.structural.R_realized.
%
% LAYOUT. opts.layout 'matrix' writes R.cal as a scalar struct with the
% states in columns (gain n x K, ev_M n x K, ev_null n x q x K); 'struct'
% writes R.cal as a 1 x K struct array of per-state outlier_null-like
% results; 'mixed' (default) alternates by seed, so a reader is tested on
% both. See eigdens_load_runs.
%
% Name-value options (defaults): n 80, T 3, fs 500, n_cal 3, n_null 3,
% networks {'reference', 'EEx1p5'}, conditions {'no_adaptation',
% 'sfa1_std1', 'sfa3_stdc2'}, seeds [1 2; 2 3], layout 'mixed'.
%
% See also: eigdens_load_runs, test_eigenvalue_density_stage

arguments
    run_dir (1,:) char
    opts.n (1,1) double = 80
    opts.T (1,1) double = 3
    opts.fs (1,1) double = 500
    opts.n_cal (1,1) double = 3
    opts.n_null (1,1) double = 3
    opts.networks = {'reference', 'EEx1p5'}
    opts.conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
    opts.seeds (:,2) double = [1 2; 2 3]
    opts.layout (1,:) char {mustBeMember(opts.layout, {'matrix', 'struct', 'mixed'})} = 'mixed'
end

preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
folder = fullfile(run_dir, 'outlier_spectra');
if ~isfolder(folder); mkdir(folder); end
k_in = max(2, round(0.2 * opts.n));
size_args = {'n', opts.n, 'indegree', k_in, 'F_ref_n', opts.n, 'F_ref_indegree', k_in};
rs = RandStream('mt19937ar', 'Seed', 4242);          % the null permutations; private stream

rows = {};
for w = 1:numel(opts.networks)
    network = opts.networks{w};
    switch network
        case 'reference'; net_args = size_args;
        case 'EEx1p5';    net_args = [size_args, {'mu_EE_relative', 7 * 1.5}];
        otherwise; error('make_synthetic_outlier_spectra:Network', 'Network ''%s'' is not built here.', network);
    end
    for c = 1:numel(opts.conditions)
        condition = opts.conditions{c};
        for s = 1:size(opts.seeds, 1)
            seed = opts.seeds(s, :);
            layout = opts.layout;
            if strcmp(layout, 'mixed'); layout = ternary(mod(s, 2) == 1, 'matrix', 'struct'); end
            R = one_run(preset, network, condition, seed, net_args, opts, layout, rs);
            file = sprintf('run_%s_%s_s%d-%d.mat', network, condition, seed(1), seed(2));
            save(fullfile(folder, file), 'R');
            rows(end + 1, :) = {network, condition, seed(1), seed(2), file, true, opts.T * 2 / 3}; %#ok<AGROW>
        end
    end
end
runs = cell2table(rows, 'VariableNames', {'network', 'condition', 'seed1', 'seed2', 'file', ...
    'has_full_series', 'window_s'});
index_file = fullfile(folder, 'ospec_index.mat');
save(index_file, 'runs');
end

%% ------------------------------------------------------------------------
function R = one_run(preset, network, condition, seed, net_args, opts, layout, rs)
model = build_from_preset(preset, condition, net_args{:}, 'rng_seeds', seed, ...
    'T_range', [0 opts.T], 'fs', opts.fs, 'store_full_state', true, 'verbose', 'near-none', ...
    'lya_method', 'none');
model.run();
p = model.get_params();
n = p.n; C = p.n_cellTypes; ti = p.type_indices; L = p.state_layout;
W = full(p.W);
type_of = zeros(n, 1);
mu = zeros(C); s2 = zeros(C); col_sq = zeros(C, n);
for a = 1:C
    type_of(ti{a}) = a;
    for b = 1:C
        blk = W(ti{a}, ti{b});
        mu(a, b) = mean(blk(:)); s2(a, b) = var(blk(:), 1);
    end
    col_sq(a, :) = sum(W(ti{a}, :).^2, 1);
end
n_per = cellfun(@numel, ti);
Wc = W;
for a = 1:C
    for b = 1:C
        Wc(ti{a}, ti{b}) = Wc(ti{a}, ti{b}) - mu(a, b);
    end
end

t_cal = linspace(opts.T / 3, opts.T, opts.n_cal);
[~, idx] = min(abs(model.t_out(:) - t_cal), [], 1);
K = numel(idx);
S = model.S_out(idx, :);                               % K x N
gain = zeros(n, K); u = gain; rate = gain;
ev_M = complex(zeros(n, K)); ev_uniform = ev_M; ev_null = complex(zeros(n, opts.n_null, K));
R_eff = zeros(1, K);
dense = struct('lam', cell(1, K), 'xfrac', cell(1, K));
for k = 1:K
    Sk = S(k, :)';
    J = SRNNCellTypePairs.compute_Jacobian_fast(Sk, p);
    M = full(J(L.x, L.x)) * p.tau_d;
    M(1:n + 1:end) = M(1:n + 1:end) + 1;
    % The gain of neuron j on its route to each postsynaptic type, then one
    % per neuron (checked: the routes of a neuron carry the same factor).
    d_route = nan(C, n);
    for a = 1:C
        ok = col_sq(a, :) > 0;
        d_route(a, ok) = sum(M(ti{a}, ok) .* W(ti{a}, ok), 1) ./ col_sq(a, ok);
    end
    d = max(d_route, [], 1, 'omitnan')';
    d(isnan(d)) = 0;
    if max(abs(M - W .* d'), [], 'all') > 1e-9 * max(1, max(abs(M(:))))
        error('make_synthetic_outlier_spectra:Gains', 'M is not W diag(d) at %s %s: the routes of a neuron differ.', ...
            network, condition);
    end
    gain(:, k) = d;
    x_eff = Sk(L.x);
    for q = 1:C
        if p.n_a(q) > 0
            a_q = reshape(Sk(L.a{q}), p.n_per_type(q), p.n_a(q));
            x_eff(ti{q}) = x_eff(ti{q}) - p.c_eff(q) .* sum(a_q, 2);
        end
    end
    u(:, k) = x_eff - p.S_c_vec(:);
    rate(:, k) = p.activation_function(x_eff);
    ev_M(:, k) = eig(M);
    ev_uniform(:, k) = eig(W * sqrt(mean(d.^2)));
    for q = 1:opts.n_null
        Wn = Wc;
        for a = 1:C
            for b = 1:C
                blk = Wc(ti{a}, ti{b});
                blk(:) = blk(randperm(rs, numel(blk)));
                Wn(ti{a}, ti{b}) = blk;
            end
        end
        ev_null(:, q, k) = eig(Wn .* d');
    end
    V = zeros(C);
    for a = 1:C
        for b = 1:C
            V(a, b) = s2(a, b) * sum(d(ti{b}).^2);
        end
    end
    R_eff(k) = sqrt(max(0, max(real(eig(V)))));
    [Vj, Dj] = eig(full(J));
    lam = diag(Dj);
    [~, order] = sort(real(lam), 'descend');
    xn2 = sum(abs(Vj(L.x, order)).^2, 1);
    dense(k).lam = single(lam(order));
    dense(k).xfrac = single(xn2(:));
end

R = struct();
R.network = network; R.condition = condition;
R.preset_name = preset; R.net_args = net_args; R.rng_seeds = seed;
R.n = n; R.N = model.N_sys_eqs; R.tau_d = p.tau_d; R.n_per_type = p.n_per_type;
R.structural = struct('R', model.R, 'R_realized', sqrt(max(0, max(real(eig(s2 .* n_per(:)'))))));
R.lya = struct();
R.cal_t = model.t_out(idx)';
switch layout
    case 'matrix'
        R.cal = struct('t', R.cal_t, 'S', S, 'gain', gain, 'u', u, 'rate', rate, 'cell_type', type_of, ...
            'ev_M', ev_M, 'ev_null', ev_null, 'ev_uniform', ev_uniform, 'R_eff', R_eff);
    case 'struct'
        cal = cell(1, K);
        for k = 1:K
            cal{k} = struct('t', R.cal_t(k), 'S', S(k, :), 'gain', gain(:, k), 'u', u(:, k), ...
                'rate', rate(:, k), 'cell_type', type_of, 'ev_M', ev_M(:, k), ...
                'ev_null', ev_null(:, :, k), 'ev_uniform', ev_uniform(:, k), 'R_eff', R_eff(k));
        end
        R.cal = [cal{:}];
end
R.dense = dense;
end

function v = ternary(c, a, b)
if c; v = a; else; v = b; end
end

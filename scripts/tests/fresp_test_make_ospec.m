function index_file = fresp_test_make_ospec(run_dir, opts)
% FRESP_TEST_MAKE_OSPEC Test fixture: files in the layout of the outlier_spectra stage, from tiny networks, for the frequency_response stage's tests.
%
%   index_file = FRESP_TEST_MAKE_OSPEC(run_dir, Name, Value, ...)
%
% NOT A STAGE AND NOT THE OUTLIER_SPECTRA STAGE (FD-040 owns that). The
% frequency_response stage (FD-041) reads calibration states that
% outlier_spectra saves, and was built in parallel with it, so its tests
% need files of that layout made some other way. This builds each network
% through the preset machinery at a small size, integrates a few seconds,
% and writes exactly the part of the contract the frequency response reads
% (docs/Plan_spectral_supplement.md, section 5; fresp_load_states):
%
%   <run_dir>/outlier_spectra/ospec_index.mat   table runs: network,
%       condition, seed1, seed2, file, has_full_series, window_s
%   <run_dir>/outlier_spectra/<file>            struct R: network,
%       condition, rng_seeds, preset_name, net_args, cal.t (1 x n_cal),
%       cal.S (n_cal x N), dense(k).cont
% and nothing else of what the real stage saves. ONCE FD-040 IS MERGED the
% tests can call the real stage at tiny size instead, and this goes.
%
% The leftover mode dense(k).cont is found by the outlier exploration's
% rule (outlier_dense_spectrum): the eigenvalue of the dense full Jacobian
% whose eigenvector, restricted to the dendritic rows, overlaps most with
% the rightmost eigenvector of the frozen-gain matrix, among eigenvectors
% with at least 10 % of their norm^2 in the dendritic rows.
%
% Name-value options (defaults)
%   networks    {'reference', 'EEx1p5'}   EEx1p5: E-to-E mean x 1.5
%   conditions  {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
%   seeds       [1 2; 2 3]                one row per seed pair
%   n           60                        neurons
%   T_settle 4, T_measure 2 (s), fs 500, n_cal 3
%   preset      the mu7chain network
%
% See also: test_frequency_response_stage, fresp_load_states

arguments
    run_dir (1,:) char
    opts.networks   (1,:) cell = {'reference', 'EEx1p5'}
    opts.conditions (1,:) cell = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
    opts.seeds      (:,2) double = [1 2; 2 3]
    opts.n          (1,1) double = 60
    opts.T_settle   (1,1) double = 4
    opts.T_measure  (1,1) double = 2
    opts.fs         (1,1) double = 500
    opts.n_cal      (1,1) double = 3
    opts.preset     (1,:) char = ['celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_' ...
        'stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain']
end

folder = fullfile(run_dir, 'outlier_spectra');
if ~isfolder(folder); mkdir(folder); end
T_end = opts.T_settle + opts.T_measure;
size_args = {'n', opts.n, 'indegree', max(2, round(0.2 * opts.n)), 'F_tracks_network', true};

rows = cell(0, 7);
for w = 1:numel(opts.networks)
    network = opts.networks{w};
    switch network
        case 'reference'
            net_args = size_args;
        case 'EEx1p5'
            net_args = [size_args, {'mu_EE_relative', 1.5 * 7}];     % the preset's E-to-E mean is 7
        otherwise
            error('fresp_test_make_ospec:Network', 'Unknown network ''%s''.', network);
    end
    for c = 1:numel(opts.conditions)
        for s = 1:size(opts.seeds, 1)
            seed = opts.seeds(s, :);
            model = build_from_preset(opts.preset, opts.conditions{c}, net_args{:}, ...
                'rng_seeds', seed, 'T_range', [0 T_end], 'fs', opts.fs, ...
                'store_full_state', true, 'verbose', 'near-none', 'lya_method', 'none');
            model.run();
            params = model.get_params();
            first = find(model.t_out >= opts.T_settle, 1, 'first');
            idx = unique(round(linspace(first, numel(model.t_out), opts.n_cal)));
            R = struct();
            R.network = network; R.condition = opts.conditions{c}; R.rng_seeds = seed;
            R.preset_name = opts.preset; R.net_args = net_args;
            R.cal = struct('t', reshape(model.t_out(idx), 1, []), 'S', model.S_out(idx, :));
            R.dense = struct('cont', cell(1, numel(idx)));
            for k = 1:numel(idx)
                R.dense(k).cont = leftover_mode(model.S_out(idx(k), :)', params);
            end
            file = sprintf('run_%s_%s_s%d-%d.mat', network, opts.conditions{c}, seed(1), seed(2));
            save(fullfile(folder, file), 'R');
            rows(end + 1, :) = {network, opts.conditions{c}, seed(1), seed(2), file, false, opts.T_measure}; %#ok<AGROW>
        end
    end
end
runs = cell2table(rows, 'VariableNames', ...
    {'network', 'condition', 'seed1', 'seed2', 'file', 'has_full_series', 'window_s'});
index_file = fullfile(folder, 'ospec_index.mat');
save(index_file, 'runs');
end

%% ------------------------------------------------------------------------
function cont = leftover_mode(S, params)
L = params.state_layout;
n = params.n;
J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
Mf = full(J(L.x, L.x)) .* params.tau_d;
Mf(1:n + 1:end) = Mf(1:n + 1:end) + 1;
[Vm, Dm] = eig(Mf);
[~, i1] = max(real(diag(Dm)));
v1 = Vm(:, i1) ./ norm(Vm(:, i1));
[V, D] = eig(full(J));
lam = diag(D);
Vx = V(L.x, :);
xn2 = sum(abs(Vx).^2, 1);
ov = abs(v1' * Vx) ./ sqrt(max(xn2, realmin));
cand = find(xn2 >= 0.1);
[~, q] = max(ov(cand));
q = cand(q);
cont = complex(real(lam(q)), abs(imag(lam(q))));
end

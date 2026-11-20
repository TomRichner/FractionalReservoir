function R = ritz_pass_at_states(model, idx, K, scalings, evals, pass_opts)
% RITZ_PASS_AT_STATES The Jacobian projected on the top-K Lyapunov subspace at stored states.
%
%   R = RITZ_PASS_AT_STATES(model, idx, K, scalings, evals)
%   R = RITZ_PASS_AT_STATES(model, idx, K, scalings, evals, pass_opts)
%
% A top-K pass (the class's projected_jacobian) over the stored trajectory
% of a model that has run with store_full_state, returning H = Q' J Q at the
% rows idx of S_out, once per state scaling ('raw', 'type_sd'), with the
% bookkeeping the Ritz overlay's figure and table read (FD-033; the fields
% are described in run_eig_heatmap's header, whose Ritz pass this is). The
% initial basis is seeded as the model's own top-K pass is.
%
% evals: the pooled eigenvalues of the full Jacobian at the same states (N
% per state, in the order of idx), for the per-state comparison.
%
% pass_opts (optional struct): fields handed to projected_jacobian beside
% the seed and the scale, 'warmup_s' and 'lya_dt'. Without it the alignment
% is all the trajectory before the first sampled state. run_local_lyapunov
% names its own warm-up here so that the pass is, argument for argument,
% the model's own top-K pass (FD-034).
%
% See also: run_eig_heatmap, run_local_lyapunov, ritz_diagnostics,
%           SRNNCellTypePairs.projected_jacobian, jacobian_spectra_at_states
if nargin < 6 || isempty(pass_opts); pass_opts = struct(); end
cls    = class(model);
params = model.get_params();
seed   = model.rng_seeds(1) + 424242;
rows   = min(idx):numel(model.t_out);                 % the sampled window, for the scales
scales = feval([cls '.state_type_scales'], model.S_out, rows, params);
for s = 1:numel(scalings)
    opts = pass_opts;
    opts.seed = seed;
    if strcmp(scalings{s}, 'type_sd'); opts.state_scale = scales.d; end
    P = feval([cls '.projected_jacobian'], model.S_out, model.t_out, model.fs, idx, K, params, opts);
    if s == 1
        n_s = numel(idx); n_sc = numel(scalings);
        R = struct('K', P.K, 'N', P.N, 'scalings', {scalings}, ...
            'H', nan(P.K, P.K, n_s, n_sc), 'values', complex(nan(P.K, n_s, n_sc)), ...
            'valid', P.valid, 'sample_idx', P.sample_idx, 't', P.t, ...
            'LE_spectrum', nan(P.K, n_sc), 'n_positive', nan(1, n_sc), 'cond_max', nan(1, n_sc), ...
            'x_fraction', nan(n_s, P.K, n_sc), 'block_fractions', nan(n_s, 4, n_sc), ...
            'scales', rmfield(scales, 'd'), 'warmup_s', P.warmup_s, 'interval', P.interval, ...
            'lya_dt', P.lya_dt, 'seed', seed, 'seconds', nan(1, n_sc), ...
            'JQ_gram', nan(P.K, P.K, n_s, n_sc), 'full_num_abscissa', nan(n_s, n_sc));
        ev = reshape(evals, P.N, n_s);
        R.full_max_re = max(real(ev), [], 1)';
        R.full_n_positive = sum(real(ev) > 0, 1)';
    end
    R.H(:, :, :, s) = P.H;
    R.JQ_gram(:, :, :, s) = P.JQ_gram;
    for k = 1:numel(idx)
        % Numerical abscissa of the full Jacobian in this norm's coordinates.
        J = feval([cls '.compute_Jacobian_fast'], model.S_out(idx(k), :)', params);
        if isfield(opts, 'state_scale'); J = (1 ./ opts.state_scale) .* J .* opts.state_scale'; end
        R.full_num_abscissa(k, s) = largest_symmetric_eigenvalue(J);
    end
    for k = find(P.valid)'
        R.values(:, k, s) = eig(P.H(:, :, k));
    end
    R.valid = R.valid & P.valid;
    R.LE_spectrum(:, s) = P.LE_spectrum;
    R.n_positive(s) = P.n_positive;
    R.cond_max(s) = P.cond_max;
    R.x_fraction(:, :, s) = P.x_fraction;
    R.block_fractions(:, :, s) = P.block_fractions;
    R.seconds(s) = P.seconds;
end
end

function w = largest_symmetric_eigenvalue(J)
% max eig((J + J')/2): sparse Lanczos, dense when that does not converge.
Js = (J + J') / 2;
w = NaN;
if size(Js, 1) > 400
    [~, w, flag] = eigs(sparse(Js), 1, 'largestreal', 'Tolerance', 1e-8, 'MaxIterations', 500, 'Display', false);
    if flag ~= 0; w = NaN; end
end
if isnan(w); w = max(eig(full(Js))); end
end

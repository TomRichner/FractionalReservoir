function R = so_run_network(condition_name, seeds, S)
% SO_RUN_NETWORK One bursting network under a sustained DC step: trajectory, frozen-gain spectra along it, bursts and crossings.
%
%   R = so_run_network(condition_name, seeds, S)
%
% condition_name  'original' (bursting_pairs sfa3_std1, the condition of
%                 fig_stim_engages_adaptation and pilot 1b) or 'chain_mts',
%                 built by bursting_condition_model (the one definition)
% seeds           1 x 2 rng_seeds ([42 42] is the stimulation figure's)
% S               settings (run_stim_outlier): T_settle, T_show, t_stim,
%                 amp, burst, win_before, win_during, spec_every,
%                 full_every, do_full, do_frozen, frozen_every
%
% THE STIMULUS is pilot 1b's dc0p4: a noise-free DC step of S.amp to every
% neuron from S.t_stim to the end (p1b_dc_step), after S.T_settle s of
% settling. The model is deterministic (rk4, no noise).
%
% AT EVERY S.spec_every-TH SAMPLE (400 Hz integration) the frozen-gain
% effective connectivity M = W diag(d) (ospec_snapshot, which takes M from
% the class's own Jacobian, M = I + tau_d J_xx) is decomposed:
%   lam1      the rightmost eigenvalue (Im >= 0)
%   coh1      fraction of its eigenvector's norm^2 in the span of the
%             population-uniform patterns (constant on E, constant on I);
%             a random vector has 2/50 = 0.04
%   lam_coh   among the 25 rightmost modes, the most coherent one, and
%             coh_max its coherence ("the population mode" when it is not
%             the rightmost)
%   lam_ref   the structural block-mean mode followed by its eigenvector
%             (ospec_mode_reference 'blockmean'), ov_ref its overlap
%   R_eff     the bulk radius predicted from the gains
%   K_eig     the eigenvalues of the 2 x 2 gain-weighted mean matrix
%   dbar      2 x 2 mean gains (post, pre), rate 1 x 2 mean rates
%   ev        all 50 eigenvalues (single), for the complex-plane panels
% At every S.full_every-th sample, the rightmost network mode of the FULL
% Jacobian (ospec_full_spectrum, dense eig of the 225 x 225 matrix).
%
% THE FROZEN CONTROL (S.do_frozen): the same run integrated on from the
% last sample before the step with every state but x held (p2_frozen_rhs,
% block B pilot 2's method), spectra every S.frozen_every-th sample.
%
% See also: run_stim_outlier, so_crossings, so_main_figure, p1b_run_case,
%           p2_run_arms, ospec_snapshot

ic = struct('amp', S.amp, 't_on', S.t_stim, 'target', [], ...
    'intrinsic_drive', [], 'generator', @p1b_dc_step);
t_build = tic;
[m, values] = bursting_condition_model(condition_name, 'rng_seeds', seeds, ...
    'T_settle', S.T_settle, 'T_show', S.T_show, 'input_config', ic, ...
    'plot_deci', 4, 'model_args', {'store_full_state', true});
m.run();
sim_s = toc(t_build);

params = m.get_params();
ws  = ospec_weight_stats(params);
ref = ospec_mode_reference(ws, 'blockmean');
t   = m.t_out(:);
Sx  = m.S_out;

R = struct();
R.condition = condition_name;
R.seeds = seeds;
R.values = rmfield(values, 'input_config');
R.stim = rmfield(ic, 'generator');
R.tau_d = params.tau_d;
R.n_per_type = params.n_per_type;
R.struct = struct('ev_W', eig(ws.Wf), 'K_struct_eig', ws.K_struct_eig, ...
    'R_struct', ws.R_struct, 'R_struct_theory', ws.R_struct_theory, ...
    'lam0_blockmean', ref.lam0, 'coh_ref', ref.coh);
R.sim_s = sim_s;

% The M check at one state during the baseline: against W diag(F phi')
% formed here from the route outputs, and against I + tau_d times a
% central finite difference of dynamics_fast.
R.check = m_check(Sx(find(t >= 20, 1), :)', t(find(t >= 20, 1)), m, params, ws);

sel = find(t >= 0 & mod((0:numel(t) - 1)', S.spec_every) == 0);
R.intact = spectral_series(Sx, t, sel, params, ws, ref, m.u_ex);
if S.do_full
    self = find(t >= 0 & mod((0:numel(t) - 1)', S.full_every) == 0);
    R.full = full_series(Sx, t, self, params, ws);
end
R.intact.bursts = p1b_detect_bursts(R.intact.t, R.intact.rate(:, 1)', S.burst);
R.intact.meas = so_crossings(R.intact, S);

if S.do_frozen
    p = m.cached_params;
    p.u_interpolant = m.u_interpolant;
    x_idx = p.state_layout.x;
    k = find(t < S.t_stim, 1, 'last');          % last sample before the step
    m.build_noise();                             % a no-op here (no noise); as p2_run_arms
    [~, S2] = m.integrate(@(tt, s) p2_frozen_rhs(tt, s, p, x_idx), t(k:end), Sx(k, :)');
    Sf = [Sx(1:k - 1, :); S2];
    R.frozen_check = struct('t_freeze', t(k), ...
        'restart_mismatch', max(abs(S2(1, :) - Sx(k, :))), ...
        'adapt_drift', max(abs(S2(:, setdiff(1:p.state_layout.N, x_idx)) - ...
                               S2(1, setdiff(1:p.state_layout.N, x_idx))), [], 'all'));
    self = find(t >= S.t_stim - 2 & mod((0:numel(t) - 1)', S.frozen_every) == 0);
    R.frozen = spectral_series(Sf, t, self, params, ws, ref, m.u_ex);
    R.frozen.bursts = p1b_detect_bursts(R.frozen.t, R.frozen.rate(:, 1)', S.burst);
    R.frozen.meas = so_crossings(R.frozen, S);
    clear S2 Sf
end
end

function Q = spectral_series(Sx, t, sel, params, ws, ref, u_ex)
ns = numel(sel); n = params.n; n_top = 25;
Q = struct('t', t(sel)', 'lam1', complex(nan(ns, 1)), 'coh1', nan(ns, 1), ...
    'lam_coh', complex(nan(ns, 1)), 'coh_max', nan(ns, 1), ...
    'lam_ref', complex(nan(ns, 1)), 'ov_ref', nan(ns, 1), ...
    'R_eff', nan(ns, 1), 'K_eig', complex(nan(ns, 2)), 'dbar', nan(ns, 2, 2), ...
    'rate', nan(ns, 2), 'frac_on', nan(ns, 2), 'n_unstable', nan(ns, 1), ...
    'ev', single(complex(nan(ns, n))), 'u_mean', mean(u_ex(:, sel), 1)');
for q = 1:ns
    s = ospec_snapshot(Sx(sel(q), :)', params, ws, n_top, n, ref);
    Q.lam1(q) = s.ev_top(1);
    Q.coh1(q) = s.coh(1);
    [Q.coh_max(q), kk] = max(s.coh);
    Q.lam_coh(q) = s.ev_top(kk);
    Q.lam_ref(q) = s.mode.lam;
    Q.ov_ref(q) = s.mode.ov;
    Q.R_eff(q) = s.R_eff;
    Q.K_eig(q, :) = sort(s.K_eig, 'descend', 'ComparisonMethod', 'real').';
    Q.dbar(q, :, :) = s.dbar;
    Q.rate(q, :) = s.rate;
    Q.frac_on(q, :) = s.frac_on;
    Q.n_unstable(q) = s.n_unstable;
    Mf = ospec_effective_matrix(Sx(sel(q), :)', params, ws);
    Q.ev(q, :) = single(eig(Mf)).';
end
end

function F = full_series(Sx, t, sel, params, ws)
ns = numel(sel);
F = struct('t', t(sel)', 'net', complex(nan(ns, 1)), 'abscissa', nan(ns, 1), ...
    'net_coh', nan(ns, 1), 'kin_edge', nan(ns, 1));
for q = 1:ns
    f = ospec_full_spectrum(Sx(sel(q), :)', params, ws, 8);
    F.net(q) = f.net;
    F.abscissa(q) = f.abscissa;
    F.kin_edge(q) = f.kin_edge;
    k = find(f.xfrac >= 0.1, 1);
    if ~isempty(k); F.net_coh(q) = f.coh(k); end
end
end

function c = m_check(S, tt, m, params, ws)
[Mf, ~, J, ~, pop] = ospec_effective_matrix(S, params, ws);
p = params; p.u_interpolant = m.u_interpolant;
L = p.state_layout; n = p.n;
route = SRNNCellTypePairs.active_routes(p);
Mi = zeros(n);
for pre = 1:p.n_cellTypes
    for post = 1:p.n_cellTypes
        Rt = route{pre, post};
        pidx = p.type_indices{pre}; qidx = p.type_indices{post};
        [~, Fo] = Rt.model.rhs(reshape(S(Rt.idx), p.n_per_type(pre), Rt.P.ns), pop.r(pidx), Rt.P);
        Mi(qidx, pidx) = full(p.W(qidx, pidx)) .* (Fo(:) .* pop.rp(pidx))';
    end
end
h = 1e-7; Jfd = zeros(n);
for j = 1:n
    Sp = S; Sp(L.x(j)) = Sp(L.x(j)) + h;
    Sm = S; Sm(L.x(j)) = Sm(L.x(j)) - h;
    fp = SRNNCellTypePairs.dynamics_fast(tt, Sp, p);
    fm = SRNNCellTypePairs.dynamics_fast(tt, Sm, p);
    Jfd(:, j) = (fp(L.x) - fm(L.x)) / (2 * h);
end
c = struct('t', tt, 'max_abs_M', max(abs(Mf), [], 'all'), ...
    'M_vs_I_plus_tau_Jxx', max(abs(Mf - (eye(n) + p.tau_d * full(J(L.x, L.x)))), [], 'all'), ...
    'M_vs_W_diag_F_rp', max(abs(Mf - Mi), [], 'all'), ...
    'M_vs_I_plus_tau_Jxx_FD', max(abs(Mf - (eye(n) + p.tau_d * Jfd)), [], 'all'));
end

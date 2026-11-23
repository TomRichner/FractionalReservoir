function so_write_sidecar(D, E, md_path, FJ)
% SO_WRITE_SIDECAR The numbers behind the stim_outlier figures, as markdown.
%
%   so_write_sidecar(D, E, md_path, FJ)
%
% D: stim_outlier_main.mat (R_orig, R_chain, S); E: stim_outlier_seeds.mat
% (R_seeds); FJ: stim_outlier_fullJ_check.mat's FJ (optional: the full
% Jacobian's rightmost network mode at five states during the step, seeds
% 42, 1, 6, 3).
%
% See also: run_stim_outlier

S = D.S; R0 = D.R_orig;
fid = fopen(md_path, 'w');
w = @(varargin) fprintf(fid, varargin{:});
w('# Bursts, stimulation and the rightmost frozen-gain eigenvalue (pilot, %s)\n\n', datestr(now, 'yyyy-mm-dd'));
w('Code: `scripts/explorations/stim_outlier/` (run_stim_outlier). Data: `data/explorations/stim_outlier/`.\n\n');
w('Network `bursting_pairs` (n = 50, 35 E / 15 I, indegree 10), condition `original` = sfa3_std1 ');
w('(the condition of fig_stim_engages_adaptation and block B pilot 1b) unless said. ');
w('%g s settling, then %g s; DC step %.2f to every neuron from t = %g s, no noise (pilot 1b case dc0p4), rk4 400 Hz. ', ...
    S.T_settle, S.T_show, S.amp, S.t_stim);
w('Baseline window %s s, step window %s s. Bursts: mean E rate > %.2f for >= %.2f s, gaps < %.2f s merged. ', ...
    mat2str(S.win_before), mat2str(S.win_during), S.burst.thresh, S.burst.min_dur, S.burst.merge);
w('M = W diag(d) from ospec_snapshot every %d samples (%.0f ms). A crossing is a run of samples with Re lambda_1(M) > 1; ', ...
    S.spec_every, 1000 * S.spec_every / 400);
w('it coincides with a burst when it overlaps [onset - %.1f s, end + %.1f s].\n\n', S.tol_pre, S.tol_post);

w('## The connectivity\n\n');
w('Seed 42: outlier of W (block-mean mode, followed) %.3f; eig of the 2x2 block-mean matrix %s; bulk radius %.3f (theory %.3f); rightmost eig(W) %.3f.\n\n', ...
    real(R0.struct.lam0_blockmean), mat2str(real(R0.struct.K_struct_eig'), 4), R0.struct.R_struct, ...
    R0.struct.R_struct_theory, max(real(R0.struct.ev_W)));

c = R0.check;
w('## Check of M (seed 42, t = %g s)\n\n', c.t);
w('max|M| = %.3f. max|M - (I + tau_d J_xx)| = %.3g (the class Jacobian; ospec reads M from it). ', c.max_abs_M, c.M_vs_I_plus_tau_Jxx);
w('max|M - W diag(F phi'')| (F the route output factor, phi'' the slope, formed independently from the route models) = %.3g. ', c.M_vs_W_diag_F_rp);
w('max|M - (I + tau_d J_xx,FD)| with a central finite difference of dynamics_fast (h = 1e-7) = %.3g.\n\n', c.M_vs_I_plus_tau_Jxx_FD);

w('## Per network\n\n');
w('fb, fd: fraction of time Re lambda_1(M) > 1, baseline and step. ff: the same, step with every adaptation and depression variable frozen at its value at %.4f s. ', R0.frozen_check.t_freeze);
w('Gains d: mean over neurons, (post <- pre). coh: coherence of the rightmost eigenvector (a random vector: 0.04).\n\n');
w('| condition | seed | W outlier | fb | fd | ff | bursts | crossings | bursts w/ crossing | crossings in bursts | median lag (s) | Re l1 mean b / d | coh b / d | R_eff b / d | d(E<-E) b / d | d(I<-E) b / d | d(.<-I) b / d | rate E b / d | rate E frozen |\n');
w('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
Rs = [{D.R_orig, D.R_chain}, E.R_seeds];
for i = 1:numel(Rs)
    R = Rs{i}; b = R.intact.meas.before; d = R.intact.meas.during; f = R.frozen.meas.during;
    w('| %s | %d | %.2f | %.3f | %.3f | %.3f | %d | %d | %d | %d | %.2f | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f | %.3f / %.3f | %.3f / %.3f | %.3f / %.3f | %.3f / %.3f | %.3f |\n', ...
        R.condition, R.seeds(1), real(R.struct.lam0_blockmean), b.frac_beyond, d.frac_beyond, f.frac_beyond, ...
        b.n_bursts, b.n_cross, b.n_bursts_with_cross, b.n_cross_in_burst, b.median_lag, ...
        b.mean_re1, d.mean_re1, b.mean_coh1, d.mean_coh1, b.mean_R_eff, d.mean_R_eff, ...
        b.dbar(1, 1), d.dbar(1, 1), b.dbar(2, 1), d.dbar(2, 1), b.dbar(1, 2), d.dbar(1, 2), ...
        b.rate(1), d.rate(1), f.rate(1));
end
w('\nDuring the step no network bursts (0 bursts in 50-100 s in every row). Lag = onset of the coinciding crossing minus burst onset.\n\n');

F = R0.full; kb = F.t >= S.win_before(1) & F.t < S.win_before(2); kd = F.t >= S.win_during(1);
w('## Full Jacobian (seed 42, every %.0f ms)\n\n', 1000 * S.full_every / 400);
w('Rightmost network mode (eigenvector >= 10%% in x): baseline, Re > 0 for %.3f of the time, max %.2f 1/s. ', ...
    mean(real(F.net(kb)) > 0), max(real(F.net(kb))));
w('Step: no network mode among the 8 rightmost (those are single-neuron slow modes, abscissa %.3f 1/s).\n\n', mean(F.abscissa(kd)));
if nargin >= 4 && ~isempty(FJ)
    w('All 223 rightmost eigenvalues at five states during the step (t | network mode Re, Im in 1/s | its coherence | frozen-gain Re lambda_1(M), coh):\n\n');
    fn = fieldnames(FJ);
    for i = 1:numel(fn)
        X = FJ.(fn{i});
        w('- seed %s: ', fn{i}(2:end));
        for q = 1:size(X, 1)
            w('t %g: %.2f%+.2fi (coh %.2f), M %.3f (coh %.2f); ', X(q, 1), X(q, 2), X(q, 3), X(q, 4), X(q, 6), X(q, 7));
        end
        w('\n');
    end
    w('\n');
end
fclose(fid);
end

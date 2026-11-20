function [ev_all, omega, alpha, t_used] = jacobian_spectra_at_states(model, idx, use_parallel)
% JACOBIAN_SPECTRA_AT_STATES Eigenvalues of the Jacobian at stored states of a run.
%
%   [ev_all, omega, alpha, t_used] = JACOBIAN_SPECTRA_AT_STATES(model, idx, use_parallel)
%
% model has run with store_full_state; idx are rows of its S_out. Shared by
% run_eig_heatmap (the baseline run) and run_local_lyapunov (the switching
% run, FD-034).
%
% Pools the Jacobian eigenvalues at the requested rows of the stored
% trajectory, and per state the
% NUMERICAL ABSCISSA omega = max eig((J_xx + J_xx')/2) beside the SPECTRAL
% ABSCISSA alpha = max real eig(J_xx) of the DENDRITIC BLOCK
% J_xx = (-I + W diag(theta'(x)))/tau_d -- the rate-network Jacobian of
% Hennequin et al., with the adaptation and depression states held fixed.
% omega bounds the instantaneous growth rate of a perturbation of x
% (d/dt ||dx|| <= omega ||dx||), so omega > alpha is the margin by which the
% non-normal recurrent coupling can amplify transiently even when every
% eigenvalue decays (Trefethen & Embree). omega >= alpha always.
%
% WHY THE x BLOCK ONLY (measured 2026-09-12): on the FULL Jacobian the
% numerical abscissa is 50-180 /s against a spectral abscissa of 0-10, and it
% ranks the regimes by how their state coordinates are scaled rather than by
% dynamics -- the symmetric part mixes rows in units of 1/tau_d = 10 (x) with
% rows in units of 1/tau_a and 1/tau_rel (a, b), and the (x,a) vs (a,x)
% blocks differ by a factor W/(c tau_d) against 1/tau_a. The numerical
% abscissa is not invariant to a diagonal rescaling of the state, so it is
% only meaningful on a block with one unit. The full-J eigenvalues are still
% pooled for the heatmap as before.
%
% Resolves compute_Jacobian_fast BY CLASS NAME, so one implementation serves
% both model classes. They do not share a state layout -- SRNNCellTypePairs
% carries b-states per ROUTE where SRNNModel2 carries them per population -- so
% indexing S_out by hand would be class-specific and fragile; the x rows are
% the LAST n of the state on both classes.
%
% See also: run_eig_heatmap, run_local_lyapunov, ritz_pass_at_states
cls    = class(model);
params = model.get_params();
S_samp = model.S_out(idx, :);          % only the sampled rows travel to the pool workers
t_out  = model.t_out;
n      = model.n;

t_used = t_out(idx);
t_used = t_used(:);

n_idx = numel(idx);
ev = cell(1, n_idx);
omega = zeros(n_idx, 1);
alpha = zeros(n_idx, 1);
if use_parallel
    parfor k = 1:n_idx
        J = full(feval([cls '.compute_Jacobian_fast'], S_samp(k, :)', params));
        ev{k} = eig(J);
        Jxx = J(end - n + 1:end, end - n + 1:end);
        alpha(k) = max(real(eig(Jxx)));
        omega(k) = max(eig((Jxx + Jxx') / 2));
    end
else
    for k = 1:n_idx
        J = full(feval([cls '.compute_Jacobian_fast'], S_samp(k, :)', params));
        ev{k} = eig(J);
        Jxx = J(end - n + 1:end, end - n + 1:end);
        alpha(k) = max(real(eig(Jxx)));
        omega(k) = max(eig((Jxx + Jxx') / 2));
    end
end
ev_all = vertcat(ev{:});
end

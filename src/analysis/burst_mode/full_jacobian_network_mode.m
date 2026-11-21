function f = full_jacobian_network_mode(J, is_x, x_thr, U, groups)
%FULL_JACOBIAN_NETWORK_MODE The rightmost network mode of the full Jacobian at one state, by a dense eigendecomposition.
%
%   f = full_jacobian_network_mode(J, is_x, x_thr, U)
%   f = full_jacobian_network_mode(J, is_x, x_thr, U, groups)
%
% J       the full Jacobian (sparse or full), N x N, rates in 1/s
% is_x    N x 1 logical, the dendritic rows (ospec_weight_stats' is_x)
% x_thr   the dendritic fraction that makes an eigenvector a network mode
%         (0.1, the rule of ospec_full_spectrum and the outlier_spectra
%         stage: at least x_thr of the eigenvector's norm^2 in the x rows)
% U       n x p, orthonormal population indicators (ws.U: the E and I
%         patterns; or the ten subnet populations of the modular network)
% groups  optional n x 1 labels (1..G); f.share then holds the mode's
%         dendritic norm^2 per group
%
% WHY DENSE. ospec_full_spectrum asks eigs for the k rightmost eigenvalues
% above 1000 state variables. In the modular IED network (350 neurons,
% 1190 state variables) eigs converged on 2, 12 and 4 of 12 eigenvalues at
% the states tried, and none of the converged ones was a network mode
% (pilot of 2026-11-21, scripts/explorations/stim_outlier/soied_*). A dense
% eig with vectors takes 0.4 s there on 14 threads (2026-11-21, the
% workstation), so this does it densely at every size; at 225 state
% variables (the 50-neuron bursting network) it takes milliseconds.
%
% Fields (rates in 1/s)
%   net        the rightmost network mode, complex with Im >= 0 (NaN if no
%              eigenvector passes x_thr, which does not happen in these
%              networks: the n dendritic directions have to go somewhere)
%   net_xfrac  its dendritic fraction
%   net_coh    fraction of its dendritic part's norm^2 in the span of U
%   share      1 x G dendritic norm^2 per group (empty without groups)
%   abscissa   the rightmost real part of any eigenvalue
%   n_net_unstable   network modes with Re > 0 (a conjugate pair counts twice)
%
% Ported from the pilots' full_series and soied_event_full (2026-11-21).
%
% See also: ospec_full_spectrum, run_burst_mode_stage, run_ied_network_mode_stage

arguments
    J
    is_x (:,1) logical
    x_thr (1,1) double
    U double
    groups (:,1) double = []
end
[V, D] = eig(full(J), 'vector');
V = V ./ vecnorm(V);
e = abs(V).^2;
xfrac = sum(e(is_x, :), 1)';
[~, o] = sort(real(D), 'descend');
net = o(xfrac(o) >= x_thr);
f = struct('net', complex(NaN, NaN), 'net_xfrac', NaN, 'net_coh', NaN, 'share', [], ...
    'abscissa', real(D(o(1))), 'n_net_unstable', nnz(real(D(net)) > 0));
if isempty(net)
    return
end
k = net(1);
vx = V(is_x, k);
f.net = complex(real(D(k)), abs(imag(D(k))));
f.net_xfrac = xfrac(k);
f.net_coh = sum(abs(U' * vx).^2) / max(sum(abs(vx).^2), realmin);
if ~isempty(groups)
    ex = abs(vx).^2;
    f.share = accumarray(groups(:), ex, [max(groups) 1])' / sum(ex);
end
end

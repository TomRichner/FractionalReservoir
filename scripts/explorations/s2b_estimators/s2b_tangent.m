function R = s2b_tangent(M, seg, varargin)
%S2B_TANGENT Local growth rates along a stored trajectory, by the linear
% (variational) equation dY/dt = J(S(t)) Y, integrated with RK4 at step 2/fs
% on the stored grid (the Jacobian at the grid states; no interpolation).
%   R = s2b_tangent(M, seg, 'K', K, 'Q0', Q0, 'V0', V0, 'ks', [1 2 3 5 10 20])
% Two kinds of vectors are carried at once, over the same segments of length
% seg (s) from t(1):
%   Q (N x K)  re-orthonormalised by QR every segment: the QR method. Column
%              order is the Gram-Schmidt order and is NEVER re-sorted here.
%              Default Q0 = the first K columns of eye(N) (the stage's 'qr'
%              method starts from eye(N)).
%   V (N x m)  tracers: each column evolved ALONE and renormalised every
%              segment. Benettin's method in the limit of an infinitesimal
%              perturbation.
% Output: t (segment starts), loc_qr (n_seg x K, log|R_jj|/seg), loc_v
% (n_seg x m), proj_v (n_seg x m x numel(ks): norm of the projection of each
% tracer on span(Q(:,1:k))), rot_q1 (angle in rad between q1 before and after
% each segment), q1_flow (|cos| between q1 and the flow direction dS/dt), and
% v_flow (same for each tracer), Q_end, V_end.
o = struct('K', M.N, 'Q0', [], 'V0', zeros(M.N, 0), 'ks', [1 2 3 5 10 20], 't_stop', M.t(end));
for k = 1:2:numel(varargin); o.(varargin{k}) = varargin{k + 1}; end
N = M.N;  K = o.K;
if isempty(o.Q0); Q = eye(N, K); else; Q = o.Q0; end
V = o.V0 ./ vecnorm(o.V0);  m = size(V, 2);
ks = o.ks(o.ks <= K);
h = 2 / M.fs;
spp = round(seg / h);                       % RK4 steps per segment
assert(abs(spp * h - seg) < 1e-12 && spp >= 1, 'seg must be a multiple of 2/fs');
nt = find(M.t <= o.t_stop + 1e-9, 1, 'last');
n_seg = floor((nt - 1) / (2 * spp));
R = struct('seg', seg, 't', M.t(1) + (0:n_seg - 1)' * seg, 'loc_qr', zeros(n_seg, K), ...
    'loc_v', zeros(n_seg, m), 'proj_v', zeros(n_seg, m, numel(ks)), 'ks', ks, ...
    'rot_q1', zeros(n_seg, 1), 'q1_flow', zeros(n_seg, 1), 'v_flow', zeros(n_seg, m));
Jf = @(i) SRNNCellTypePairs.compute_Jacobian_fast(M.S(i, :)', M.params);
i = 1;
J2 = Jf(1);
for s = 1:n_seg
    Y = [Q, V];
    for q = 1:spp
        J0 = J2;  J1 = Jf(i + 1);  J2 = Jf(i + 2);
        k1 = J0 * Y;
        k2 = J1 * (Y + (h / 2) * k1);
        k3 = J1 * (Y + (h / 2) * k2);
        k4 = J2 * (Y + h * k3);
        Y = Y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        i = i + 2;
    end
    q1_old = Q(:, 1);
    [Q, Rr] = qr(Y(:, 1:K), 0);
    sg = sign(diag(Rr));  sg(sg == 0) = 1;
    Q = Q .* sg';                           % keep R's diagonal positive
    R.loc_qr(s, :) = log(abs(diag(Rr)))' / seg;
    R.rot_q1(s) = acos(min(1, abs(q1_old' * Q(:, 1))));
    f = M.rhs(M.t(i), M.S(i, :)');  f = f / norm(f);
    R.q1_flow(s) = abs(f' * Q(:, 1));
    if m > 0
        nv = vecnorm(Y(:, K + 1:end));
        V = Y(:, K + 1:end) ./ nv;
        R.loc_v(s, :) = log(nv) / seg;
        R.v_flow(s, :) = abs(f' * V);
        for a = 1:numel(ks)
            R.proj_v(s, :, a) = vecnorm(Q(:, 1:ks(a))' * V, 2, 1);
        end
    end
end
R.Q_end = Q;  R.V_end = V;
end

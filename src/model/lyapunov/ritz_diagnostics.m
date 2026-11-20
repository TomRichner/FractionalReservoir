function D = ritz_diagnostics(H, G, k)
%RITZ_DIAGNOSTICS What the Ritz values of a projected Jacobian are, and how far to trust them.
%   D = ritz_diagnostics(H, G)
%   D = ritz_diagnostics(H, G, k)
%
%   H = Q' J Q is the Jacobian projected on an orthonormal basis Q (N x K)
%   and G = (J Q)' (J Q) the Gram matrix of J Q, both K x K, as
%   lyapunov_topk returns them for one state (H_t, H_gram_t). k (default K)
%   takes the leading k x k blocks: the top-k subspace, since the columns
%   are in Gram-Schmidt order. Everything below follows from H and G alone,
%   because Q' (J Q) = H: for any vector y in the subspace,
%       ||J Q y - Q H y||^2 = y' G y - ||H y||^2 .
%
%   The eigenvalues theta of H (Ritz values) are eigenvalues of J only if
%   span(Q) is an invariant subspace of J. The subspace the Lyapunov
%   computation tracks is carried by the DYNAMICS, not by the instantaneous
%   Jacobian, so in general it is not, and then the Ritz values are points
%   of the numerical range of J: what J does to vectors of that subspace,
%   seen from inside it. The RESIDUALS say which case holds.
%
%   OUTPUT struct D
%     values           k x 1 Ritz values (eigenvalues of the block)
%     residual         ||J Q - Q H||_F / ||J Q||_F in [0, 1]: 0 for an
%                      invariant subspace; near 1 when J throws the
%                      subspace's vectors almost entirely out of it
%     pair_residual    k x 1, per Ritz pair (theta, y): ||J Q y - theta Q y||
%                      / ||J Q y||; small means (theta, Q y) is nearly an
%                      eigenpair of J (floor about 1e-8: a difference of
%                      squares)
%     diag             k x 1 diagonal of the block: the instantaneous local
%                      Lyapunov rates of the k directions
%     share_diag_pos, share_ritz_pos   shares of local rates and of Ritz
%                      values with positive (real part)
%     trace_err        |sum of the Ritz real parts - sum of the diagonal|
%     THE FORM THE THEORY GIVES. With Phi = Q R, H = T + S, T = (dR/dt) R^-1
%     upper triangular (the local rates on its diagonal) and S = Q' dQ/dt
%     skew: the turning of the frame. The strictly lower triangle of H is
%     that of S, so S and T are read off H:
%     norm_diag        ||diag(H)||_2
%     norm_lower       ||strictly lower triangle of H||_F = ||S||_F / sqrt(2)
%     norm_upper_T     ||strictly upper triangle of T||_F (the feed-forward,
%                      non-normal coupling once the frame's turning is removed)
%     coupling_TS      max over i < j of sqrt(|T_ij S_ji|) (1/s): in a 2 x 2
%                      block the eigenvalues leave the diagonal by about the
%                      square root of the product of the two off-diagonal
%                      entries, so this is the scale by which feed-forward
%                      weight times frame rotation moves a Ritz value
%     coupling_H       the same for the entries of H itself, sqrt(|H_ij H_ji|)
%     num_abscissa     largest eigenvalue of (H + H')/2: the largest
%                      instantaneous growth rate within the subspace, an
%                      upper bound on every Ritz real part (and itself
%                      bounded by the numerical abscissa of J)
%     max_re           largest Ritz real part
%
%   See also lyapunov_topk, SRNNCellTypePairs.projected_jacobian, fig_eig_heatmap
    if nargin < 3 || isempty(k); k = size(H, 1); end
    Hk = H(1:k, 1:k);
    Gk = G(1:k, 1:k);
    [Y, Th] = eig(Hk);                              % columns of unit 2-norm
    theta = diag(Th);
    jq2 = real(sum(conj(Y) .* (Gk * Y), 1)).';      % ||J Q y||^2
    hy2 = real(sum(abs(Hk * Y) .^ 2, 1)).';         % ||H y||^2 = |theta|^2 for unit y
    d = diag(Hk);
    L = tril(Hk, -1);
    U = triu(Hk, 1);
    T_upper = U + L.';                              % strict upper triangle of T = H - (L - L')
    D = struct();
    D.values         = theta;
    D.residual       = sqrt(max(trace(Gk) - norm(Hk, 'fro') ^ 2, 0) / trace(Gk));
    D.pair_residual  = sqrt(max(jq2 - hy2, 0) ./ jq2);
    D.diag           = d;
    D.share_diag_pos = mean(d > 0);
    D.share_ritz_pos = mean(real(theta) > 0);
    D.trace_err      = abs(sum(real(theta)) - sum(d));
    D.norm_diag      = norm(d);
    D.norm_lower     = norm(L, 'fro');
    D.norm_upper_T   = norm(T_upper, 'fro');
    D.coupling_TS    = sqrt(max(abs(T_upper .* L.'), [], 'all'));
    D.coupling_H     = sqrt(max(abs(U .* L.'), [], 'all'));
    D.num_abscissa   = max(eig((Hk + Hk') / 2));
    D.max_re         = max(real(theta));
    if k == 1; D.coupling_TS = 0; D.coupling_H = 0; end
end

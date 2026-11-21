function ref = ospec_mode_reference(ws, target, target_value, top_frac)
%OSPEC_MODE_REFERENCE The intended mode of a network: one eigenvalue of the structural W, its two eigenvectors and which neurons carry it.
%
%   ref = ospec_mode_reference(ws, target, target_value)
%
% ws is ospec_weight_stats(params): ws.Wf is the matrix the dynamics use,
% and with every gain equal to d the effective connectivity is d ws.Wf, so
% the eigenvectors found here are those of the uniform-gain picture at any
% d. target says which eigenvalue is meant (ospec_network):
%   'lowrank'    the one nearest target_value (the eigenvalue the added
%                structure was built to have, units of W)
%   'blockmean'  the one nearest the rightmost eigenvalue of the realized
%                block-mean matrix ws.K_struct
%   'rightmost'  the rightmost
% A conjugate pair is one mode; the member with Im >= 0 represents it.
%
% WHY THE LEFT EIGENVECTOR. The gains scale the columns of W, M = W D, and
% to first order in the departure of D from uniform an eigenvalue moves by
%   delta lambda / lambda = sum_j p_j (d_j / d_ref - 1),
%   p_j = conj(l_j) r_j / (l' r),      sum_j p_j = 1,
% with r and l the right and left eigenvectors (l' W = lambda l'). p_j is
% how much of the mode neuron j carries: its gain matters in proportion to
% both how strongly the mode drives the network through it (r_j) and how
% strongly the network reads the mode back from it (l_j). The p_j are
% complex in general; |p_j| ranks the carriers.
%
% Fields
%   lam0       the eigenvalue (units of W), Im >= 0
%   r0, l0     unit right and left eigenvectors (n x 1 complex)
%   p          n x 1 complex participation, sums to 1
%   top        the top_frac (default 0.1) of neurons by |p_j| (indices)
%   pr         participation ratio of |p| (how many neurons carry the mode)
%   coh        fraction of r0's norm^2 in the population-uniform span
%   cond       1 / |l0' r0| (eigenvalue condition number: 1 for a normal
%              mode, large when the first-order formula is fragile)
%   p_type     1 x C: sum of |p_j| over each type / sum of |p_j|
%   target, target_value, R_struct, bulk_ratio (|lam0| / R_struct),
%   re_ratio (Re lam0 / R_struct), is_rightmost
%   ev_W       n x 1: the whole structural spectrum
%
% Ported from scripts/explorations/outliers/outlier_mode_reference.m (FD-040), unchanged.
%
% See also: ospec_network, ospec_snapshot, ospec_run

if nargin < 3; target_value = []; end
if nargin < 4 || isempty(top_frac); top_frac = 0.1; end

[V, D, L] = eig(ws.Wf);
ev = diag(D);
[~, i_right] = max(real(ev));
switch target
    case 'lowrank'
        [~, i0] = min(abs(ev - target_value));
    case 'blockmean'
        kev = ws.K_struct_eig;
        [~, q] = max(real(kev));
        target_value = complex(real(kev(q)), abs(imag(kev(q))));
        [~, i0] = min(abs(ev - target_value));
    case 'rightmost'
        i0 = i_right;
    otherwise
        error('ospec_mode_reference:Target', 'Unknown target ''%s''.', target);
end
if imag(ev(i0)) < 0
    [~, i0] = min(abs(ev - conj(ev(i0))));
end
r0 = V(:, i0) ./ norm(V(:, i0));
l0 = L(:, i0) ./ norm(L(:, i0));
lr = l0' * r0;
p = conj(l0) .* r0 ./ lr;
ap = abs(p);
[~, order] = sort(ap, 'descend');
n = numel(p);
top = sort(order(1:max(1, round(top_frac * n))));
C = numel(ws.ti);
p_type = zeros(1, C);
for q = 1:C
    p_type(q) = sum(ap(ws.ti{q})) / sum(ap);
end

ref = struct('lam0', ev(i0), 'r0', r0, 'l0', l0, 'p', p, 'top', top, ...
    'pr', sum(ap)^2 / sum(ap.^2), 'coh', sum(abs(ws.U' * r0).^2), ...
    'cond', 1 / abs(lr), 'p_type', p_type, ...
    'target', target, 'target_value', target_value, 'R_struct', ws.R_struct, ...
    'bulk_ratio', abs(ev(i0)) / ws.R_struct, 're_ratio', real(ev(i0)) / ws.R_struct, ...
    'is_rightmost', abs(real(ev(i0)) - real(ev(i_right))) < 1e-9 * max(1, abs(ev(i_right))), ...
    'ev_W', ev);
end

function f = ospec_dense_spectrum(S, params, ws, ref, x_thr)
%OSPEC_DENSE_SPECTRUM The WHOLE spectrum of the full Jacobian at one state, and where the frozen-gain outlier went in it.
%
%   f = ospec_dense_spectrum(S, params, ws)
%   f = ospec_dense_spectrum(S, params, ws, ref)    also follows the INTENDED mode
%   f = ospec_dense_spectrum(S, params, ws, ref, x_thr)   the dendritic fraction that makes a network mode (0.1)
%
% WHY (found on the first real run, 2026-11-19). The rightmost eigenvalues
% of the full Jacobian (ospec_full_spectrum) are, with adaptation, the
% edge of the band of slow single-neuron modes: they say whether anything
% is unstable, not what became of the outlier. At E-to-E x 1.5 with
% multiple timescales the frozen-gain matrix M has a real outlier at about
% +1 1/s and the full Jacobian has no eigenvalue right of zero at all: with
% the slow variables free the same population pattern is a damped
% oscillation, a complex pair at about -1.6 +- 5.8i 1/s. To follow the
% outlier one has to look for its pattern, not for the rightmost
% eigenvalue. A dense eig of the 4000 x 4000 Jacobian takes about 15 s, so
% this is done at a few states per run.
%
% THE CONTINUATION of the frozen-gain outlier is the eigenvalue of the full
% Jacobian whose eigenvector, restricted to the dendritic rows and
% normalized there, overlaps most with the rightmost eigenvector v_1 of M,
% among eigenvectors with at least 10% of their norm^2 in the dendritic
% rows (a slow single-neuron mode has almost none and its x part is noise).
% Everything needed to choose differently is returned.
%
% THE RIGHTMOST NETWORK MODE (added for the stage, FD-040) is the rightmost
% eigenvalue whose eigenvector has at least x_thr (0.1) of its norm^2 in
% the dendritic rows: the rule of ospec_full_spectrum, here over the WHOLE
% spectrum, so it is found also in a stable network, where it lies to the
% left of hundreds of single-neuron slow eigenvalues. Inside that band the
% two kinds mix and the answer depends on the threshold: at the first
% calibration state of the exploration's E-to-E x 1.5 run with multiple
% timescales it is -0.30 1/s at 0.1 (eigenvalue number 679 from the
% right), -0.18 at 0.02 and -1.32 at 0.3; with one timescale -0.65, -0.64
% and -1.2. Where the rightmost eigenvalue is positive the three agree.
% lam and xfrac are saved whole, so another threshold needs no new run.
%
% Fields
%   lam       N x 1 complex single: all eigenvalues, by real part, descending
%   ov        N x 1 single: overlap |v_1' u| of each eigenvector's normalized
%             x part u with v_1
%   xfrac     N x 1 single: norm^2 fraction of each eigenvector in the x rows
%   n_pos     eigenvalues with positive real part
%   abscissa  rightmost real part
%   cont      the continuation (complex, Im >= 0), cont_ov, cont_xfrac
%   net       the rightmost network mode (complex, Im >= 0; NaN if no
%             eigenvector reaches x_thr), net_rank its number from the right
%   kin_edge  the kinetic edge (ospec_kinetic_edge), 1/s
%   lam1_M    the rightmost eigenvalue of M at this state (units of W)
%   rate1_M   its rate (lam1_M - 1)/tau_d
%   seconds
% With ref (ospec_mode_reference; added 2026-11-20), the same for the
% INTENDED mode, which need not be the rightmost: the eigenvector of M
% that overlaps most with the structural right eigenvector ref.r0
% (ospec_snapshot's rule), and its continuation in the full
% Jacobian by the rule above:
%   mode_lam_M   its eigenvalue in M (units of W, Im >= 0), mode_ov_M its
%                overlap with ref.r0, mode_is_rightmost
%   mode_cont    the continuation (1/s, Im >= 0), mode_cont_ov,
%                mode_cont_xfrac
%
% Ported from scripts/explorations/outliers/outlier_dense_spectrum.m
% (FD-040); net, net_rank and kin_edge are new, the rest is unchanged.
%
% See also: ospec_full_spectrum, ospec_kinetic_edge, ospec_effective_matrix, ospec_run

t0 = tic;
S = S(:);
[Mf, ~, J] = ospec_effective_matrix(S, params, ws);
[Vm, Dm] = eig(Mf);
lm = diag(Dm);
[~, i1] = max(real(lm));
v1 = Vm(:, i1) ./ norm(Vm(:, i1));

[V, D] = eig(full(J));
lam = diag(D);
[~, order] = sort(real(lam), 'descend');
lam = lam(order);
Vx = V(ws.is_x, order);
V = [];                                        %#ok<NASGU>  256 MB at N = 4000
xn2 = sum(abs(Vx).^2, 1);                      % eig returns unit columns
ov = abs(v1' * Vx) ./ sqrt(max(xn2, realmin));

cand = find(xn2 >= 0.1);
[~, q] = max(ov(cand));
q = cand(q);
cont = complex(real(lam(q)), abs(imag(lam(q))));

if nargin < 5 || isempty(x_thr); x_thr = 0.1; end
net = complex(NaN, NaN);
net_rank = find(xn2 >= x_thr, 1, 'first');     % lam is sorted by real part
if isempty(net_rank)
    net_rank = NaN;
else
    net = complex(real(lam(net_rank)), abs(imag(lam(net_rank))));
end

f = struct('lam', single(lam), 'ov', single(ov(:)), 'xfrac', single(xn2(:)), ...
    'n_pos', sum(real(lam) > 0), 'abscissa', real(lam(1)), ...
    'net', net, 'net_rank', net_rank, 'kin_edge', ospec_kinetic_edge(J, ws), ...
    'cont', cont, 'cont_ov', ov(q), 'cont_xfrac', xn2(q), ...
    'lam1_M', lm(i1), 'rate1_M', (real(lm(i1)) - 1) / params.tau_d);

if nargin >= 4 && ~isempty(ref)
    [ov_M, im] = max(abs(ref.r0' * Vm));       % eig returns unit columns
    vm = Vm(:, im);
    ov_mode = abs(vm' * Vx) ./ sqrt(max(xn2, realmin));
    [~, qm] = max(ov_mode(cand));
    qm = cand(qm);
    f.mode_lam_M = complex(real(lm(im)), abs(imag(lm(im))));
    f.mode_ov_M = ov_M;
    f.mode_is_rightmost = abs(real(lm(im)) - real(lm(i1))) < 1e-9 * max(1, abs(lm(i1)));
    f.mode_cont = complex(real(lam(qm)), abs(imag(lam(qm))));
    f.mode_cont_ov = ov_mode(qm);
    f.mode_cont_xfrac = xn2(qm);
end
f.seconds = toc(t0);
end

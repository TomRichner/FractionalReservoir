function z = outlier_with_conjugates(z)
%OUTLIER_WITH_CONJUGATES Eigenvalues with the conjugate of each complex one appended, so that both members of a pair are marked.
%
%   z = outlier_with_conjugates(z)
%
% A real one (|Im| below 1e-9 of the larger of its modulus and 1) is kept
% once, on the real axis. Returns a row.
%
% See also: fig_outliers, fig_outlier_structures
z = z(:).';
is_real = abs(imag(z)) <= 1e-9 * max(abs(z), 1);
z(is_real) = real(z(is_real));
z = [z, conj(z(~is_real))];
end

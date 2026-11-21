function [R, K_eig, rms_gain] = fresp_radius(d, ws)
%FRESP_RADIUS Bulk radius and structural outliers of M(z) = W .* d, predicted from the gains alone.
%
%   [R, K_eig, rms_gain] = fresp_radius(d, ws)
%
% Predicted from the gains (the frozen-gain predictions of the outlier
% stage, for complex gains):
%   R        1 x nz: sqrt of the Perron root of
%            V(a, b) = s2(a, b) sum_{j in b} |d(a, j)|^2
%            (the bulk of a random matrix with a variance profile depends
%            on the moduli of the column gains only, not on their phases)
%   K_eig    C x nz complex: the eigenvalues of the mean part,
%            K(a, b) = mu(a, b) sum_{j in b} d(a, j), sorted by real part,
%            descending. Here the phases matter.
%   rms_gain 1 x nz: sqrt(mean |d|^2) over all neurons and routes
%
%   d   C x n x nz gains (fresp_gain), ws = fresp_weight_stats(params)
%
% See also: fresp_gain, fresp_weight_stats

[C, ~, nz] = size(d);
R = zeros(1, nz);
K_eig = complex(zeros(C, nz));
for q = 1:nz
    V = zeros(C); K = complex(zeros(C));
    for a = 1:C
        for b = 1:C
            g = d(a, ws.ti{b}, q);
            V(a, b) = ws.s2(a, b) * sum(abs(g).^2);
            K(a, b) = ws.mu(a, b) * sum(g);
        end
    end
    R(q) = sqrt(max(0, max(real(eig(V)))));
    e = eig(K);
    [~, order] = sort(real(e), 'descend');
    K_eig(:, q) = e(order);
end
rms_gain = reshape(sqrt(mean(abs(d).^2, [1 2])), 1, nz);
end

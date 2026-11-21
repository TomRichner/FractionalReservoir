function Mz = freqresp_matrix(d, ws)
%FREQRESP_MATRIX M(z) = W .* d: the weights, each column scaled by its neuron's gain on the route to the row's type.
%
%   Mz = freqresp_matrix(d, ws)      d is C x n (one frequency)
%
% See also: freqresp_gain, freqresp_schur
Mz = complex(zeros(size(ws.Wf)));
for a = 1:size(d, 1)
    Mz(ws.ti{a}, :) = ws.Wf(ws.ti{a}, :) .* d(a, :);
end
end

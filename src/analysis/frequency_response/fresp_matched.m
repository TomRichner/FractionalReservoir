function mt = fresp_matched(params, f, r0)
%FRESP_MATCHED Single-neuron responses of every excitatory neuron at a matched operating point.
%
%   mt = fresp_matched(params, f, r0)
%
% Every neuron is put at rate r0 with its slow variables at their steady
% state for that rate (fresp_matched_state; slope 1 on the linear part of
% the nonlinearity), so that conditions are compared at one operating point
% and what is left is the timescale structure. On one network seed a neuron
% has the same fast adaptation constant and the same releasable-pool times
% in the single- and multiple-timescale conditions, so the comparison is
% neuron by neuron.
%
%   params  model.get_params() of the network (any state: none is used)
%   f       1 x nf frequencies in Hz
%   r0      the rate
%
% mt fields
%   r0, f
%   F        n_E x 1: the frozen gain on the E-to-E route (slope 1, so the
%            route's output factor at this rate)
%   G        struct: .both, .sfa, .std, each n_E x nf single complex: the
%            gain on the E-to-E route over the frozen gain, with adaptation
%            and depression free, adaptation alone, depression alone
%   G0       struct of n_E x 1: the same at zero frequency
%
% See also: fresp_matched_state, fresp_gain, fig_freqresp_phase

modes = struct('both', struct('sfa', true, 'std', true), ...
    'sfa', struct('sfa', true, 'std', false), 'std', struct('sfa', false, 'std', true));
names = fieldnames(modes);
z = [2i * pi * f, 0];
iE = params.type_indices{1};
S = fresp_matched_state(params, r0);
[~, d_frozen] = fresp_gain(S, params, 0);
mt = struct('r0', r0, 'f', f, 'F', d_frozen(1, iE)', 'G', struct(), 'G0', struct());
for m = 1:numel(names)
    [d, d_inf] = fresp_gain(S, params, z, modes.(names{m}));
    G = reshape(d(1, iE, :), numel(iE), []) ./ d_inf(1, iE)';
    mt.G.(names{m}) = single(G(:, 1:end - 1));
    mt.G0.(names{m}) = real(G(:, end));
end
end

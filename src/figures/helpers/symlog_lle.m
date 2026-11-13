function [y, tick_pos, tick_lab, lim] = symlog_lle(lambda, lambda0, clip)
% SYMLOG_LLE Symmetric-log transform of lambda_1, for block C's colour scales and axes.
%
%   [y, tick_pos, tick_lab, lim] = SYMLOG_LLE(lambda, lambda0, clip)
%
%   y = sign(lambda) * log10(1 + min(|lambda|, clip) / lambda0)
%
% Linear near zero (|lambda| << lambda0), logarithmic beyond, and the same
% on both sides of zero, so the sign of lambda_1 is readable whether it
% spans -10 to +8 (no adaptation) or -0.1 to +3 (chained MTS, whose stable
% region sits at the slowest adaptation mode, about -0.10). Values beyond
% +-clip are clipped to it; the caller states the clipping.
%
%   lambda0   the linear scale (1/s), default 0.05
%   clip      |lambda| is clipped here (1/s), default 10
%   tick_pos  tick positions, in y: zero and both signs of each decade
%             (0.1, 1, 10, 100) that lies between 1.5 lambda0 and clip
%   tick_lab  their labels, in lambda
%   lim       [-1 1] * y(clip), the axis or colour limits
%
% See also: fig_mc_gain_bias_maps, fig_mc_gain_curves

arguments
    lambda  double
    lambda0 (1,1) double {mustBePositive} = 0.05
    clip    (1,1) double {mustBePositive} = 10
end

DECADES = [0.1 1 10 100];
% A tick closer to zero than this many lambda0 would crowd the zero tick.
MIN_TICK_IN_LAMBDA0 = 1.5;

transform = @(v) sign(v) .* log10(1 + min(abs(v), clip) / lambda0);
y = transform(lambda);
lim = [-1 1] * transform(clip);
% The small factor keeps a decade equal to clip despite round-off.
magnitudes = DECADES(DECADES >= MIN_TICK_IN_LAMBDA0 * lambda0 & DECADES <= clip * 1.001);
values = [-fliplr(magnitudes), 0, magnitudes];
tick_pos = transform(values);
tick_lab = arrayfun(@(v) sprintf('%g', v), values, 'UniformOutput', false);
end

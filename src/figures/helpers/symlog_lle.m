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
% tick_pos / tick_lab: tick positions (in y) and labels (in lambda) at
% 0, +-lambda0 * {1, 2, 5} * 10^k up to clip; lim = [-1 1] * y(clip).
%
% Defaults: lambda0 = 0.05 (1/s), clip = 10 (1/s).
%
% See also: fig_mc_gain_bias_maps, fig_mc_gain_curves

if nargin < 2 || isempty(lambda0); lambda0 = 0.05; end
if nargin < 3 || isempty(clip); clip = 10; end
f = @(v) sign(v) .* log10(1 + min(abs(v), clip) / lambda0);
y = f(lambda);
lim = [-1 1] * f(clip);
mags = [0.1 1 10 100];
mags = mags(mags >= lambda0 * 1.5 & mags <= clip * 1.001);
vals = [-fliplr(mags), 0, mags];
tick_pos = f(vals);
tick_lab = arrayfun(@(v) sprintf('%g', v), vals, 'UniformOutput', false);
end

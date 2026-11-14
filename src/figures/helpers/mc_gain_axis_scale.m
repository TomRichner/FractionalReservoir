function scale = mc_gain_axis_scale(gains)
% MC_GAIN_AXIS_SCALE How a gain sweep's gain axis is drawn: 'linear' for evenly spaced levels, else 'log'.
%
%   scale = MC_GAIN_AXIS_SCALE(D.gains)     % 'log' or 'linear'
%
% A block C grid spaces its gains logarithmically (gain against the bias) or
% evenly (gain against the self-coupling), and every figure draws the gain
% the way the grid was laid out: a linspace grid on a log axis crowds its
% top levels and their tick labels together, and a log grid on a linear axis
% crowds its low ones. The rule is here once, for the maps' x axis, the gain
% curves' x axis, and the gain as an axis or a color in the scatters
% (mc_sweep_quantity), so a color bar does not say "log scale" over a
% linear grid.
%
% Evenly spaced means at least three levels whose steps agree to one part
% in a million. Fewer than three levels do not show a spacing; they are
% drawn on a log axis, as the gain always was.
%
% See also: mc_sweep_quantity, fig_mc_gain_bias_maps, fig_mc_gain_curves

arguments
    gains double
end
RELATIVE_TOLERANCE = 1e-6;
steps = diff(sort(gains(:)));
scale = 'log';
if numel(steps) >= 2 && max(abs(steps - mean(steps))) <= RELATIVE_TOLERANCE * mean(steps)
    scale = 'linear';
end
end

function [spacing, to_position] = mc_axis_spacing(levels)
% MC_AXIS_SPACING How the levels of a block C sweep axis are laid out, and where a map draws a value along that axis.
%
%   [spacing, to_position] = MC_AXIS_SPACING(D.second_levels)
%
%   spacing      'zero_log'  a zero level, then positive levels evenly
%                            spaced in log (mc_gain_sweep_grid's spacing of
%                            that name: a noise amplitude from "off")
%                otherwise   mc_gain_axis_scale's answer: 'log' for positive
%                            levels that are not evenly spaced, 'linear' for
%                            evenly spaced ones or ones that are not all
%                            positive
%   to_position  @(v) where the value v sits along the drawn axis: v itself
%                ('linear'), log10(v) ('log'), or the zero_log rule below
%
% A ZERO_LOG AXIS IS DRAWN AS ITS LOG LADDER, CELL FOR CELL EVENLY SPACED:
% each positive level at log10 of its value, and the zero ONE LOG STEP
% BELOW THE LOWEST POSITIVE LEVEL, where the grid table had the ladder's
% lowest rung before it replaced that rung with 0. The tick labels carry the
% levels' values (the zero labelled 0), and the stretch between 0 and the
% first positive level is a break in the axis, which the maps state in
% their caption. A value between levels is placed by the same rule (a
% positive value at log10 of itself, zero or anything below the zero's place
% at the zero's place), so a preset's star lands where it belongs.
%
% The test: the first level (after sorting) is exactly 0, at least two
% levels are positive, and with three or more positive levels their log
% steps agree to one part in a million. Two positive levels are always a
% log ladder, so three levels {0, a, b} with 0 < a < b are zero_log unless
% they are evenly spaced (b = 2a), which stays linear.
%
% mc_gain_axis_scale keeps the two-way answer ('log' or 'linear') that the
% line and scatter figures use; for a zero_log axis it says 'linear'.
%
% See also: mc_gain_sweep_grid, fig_mc_gain_bias_maps, mc_gain_axis_scale

arguments
    levels double
end
RELATIVE_TOLERANCE = 1e-6;
v = sort(levels(:))';
spacing = mc_gain_axis_scale(v);
if numel(v) >= 3 && v(1) == 0 && all(v(2:end) > 0) && strcmp(spacing, 'linear') && ...
        ~evenly_spaced(v, RELATIVE_TOLERANCE) && ...
        (numel(v) == 3 || evenly_spaced(log(v(2:end)), RELATIVE_TOLERANCE))
    spacing = 'zero_log';
end
switch spacing
    case 'linear'
        to_position = @(x) x;
    case 'log'
        to_position = @log10;
    case 'zero_log'
        lp = log10(v(2:end));
        step = (lp(end) - lp(1)) / (numel(lp) - 1);
        zero_at = lp(1) - step;
        to_position = @(x) max(log10(max(x, realmin)), zero_at);
end
end

function tf = evenly_spaced(v, tol)
% At least three values whose steps agree to tol of their mean.
steps = diff(v);
tf = numel(steps) >= 2 && max(abs(steps - mean(steps))) <= tol * abs(mean(steps));
end

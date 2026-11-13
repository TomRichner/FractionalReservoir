function Q = mc_sweep_quantity(name, D, opts)
% MC_SWEEP_QUANTITY The one registry of quantities of a block C gain sweep: values, labels and units, read off a result.
%
%   Q     = MC_SWEEP_QUANTITY(name, D)              % D from read_mc_gain_sweep
%   Q     = MC_SWEEP_QUANTITY(name, D, 'clip', b)   % clip |value| at b instead
%   names = MC_SWEEP_QUANTITY()                     % every name the registry knows
%
% A quantity is either a column of the result's table or DERIVED from its
% columns, and either way it is defined here once: its values (one per
% table row), how it is labeled and how its axis is drawn. fig_mc_sweep_scatter
% plots any two against each other, colored by a third, and every block C
% figure and table that names the sweep's second axis takes its wording from
% here (Q.label, Q.symbol, Q.variable, Q.unit, Q.plural, Q.at).
%
% COLUMNS: gain, mc_total, horizon, lle_input, lle_noinput, rate_E, rate_I,
% frac_pinned, frac_silent (see read_mc_gain_sweep), and THE SECOND AXIS of
% the result, under its column name (D.second_axis.column):
%   bias           S_c, the center of the per-neuron setpoint draw
%   self_feedback  k, the negative self-feedback on the diagonal of W, in
%                  multiples of R_W, the bulk radius of W at gain 1
%   self_coupling  s = -k, the same diagonal term signed as a connection:
%                  s < 0 negative, s > 0 positive feedback (the model's
%                  alias SRNNCellTypePairs.self_coupling)
% A result has one of them; asking for the other is an error. A second axis
% on a property with no entry here is still a quantity, labeled by its
% property name: give it an entry (the local function second_axis_wording)
% when it should read better than that.
% DERIVED:
%   rate_mean      population mean rate, f_E rate_E + f_I rate_I, with the
%                  network's E and I fractions from mc_sweep_population_fractions
%                  (the run's settings or its preset; never assumed)
%   frac_on_slope  1 - frac_pinned - frac_silent: the share of neurons whose
%                  time-mean rate is neither pinned (> 0.9) nor silent
%                  (< 0.01), a proxy for the share on the sloped part of
%                  the nonlinearity (section 10.7 of Equations_chained_MTS.md)
%
% Q has the fields
%   name        the quantity's name in the registry
%   values      numeric column, one per row of D.T (NaN where unknown)
%   symbol      short name for titles ('\lambda_1', 'Memory capacity')
%   label       axis label, with units where it has any
%   variable    the quantity's mathematical symbol where it has one ('S_c',
%               'k'), else ''
%   unit        its unit as text ('R_W', '1/s'), '' where it has none
%   note        a short gloss of how to read the values, for a caption or
%               a legend ('< 0 negative, > 0 positive feedback'); '' where
%               there is none
%   plural      what several levels of it are called in a caption ('biases')
%   at          @(v) text stating one value ('bias S_c = 0.2')
%   scale       'linear' or 'log' (gain is log, as in every block C figure)
%   clip        [] or a positive bound: |value| beyond it is drawn at the
%               axis edge as an open marker (lambda_1 at 5 1/s, as the maps);
%               the option 'clip' replaces it, and the axis then spans
%               +-1.05 times the clip
%   zero_split  true when the sign of the quantity splits the statistics
%               (lambda_1: stable and unstable sides)
%   zero_line   true when a green dashed line marks zero (lambda_1)
%   levels      the grid's levels when the quantity is a grid axis (gain,
%               the second axis), else []: used for colorbar ticks
%   limits      @(v) [lo hi], the axis limits given the finite values drawn
%
% See also: fig_mc_sweep_scatter, mc_sweep_population_fractions, read_mc_gain_sweep,
%           mc_sweep_second_axis

arguments
    name (1,:) char = ''
    D    struct = struct()
    opts.clip double {mustBeScalarOrEmpty, mustBePositive} = []
end

LLE_CLIP = 5;   % 1/s, as the linear lambda_1 map and fig_mc_vs_lambda1

SECOND_AXES = {'bias', 'self_feedback', 'self_coupling'};   % the second axes with wording of their own
names = [{'gain'}, SECOND_AXES, {'mc_total', 'horizon', 'lle_input', 'lle_noinput', ...
    'rate_E', 'rate_I', 'rate_mean', 'frac_pinned', 'frac_silent', 'frac_on_slope'}];
second_column = '';
if isfield(D, 'second_axis'); second_column = D.second_axis.column; end
if ~isempty(second_column) && ~any(strcmp(second_column, names))
    names = [names, {second_column}];      % a second axis with no wording here
end
if isempty(name)
    Q = names;
    return
end
if ~any(strcmp(name, names))
    error('mc_sweep_quantity:UnknownQuantity', 'Unknown quantity ''%s''. Known: %s.', ...
        name, strjoin(names, ', '));
end
T = D.T;

unit_range = @(v) [-0.03, 1.03];               % a fraction or a rate, bounded by [0, 1]
from_zero  = @(v) [0, 1.05 * max([v(:); eps])];
padded     = @(v) pad_range(v);

% Defaults, overridden per quantity below
Q = struct('name', name, 'values', [], 'symbol', name, 'label', name, 'variable', '', ...
    'unit', '', 'note', '', 'plural', '', 'at', [], 'scale', 'linear', ...
    'clip', [], 'zero_split', false, 'zero_line', false, 'levels', [], 'limits', padded);

is_second_axis = any(strcmp(name, SECOND_AXES)) || strcmp(name, second_column);
if is_second_axis
    if ~any(strcmp(name, T.Properties.VariableNames))
        error('mc_sweep_quantity:NotInResult', ...
            'This result has no column ''%s''; its second axis is ''%s''.', name, second_column);
    end
    [Q.symbol, Q.variable, Q.unit, Q.plural, Q.note] = second_axis_wording(name);
    Q.label = strtrim(sprintf('%s %s', Q.symbol, Q.variable));
    if ~isempty(Q.unit); Q.label = sprintf('%s (%s)', Q.label, Q.unit); end
    Q.values = T.(name);
    if strcmp(name, second_column); Q.levels = D.second_levels; end
end

switch name
    case 'gain'
        Q.values = T.gain;  Q.symbol = 'gain';  Q.label = 'gain';
        Q.scale = 'log';  Q.levels = D.gains;  Q.plural = 'gains';
        Q.limits = @(v) [min(v) / 1.25, max(v) * 1.25];
    case 'mc_total'
        Q.values = T.mc_total;  Q.symbol = 'Memory capacity';  Q.label = 'Memory capacity';
        Q.limits = @(v) [0, ceil(max([v(:); 0]) + 0.5)];   % just above the largest
    case 'horizon'
        Q.values = T.horizon;  Q.symbol = 'memory horizon';  Q.label = 'memory horizon (s)';
        Q.unit = 's';
        Q.limits = from_zero;
    case {'lle_input', 'lle_noinput'}
        Q.values = T.(name);  Q.symbol = '\lambda_1';
        if strcmp(name, 'lle_input')
            Q.label = '\lambda_1 with input (1/s)';
        else
            Q.label = '\lambda_1 without input (1/s)';
        end
        Q.clip = LLE_CLIP;  Q.zero_split = true;  Q.zero_line = true;  Q.unit = '1/s';
        Q.limits = @(v) LLE_CLIP * [-1.05 1.05];
    case {'rate_E', 'rate_I'}
        Q.values = T.(name);  Q.symbol = sprintf('mean rate %s', name(end));
        Q.label = sprintf('mean rate, %s', name(end));  Q.limits = from_zero;
    case 'rate_mean'
        fr = mc_sweep_population_fractions(D);
        Q.values = fr(1) * T.rate_E + fr(2) * T.rate_I;
        Q.symbol = 'mean rate';
        Q.label = sprintf('population mean rate (%.2g E + %.2g I)', fr(1), fr(2));
        Q.limits = from_zero;
    case 'frac_pinned'
        Q.values = T.frac_pinned;  Q.symbol = 'share pinned';
        Q.label = 'share of neurons pinned';  Q.limits = unit_range;
    case 'frac_silent'
        Q.values = T.frac_silent;  Q.symbol = 'share silent';
        Q.label = 'share of neurons silent';  Q.limits = unit_range;
    case 'frac_on_slope'
        Q.values = 1 - T.frac_pinned - T.frac_silent;  Q.symbol = 'share on the slope';
        Q.label = 'share neither pinned nor silent';  Q.limits = unit_range;
end
Q.values = double(Q.values(:));
if isempty(Q.plural); Q.plural = sprintf('%s levels', Q.symbol); end
named = strtrim(sprintf('%s %s', Q.symbol, Q.variable));
unit = Q.unit;
Q.at = @(v) strtrim(sprintf('%s = %.3g %s', named, v, unit));
if ~isempty(opts.clip)
    % A clipped axis spans the clip on both sides, as lambda_1's does.
    Q.clip = opts.clip;
    Q.limits = @(v) opts.clip * [-1.05 1.05];
end
end

%% ------------------------------------------------------------------------
function [symbol, variable, unit, plural, note] = second_axis_wording(column)
% How a second axis of the sweep reads: its name in words, its mathematical
% symbol, its unit, what several of its levels are called, and a gloss of
% how to read its sign where that is not obvious. THE ONE PLACE a new
% second axis gets its wording; without an entry it goes by its column
% name (underscores escaped for the tex interpreter).
note = '';
switch column
    case 'bias'
        symbol = 'bias';  variable = 'S_c';  unit = '';  plural = 'biases';
    case 'self_feedback'
        symbol = 'self-feedback';  variable = 'k';  unit = 'R_W';
        plural = 'self-feedback levels';
        note = 'k > 0 negative feedback';
    case 'self_coupling'
        symbol = 'self-coupling';  variable = 's';  unit = 'R_W';
        plural = 'self-coupling levels';
        note = '< 0 negative, > 0 positive feedback';
    otherwise
        symbol = strrep(column, '_', '\_');  variable = '';  unit = '';
        plural = sprintf('%s levels', symbol);
end
end

function lim = pad_range(v)
% The data range, padded by 5 % of its span (or by 0.5 for a single value).
lo = min(v); hi = max(v);
if isempty(lo) || ~isfinite(lo); lim = [0 1]; return; end
pad = 0.05 * (hi - lo);
if pad == 0; pad = 0.5; end
lim = [lo - pad, hi + pad];
end

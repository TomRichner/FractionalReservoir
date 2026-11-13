function Q = mc_sweep_quantity(name, D, opts)
% MC_SWEEP_QUANTITY The one registry of quantities a block C scatter can plot, read off a gain-by-bias result.
%
%   Q     = MC_SWEEP_QUANTITY(name, D)              % D from read_mc_gain_sweep
%   Q     = MC_SWEEP_QUANTITY(name, D, 'clip', b)   % clip |value| at b instead
%   names = MC_SWEEP_QUANTITY()                     % every name the registry knows
%
% A quantity is either a column of the result's table or DERIVED from its
% columns, and either way it is defined here once: its values (one per
% table row), how it is labeled and how its axis is drawn. fig_mc_sweep_scatter
% plots any two against each other, colored by a third.
%
% COLUMNS: gain, bias, mc_total, horizon, lle_input, lle_noinput, rate_E,
% rate_I, frac_pinned, frac_silent (see read_mc_gain_sweep).
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
%   scale       'linear' or 'log' (gain is log, as in every block C figure)
%   clip        [] or a positive bound: |value| beyond it is drawn at the
%               axis edge as an open marker (lambda_1 at 5 1/s, as the maps);
%               the option 'clip' replaces it, and the axis then spans
%               +-1.05 times the clip
%   zero_split  true when the sign of the quantity splits the statistics
%               (lambda_1: stable and unstable sides)
%   zero_line   true when a green dashed line marks zero (lambda_1)
%   levels      the grid's levels when the quantity is a grid axis (gain,
%               bias), else []: used for colorbar ticks
%   limits      @(v) [lo hi], the axis limits given the finite values drawn
%
% See also: fig_mc_sweep_scatter, mc_sweep_population_fractions, read_mc_gain_sweep

arguments
    name (1,:) char = ''
    D    struct = struct()
    opts.clip double {mustBeScalarOrEmpty, mustBePositive} = []
end

LLE_CLIP = 5;   % 1/s, as the linear lambda_1 map and fig_mc_vs_lambda1

names = {'gain', 'bias', 'mc_total', 'horizon', 'lle_input', 'lle_noinput', ...
    'rate_E', 'rate_I', 'rate_mean', 'frac_pinned', 'frac_silent', 'frac_on_slope'};
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
Q = struct('name', name, 'values', [], 'symbol', name, 'label', name, 'scale', 'linear', ...
    'clip', [], 'zero_split', false, 'zero_line', false, 'levels', [], 'limits', padded);

switch name
    case 'gain'
        Q.values = T.gain;  Q.symbol = 'gain';  Q.label = 'gain';
        Q.scale = 'log';  Q.levels = D.gains;
        Q.limits = @(v) [min(v) / 1.25, max(v) * 1.25];
    case 'bias'
        Q.values = T.bias;  Q.symbol = 'bias';  Q.label = 'bias S_c';
        Q.levels = D.biases;
    case 'mc_total'
        Q.values = T.mc_total;  Q.symbol = 'Memory capacity';  Q.label = 'Memory capacity';
        Q.limits = @(v) [0, ceil(max([v(:); 0]) + 0.5)];   % just above the largest
    case 'horizon'
        Q.values = T.horizon;  Q.symbol = 'memory horizon';  Q.label = 'memory horizon (s)';
        Q.limits = from_zero;
    case {'lle_input', 'lle_noinput'}
        Q.values = T.(name);  Q.symbol = '\lambda_1';
        if strcmp(name, 'lle_input')
            Q.label = '\lambda_1 with input (1/s)';
        else
            Q.label = '\lambda_1 without input (1/s)';
        end
        Q.clip = LLE_CLIP;  Q.zero_split = true;  Q.zero_line = true;
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
if ~isempty(opts.clip)
    % A clipped axis spans the clip on both sides, as lambda_1's does.
    Q.clip = opts.clip;
    Q.limits = @(v) opts.clip * [-1.05 1.05];
end
end

%% ------------------------------------------------------------------------
function lim = pad_range(v)
% The data range, padded by 5 % of its span (or by 0.5 for a single value).
lo = min(v); hi = max(v);
if isempty(lo) || ~isfinite(lo); lim = [0 1]; return; end
pad = 0.05 * (hi - lo);
if pad == 0; pad = 0.5; end
lim = [lo - pad, hi + pad];
end

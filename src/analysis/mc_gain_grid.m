function [G, n_ok] = mc_gain_grid(D, condition, field, stat)
% MC_GAIN_GRID One quantity of the block C table on the bias-by-gain grid.
%
%   [G, n_ok] = MC_GAIN_GRID(D, condition, field)          % median over seeds
%   [G, n_ok] = MC_GAIN_GRID(D, condition, field, stat)    % 'median' | 'iqr' |
%                                                          % 'p25' | 'p75' | function handle
%
%   D          from read_mc_gain_sweep
%   G          numel(D.biases) x numel(D.gains): row = bias level, column =
%              gain level. A cell with no successful job is NaN.
%   n_ok       the number of successful seeds behind each cell
%
% Failed jobs are skipped (read_mc_gain_sweep sets their values to NaN).
%
% See also: read_mc_gain_sweep, mc_gain_bias_statistics

if nargin < 4 || isempty(stat); stat = 'median'; end
if ischar(stat)
    switch stat
        case 'median', f = @(v) median(v);
        case 'iqr',    f = @(v) prctile(v, 75) - prctile(v, 25);
        case 'p25',    f = @(v) prctile(v, 25);
        case 'p75',    f = @(v) prctile(v, 75);
        otherwise
            error('mc_gain_grid:BadStat', 'Unknown statistic ''%s''.', stat);
    end
else
    f = stat;
end

T = D.T(strcmp(D.T.condition, condition), :);
nb = numel(D.biases); ng = numel(D.gains);
G = nan(nb, ng);
n_ok = zeros(nb, ng);
[~, ib] = ismember(T.bias, D.biases);
[~, ig] = ismember(T.gain, D.gains);
v = T.(field);
for b = 1:nb
    for g = 1:ng
        x = v(ib == b & ig == g);
        x = x(isfinite(x));
        n_ok(b, g) = numel(x);
        if ~isempty(x); G(b, g) = f(x); end
    end
end
end

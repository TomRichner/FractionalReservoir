function cmap = excitatory_colormap(n_colors)
% excitatory_colormap - Custom colormap for excitatory neurons
%
% Syntax:
%   cmap = excitatory_colormap()
%   cmap = excitatory_colormap(n_colors)
%
% Description:
%   Returns a colormap with warm reds, ambers, and oranges suitable
%   for plotting excitatory neurons. The colormap can be used with
%   MATLAB's colororder() function or set() for line plots.
%
% Inputs:
%   n_colors - (optional) Number of discrete colors to return (default: 8)
%
% Outputs:
%   cmap - n_colors x 3 RGB matrix with values in [0, 1]
%
% Example:
%   % Set colororder for current axes
%   set(gca, 'ColorOrder', inhibitory_colormap());
%
%   % Get 12 colors instead of default 8
%   cmap = inhibitory_colormap(12);

% Default to 8 colors if not specified
if nargin < 1
    n_colors = 8;
end

% Base palette of 8 carefully selected warm colors
% Shuffled order for maximum visual distinction between adjacent colors
base_palette = [
    1.00, 0.00, 0.00;  % Pure red
    1.00, 0.75, 0.00;  % Golden yellow
    0.85, 0.20, 0.45;  % Rose/magenta
    0.90, 0.10, 0.60;  % Bright magenta-pink
    0.90, 0.55, 0.00;  % Amber/gold
    0.55, 0.27, 0.27;  % Grayish red (desaturated)
    0.86, 0.08, 0.24;  % Crimson
    0.60, 0.15, 0.45;  % Dark rose
    ];

n_base = size(base_palette, 1);

if n_colors == n_base
    % Return base palette directly
    cmap = base_palette;
elseif n_colors < n_base
    % Sample evenly from base palette
    indices = round(linspace(1, n_base, n_colors));
    cmap = base_palette(indices, :);
else
    % Interpolate to get more colors
    % Create interpolation points
    x_base = linspace(1, n_colors, n_base);
    x_new = 1:n_colors;

    % Interpolate each RGB channel
    cmap = zeros(n_colors, 3);
    for i = 1:3
        cmap(:, i) = interp1(x_base, base_palette(:, i), x_new, 'pchip');
    end

    % Clamp values to [0, 1] range
    cmap = max(0, min(1, cmap));
end
end


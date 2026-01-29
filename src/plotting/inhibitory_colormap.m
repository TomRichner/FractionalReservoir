function cmap = inhibitory_colormap(n_colors)
% inhibitory_colormap - Custom colormap for inhibitory neurons
%
% Syntax:
%   cmap = inhibitory_colormap()
%   cmap = inhibitory_colormap(n_colors)
%
% Description:
%   Returns a colormap with discrete blues and cyans with varying saturation
%   suitable for plotting excitatory neurons. The colormap can be used with
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
%   set(gca, 'ColorOrder', excitatory_colormap());
%
%   % Get 12 colors instead of default 8
%   cmap = excitatory_colormap(12);

% Default to 8 colors if not specified
if nargin < 1
    n_colors = 8;
end

% Base palette of 8 carefully selected blues and cyans with varying saturation
% Includes saturated colors and desaturated (grayish) variants
% Ordered to provide good visual distinction
base_palette = [
    0.00, 0.45, 0.74;  % Deep blue (saturated)
    0.00, 0.75, 1.00;  % Sky blue / bright cyan (saturated)
    0.20, 0.47, 0.62;  % Desaturated deep blue (grayish)
    0.00, 0.50, 0.50;  % Teal (saturated)
    0.30, 0.75, 0.93;  % Light cyan (saturated)
    0.25, 0.62, 0.75;  % Desaturated light cyan (grayish)
    0.00, 0.80, 0.80;  % Turquoise (saturated)
    0.15, 0.55, 0.65;  % Desaturated teal (grayish)
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


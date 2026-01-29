function cmap = blue_gray_red_colormap(n)
% BLUE_GRAY_RED_COLORMAP Create a diverging colormap from blue to red with gray in the middle
%
% Usage:
%   cmap = blue_gray_red_colormap()       % Returns 256-level colormap
%   cmap = blue_gray_red_colormap(n)      % Returns n-level colormap
%   colormap(blue_gray_red_colormap(64))  % Apply to current axes
%
% Description:
%   Creates a perceptually smooth colormap transitioning from deep blue
%   through a magenta-gray midpoint to rich red. Uses pchip interpolation
%   for organic transitions.
%
% Input:
%   n - Number of colormap levels (default: 256)
%
% Output:
%   cmap - n x 3 matrix of RGB values in [0, 1]
%
% See also: parula, turbo, coolwarm

arguments
    n (1,1) double {mustBePositive, mustBeInteger} = 256
end

% Define anchor colors
% Blue: Deep but not black
% Mid: Magenta-Grey (elevated R and B for "overlap" feel)
% Red: Rich but not blown out
blue_end    = [0.10, 0.30, 0.90];
magenta_mid = [0.55, 0.45, 0.60];
red_end     = [0.90, 0.20, 0.15];

% Create control point array
nodes = [blue_end; magenta_mid; red_end];
x = [1, n/2, n];  % Map nodes to indices

% Use pchip (Piecewise Cubic Hermite Interpolating Polynomial)
% for a smoother, more "organic" transition than linear interp
map_idx = 1:n;
r = pchip(x, nodes(:,1), map_idx);
g = pchip(x, nodes(:,2), map_idx);
b = pchip(x, nodes(:,3), map_idx);

cmap = [r', g', b'];

% Constrain to [0,1] to avoid clipping artifacts
cmap = max(0, min(1, cmap));

end

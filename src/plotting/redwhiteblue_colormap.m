function cmap = redwhiteblue_colormap(n_colors)
% redwhiteblue_colormap - Diverging colormap from blue through white to red
%
% Syntax:
%   cmap = redwhiteblue_colormap()
%   cmap = redwhiteblue_colormap(n_colors)
%
% Description:
%   Returns a diverging colormap that transitions from blue (negative values)
%   through white (zero) to red (positive values). Useful for visualizing
%   weight matrices where positive (excitatory) connections are red and
%   negative (inhibitory) connections are blue.
%
% Inputs:
%   n_colors - (optional) Number of colors in the colormap (default: 256)
%              Must be an even number for symmetry around zero
%
% Outputs:
%   cmap - n_colors x 3 RGB matrix with values in [0, 1]
%
% Example:
%   % Apply to current figure
%   colormap(redwhiteblue_colormap(256));
%
%   % Use with imagesc
%   imagesc(W);
%   colormap(redwhiteblue_colormap());
%   colorbar;
%
% See also: bluewhitered_colormap

% Default to 256 colors if not specified
if nargin < 1
    n_colors = 256;
end

% Ensure even number of colors for symmetry
if mod(n_colors, 2) ~= 0
    n_colors = n_colors + 1;
end

half = n_colors / 2;

% Blue to white for negative values
% RGB: (0,0,1) -> (1,1,1)
blues = [linspace(0, 1, half)', linspace(0, 1, half)', linspace(1, 1, half)'];

% White to red for positive values
% RGB: (1,1,1) -> (1,0,0)
reds = [linspace(1, 1, half)', linspace(1, 0, half)', linspace(1, 0, half)'];

% Concatenate blue and red halves
cmap = [blues; reds];

end

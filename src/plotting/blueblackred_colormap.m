function cmap = blueblackred_colormap(n_colors)
% blueblackred_colormap - Diverging colormap from blue through black to red
%
% Syntax:
%   cmap = blueblackred_colormap()
%   cmap = blueblackred_colormap(n_colors)
%
% Description:
%   Returns a diverging colormap that transitions from red (negative values)
%   through black (zero) to blue (positive values). Useful for visualizing
%   matrices with positive and negative values, such as Jacobian matrices.
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
%   colormap(blueblackred_colormap(256));
%   
%   % Use with imagesc
%   imagesc(J_eff);
%   colormap(blueblackred_colormap());
%   colorbar;
%
% See also: blueblackred_whitecenter_colormap

    % Default to 256 colors if not specified
    if nargin < 1
        n_colors = 256;
    end
    
    % Ensure even number of colors for symmetry
    if mod(n_colors, 2) ~= 0
        n_colors = n_colors + 1;
    end
    
    half = n_colors / 2;
    
    % Red to black for negative values
    % RGB: (1,0,0) -> (0,0,0)
    reds = [linspace(1, 0, half)', linspace(0, 0, half)', linspace(0, 0, half)'];
    
    % Black to blue for positive values
    % RGB: (0,0,0) -> (0,0,1)
    blues = [linspace(0, 0, half)', linspace(0, 0, half)', linspace(0, 1, half)'];
    
    % Concatenate red and blue halves
    cmap = [reds; blues];
    
end


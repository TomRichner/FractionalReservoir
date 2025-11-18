function cmap = bluewhitered_colormap(n_colors)
% bluewhitered_colormap - Diverging colormap from blue through white to red
%
% Syntax:
%   cmap = bluewhitered_colormap()
%   cmap = bluewhitered_colormap(n_colors)
%
% Description:
%   Returns a diverging colormap that transitions from red (negative values)
%   through white (zero) to blue (positive values). Useful for visualizing
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
%   colormap(bluewhitered_colormap(256));
%   
%   % Use with imagesc
%   imagesc(J_eff);
%   colormap(bluewhitered_colormap());
%   colorbar;
%
% See also: blueblackred_colormap, blueblackred_whitecenter_colormap

    % Default to 256 colors if not specified
    if nargin < 1
        n_colors = 256;
    end
    
    % Ensure even number of colors for symmetry
    if mod(n_colors, 2) ~= 0
        n_colors = n_colors + 1;
    end
    
    half = n_colors / 2;
    
    % Red to white for negative values
    % RGB: (1,0,0) -> (1,1,1)
    reds = [linspace(1, 1, half)', linspace(0, 1, half)', linspace(0, 1, half)'];
    
    % White to blue for positive values
    % RGB: (1,1,1) -> (0,0,1)
    blues = [linspace(1, 0, half)', linspace(1, 0, half)', linspace(1, 1, half)'];
    
    % Concatenate red and blue halves
    cmap = [reds; blues];
    
end


function cmap = bluewhitered_blackcenter_colormap(n_colors, black_line_width)
% bluewhitered_blackcenter_colormap - Diverging colormap with black center line
%
% Syntax:
%   cmap = bluewhitered_blackcenter_colormap()
%   cmap = bluewhitered_blackcenter_colormap(n_colors)
%   cmap = bluewhitered_blackcenter_colormap(n_colors, black_line_width)
%
% Description:
%   Returns a diverging colormap that transitions from red (negative values)
%   through white (zero) to blue (positive values), with a thin black line
%   at the center to emphasize the zero point. Useful for visualizing
%   matrices where the zero point is particularly important.
%
% Inputs:
%   n_colors         - (optional) Number of colors in the colormap (default: 256)
%                      Must be an even number for symmetry around zero
%   black_line_width - (optional) Width of black line in pixels (default: 2)
%
% Outputs:
%   cmap - n_colors x 3 RGB matrix with values in [0, 1]
%
% Example:
%   % Apply to current figure with default 2-pixel black line
%   colormap(bluewhitered_blackcenter_colormap(256));
%   
%   % Use wider black line (4 pixels)
%   colormap(bluewhitered_blackcenter_colormap(256, 4));
%   
%   % Use with imagesc and colorbar
%   imagesc(J_eff);
%   colormap(bluewhitered_blackcenter_colormap());
%   colorbar;
%
% See also: bluewhitered_colormap, blueblackred_whitecenter_colormap

    % Default to 256 colors if not specified
    if nargin < 1
        n_colors = 256;
    end
    
    % Default black line width is 2 pixels
    if nargin < 2
        black_line_width = 2;
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
    
    % Replace center pixels with black line
    center_idx = round(size(cmap, 1) / 2);
    black_indices = center_idx + (-floor(black_line_width/2):floor(black_line_width/2));
    
    % Clamp indices to valid range
    black_indices = black_indices(black_indices >= 1 & black_indices <= n_colors);
    
    % Set to black (0, 0, 0)
    cmap(black_indices, :) = repmat([0 0 0], length(black_indices), 1);
    
end


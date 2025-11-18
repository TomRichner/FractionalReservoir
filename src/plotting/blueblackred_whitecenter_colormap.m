function cmap = blueblackred_whitecenter_colormap(n_colors, white_line_width)
% blueblackred_whitecenter_colormap - Diverging colormap with white center line
%
% Syntax:
%   cmap = blueblackred_whitecenter_colormap()
%   cmap = blueblackred_whitecenter_colormap(n_colors)
%   cmap = blueblackred_whitecenter_colormap(n_colors, white_line_width)
%
% Description:
%   Returns a diverging colormap that transitions from red (negative values)
%   through black (zero) to blue (positive values), with a thin white line
%   at the center to emphasize the zero point. Useful for visualizing
%   matrices where the zero point is particularly important.
%
% Inputs:
%   n_colors         - (optional) Number of colors in the colormap (default: 256)
%                      Must be an even number for symmetry around zero
%   white_line_width - (optional) Width of white line in pixels (default: 2)
%
% Outputs:
%   cmap - n_colors x 3 RGB matrix with values in [0, 1]
%
% Example:
%   % Apply to current figure with default 2-pixel white line
%   colormap(blueblackred_whitecenter_colormap(256));
%   
%   % Use wider white line (4 pixels)
%   colormap(blueblackred_whitecenter_colormap(256, 4));
%   
%   % Use with imagesc and colorbar
%   imagesc(J_eff);
%   colormap(blueblackred_whitecenter_colormap());
%   colorbar;
%
% See also: blueblackred_colormap

    % Default to 256 colors if not specified
    if nargin < 1
        n_colors = 256;
    end
    
    % Default white line width is 2 pixels
    if nargin < 2
        white_line_width = 2;
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
    
    % Replace center pixels with white line
    center_idx = round(size(cmap, 1) / 2);
    white_indices = center_idx + (-floor(white_line_width/2):floor(white_line_width/2));
    
    % Clamp indices to valid range
    white_indices = white_indices(white_indices >= 1 & white_indices <= n_colors);
    
    % Set to white (1, 1, 1)
    cmap(white_indices, :) = repmat([1 1 1], length(white_indices), 1);
    
end


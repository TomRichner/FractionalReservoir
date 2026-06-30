function h = confplot(X, Y, STD1, STD2, C)
% confplot - Plot a mean line with a shaded confidence/std band
%
% Syntax:
%   confplot(X, Y, STD)
%   confplot(X, Y, STD1, STD2)
%   confplot(X, Y, STD1, STD2, C)
%   h = confplot(...)
%
% Description:
%   Plots the mean curve Y(X) as a line and shades the band
%   [Y-STD2, Y+STD1] beneath it as a filled patch (no edges). The band can be
%   symmetric (pass a single STD) or asymmetric (pass STD1 for the upper offset
%   and STD2 for the lower). Draws on the current axes; leaves hold state off.
%
% Inputs:
%   X     - x values (vector)
%   Y     - mean values (vector, same length as X)
%   STD1  - upper band offset added to Y (vector, same length as X)
%   STD2  - lower band offset subtracted from Y (optional; default STD2 = STD1)
%   C     - 2x3 RGB matrix: row 1 = mean line color, row 2 = band fill color
%           (optional; default = blue line, light-blue fill)
%
% Output:
%   h - handle to the mean line
%
% Example:
%   confplot(1:nL, mean_lle, std_lle, std_lle, [0 0 0.8; 0.7 0.8 1]);

    if nargin < 4 || isempty(STD2)
        STD2 = STD1;
    end
    if nargin < 5 || isempty(C)
        C = [0 0 0.8; 0.7 0.8 1.0];   % mean line color; band fill color
    end

    % Row vectors for the patch construction
    X    = X(:)';
    Y    = Y(:)';
    STD1 = STD1(:)';
    STD2 = STD2(:)';

    % Closed polygon: upper edge left->right, then lower edge right->left
    p_X = [X,        fliplr(X)];
    p_Y = [Y + STD1, fliplr(Y) - fliplr(STD2)];

    was_held = ishold;

    f = fill(p_X, p_Y, C(2, :));
    set(f, 'EdgeColor', 'none');
    hold on;
    h = plot(X, Y, 'Color', C(1, :), 'LineWidth', 1.5);

    if ~was_held
        hold off;
    end
end

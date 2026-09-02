function y = piecewiseSigmoid(x, a, c)
% PIECEWISESIGMOID A piecewise linear/quadratic sigmoid activation function.
%
%   y = piecewiseSigmoid(x, a, c) computes the sigmoid function with
%   domain (-inf, inf) and range [0, 1].
%
%   Inputs:
%     x - Input array or scalar.
%     a - Fraction of the domain from -0.5 to 0.5 that is linear (slope=1).
%         Must be in the range [0, 1].
%         a=0 gives a purely quadratic sigmoid.
%         a=1 gives a hard sigmoid (piecewise linear).
%     c - Horizontal shift (center) of the sigmoid.
%
%   Output:
%     y - The computed sigmoid output, same size as x.
%

% --- Input Validation ---
if a < 0 || a > 1
    error('Parameter "a" must be between 0 and 1.');
end

% Convert a to half-width of the linear segment
a = a / 2;

% --- Case 1: Hard Sigmoid (a = 0.5) ---
% This is a special case to avoid division by zero when a=0.5.
% The function simplifies to a clipped linear function.
if a == 0.5
    y_linear = (x - c) + 0.5;
    y = min(max(y_linear, 0), 1); % Clip between 0 and 1

% --- Case 2: General Case (0 <= a < 0.5) ---
else
    % Initialize the output array
    y = zeros(size(x));
    
    % Calculate the scaling constant for the quadratic parts
    % k = 1 / (2 * (1 - 2*a))
    % We handle a=0 separately to avoid 1/inf if 1-2a is exactly 0
    % but in practice, the 'if a == 0.5' check already handles this.
    k = 0.5 / (1 - 2*a);
    
    % Define the five breakpoints
    x1 = c + a - 1;  % Start of left quadratic
    x2 = c - a;      % Start of linear segment
    x3 = c + a;      % End of linear segment
    x4 = c + 1 - a;  % End of right quadratic
    
    % Define logical masks for the 5 regions
    mask_left_sat = (x < x1);
    mask_left_quad = (x >= x1) & (x < x2);
    mask_linear = (x >= x2) & (x <= x3);
    mask_right_quad = (x > x3) & (x <= x4);
    mask_right_sat = (x > x4);
    
    % Region 1: y = 0 (Left saturation)
    % y(mask_left_sat) is already 0 from initialization
    
    % Region 2: y = k * (x - x1)^2 (Left quadratic)
    % Note: x(mask_left_quad) is used for efficient vectorization
    if any(mask_left_quad, 'all')
        y(mask_left_quad) = k * (x(mask_left_quad) - x1).^2;
    end
    
    % Region 3: y = (x - c) + 0.5 (Linear segment)
    if any(mask_linear, 'all')
        y(mask_linear) = (x(mask_linear) - c) + 0.5;
    end
    
    % Region 4: y = 1 - k * (x - x4)^2 (Right quadratic)
    if any(mask_right_quad, 'all')
        y(mask_right_quad) = 1 - k * (x(mask_right_quad) - x4).^2;
    end
    
    % Region 5: y = 1 (Right saturation)
    if any(mask_right_sat, 'all')
        y(mask_right_sat) = 1;
    end
end

end
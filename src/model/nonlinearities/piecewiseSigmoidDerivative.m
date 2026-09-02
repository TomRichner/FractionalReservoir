function dy = piecewiseSigmoidDerivative(x, a, c)
% PIECEWISESIGMOIDDERIVATIVE First derivative of the piecewise sigmoid.
%
%   dy = piecewiseSigmoidDerivative(x, a, c) computes the first
%   derivative of the piecewiseSigmoid function.
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
%     dy - The computed derivative, same size as x.
%

% --- Input Validation ---
if a < 0 || a > 1
    error('Parameter "a" must be between 0 and 1.');
end

% Convert a to half-width of the linear segment
a = a / 2;

% Initialize output array with the same type and size as x
dy = zeros(size(x), 'like', x);

% --- Case 1: Hard Sigmoid (a = 0.5) ---
% The derivative is a simple boxcar function (0 or 1).
if a == 0.5
    breakpoint1 = c - 0.5;
    breakpoint2 = c + 0.5;
    
    % Define the central mask where the derivative is 1
    % We use >= and <= to include the breakpoints in the '1' region
    mask_linear = (x >= breakpoint1) & (x <= breakpoint2);
    
    dy(mask_linear) = 1;

% --- Case 2: General Case (0 <= a < 0.5) ---
% The derivative is a continuous "bump".
else
    % Calculate the scaling constant for the quadratic parts
    k = 0.5 / (1 - 2*a);
    
    % Define the four breakpoints
    x1 = c + a - 1;  % Start of left quadratic ramp
    x2 = c - a;      % Start of linear segment (slope=1)
    x3 = c + a;      % End of linear segment (slope=1)
    x4 = c + 1 - a;  % End of right quadratic ramp
    
    % Define logical masks for the 3 active regions
    % (Regions 1 and 5 are already 0)
    mask_left_quad = (x >= x1) & (x < x2);
    mask_linear = (x >= x2) & (x <= x3);
    mask_right_quad = (x > x3) & (x <= x4);
    
    % Region 2: dy = 2*k*(x - x1) (Left ramp)
    if any(mask_left_quad)
        dy(mask_left_quad) = 2 * k * (x(mask_left_quad) - x1);
    end
    
    % Region 3: dy = 1 (Central constant)
    if any(mask_linear)
        dy(mask_linear) = 1;
    end
    
    % Region 4: dy = -2*k*(x - x4) (Right ramp)
    if any(mask_right_quad)
        dy(mask_right_quad) = -2 * k * (x(mask_right_quad) - x4);
    end
end

end
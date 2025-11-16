function y = tanhActivation(x)
% TANHACTIVATION Hyperbolic tangent activation function.
%
%   y = tanhActivation(x) computes the hyperbolic tangent function
%   with domain (-inf, inf) and range (-1, 1).
%   The slope at x=0 is exactly 1.
%
%   Inputs:
%     x - Input array or scalar.
%
%   Output:
%     y - The computed tanh output, same size as x.
%
%   Example:
%     x = -5:0.1:5;
%     y = tanhActivation(x);
%     plot(x, y);
%
%   See also: tanhActivationDerivative, logisticSigmoid

% --- Compute hyperbolic tangent ---
y = tanh(x);

end


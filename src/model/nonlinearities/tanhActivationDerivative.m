function dy = tanhActivationDerivative(x)
% TANHACTIVATIONDERIVATIVE First derivative of the hyperbolic tangent.
%
%   dy = tanhActivationDerivative(x) computes the first
%   derivative of the tanhActivation function.
%   The slope at x=0 is exactly 1.
%
%   Inputs:
%     x - Input array or scalar.
%
%   Output:
%     dy - The computed derivative, same size as x.
%
%   Example:
%     x = -5:0.1:5;
%     dy = tanhActivationDerivative(x);
%     plot(x, dy);
%
%   See also: tanhActivation, logisticSigmoidDerivative

% --- Compute derivative ---
% Formula: dy/dx = 1 - tanh(x)^2
% This is equivalent to: dy/dx = sech(x)^2
dy = 1 - tanh(x).^2;

end


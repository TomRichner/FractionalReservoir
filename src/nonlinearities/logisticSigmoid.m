function y = logisticSigmoid(x)
% LOGISTICSIGMOID Standard logistic sigmoid activation function.
%
%   y = logisticSigmoid(x) computes the logistic sigmoid function
%   with domain (-inf, inf) and range (0, 1).
%   The slope at x=0 is exactly 1 (beta is fixed at 4).
%
%   Inputs:
%     x - Input array or scalar.
%
%   Output:
%     y - The computed sigmoid output, same size as x.
%
%   Example:
%     x = -5:0.1:5;
%     y = logisticSigmoid(x);
%     plot(x, y);
%
%   See also: logisticSigmoidDerivative, piecewiseSigmoid

% --- Compute logistic sigmoid ---
% Formula: y = 1 / (1 + exp(-4 * x))
% beta = 4 is hardcoded to ensure slope at x=0 is 1
y = 1 ./ (1 + exp(-4 * x));

end


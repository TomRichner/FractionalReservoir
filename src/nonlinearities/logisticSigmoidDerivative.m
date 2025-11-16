function dy = logisticSigmoidDerivative(x)
% LOGISTICSIGMOIDDERIVATIVE First derivative of the logistic sigmoid.
%
%   dy = logisticSigmoidDerivative(x) computes the first
%   derivative of the logisticSigmoid function.
%   The slope at x=0 is exactly 1 (beta is fixed at 4).
%
%   Inputs:
%     x - Input array or scalar.
%
%   Output:
%     dy - The computed derivative, same size as x.
%
%   Example:
%     x = -5:0.1:5;
%     dy = logisticSigmoidDerivative(x);
%     plot(x, dy);
%
%   See also: logisticSigmoid, piecewiseSigmoidDerivative

% --- Compute derivative ---
% Formula: dy/dx = 4 * sigmoid(x) * (1 - sigmoid(x))
% beta = 4 is hardcoded to ensure slope at x=0 is 1
% We compute the sigmoid value to reuse it in the derivative calculation
sigmoid_val = 1 ./ (1 + exp(-4 * x));
dy = 4 * sigmoid_val .* (1 - sigmoid_val);

end


function metrics = compute_metrics(predictions, targets)
% compute_metrics: Calculate regression performance metrics
%
% Computes Mean Squared Error (MSE) and Normalized Root Mean Squared Error
% (NRMSE) for time series prediction tasks.
%
% Inputs:
%   predictions - Predicted values (n_samples x n_outputs)
%   targets     - Target/ground truth values (n_samples x n_outputs)
%
% Output:
%   metrics - struct with fields:
%     mse   - Mean Squared Error
%     rmse  - Root Mean Squared Error
%     nrmse - Normalized RMSE (divided by standard deviation of targets)
%
% Example:
%   metrics = compute_metrics(Y_pred, Y_true);
%   fprintf('MSE: %.6f, NRMSE: %.4f\n', metrics.mse, metrics.nrmse);

    % Ensure inputs are the same size
    if ~isequal(size(predictions), size(targets))
        error('compute_metrics:SizeMismatch', ...
              'Predictions and targets must have the same dimensions');
    end
    
    % Compute errors
    errors = predictions - targets;
    
    % Mean Squared Error
    mse = mean(errors(:).^2);
    
    % Root Mean Squared Error
    rmse = sqrt(mse);
    
    % Normalized RMSE (normalized by standard deviation of targets)
    std_targets = std(targets(:));
    
    % Handle edge case where targets have zero variance
    if std_targets == 0
        nrmse = inf;
        warning('compute_metrics:ZeroVariance', ...
                'Target data has zero variance. NRMSE set to inf.');
    else
        nrmse = rmse / std_targets;
    end
    
    % Package results
    metrics = struct();
    metrics.mse = mse;
    metrics.rmse = rmse;
    metrics.nrmse = nrmse;
end


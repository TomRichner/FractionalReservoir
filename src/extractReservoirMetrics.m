function metrics = extractReservoirMetrics(states_r, S_history, params, silence_start_idx)
%EXTRACTRESERVOIRMETRICS Extract comprehensive metrics from reservoir simulation
%
% This function computes various metrics from reservoir state history:
%   - Mean/std firing rates
%   - Participation ratio (dimensionality)
%   - Spectral abscissa (from Jacobian analysis)
%   - Population autocorrelation
%   - E/I balance metrics
%
% Inputs:
%   states_r - Firing rates (n_timesteps x n_neurons)
%   S_history - Full state history (n_timesteps x n_states)
%   params - Parameter struct with n, n_E, n_I, etc.
%   silence_start_idx - Index where silence period begins
%
% Output:
%   metrics - struct with all computed metrics

    metrics = struct();
    
    n = params.n;
    n_E = params.n_E;
    n_I = params.n_I;
    n_timesteps = size(states_r, 1);
    
    % Basic firing rate statistics
    mean_activity = mean(states_r, 1);
    std_activity = std(states_r, 1);
    max_activity = max(states_r, [], 1);
    min_activity = min(states_r, [], 1);
    
    metrics.mean_firing_rate = mean(mean_activity);
    metrics.std_firing_rate = mean(std_activity);
    metrics.max_firing_rate = max(max_activity);
    metrics.min_firing_rate = min(min_activity);
    
    % E/I population statistics
    if n_E > 0
        E_activity = states_r(:, 1:n_E);
        metrics.mean_firing_rate_E = mean(E_activity(:));
        metrics.std_firing_rate_E = std(E_activity(:));
    else
        metrics.mean_firing_rate_E = 0;
        metrics.std_firing_rate_E = 0;
    end
    
    if n_I > 0
        I_activity = states_r(:, n_E+1:end);
        metrics.mean_firing_rate_I = mean(I_activity(:));
        metrics.std_firing_rate_I = std(I_activity(:));
    else
        metrics.mean_firing_rate_I = 0;
        metrics.std_firing_rate_I = 0;
    end
    
    % E/I balance ratio
    if metrics.mean_firing_rate_I > 0
        metrics.EI_balance = metrics.mean_firing_rate_E / metrics.mean_firing_rate_I;
    else
        metrics.EI_balance = Inf;
    end
    
    % Participation ratio (effective dimensionality)
    try
        C = cov(states_r);
        eigvals = eig(C);
        eigvals = eigvals(eigvals > 1e-10);  % Remove numerical zeros
        
        if ~isempty(eigvals)
            participation_ratio = sum(eigvals)^2 / sum(eigvals.^2);
            metrics.participation_ratio = participation_ratio;
            metrics.effective_dimensionality = participation_ratio;
        else
            metrics.participation_ratio = 0;
            metrics.effective_dimensionality = 0;
        end
    catch
        metrics.participation_ratio = NaN;
        metrics.effective_dimensionality = NaN;
    end
    
    % Population autocorrelation (during silence period)
    if silence_start_idx < n_timesteps
        recovery_activity = mean(states_r(silence_start_idx:end, :), 2);
        
        if length(recovery_activity) > 10
            try
                [acf, corr_lags] = autocorr(recovery_activity, ...
                    'NumLags', min(200, length(recovery_activity)-1));
                acf = acf(corr_lags >= 0);
                
                % Autocorrelation time (first zero crossing or decay to 1/e)
                if length(acf) > 1
                    decay_threshold = 1/exp(1);
                    decay_idx = find(acf < decay_threshold, 1, 'first');
                    if ~isempty(decay_idx)
                        metrics.autocorr_time = corr_lags(decay_idx);
                    else
                        metrics.autocorr_time = corr_lags(end);
                    end
                    
                    % First zero crossing
                    zero_crossing = find(acf < 0, 1, 'first');
                    if ~isempty(zero_crossing)
                        metrics.autocorr_zero_crossing = corr_lags(zero_crossing);
                    else
                        metrics.autocorr_zero_crossing = NaN;
                    end
                else
                    metrics.autocorr_time = NaN;
                    metrics.autocorr_zero_crossing = NaN;
                end
            catch
                metrics.autocorr_time = NaN;
                metrics.autocorr_zero_crossing = NaN;
            end
        else
            metrics.autocorr_time = NaN;
            metrics.autocorr_zero_crossing = NaN;
        end
    else
        metrics.autocorr_time = NaN;
        metrics.autocorr_zero_crossing = NaN;
    end
    
    % Spectral abscissa (if Jacobian computation is available)
    % This would require computing Jacobian at multiple timepoints
    % For now, we'll skip this or compute it separately if needed
    metrics.spectral_abscissa = NaN;  % Placeholder
    
    % Activity range during silence
    if silence_start_idx < n_timesteps
        silence_activity = states_r(silence_start_idx:end, :);
        metrics.silence_mean_rate = mean(silence_activity(:));
        metrics.silence_std_rate = std(silence_activity(:));
        metrics.silence_max_rate = max(silence_activity(:));
        metrics.silence_min_rate = min(silence_activity(:));
    else
        metrics.silence_mean_rate = NaN;
        metrics.silence_std_rate = NaN;
        metrics.silence_max_rate = NaN;
        metrics.silence_min_rate = NaN;
    end
    
    % Coefficient of variation (CV) as measure of variability
    if metrics.mean_firing_rate > 0
        metrics.coefficient_of_variation = metrics.std_firing_rate / metrics.mean_firing_rate;
    else
        metrics.coefficient_of_variation = Inf;
    end
end


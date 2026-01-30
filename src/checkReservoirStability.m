function stability = checkReservoirStability(states_r, silence_start_idx, options)
%CHECKRESERVOIRSTABILITY Detect unstable and dead reservoir simulations
%
% This function analyzes reservoir firing rates to detect:
%   - Unstable: exponential growth during silence OR max rate > 100 Hz
%   - Dead: std of activity < 1e-6 (no variation)
%
% Inputs:
%   states_r - Firing rates (n_timesteps x n_neurons)
%   silence_start_idx - Index where silence period begins (after stimulus)
%   options - struct with optional fields:
%     max_rate_threshold - Threshold for max rate (default: 100)
%     growth_rate_threshold - Threshold for exponential growth slope (default: 0.01)
%     dead_std_threshold - Threshold for dead detection (default: 1e-6)
%
% Output:
%   stability - struct with fields:
%     is_unstable - true if unstable
%     is_dead - true if dead
%     reason - string describing the issue ('stable', 'unstable_max_rate', 
%              'unstable_exponential_growth', 'dead_no_variation', 'dead_zero_state')
%     metrics - struct with diagnostic metrics:
%       max_rate - Maximum firing rate during silence
%       growth_rate - Exponential growth rate (slope of log-linear fit)
%       activity_std - Standard deviation of activity during silence
%       mean_rate - Mean firing rate during silence
%       has_nan_inf - true if NaN/Inf detected in states

    % Default options
    if nargin < 3
        options = struct();
    end
    
    max_rate_threshold = getFieldOrDefault(options, 'max_rate_threshold', 100);
    growth_rate_threshold = getFieldOrDefault(options, 'growth_rate_threshold', 0.01);
    dead_std_threshold = getFieldOrDefault(options, 'dead_std_threshold', 1e-6);
    
    % Initialize output
    stability = struct();
    stability.is_unstable = false;
    stability.is_dead = false;
    stability.reason = 'stable';
    stability.metrics = struct();
    
    % Extract silence period (after stimulus ends)
    if silence_start_idx > size(states_r, 1)
        silence_start_idx = size(states_r, 1);
    end
    
    silence_activity = states_r(silence_start_idx:end, :);
    
    % Check for NaN/Inf
    has_nan_inf = any(~isfinite(silence_activity(:)));
    stability.metrics.has_nan_inf = has_nan_inf;
    
    if has_nan_inf
        stability.is_unstable = true;
        stability.reason = 'unstable_nan_inf';
        return;
    end
    
    % Compute basic metrics
    max_rate = max(silence_activity(:));
    mean_rate = mean(silence_activity(:));
    activity_std = std(silence_activity(:));
    
    stability.metrics.max_rate = max_rate;
    stability.metrics.mean_rate = mean_rate;
    stability.metrics.activity_std = activity_std;
    
    % Check for unstable: max rate > threshold
    if max_rate > max_rate_threshold
        stability.is_unstable = true;
        stability.reason = 'unstable_max_rate';
        return;
    end
    
    % Check for unstable: exponential growth
    % Fit log-linear model to detect exponential growth
    if size(silence_activity, 1) > 10  % Need enough points for fitting
        % Use population mean activity for growth detection
        pop_activity = mean(silence_activity, 2);
        
        % Avoid log(0) by adding small offset
        log_activity = log(max(pop_activity, 1e-10));
        
        % Fit linear model: log(activity) = a*t + b
        time_vec = (1:length(log_activity))';
        if length(time_vec) > 1
            p = polyfit(time_vec, log_activity, 1);
            growth_rate = p(1);  % Slope = growth rate
            
            stability.metrics.growth_rate = growth_rate;
            
            % Check if growth rate exceeds threshold (positive = exponential growth)
            if growth_rate > growth_rate_threshold
                stability.is_unstable = true;
                stability.reason = 'unstable_exponential_growth';
                return;
            end
        end
    else
        stability.metrics.growth_rate = NaN;
    end
    
    % Check for dead: std < threshold (no variation)
    if activity_std < dead_std_threshold
        stability.is_dead = true;
        stability.reason = 'dead_no_variation';
        return;
    end
    
    % Check for dead: all states near zero
    if mean_rate < 1e-6 && max_rate < 1e-5
        stability.is_dead = true;
        stability.reason = 'dead_zero_state';
        return;
    end
    
    % If we get here, reservoir is stable and alive
    stability.reason = 'stable';
end

function val = getFieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end


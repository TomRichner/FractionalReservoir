function [t, x] = generate_mackey_glass(varargin)
% generate_mackey_glass: Generate Mackey-Glass time series
%
% The Mackey-Glass equation is a classic chaotic time-delay differential equation:
%   dx/dt = beta * x(t-tau) / (1 + x(t-tau)^n) - gamma * x(t)
%
% This implementation uses the 4th-order Runge-Kutta method with history
% for the delay term.
%
% Usage:
%   [t, x] = generate_mackey_glass()  % Use default parameters
%   [t, x] = generate_mackey_glass('tau', 17, 'n_samples', 5000)
%
% Parameters (Name-Value pairs):
%   tau        - Time delay (default: 17, chaotic regime)
%   beta       - Production rate constant (default: 0.2)
%   gamma      - Decay rate constant (default: 0.1)
%   n          - Nonlinearity exponent (default: 10)
%   dt         - Time step (default: 1.0)
%   n_samples  - Number of samples to generate (default: 10000)
%   x0         - Initial condition (default: 1.2)
%   discard    - Number of initial samples to discard (default: 1000)
%
% Outputs:
%   t - Time vector (n_samples x 1)
%   x - Mackey-Glass time series (n_samples x 1)
%
% Example:
%   [t, x] = generate_mackey_glass('tau', 17, 'n_samples', 5000);
%   plot(t, x);
%   xlabel('Time'); ylabel('x(t)'); title('Mackey-Glass Time Series');

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'tau', 17, @isnumeric);
    addParameter(p, 'beta', 0.2, @isnumeric);
    addParameter(p, 'gamma', 0.1, @isnumeric);
    addParameter(p, 'n', 10, @isnumeric);
    addParameter(p, 'dt', 1.0, @isnumeric);
    addParameter(p, 'n_samples', 10000, @isnumeric);
    addParameter(p, 'x0', 1.2, @isnumeric);
    addParameter(p, 'discard', 1000, @isnumeric);
    parse(p, varargin{:});
    
    tau = p.Results.tau;
    beta = p.Results.beta;
    gamma = p.Results.gamma;
    n_exp = p.Results.n;
    dt = p.Results.dt;
    n_samples = p.Results.n_samples;
    x0 = p.Results.x0;
    discard = p.Results.discard;
    
    % Total samples including discard period
    n_total = n_samples + discard;
    
    % Number of history steps needed for delay
    delay_steps = ceil(tau / dt);
    
    % Initialize history with constant initial value
    history_length = delay_steps + 1;
    x_history = x0 * ones(history_length, 1);
    
    % Allocate output
    x = zeros(n_total, 1);
    x(1) = x0;
    
    % 4th-order Runge-Kutta integration with delay
    for i = 2:n_total
        % Current state
        x_current = x_history(end);
        
        % Delayed state (from tau time units ago)
        x_delayed = x_history(1);
        
        % RK4 steps
        k1 = dt * mackey_glass_derivative(x_current, x_delayed, beta, gamma, n_exp);
        
        x_mid1 = x_current + 0.5 * k1;
        k2 = dt * mackey_glass_derivative(x_mid1, x_delayed, beta, gamma, n_exp);
        
        x_mid2 = x_current + 0.5 * k2;
        k3 = dt * mackey_glass_derivative(x_mid2, x_delayed, beta, gamma, n_exp);
        
        x_next_est = x_current + k3;
        k4 = dt * mackey_glass_derivative(x_next_est, x_delayed, beta, gamma, n_exp);
        
        % Update
        x_next = x_current + (k1 + 2*k2 + 2*k3 + k4) / 6;
        
        % Store
        x(i) = x_next;
        
        % Update history (shift and append)
        x_history = [x_history(2:end); x_next];
    end
    
    % Discard initial transient
    x = x(discard+1:end);
    
    % Create time vector
    t = (0:n_samples-1)' * dt;
end

function dxdt = mackey_glass_derivative(x_current, x_delayed, beta, gamma, n_exp)
    % Compute derivative of Mackey-Glass equation
    % dx/dt = beta * x(t-tau) / (1 + x(t-tau)^n) - gamma * x(t)
    
    production = beta * x_delayed / (1 + x_delayed^n_exp);
    decay = gamma * x_current;
    dxdt = production - decay;
end


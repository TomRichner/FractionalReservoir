%% Comprehensive Test Suite for piecewiseSigmoid and piecewiseSigmoidDerivative
% This script tests the correctness of both the sigmoid function and its derivative

clear; close all; clc;

%% Test 1: Visual Inspection of Function and Derivative
fprintf('=== Test 1: Visual Inspection ===\n');

x = linspace(-2.5, 3.5, 1000);

% Test different parameter combinations
test_params = [
    0.0, 0.0;   % Purely quadratic, centered at 0
    0.9, 0.4;   % Half linear, centered at 0
    0.8, 0.0;   % Mostly linear, centered at 0
    1.0, 0.0;   % Hard sigmoid, centered at 0
    1.0, 1.0;   % Hard sigmoid, shifted right
];

figure();

for i = 1:size(test_params, 1)
    a = test_params(i, 1);
    c = test_params(i, 2);
    
    % Compute function and analytical derivative
    y = piecewiseSigmoid(x, a, c);
    dy_analytical = piecewiseSigmoidDerivative(x, a, c);
    
    % Compute numerical derivative for comparison
    dy_numerical = gradient(y) ./ gradient(x);
    
    % Plot function and derivatives on same axes
    subplot(2, 3, i);
    plot(x, y, 'b-', 'LineWidth', 2, 'DisplayName', 'f(x)');
    hold on;
    plot(x, dy_analytical, 'r-', 'LineWidth', 2, 'DisplayName', 'Analytical dy/dx');
    plot(x, dy_numerical, 'k--', 'LineWidth', 1, 'DisplayName', 'Numerical dy/dx');
    hold off;
    
    xlabel('x');
    ylabel('y');
    ylim([-0.1, 1.5]);
    title(sprintf('a=%.1f, c=%.1f', a, c));
    legend('Location', 'best');
    grid on;
end

sgtitle('Piecewise Sigmoid Function and Derivative Tests');

%% Test 2: Continuity at Breakpoints
fprintf('\n=== Test 2: Continuity at Breakpoints ===\n');

test_cases = [0.0, 0.0; 0.5, 0.0; 0.8, 0.0; 1.0, 1.0];
tolerance = 1e-6;  % Relaxed tolerance for numerical precision
all_continuous = true;

for i = 1:size(test_cases, 1)
    a = test_cases(i, 1);
    c = test_cases(i, 2);
    
    % Calculate breakpoints (need to replicate logic from function)
    a_half = a / 2;
    
    if a_half == 0.5
        breakpoints = [c - 0.5, c + 0.5];
    else
        x1 = c + a_half - 1;
        x2 = c - a_half;
        x3 = c + a_half;
        x4 = c + 1 - a_half;
        breakpoints = [x1, x2, x3, x4];
    end
    
    % Test continuity at each breakpoint
    for bp = breakpoints
        x_left = bp - 1e-8;
        x_right = bp + 1e-8;
        x_mid = bp;
        
        y_left = piecewiseSigmoid(x_left, a, c);
        y_mid = piecewiseSigmoid(x_mid, a, c);
        y_right = piecewiseSigmoid(x_right, a, c);
        
        diff_left = abs(y_mid - y_left);
        diff_right = abs(y_right - y_mid);
        
        if diff_left > tolerance || diff_right > tolerance
            fprintf('  FAIL: a=%.1f, c=%.1f at x=%.4f: left_diff=%.2e, right_diff=%.2e\n', ...
                    a, c, bp, diff_left, diff_right);
            all_continuous = false;
        end
    end
end

if all_continuous
    fprintf('  PASS: All breakpoints are continuous (within tolerance %.2e)\n', tolerance);
end

%% Test 3: Derivative Correctness (Analytical vs Numerical)
fprintf('\n=== Test 3: Derivative Accuracy ===\n');

x = linspace(-2, 3, 500);
test_cases = [0.0, 0.0; 0.5, 0.0; 0.8, 0.0];  % Exclude a=1.0 (hard sigmoid has discontinuous derivative)
tolerance = 1e-2;  % Numerical derivatives are less accurate
all_derivatives_correct = true;

for i = 1:size(test_cases, 1)
    a = test_cases(i, 1);
    c = test_cases(i, 2);
    
    y = piecewiseSigmoid(x, a, c);
    dy_analytical = piecewiseSigmoidDerivative(x, a, c);
    dy_numerical = gradient(y) ./ gradient(x);
    
    % Compute max absolute error
    max_error = max(abs(dy_analytical - dy_numerical));
    mean_error = mean(abs(dy_analytical - dy_numerical));
    
    if max_error > tolerance
        fprintf('  FAIL: a=%.1f, c=%.1f: max_error=%.4e, mean_error=%.4e\n', ...
                a, c, max_error, mean_error);
        all_derivatives_correct = false;
    else
        fprintf('  PASS: a=%.1f, c=%.1f: max_error=%.4e, mean_error=%.4e\n', ...
                a, c, max_error, mean_error);
    end
end

%% Test 4: Boundary Conditions
fprintf('\n=== Test 4: Boundary Conditions ===\n');

test_cases = [0.0, 0.0; 0.5, 0.0; 1.0, 0.0];
all_boundaries_correct = true;

for i = 1:size(test_cases, 1)
    a = test_cases(i, 1);
    c = test_cases(i, 2);
    
    % Test saturation at extremes
    y_left = piecewiseSigmoid(-1000, a, c);
    y_right = piecewiseSigmoid(1000, a, c);
    
    % Test that y is in [0, 1]
    x_test = linspace(-5, 5, 1000);
    y_test = piecewiseSigmoid(x_test, a, c);
    
    if abs(y_left) > 1e-10
        fprintf('  FAIL: a=%.1f, c=%.1f: Left saturation = %.4e (expected 0)\n', a, c, y_left);
        all_boundaries_correct = false;
    end
    
    if abs(y_right - 1) > 1e-10
        fprintf('  FAIL: a=%.1f, c=%.1f: Right saturation = %.4e (expected 1)\n', a, c, y_right);
        all_boundaries_correct = false;
    end
    
    if any(y_test < -1e-10) || any(y_test > 1 + 1e-10)
        fprintf('  FAIL: a=%.1f, c=%.1f: Output outside [0,1] range\n', a, c);
        all_boundaries_correct = false;
    end
end

if all_boundaries_correct
    fprintf('  PASS: All boundary conditions satisfied\n');
end

%% Test 5: Derivative Continuity (C1 Smoothness)
fprintf('\n=== Test 5: Derivative Continuity (C1 Smoothness) ===\n');

test_cases = [0.0, 0.0; 0.5, 0.0; 0.8, 0.0];
tolerance = 1e-6;  % Relaxed tolerance for numerical precision
all_smooth = true;

for i = 1:size(test_cases, 1)
    a = test_cases(i, 1);
    c = test_cases(i, 2);
    
    % Calculate breakpoints
    a_half = a / 2;
    
    if a_half == 0.5
        % Hard sigmoid: derivative is discontinuous (expected)
        continue;
    else
        x1 = c + a_half - 1;
        x2 = c - a_half;
        x3 = c + a_half;
        x4 = c + 1 - a_half;
        breakpoints = [x1, x2, x3, x4];
    end
    
    % Test derivative continuity at each breakpoint
    for bp = breakpoints
        x_left = bp - 1e-8;
        x_right = bp + 1e-8;
        x_mid = bp;
        
        dy_left = piecewiseSigmoidDerivative(x_left, a, c);
        dy_mid = piecewiseSigmoidDerivative(x_mid, a, c);
        dy_right = piecewiseSigmoidDerivative(x_right, a, c);
        
        diff_left = abs(dy_mid - dy_left);
        diff_right = abs(dy_right - dy_mid);
        
        if diff_left > tolerance || diff_right > tolerance
            fprintf('  FAIL: a=%.1f, c=%.1f at x=%.4f: left_diff=%.2e, right_diff=%.2e\n', ...
                    a, c, bp, diff_left, diff_right);
            all_smooth = false;
        end
    end
end

if all_smooth
    fprintf('  PASS: Derivative is continuous at all breakpoints (C1 smooth)\n');
end

%% Test 6: Symmetry Check (for c=0)
fprintf('\n=== Test 6: Symmetry Check (for c=0) ===\n');

test_cases = [0.0; 0.5; 0.8; 1.0];
tolerance = 1e-6;  % Relaxed tolerance for numerical precision
all_symmetric = true;

for i = 1:length(test_cases)
    a = test_cases(i);
    c = 0.0;
    
    x_test = linspace(-3, 3, 100);
    y_pos = piecewiseSigmoid(x_test, a, c);
    y_neg = piecewiseSigmoid(-x_test, a, c);
    
    % Check if y(x) + y(-x) = 1 (point symmetry about (0, 0.5))
    symmetry_check = abs(y_pos + y_neg - 1);
    max_symmetry_error = max(symmetry_check);
    
    if max_symmetry_error > tolerance
        fprintf('  FAIL: a=%.1f: max symmetry error = %.2e\n', a, max_symmetry_error);
        all_symmetric = false;
    else
        fprintf('  PASS: a=%.1f: symmetric about (c, 0.5)\n', a);
    end
end

%% Test 7: Monotonicity
fprintf('\n=== Test 7: Monotonicity Check ===\n');

test_cases = [0.0, 0.0; 0.5, 0.0; 0.8, 0.0; 1.0, 1.0];
all_monotonic = true;

for i = 1:size(test_cases, 1)
    a = test_cases(i, 1);
    c = test_cases(i, 2);
    
    x_test = linspace(-5, 5, 1000);
    y_test = piecewiseSigmoid(x_test, a, c);
    dy_test = piecewiseSigmoidDerivative(x_test, a, c);
    
    % Check that derivative is non-negative (monotonically increasing)
    if any(dy_test < -1e-10)
        fprintf('  FAIL: a=%.1f, c=%.1f: Negative derivative found\n', a, c);
        all_monotonic = false;
    end
    
    % Check that y is monotonically increasing
    dy_numerical = diff(y_test);
    if any(dy_numerical < -1e-10)
        fprintf('  FAIL: a=%.1f, c=%.1f: Function is not monotonically increasing\n', a, c);
        all_monotonic = false;
    end
end

if all_monotonic
    fprintf('  PASS: All functions are monotonically increasing\n');
end

%% Summary
fprintf('\n=== TEST SUMMARY ===\n');
fprintf('All tests completed. Review results above.\n');
fprintf('Visual inspection plots are displayed in Figure 1.\n');

% Define the input range
x = -2.5:0.01:3.5;

% --- Calculate different sigmoids ---
% 1. General case: a=0.8 (80% linear), centered at c=0
y1 = piecewiseSigmoid(x, 0.8, 0);

% 2. Hard sigmoid: a=1.0 (100% linear), shifted by c=1.0
y2 = piecewiseSigmoid(x, 1.0, 1.0);

% 3. Purely quadratic sigmoid: a=0 (0% linear), at center c=0
y3 = piecewiseSigmoid(x, 0, 0);

% --- Plot the results ---
figure;
plot(x, y1, 'b-', 'LineWidth', 2, 'DisplayName', 'a=0.8, c=0');
hold on;
plot(x, y2, 'r--', 'LineWidth', 2, 'DisplayName', 'a=1.0, c=1.0 (Hard)');
plot(x, y3, 'g-.', 'LineWidth', 2, 'DisplayName', 'a=0, c=0 (Quadratic)');

% --- Add derivative plot for a=0.8 to show slope ---
% We can compute the derivative numerically
dy1_dx = gradient(y1) ./ gradient(x);
plot(x, dy1_dx, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Derivative (a=0.8)');

% --- Formatting ---
grid on;
ylim([-0.1, 1.2]);
xlabel('Input (x)');
ylabel('Output (y)');
title('Piecewise Sigmoid Activation Function');
legend('show', 'Location', 'northwest');
hold off;
axis equal
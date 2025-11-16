% Define the input range
x = -2.5:0.01:3.5;

% --- Calculate different sigmoids ---
% 1. General case: a=0.25, shifted by c=0.5
y1 = piecewiseSigmoid(x, 0.4, 0);

% 2. Hard sigmoid: a=0.5, shifted by c=1.0
y2 = piecewiseSigmoid(x, 0.5, 1.0);

% 3. Purely quadratic sigmoid: a=0, at center c=0
y3 = piecewiseSigmoid(x, 0, 0);

% --- Plot the results ---
figure;
plot(x, y1, 'b-', 'LineWidth', 2, 'DisplayName', 'a=0.4, c=0');
hold on;
plot(x, y2, 'r--', 'LineWidth', 2, 'DisplayName', 'a=0.5, c=1.0 (Hard)');
plot(x, y3, 'g-.', 'LineWidth', 2, 'DisplayName', 'a=0, c=0 (Quadratic)');

% --- Add derivative plot for a=0.25 to show slope ---
% We can compute the derivative numerically
dy1_dx = gradient(y1) ./ gradient(x);
plot(x, dy1_dx, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Derivative (a=0.25)');

% --- Formatting ---
grid on;
ylim([-0.1, 1.2]);
xlabel('Input (x)');
ylabel('Output (y)');
title('Piecewise Sigmoid Activation Function');
legend('show', 'Location', 'northwest');
hold off;
axis equal
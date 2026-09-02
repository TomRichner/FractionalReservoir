% plot_unit_slope_sigmoids.m
%
% Overlay two smooth sigmoids with range (0, 1), midpoint 0.5 at x = 0,
% and unit slope at x = 0.

close all; clear; clc;

x = linspace(-3, 3, 2001);

% The factor 4 gives the logistic sigmoid unit slope at its midpoint.
logistic_sigmoid = 1 ./ (1 + exp(-4 .* x));

% Scaling tanh by 1/2 sets its range to (0, 1); the factor 2 inside tanh
% gives the transformed function unit slope at its midpoint.
scaled_tanh = 0.5 .* (tanh(2 .* x) + 1);

colors = lines(2);
figure('Color', 'w', 'Name', 'Unit-slope sigmoid comparison');
hold on;
plot(x, logistic_sigmoid, ...
    'Color', colors(1, :), 'LineWidth', 2);
plot(x, scaled_tanh, ...
    'Color', colors(2, :), 'LineWidth', 2);
xline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
yline(0.5, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
plot(0, 0.5, 'o', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', [0.25 0.25 0.25], ...
    'MarkerEdgeColor', [0.25 0.25 0.25], ...
    'HandleVisibility', 'off');
hold off;

xlim([-3 3]);
ylim([0 1]);
xlabel('Input, x');
ylabel('Activation');
title('Range-(0,1) sigmoids with unit slope at x = 0');
legend({ ...
    '$[1+\exp(-4x)]^{-1}$', ...
    '$[\tanh(2x)+1]/2$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast');
grid on;
box off;

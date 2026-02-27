% plot_piecewiseSigmoid_example.m
% Plots the piecewiseSigmoid nonlinearity with a=0.9, c=0.4
% and saves the result as a PNG.

%% Parameters
a = 0.9;
c = 0.4;

%% Generate data
x = linspace(-1.5, 1.5, 1000);
y = piecewiseSigmoid(x, a, c);

%% Plot
fig = figure;
plot(x, y, 'k', 'LineWidth', 2);
xlabel('x');
ylabel('y');
box off;
axis equal;
xlim([-1.5, 1.5]);
ylim([-0.1, 1.1]);

%% Save
saveas(fig, fullfile(fileparts(mfilename('fullpath')), ...
    'hard_sigmoid_rounded_corners.png'));

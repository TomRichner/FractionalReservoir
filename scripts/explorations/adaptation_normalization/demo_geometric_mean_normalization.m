% demo_geometric_mean_normalization.m
%
% How should M short-term-depression variables b_m(t), each a fraction in
% [0, 1], be combined into one synaptic gain? The model currently uses their
% PRODUCT. This script compares the product with the GEOMETRIC MEAN
% prod(b)^(1/M) and the ARITHMETIC MEAN mean(b) on three simple sinusoids,
% b_m(t) = 0.5 + 0.5*sin(2*pi*f_m*t). It is a visual aid, not model code.
%
% Points worth noticing:
%   - The product is always the smallest of the three combinations, and it
%     gets deeper as more variables are added (M fractions multiplied).
%   - The geometric mean stays on the scale of a single variable: if every
%     b_m equals the same value b, the geometric mean is exactly b.
%   - The geometric mean still goes to zero whenever ANY one variable is
%     zero (it is a normalised product); the arithmetic mean does not.
%   - Arithmetic mean >= geometric mean >= product at every instant
%     (AM-GM inequality), with min(b) shown as a reference.
%
% Standalone, base MATLAB only. Leaves its figure open.

%% Parameters
f_Hz   = [0.1 0.16 0.5];   % frequency of each b_m (Hz); add entries to add variables
T_sec  = 20;               % duration (s); 20 s gives every signal >= 2 cycles. Set to 2 for a 0-2 s view.
fs_Hz  = 200;              % sample rate (Hz)

%% Signals
t = (0:1/fs_Hz:T_sec).';                 % column time vector
M = numel(f_Hz);
b = 0.5 + 0.5*sin(2*pi*t*f_Hz);          % numel(t) x M, each column one b_m

b_prod  = prod(b, 2);                    % current model
b_geo   = b_prod.^(1/M);                 % geometric mean
b_arith = mean(b, 2);                    % arithmetic mean
b_min   = min(b, [], 2);                 % reference: the deepest single variable

%% Numbers making the same points
fprintf('\nCombining M = %d depression variables over %g s\n', M, T_sec);
fprintf('  time-averaged value:  product %.3f | geometric mean %.3f | arithmetic mean %.3f | single b_m %.3f\n', ...
    mean(b_prod), mean(b_geo), mean(b_arith), mean(b(:)));
fprintf('  minimum value:        product %.3f | geometric mean %.3f | arithmetic mean %.3f\n', ...
    min(b_prod), min(b_geo), min(b_arith));
fprintf('  fraction of time AM >= GM >= product: %.3f\n', ...
    mean(b_arith >= b_geo - 1e-12 & b_geo >= b_prod - 1e-12));
b_eq = 0.6;   % all variables equal: product deepens with M, geometric mean does not
fprintf('  all b_m = %.1f:  product for M = 1,2,3,5 is %s; geometric mean is %.1f for every M\n\n', ...
    b_eq, mat2str(b_eq.^[1 2 3 5], 3), b_eq);

%% Figure
fig_gm = figure('Name', 'Geometric-mean normalization of STD variables', 'Color', 'w');
tl = tiledlayout(fig_gm, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(5, 1);
c_prod  = [0.85 0.33 0.10];
c_geo   = [0.00 0.45 0.74];
c_arith = [0.47 0.67 0.19];

ax(1) = nexttile(tl);
plot(ax(1), t, b, 'LineWidth', 1);
legend(ax(1), compose('b_%d (%.2g Hz)', (1:M).', f_Hz.'), 'Location', 'eastoutside');
title(ax(1), 'Individual depression variables b_m(t)');

ax(2) = nexttile(tl);
plot(ax(2), t, b_prod, 'Color', c_prod, 'LineWidth', 1.5);
title(ax(2), 'Product  \Pi_m b_m  (current model)');

ax(3) = nexttile(tl);
plot(ax(3), t, b_geo, 'Color', c_geo, 'LineWidth', 1.5);
title(ax(3), sprintf('Geometric mean  (\\Pi_m b_m)^{1/%d}', M));

ax(4) = nexttile(tl);
plot(ax(4), t, b_arith, 'Color', c_arith, 'LineWidth', 1.5);
title(ax(4), 'Arithmetic mean  (1/M) \Sigma_m b_m');

ax(5) = nexttile(tl);
hold(ax(5), 'on');
plot(ax(5), t, b_prod,  'Color', c_prod,  'LineWidth', 1.5);
plot(ax(5), t, b_geo,   'Color', c_geo,   'LineWidth', 1.5);
plot(ax(5), t, b_arith, 'Color', c_arith, 'LineWidth', 1.5);
plot(ax(5), t, b_min,   'k:', 'LineWidth', 0.75);
hold(ax(5), 'off');
legend(ax(5), {'product', 'geometric mean', 'arithmetic mean', 'min_m b_m'}, 'Location', 'eastoutside');
title(ax(5), 'Overlay');

for k = 1:5
    ylim(ax(k), [0 1]);
    ylabel(ax(k), 'gain');
    grid(ax(k), 'on');
    box(ax(k), 'off');
    if k < 5
        ax(k).XTickLabel = [];
    end
end
xlabel(ax(5), 'Time (s)');
linkaxes(ax, 'x');
xlim(ax(5), [t(1) t(end)]);

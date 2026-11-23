function fig = so_main_figure(R, S, out_base)
% SO_MAIN_FIGURE The rightmost frozen-gain eigenvalue of the bursting network through bursts and a sustained step.
%
%   fig = so_main_figure(R, S, out_base)
%
% R is so_run_network's result for the network of the stimulation figure
% (original, seed [42 42], with its frozen arm); S the settings. Writes
% out_base.png (300 dpi) and out_base.fig.
%
% Rows, on one time axis (stimulus onset dashed):
%   (A) the input to every neuron
%   (B) population mean rate of E and of I; dashed, E with every
%       adaptation and depression variable held at its pre-step value
%   (C) Re of the rightmost eigenvalue of M = W diag(d), the boundary
%       Re = 1 and the bulk radius predicted from the gains; dashed, the
%       same for the frozen arm. Markers: the states of row (F)
%   (D) coherence of the rightmost eigenvector: fraction of its norm^2
%       constant within E and within I (a random vector has 0.04)
% (E) one burst cycle of (B) and (C), enlarged.
% (F) the eigenvalues of W (every gain 1) and of M at the four marked
%     states, on one scale, the circle being the bulk radius (for W the
%     structural one), the rightmost eigenvalue filled.
%
% See also: run_stim_outlier, so_run_network, complex_plane_axes

Q = R.intact; F = R.frozen;
t = Q.t(:); re1 = real(Q.lam1(:));
t_win = [10 70];
fs = 14;
cE = [0.80 0.15 0.15]; cI = [0.15 0.30 0.80]; cF = [0.35 0.35 0.35];
cB = [0.55 0.55 0.55];
st_col = [0.30 0.55 0.75; 0.95 0.55 0.10; 0.60 0.05 0.10; 0.10 0.60 0.35];
st_name = {'quiet', 'before burst', 'burst', 'stimulation'};

% The burst cycle enlarged and the four states
bb = Q.bursts;
[~, ib] = min(abs(bb.t_on - 25));
ton = bb.t_on(ib); toff = bb.t_off(ib);
zoom = [ton - 2, toff + 1];
k_quiet = find(t >= ton - 1.5, 1);
kk = find(t >= ton - 0.6 & t <= ton + 0.05);
[~, j] = max(re1(kk)); k_pre = kk(j);
kk = find(t >= ton & t <= toff);
[~, j] = max(Q.rate(kk, 1)); k_burst = kk(j);
k_stim = find(t >= 60, 1);
k_st = [k_quiet, k_pre, k_burst, k_stim];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 1500 1120]);
set(fig, 'DefaultAxesFontSize', fs, 'DefaultAxesLineWidth', 1, 'DefaultAxesTickDir', 'out', ...
    'DefaultAxesBox', 'off', 'DefaultAxesTitleFontSizeMultiplier', 1, ...
    'DefaultAxesLabelFontSizeMultiplier', 1, 'DefaultTextFontSize', fs);

xl = 0.08; xw = 0.55;
rows = [0.890 0.07; 0.700 0.16; 0.465 0.21; 0.335 0.10];
ax = gobjects(4, 1);
for r = 1:4
    ax(r) = axes(fig, 'Position', [xl rows(r, 1) xw rows(r, 2)]);
    hold(ax(r), 'on');
end

% (A) input
plot(ax(1), t, Q.u_mean, 'k', 'LineWidth', 1.5);
ylim(ax(1), [-0.05 0.5]); yticks(ax(1), [0 0.4]); ylabel(ax(1), 'u');

% (B) rates
hE = plot(ax(2), t, Q.rate(:, 1), 'Color', cE, 'LineWidth', 1.2);
hI = plot(ax(2), t, Q.rate(:, 2), 'Color', cI, 'LineWidth', 1.2);
kf = F.t(:) >= S.t_stim;
hF = plot(ax(2), F.t(kf), F.rate(kf, 1), '--', 'Color', cF, 'LineWidth', 1.5);
ylim(ax(2), [0 1]); yticks(ax(2), [0 0.5 1]); ylabel(ax(2), 'mean rate');
legend(ax(2), [hE hI hF], {'E', 'I', 'E, adaptation frozen'}, 'Location', 'northwest', ...
    'Box', 'off', 'FontSize', fs - 1, 'NumColumns', 3, 'AutoUpdate', 'off');

% (C) the rightmost eigenvalue
hR = plot(ax(3), t, Q.R_eff, 'Color', cB, 'LineWidth', 1.5);
h1 = plot(ax(3), t, re1, 'k', 'LineWidth', 1.2);
hf = plot(ax(3), F.t(kf), real(F.lam1(kf)), '--', 'Color', cF, 'LineWidth', 1.5);
hb = yline(ax(3), 1, ':', 'Color', [0.8 0 0], 'LineWidth', 1.5);
for s = 1:3
    plot(ax(3), t(k_st(s)), 3.15, 'v', 'MarkerFaceColor', st_col(s, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 9);
end
plot(ax(3), t(k_st(4)), 3.15, 'v', 'MarkerFaceColor', st_col(4, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 9);
ylim(ax(3), [-0.1 3.3]); yticks(ax(3), [0 1 2 3]);
ylabel(ax(3), 'Re \lambda_1(M)');
legend(ax(3), [h1 hR hb hf], {'rightmost eigenvalue', 'bulk radius', 'Re = 1', 'adaptation frozen'}, ...
    'Location', 'east', 'Box', 'off', 'FontSize', fs - 1, 'NumColumns', 2, 'AutoUpdate', 'off');

% (D) coherence
plot(ax(4), t, Q.coh1, 'k', 'LineWidth', 1);
ylim(ax(4), [0 1]); yticks(ax(4), [0 1]); ylabel(ax(4), 'coherence');
xlabel(ax(4), 'time (s)');

for r = 1:4
    xlim(ax(r), t_win);
    xline(ax(r), S.t_stim, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    if r < 4; set(ax(r), 'XTickLabel', []); end
    patch_zoom(ax(r), zoom);
end

% (E) one burst cycle
axz1 = axes(fig, 'Position', [0.72 0.700 0.25 0.26]); hold(axz1, 'on');
plot(axz1, t, Q.rate(:, 1), 'Color', cE, 'LineWidth', 1.5);
plot(axz1, t, Q.rate(:, 2), 'Color', cI, 'LineWidth', 1.5);
yline(axz1, S.burst.thresh, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
xlim(axz1, zoom); ylim(axz1, [0 0.7]); yticks(axz1, [0 0.5]); set(axz1, 'XTickLabel', []);
ylabel(axz1, 'mean rate');
axz2 = axes(fig, 'Position', [0.72 0.465 0.25 0.21]); hold(axz2, 'on');
plot(axz2, t, Q.R_eff, 'Color', cB, 'LineWidth', 1.5);
plot(axz2, t, re1, 'k', 'LineWidth', 1.5);
yline(axz2, 1, ':', 'Color', [0.8 0 0], 'LineWidth', 1.5);
for s = 1:3
    plot(axz2, t(k_st(s)), re1(k_st(s)), 'o', 'MarkerFaceColor', st_col(s, :), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 8);
end
xlim(axz2, zoom); ylim(axz2, [-0.1 3.3]); yticks(axz2, [0 1 2 3]);
ylabel(axz2, 'Re \lambda_1(M)'); xlabel(axz2, 'time (s)');
for a = [axz1 axz2]
    patch(a, [ton toff toff ton], [a.YLim(1) a.YLim(1) a.YLim(2) a.YLim(2)], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    a.Children = a.Children([2:end 1]);
end

% (F) complex planes
opts = {'re_lim', [-2.3 3.8], 'im_lim', [-2.3 2.3], 're_ticks', [-2 2 3], ...
    'im_ticks', [-2 2], 're_clear', [-1.99 1.99], 'im_clear', [-1.99 1.99], ...
    're_label', 'Re', 'im_label', 'Im', 'font_size', fs - 1};
evs = {R.struct.ev_W, Q.ev(k_quiet, :), Q.ev(k_pre, :), Q.ev(k_burst, :), Q.ev(k_stim, :)};
rad = [R.struct.R_struct, Q.R_eff(k_st)'];
ttl = [{'W'}, st_name];
col = [0 0 0; st_col];
pw = 0.17; gap = 0.018; x0 = 0.06;
th = linspace(0, 2 * pi, 200);
axp = gobjects(1, 5);
for p = 1:5
    axp(p) = axes(fig, 'Position', [x0 + (p - 1) * (pw + gap), 0.035, pw, 0.22]);
    hold(axp(p), 'on');
    cp = complex_plane_axes(axp(p), opts{:});
    z = double(evs{p}(:));
    plot(axp(p), cp.x(rad(p) * cos(th)), rad(p) * sin(th), '-', 'Color', cB, 'LineWidth', 1.2);
    plot(axp(p), cp.x([1 1]), [-2.2 2.2], ':', 'Color', [0.8 0 0], 'LineWidth', 1.5);
    plot(axp(p), cp.x(real(z)), imag(z), '.', 'Color', [0.25 0.25 0.25], 'MarkerSize', 9);
    [~, j] = max(real(z));
    plot(axp(p), cp.x(real(z(j))), imag(z(j)), 'o', 'MarkerFaceColor', col(p, :), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 9);
    title(axp(p), ttl{p}, 'FontWeight', 'normal', 'FontSize', fs, 'Color', col(p, :));
end

% Panel letters
lab = {'(A)', ax(1); '(B)', ax(2); '(C)', ax(3); '(D)', ax(4); '(E)', axz1; '(F)', axp(1)};
for i = 1:size(lab, 1)
    p = lab{i, 2}.Position;
    annotation(fig, 'textbox', [p(1) - 0.06, p(2) + p(4) - 0.005, 0.05, 0.03], 'String', lab{i, 1}, ...
        'EdgeColor', 'none', 'FontSize', fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

if nargin >= 3 && ~isempty(out_base)
    exportgraphics(fig, [out_base '.png'], 'Resolution', 300);
    set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
    savefig(fig, [out_base '.fig']);
    close(fig);
end
end

function patch_zoom(a, zoom)
yl = a.YLim;
patch(a, [zoom(1) zoom(2) zoom(2) zoom(1)], [yl(1) yl(1) yl(2) yl(2)], [0.93 0.93 0.93], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
a.Children = a.Children([2:end 1]);
end

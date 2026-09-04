% Multi-timescale STD: steady state and step response. Edit the block below.
%
%   b_m(r) = 1/(1 + r/rho_m),  rho_m = tau_rel_m/tau_rec_m,  theta = r*prod(b)
%
% Only the RATIO rho sets the steady state; the absolute taus set the speed.
% K equal ratios: r_peak = rho/(K-1), theta_peak = r_peak*(1-1/K)^K. K=1 has no
% peak. Two timescales: r_peak = sqrt(rho1*rho2), independent of their spread.

%% ---- parameters -----------------------------------------------------------
tau_rel    = [0.25 1];
tau_rec    = 4 * tau_rel;
step_rates = [0.125 0.52 0.5 1.0];
on_s       = 5;
off_s      = 10;      % < ~4*max(tau_rec) and steps start already depressed
settle_s   = 5;
fs         = 1000;

%% ---- steady state ---------------------------------------------------------
rho    = tau_rel ./ tau_rec;
K      = numel(rho);
r      = linspace(0, 1, 4000);
b_each = 1 ./ (1 + (1 ./ rho(:)) * r);
b_prod = prod(b_each, 1);
theta  = b_prod .* r;

rf  = linspace(0, 1, 4e5);
thf = rf .* prod(1 ./ (1 + (1 ./ rho(:)) * rf), 1);
[theta_peak, ip] = max(thf);
r_peak   = rf(ip);
has_peak = ip < numel(rf);

%% ---- step response (exact: r is piecewise constant, so each segment is one
%%      exponential with tau_eff = 1/(1/tau_rec + r/tau_rel)) -----------------
seg_rate = [0, reshape([step_rates; zeros(1, numel(step_rates))], 1, [])];
seg_dur  = [settle_s, repmat([on_s, off_s], 1, numel(step_rates))];
dt = 1 / fs;
t = []; r_t = []; b_t = zeros(K, 0); b0 = ones(K, 1); t0 = 0;
for s = 1:numel(seg_rate)
    ts   = dt : dt : seg_dur(s);
    teff = 1 ./ (1 ./ tau_rec(:) + seg_rate(s) ./ tau_rel(:));
    binf = teff ./ tau_rec(:);
    bs   = binf + (b0 - binf) .* exp(-ts ./ teff);
    t    = [t, t0 + ts];                            %#ok<AGROW>
    r_t  = [r_t, seg_rate(s) * ones(1, numel(ts))]; %#ok<AGROW>
    b_t  = [b_t, bs];                               %#ok<AGROW>
    b0   = bs(:, end);  t0 = t0 + seg_dur(s);
end
theta_t   = prod(b_t, 1) .* r_t;
b_prod_t  = prod(b_t, 1);
seg_start = cumsum([0, seg_dur]);
ss_theta  = arrayfun(@(x) x * prod(1 ./ (1 + x ./ rho)), step_rates);
ss_b      = arrayfun(@(x)     prod(1 ./ (1 + x ./ rho)), step_rates);

%% ---- report ---------------------------------------------------------------
fprintf('\ntau_rec = [%s]   tau_rel = [%s]\n', ...
    strjoin(compose('%g', tau_rec), ' '), strjoin(compose('%g', tau_rel), ' '));
fprintf('rho = [%s]   K = %d\n', strjoin(compose('%.4g', rho), ' '), K);
if has_peak
    fprintf('r_peak = %.4g   theta_peak = %.4g   ratio = %.4g\n', ...
        r_peak, theta_peak, theta_peak / r_peak);
    fprintf('theta_peak/theta(1) = %.3f\n', theta_peak / theta(end));
else
    fprintf('no peak (K = 1): theta rises monotonically to rho = %.4g\n', rho(1));
end
fprintf('recovery between steps = %.1f%%\n', 100 * (1 - exp(-off_s / max(tau_rec))));

%% ---- figure ---------------------------------------------------------------
c_prod = [0.85 0.325 0.098];
c_one  = [0.5 0.5 0.5];
c_id   = [0.55 0.80 0.55];
fs_t = 12; fs_l = 14; lw = 2;

sz  = [900 900];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', 'Position', [scr(1:2) + max((scr(3:4)-sz)/2, 0), sz]);
tl  = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('\\tau_{rec} = [%s],  \\tau_{rel} = [%s],  \\rho = [%s]', ...
    strjoin(compose('%g', tau_rec), ' '), strjoin(compose('%g', tau_rel), ' '), ...
    strjoin(compose('%.3g', rho), ' ')), 'FontSize', 15);

% 1: prod(b) vs r
ax1 = nexttile(tl); hold(ax1, 'on');
plot(ax1, r, b_each(1,:), '--', 'LineWidth', 1, 'Color', c_one);
plot(ax1, r, b_prod, 'LineWidth', lw, 'Color', c_prod);
box(ax1,'off'); set(ax1,'FontSize',fs_t,'XTick',[0 1],'YTick',[0 1]);
xlim(ax1,[0 1]); ylim(ax1,[0 1]);
xlabel(ax1,'firing rate  r','FontSize',fs_l);
ylabel(ax1,'depression  $\prod_k b_k$','Interpreter','latex','FontSize',fs_l);
title(ax1,'Steady-state depression','FontWeight','normal','FontSize',15);
legend(ax1,{'single $b_k$','$\prod_k b_k$'},'Interpreter','latex','Box','off', ...
    'FontSize',10,'Location','northeast');

% 2: prod(b)*r vs r
ymax2 = max([b_each(1,:).*r, theta]) * 1.15;
ax2 = nexttile(tl); hold(ax2, 'on');
h_id  = plot(ax2, [0 1], [0 1], '-', 'LineWidth', 1, 'Color', c_id);
h_one = plot(ax2, r, b_each(1,:).*r, '--', 'LineWidth', 1, 'Color', c_one);
h_pr  = plot(ax2, r, theta, 'LineWidth', lw, 'Color', c_prod);
if has_peak
    plot(ax2, r_peak, theta_peak, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', c_prod, 'MarkerEdgeColor', 'none');
    text(ax2, r_peak, theta_peak, sprintf('  (%.2f, %.2f)', r_peak, theta_peak), ...
        'FontSize', 11, 'VerticalAlignment','bottom','HorizontalAlignment','left');
    plot(ax2, 1, theta(end), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', c_prod, 'MarkerEdgeColor', 'none');
    text(ax2, 1, theta(end), sprintf('R = %.2f  ', theta_peak / theta(end)), ...
        'FontSize', 11, 'VerticalAlignment','bottom','HorizontalAlignment','right');
end
box(ax2,'off'); set(ax2,'FontSize',fs_t,'XTick',[0 1],'YTick',[0 round(ymax2,3)]);
xlim(ax2,[0 1]); ylim(ax2,[0 ymax2]);
xlabel(ax2,'firing rate  r','FontSize',fs_l);
ylabel(ax2,'synaptic output  $\prod_k b_k \cdot r$','Interpreter','latex','FontSize',fs_l);
title(ax2,'Steady-state synaptic output','FontWeight','normal','FontSize',15);
legend(ax2,[h_id h_one h_pr], ...
    {'undepressed  $y = r$','single $b_k \cdot r$','$\prod_k b_k \cdot r$'}, ...
    'Interpreter','latex','Box','off','FontSize',10,'Location','southeast');

% 3 (tiles 3-4): depression step response
ax3 = nexttile(tl, 3, [1 2]); hold(ax3, 'on');
h3r = plot(ax3, t, r_t, '-', 'LineWidth', 1.5, 'Color', c_id);
for k = 1:numel(ss_b)
    h3s = plot(ax3, seg_start([2*k 2*k+1]), ss_b([k k]), '-', ...
        'LineWidth', 1.5, 'Color', c_one);
end
h3o = plot(ax3, t, b_t(1,:), '--', 'LineWidth', 1, 'Color', c_one);
h3p = plot(ax3, t, b_prod_t, 'LineWidth', lw, 'Color', c_prod);
box(ax3,'off'); set(ax3,'FontSize',fs_t); xlim(ax3,[0 t(end)]); ylim(ax3,[0 1.05]);
xlabel(ax3,'time (s)','FontSize',fs_l);
ylabel(ax3,'depression  $\prod_k b_k$','Interpreter','latex','FontSize',fs_l);
title(ax3,'Step response: depression','FontWeight','normal','FontSize',15);
legend(ax3,[h3r h3s h3o h3p], ...
    {'rate  $r(t)$','steady state','single $b_k$','depression  $\prod_k b_k$'}, ...
    'Interpreter','latex','Box','off','FontSize',10,'Location','southeast');

% 4 (tiles 5-6): synaptic output step response
ax4 = nexttile(tl, 5, [1 2]); hold(ax4, 'on');
h4r = plot(ax4, t, r_t, '-', 'LineWidth', 1.5, 'Color', c_id);
for k = 1:numel(ss_theta)
    h4s = plot(ax4, seg_start([2*k 2*k+1]), ss_theta([k k]), '-', ...
        'LineWidth', 1.5, 'Color', c_one);
end
h4o = plot(ax4, t, b_t(1,:).*r_t, '--', 'LineWidth', 1, 'Color', c_one);
h4p = plot(ax4, t, theta_t, 'LineWidth', lw, 'Color', c_prod);
box(ax4,'off'); set(ax4,'FontSize',fs_t);
xlim(ax4,[0 t(end)]); ylim(ax4,[0 max([theta_t r_t])*1.12]);
xlabel(ax4,'time (s)','FontSize',fs_l);
ylabel(ax4,'synaptic output  $\prod_k b_k \cdot r$','Interpreter','latex','FontSize',fs_l);
title(ax4,'Step response: synaptic output','FontWeight','normal','FontSize',15);
legend(ax4,[h4r h4s h4o h4p], ...
    {'rate  $r(t)$','steady state','single $b_k \cdot r$', ...
    'synaptic output  $\prod_k b_k \cdot r$'}, ...
    'Interpreter','latex','Box','off','FontSize',10,'Location','northeast');

linkaxes([ax3 ax4], 'x');

function out = explore_std_steady_state(opts)
% EXPLORE_STD_STEADY_STATE Multi-timescale STD for a chosen (tau_rec, tau_rel).
%
%   out = EXPLORE_STD_STEADY_STATE()                        % the preset's values
%   out = EXPLORE_STD_STEADY_STATE('tau_rec', [2 4], 'tau_rel', [0.5 1.0])
%   out = EXPLORE_STD_STEADY_STATE('tau_rel', 1, 'tau_rec', 2)   % single timescale
%
% Three panels for one set of depression timescales:
%
%   1. prod(b) against rate  -- the steady-state gain.
%   2. prod(b)*r against rate -- what the recurrent sum actually receives, with
%      the peak marked and the identity line y = r for the undepressed synapse.
%   3. The STEP RESPONSE: rate held at each of step_rates for on_s seconds,
%      separated by off_s seconds at zero, with the synaptic output plotted
%      against TIME.
%
% Panels 2 and 3 are the same quantity seen two ways, and reading them together
% is the point: the dashed line each plateau in panel 3 settles onto is exactly
% the value panel 2 reports at that rate. Panel 3 adds what panel 2 cannot show
% -- the TRANSIENT. At the onset of a step from a rested synapse b = 1, so the
% output jumps to the full r before depression pulls it down, which means the
% peak transient grows with r even where the steady state SHRINKS with r.
%
% WHY THIS IS AN EXPLORATION AND NOT A FIGURE. src/figures/fig_STD_steady_state
% draws panels 1 and 2 from whatever the PRESET carries and is what the paper
% ships. This takes the timescales as arguments so they can be varied without
% touching a preset -- editing the preset is a physics change that invalidates
% the frozen fixtures (test_preset_golden, test_c_over_K) and should be done
% deliberately, not while exploring.
%
% THE STEADY STATE DEPENDS ONLY ON THE RATIO rho = tau_rel/tau_rec:
%
%   b_m(r) = 1/(1 + r/rho_m),   theta(r) = r * prod_m b_m(r)
%
% The absolute timescales set how FAST b gets there -- panel 3 -- and nothing
% about where. For two timescales the peak has a closed form:
%
%   r_peak     = sqrt(rho_1 * rho_2)                  (the GEOMETRIC MEAN)
%   theta_peak = rho_1*rho_2 / (sqrt(rho_1) + sqrt(rho_2))^2
%
% so separating the ratios lowers the peak WITHOUT MOVING IT, and
% theta_peak <= r_peak/4 always, with equality only when rho_1 == rho_2. More
% generally, K equal ratios give theta_peak = r_peak*(1 - 1/K)^K, which is 1/4
% at K = 2 and approaches 1/e. A SINGLE timescale has no peak at all: theta
% rises monotonically to rho. The turnover is a multi-timescale phenomenon.
%
% Note the defaults have rho = [0.125 0.125] -- the paper preset's two
% timescales differ in SPEED (tau_rec 2 s vs 4 s) but share a ratio, so they are
% indistinguishable in panels 1 and 2 and differ only in panel 3.
%
% Depends on nothing in src/, so it runs in a cold session.
%
% See also: fig_STD_steady_state, srnn_param_preset

arguments
    % Defaults are celltype_pairs_..._dualStd_*'s E->E route.
    opts.tau_rec     (1,:) double  = [2 4]
    opts.tau_rel     (1,:) double  = [0.25 0.5]
    opts.step_rates  (1,:) double  = [0.25 0.5 1]
    opts.on_s        (1,1) double  = 5
    % Long enough for the SLOWEST tau_rec to recover between steps: at the
    % default 15 s and tau_rec = 4 s, b returns to 1 - exp(-15/4) = 97.6%. Too
    % short and each step starts more depressed than the last, which reads as a
    % property of the rate rather than of the protocol. The realised recovery is
    % returned in out.recovery_frac -- check it after changing tau_rec.
    opts.off_s       (1,1) double  = 15
    opts.settle_s    (1,1) double  = 2
    opts.fs          (1,1) double  = 1000
    opts.save        (1,1) logical = false
    opts.out_dir     (1,:) char    = ''
end

tau_rec = opts.tau_rec(:)';
tau_rel = opts.tau_rel(:)';
if numel(tau_rec) ~= numel(tau_rel)
    error('explore_std_steady_state:LengthMismatch', ...
        'tau_rec and tau_rel must have the same number of timescales (got %d and %d).', ...
        numel(tau_rec), numel(tau_rel));
end
if any(tau_rec <= 0) || any(tau_rel <= 0)
    error('explore_std_steady_state:NonPositiveTau', ...
        'All tau_rec and tau_rel must be positive.');
end
K   = numel(tau_rec);
rho = tau_rel ./ tau_rec;          % the ONLY combination the steady state sees

%% ---- Steady state ---------------------------------------------------------
r      = linspace(0, 1, 4000);
b_each = 1 ./ (1 + (1 ./ rho(:)) * r);      % rows = timescale
b_prod = prod(b_each, 1);
theta  = b_prod .* r;
theta_single = b_each(1, :) .* r;

% Peak found on a fine grid so it is right for any K, then cross-checked against
% the closed form when K == 2. A mismatch there means the algebra in the header
% has stopped describing the code.
r_fine     = linspace(0, 1, 400000);
th_fine    = r_fine .* prod(1 ./ (1 + (1 ./ rho(:)) * r_fine), 1);
[theta_peak, ip] = max(th_fine);
r_peak     = r_fine(ip);
has_peak   = ip < numel(r_fine);            % K == 1 is monotone: no turnover
if K == 2
    r_peak_exact     = sqrt(prod(rho));
    theta_peak_exact = prod(rho) / sum(sqrt(rho))^2;
else
    r_peak_exact = NaN; theta_peak_exact = NaN;
end

%% ---- Step response, integrated EXACTLY ------------------------------------
% r(t) is piecewise constant, and on a segment of constant r the ODE
%   db/dt = (1 - b)/tau_rec - b*r/tau_rel
% is linear with time constant tau_eff = 1/(1/tau_rec + r/tau_rel) and fixed
% point b_inf = tau_eff/tau_rec. So each segment is one exponential, evaluated
% in closed form -- no integrator, no tolerance, and the plateaus are exactly
% the steady state panel 2 plots rather than approximately it.
seg_rate = [0, reshape([opts.step_rates; zeros(1, numel(opts.step_rates))], 1, [])];
seg_dur  = [opts.settle_s, repmat([opts.on_s, opts.off_s], 1, numel(opts.step_rates))];

dt = 1 / opts.fs;
t  = []; r_t = []; b_t = zeros(K, 0);
b0 = ones(K, 1);                            % rested synapse
t0 = 0;
for s = 1:numel(seg_rate)
    ts   = (dt : dt : seg_dur(s));
    rs   = seg_rate(s);
    teff = 1 ./ (1 ./ tau_rec(:) + rs ./ tau_rel(:));
    binf = teff ./ tau_rec(:);
    bs   = binf + (b0 - binf) .* exp(-ts ./ teff);
    t    = [t,   t0 + ts];                  %#ok<AGROW>
    r_t  = [r_t, rs * ones(1, numel(ts))];  %#ok<AGROW>
    b_t  = [b_t, bs];                       %#ok<AGROW>
    b0   = bs(:, end);
    t0   = t0 + seg_dur(s);
end
theta_t = prod(b_t, 1) .* r_t;

% How far b actually recovers in one off period, from fully depressed. The
% slowest timescale is the binding one.
recovery_frac = 1 - exp(-opts.off_s / max(tau_rec));

%% ---- Report --------------------------------------------------------------
fprintf('\n--- STD steady state, %d timescale(s) ---\n', K);
fprintf('  tau_rec = [%s]\n', strjoin(compose('%g', tau_rec), ' '));
fprintf('  tau_rel = [%s]\n', strjoin(compose('%g', tau_rel), ' '));
fprintf('  rho = tau_rel/tau_rec = [%s]\n', strjoin(compose('%.4f', rho), ' '));
if K == 2 && abs(rho(1) - rho(2)) < 1e-12
    fprintf('  NOTE: both ratios equal -- the timescales are indistinguishable\n');
    fprintf('        in panels 1-2 and differ only in the transients of panel 3.\n');
end
if has_peak
    fprintf('  r_peak     = %.4f', r_peak);
    if K == 2; fprintf('   (closed form %.4f)', r_peak_exact); end
    fprintf('\n  theta_peak = %.5f', theta_peak);
    if K == 2; fprintf('   (closed form %.5f)', theta_peak_exact); end
    fprintf('\n  theta_peak/r_peak = %.4f   ((1-1/K)^K = %.4f)\n', ...
        theta_peak / r_peak, (1 - 1/K)^K);
else
    fprintf('  NO PEAK: a single timescale is monotone, saturating at rho = %.4f\n', rho(1));
end
fprintf('  recovery between steps: %.1f%% (off_s = %g s, slowest tau_rec = %g s)\n', ...
    100 * recovery_frac, opts.off_s, max(tau_rec));
fprintf('  steady-state theta at each step rate:\n');
for k = 1:numel(opts.step_rates)
    rk = opts.step_rates(k);
    fprintf('     r = %.3f  ->  theta = %.5f   (onset transient reaches %.3f)\n', ...
        rk, rk * prod(1 ./ (1 + rk ./ rho)), rk);
end
fprintf('\n');

%% ---- Figure --------------------------------------------------------------
prod_color     = [0.85 0.325 0.098];
single_color   = [0.5 0.5 0.5];
identity_color = [0.55 0.80 0.55];
tick_fs = 12; label_fs = 14; title_fs = 15; lw = 2;

fig_size = [1180, 340];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', ...
    'Position', [scr(1:2) + max((scr(3:4) - fig_size)/2, 0), fig_size]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('\\tau_{rec} = [%s] s,   \\tau_{rel} = [%s] s,   \\rho = [%s]', ...
    strjoin(compose('%g', tau_rec), ' '), strjoin(compose('%g', tau_rel), ' '), ...
    strjoin(compose('%.3g', rho), ' ')), 'FontSize', title_fs);

% --- 1. prod(b) vs r ---
ax1 = nexttile(tl); hold(ax1, 'on');
plot(ax1, r, b_each(1, :), '--', 'LineWidth', 1, 'Color', single_color);
plot(ax1, r, b_prod, 'LineWidth', lw, 'Color', prod_color);
hold(ax1, 'off'); box(ax1, 'off'); set(ax1, 'FontSize', tick_fs);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'depression  $\prod_k b_k$', 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax1, 'Steady-state depression', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax1, [0 1]); ylim(ax1, [0 1]); set(ax1, 'XTick', [0 1], 'YTick', [0 1]);
legend(ax1, {'single $b_k$', '$\prod_k b_k$'}, 'Interpreter', 'latex', ...
    'Box', 'off', 'FontSize', 11, 'Location', 'northeast');

% --- 2. prod(b)*r vs r ---
% Zoomed like the paper figure's 'zoom' variant: at full scale the curves sit in
% the bottom tenth and the turnover -- the point of the panel -- is invisible.
% The 1:1 aspect is therefore dropped, so no angle here means anything.
zoom_ymax = max([theta_single, theta]) * 1.15;
ax2 = nexttile(tl); hold(ax2, 'on');
plot(ax2, [0 1], [0 1], '-', 'LineWidth', 1, 'Color', identity_color);
plot(ax2, r, theta_single, '--', 'LineWidth', 1, 'Color', single_color);
plot(ax2, r, theta, 'LineWidth', lw, 'Color', prod_color);
if has_peak
    plot(ax2, r_peak, theta_peak, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', prod_color, 'MarkerEdgeColor', 'none');
    text(ax2, r_peak, theta_peak, sprintf('  r = %.3f', r_peak), 'FontSize', 11, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
hold(ax2, 'off'); box(ax2, 'off'); set(ax2, 'FontSize', tick_fs);
xlabel(ax2, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax2, 'synaptic output  $\prod_k b_k \cdot r$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax2, 'Delivered output', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax2, [0 1]); ylim(ax2, [0 zoom_ymax]);
set(ax2, 'XTick', [0 1], 'YTick', [0 round(zoom_ymax, 3)]);

% --- 3. step response ---
ax3 = nexttile(tl); hold(ax3, 'on');
% The drive at TRUE SCALE, sharing the output's axis. This was briefly plotted
% rescaled to fill the panel, which made the steps read as 0.3/0.6/1.2 -- the
% rates were right and only the reference trace was stretched, but it invited
% exactly the misreading it got, and it threw away the reason to draw r(t) at
% all: because b = 1 at rest, the onset of each step delivers very nearly the
% FULL r, so theta rises to meet the drive and then collapses away from it. The
% small visible shortfall at each onset is incomplete recovery during off_s.
% Same units, same axis, no scale factor.
ss = arrayfun(@(rk) rk * prod(1 ./ (1 + rk ./ rho)), opts.step_rates);
y_max = max([theta_t, r_t, ss]) * 1.12;
plot(ax3, t, r_t, '-', 'LineWidth', 1, 'Color', [0.82 0.82 0.82]);
% Each plateau's asymptote -- the SAME number panel 2 plots at that rate.
for k = 1:numel(ss)
    plot(ax3, [0 t(end)], [ss(k) ss(k)], ':', 'LineWidth', 1, 'Color', identity_color);
    text(ax3, t(end), ss(k), sprintf(' r=%.2g ', opts.step_rates(k)), ...
        'FontSize', 10, 'Color', [0.3 0.5 0.3], ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
end
plot(ax3, t, theta_t, 'LineWidth', lw, 'Color', prod_color);
hold(ax3, 'off'); box(ax3, 'off'); set(ax3, 'FontSize', tick_fs);
xlabel(ax3, 'time (s)', 'FontSize', label_fs);
ylabel(ax3, 'synaptic output  $\prod_k b_k \cdot r$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax3, 'Step response', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax3, [0 t(end)]); ylim(ax3, [0 y_max]);
legend(ax3, {'rate  r(t)', 'steady state'}, 'Box', 'off', 'FontSize', 10, ...
    'Location', 'northwest');

%% ---- Outputs -------------------------------------------------------------
out = struct('fig', fig, 'rho', rho, 'r_peak', r_peak, 'theta_peak', theta_peak, ...
    'has_peak', has_peak, 'r_peak_exact', r_peak_exact, ...
    'theta_peak_exact', theta_peak_exact, 'recovery_frac', recovery_frac, ...
    't', t, 'r_t', r_t, 'b_t', b_t, 'theta_t', theta_t, ...
    'r', r, 'b_prod', b_prod, 'theta', theta, 'files', {{}});

if opts.save
    out_dir = opts.out_dir;
    if isempty(out_dir)
        out_dir = fullfile(fileparts(mfilename('fullpath')), 'output');
    end
    if ~isfolder(out_dir); mkdir(out_dir); end
    tag = sprintf('STD_rec_%s_rel_%s', ...
        strjoin(compose('%g', tau_rec), '_'), strjoin(compose('%g', tau_rel), '_'));
    tag = strrep(tag, '.', 'p');
    for ext = {'png', 'fig'}
        f = fullfile(out_dir, [tag '.' ext{1}]);
        if strcmp(ext{1}, 'png')
            exportgraphics(fig, f, 'Resolution', 300);
        else
            savefig(fig, f);
        end
        out.files{end+1} = f;
    end
    fprintf('Saved: %s\n', strjoin(out.files, '\n       '));
end
end

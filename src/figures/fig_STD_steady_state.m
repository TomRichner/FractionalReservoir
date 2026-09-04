function out = fig_STD_steady_state(cfg)
% FIG_STD_STEADY_STATE Multi-timescale STD: steady state and step response.
%
%   out = FIG_STD_STEADY_STATE()
%   out = FIG_STD_STEADY_STATE('preset_name', p)
%
% Conceptual, ANALYTIC apart from the step response, which is integrated in
% closed form. Setting db/dt = 0 in
%     db_k/dt = (1 - b_k)/tau_rec_k - b_k*r/tau_rel_k
% gives b_k(r) = 1/(1 + r/rho_k) with rho_k = tau_rel_k/tau_rec_k, and the
% synapse multiplies the timescales, so it sees prod_k b_k(r).
%
% Four panels, in a 3x2 grid:
%
%   1. prod(b) against rate           -- steady-state gain.
%   2. prod(b)*r against rate         -- what the recurrent sum receives, with
%      the peak marked and R = theta_peak/theta(1) at the right-hand end.
%   3. prod(b) against TIME           -- step response of depression.
%   4. prod(b)*r against TIME         -- step response of synaptic output.
%
% Reading the columns together is the point. The grey level each plateau in
% rows 2-3 settles onto is exactly the value the panel above it reports at that
% rate; what the step panels add is the TRANSIENT. From a rested synapse b = 1,
% so each onset delivers very nearly the full r before depression pulls it down
% -- the onset grows with r even where the STEADY STATE SHRINKS with r.
%
% THE STEADY STATE DEPENDS ONLY ON THE RATIO rho. The absolute timescales set
% how fast b gets there (rows 2-3) and nothing about where (row 1). For K equal
% ratios, r_peak = rho/(K-1) and theta_peak = r_peak*(1 - 1/K)^K; for two
% timescales r_peak = sqrt(rho_1*rho_2), independent of their separation. A
% SINGLE timescale has no peak at all -- theta rises monotonically to rho -- so
% the turnover is a multi-timescale phenomenon, not a stronger version of the
% same thing.
%
% R = theta_peak/theta(1) is (1+rho_1)(1+rho_2)/(sqrt(rho_1)+sqrt(rho_2))^2 at
% K = 2, in which the product rho_1*rho_2 cancels. It is unbounded, grows as rho
% falls, and at fixed geometric mean is largest when the ratios are EQUAL.
%
% CONTRAST WITH SFA (fig_SFA_steady_state): SFA enters as a SUM, so splitting c
% as a budget makes the timescale count invisible at steady state. STD enters as
% a PRODUCT, so two timescales SQUARE the depression however the taus are
% chosen. There is no budget-split of tau that would make dual STD match single.
%
% THE TIMESCALES COME FROM THE PRESET, and specifically from its CONDITIONS:
% synapse_config can only reach the model through a condition, so that is where
% a preset puts its depression routes -- never in the model_defaults struct.
%
% Was two figures (a full-scale 'square' variant with daspect [1 1 1] and a
% zoomed one) showing panels 1-2 only. It is now ONE figure carrying the step
% responses as well, so it writes 3 files rather than 6.
%
% See also: fig_SFA_steady_state, srnn_param_preset,
%           scripts/explorations/explore_std_steady_state

arguments
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.route_pre   (1,:) char    = 'E'
    cfg.route_post  (1,:) char    = 'E'
    cfg.step_rates  (1,:) double  = [0.25 0.5 1]
    cfg.on_s        (1,1) double  = 5
    % Long enough for the slowest tau_rec to recover between steps, so the steps
    % read as independent responses rather than as accumulating depression. At
    % tau_rec = 4 s this gives 1 - exp(-15/4) = 98%; the realised figure is
    % reported in out.recovery_frac.
    cfg.off_s       (1,1) double  = 15
    cfg.settle_s    (1,1) double  = 5
    cfg.fs          (1,1) double  = 1000
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

% Pull the depression timescales out of the preset's own STD routes. They live
% on the conditions, not on the model_defaults struct: synapse_config can only
% reach the model through a condition, so that is where a preset puts them.
% The most-adapted regime, resolved rather than named: which condition carries
% the full route set differs by preset (sfa3_std2 here, sfa3_std1 for a
% single-timescale network).
[~, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
sc = conditions{strcmp(cond_names, full_adaptation_condition(conditions))}.synapse_config;
% route_pre/route_post were declared but ignored -- the route was hardwired to
% E.E. Honoured now; the defaults reproduce the previous behaviour exactly.
if ~isfield(sc, cfg.route_pre) || ~isfield(sc.(cfg.route_pre), cfg.route_post) ...
        || ~isfield(sc.(cfg.route_pre).(cfg.route_post), 'std')
    error('fig_STD_steady_state:NoSuchRoute', ...
        'Preset ''%s'' has no STD on route %s->%s.', ...
        cfg.preset_name, cfg.route_pre, cfg.route_post);
end
route   = sc.(cfg.route_pre).(cfg.route_post).std;
tau_rec = route.tau_rec(:)';
tau_rel = route.tau_rel(:)';
rho     = tau_rel ./ tau_rec;       % the only combination that sets the steady state
K       = numel(rho);

%% ---- Steady state ---------------------------------------------------------
r      = linspace(0, 1, 4000);      % rate, over the full range of the nonlinearity
b_each = 1 ./ (1 + (1 ./ rho(:)) * r);
b_prod = prod(b_each, 1);
theta  = b_prod .* r;
theta_single = b_each(1, :) .* r;

rf  = linspace(0, 1, 4e5);
thf = rf .* prod(1 ./ (1 + (1 ./ rho(:)) * rf), 1);
[theta_peak, ip] = max(thf);
r_peak   = rf(ip);
has_peak = ip < numel(rf);          % K == 1 is monotone: no turnover

%% ---- Step response, integrated exactly ------------------------------------
% r(t) is piecewise constant, and on a segment of constant r the ODE is linear
% with time constant tau_eff = 1/(1/tau_rec + r/tau_rel) and fixed point
% tau_eff/tau_rec. Each segment is therefore one exponential in closed form, so
% the plateaus are exactly the steady state of the panels above rather than
% approximately it.
seg_rate = [0, reshape([cfg.step_rates; zeros(1, numel(cfg.step_rates))], 1, [])];
seg_dur  = [cfg.settle_s, repmat([cfg.on_s, cfg.off_s], 1, numel(cfg.step_rates))];
dt = 1 / cfg.fs;
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
b_prod_t  = prod(b_t, 1);
theta_t   = b_prod_t .* r_t;
seg_start = cumsum([0, seg_dur]);
ss_b      = arrayfun(@(x)     prod(1 ./ (1 + x ./ rho)), cfg.step_rates);
ss_theta  = arrayfun(@(x) x * prod(1 ./ (1 + x ./ rho)), cfg.step_rates);
recovery_frac = 1 - exp(-cfg.off_s / max(tau_rec));

%% ---- Figure ---------------------------------------------------------------
prod_color     = [0.85 0.325 0.098];   % warm, matching the E colour used elsewhere
single_color   = [0.5 0.5 0.5];
% Green is the UNDEPRESSED reference in every panel: y = r in row 1, and r(t)
% itself in the step panels, which is what this synapse would deliver at b = 1.
identity_color = [0.55 0.80 0.55];
tick_fs = st.tick_fs; label_fs = st.label_fs; title_fs = 15; lw = st.line_lw;

% Size computed rather than copied: the older version hardcoded a second-monitor
% x coordinate that landed the window off-screen on a single display.
fig_size = [900, 900];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', ...
    'Position', [scr(1:2) + max((scr(3:4) - fig_size)/2, 0), fig_size]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('\\tau_{rec} = [%s] s,   \\tau_{rel} = [%s] s,   \\rho = [%s]', ...
    strjoin(compose('%g', tau_rec), ' '), strjoin(compose('%g', tau_rel), ' '), ...
    strjoin(compose('%.3g', rho), ' ')), 'FontSize', title_fs);

% --- 1: steady-state depression ---
ax1 = nexttile(tl); hold(ax1, 'on');
plot(ax1, r, b_each(1, :), '--', 'LineWidth', 1, 'Color', single_color);
plot(ax1, r, b_prod, 'LineWidth', lw, 'Color', prod_color);
box(ax1, 'off'); set(ax1, 'FontSize', tick_fs, 'XTick', [0 1], 'YTick', [0 1]);
xlim(ax1, [0 1]); ylim(ax1, [0 1]);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'depression  $\prod_k b_k$', 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax1, 'Steady-state depression', 'FontWeight', 'normal', 'FontSize', title_fs);
legend(ax1, {'single $b_k$', '$\prod_k b_k$'}, 'Interpreter', 'latex', ...
    'Box', 'off', 'FontSize', 10, 'Location', 'northeast');

% --- 2: steady-state synaptic output ---
% Zoomed: at full scale both curves sit in the bottom tenth of the panel and the
% turnover, which is the point, is invisible. The 1:1 aspect the old 'square'
% variant used is therefore dropped, so no angle here means anything.
ymax2 = max([theta_single, theta]) * 1.15;
ax2 = nexttile(tl); hold(ax2, 'on');
h_id  = plot(ax2, [0 1], [0 1], '-', 'LineWidth', 1, 'Color', identity_color);
h_one = plot(ax2, r, theta_single, '--', 'LineWidth', 1, 'Color', single_color);
h_pr  = plot(ax2, r, theta, 'LineWidth', lw, 'Color', prod_color);
if has_peak
    plot(ax2, r_peak, theta_peak, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', prod_color, 'MarkerEdgeColor', 'none');
    text(ax2, r_peak, theta_peak, sprintf('  (%.2f, %.2f)', r_peak, theta_peak), ...
        'FontSize', 11, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    plot(ax2, 1, theta(end), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', prod_color, 'MarkerEdgeColor', 'none');
    text(ax2, 1, theta(end), sprintf('R = %.2f  ', theta_peak / theta(end)), ...
        'FontSize', 11, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs, 'XTick', [0 1], 'YTick', [0 round(ymax2, 3)]);
xlim(ax2, [0 1]); ylim(ax2, [0 ymax2]);
xlabel(ax2, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax2, 'synaptic output  $\prod_k b_k \cdot r$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax2, 'Steady-state synaptic output', 'FontWeight', 'normal', 'FontSize', title_fs);
legend(ax2, [h_id h_one h_pr], ...
    {'undepressed  $y = r$', 'single $b_k \cdot r$', '$\prod_k b_k \cdot r$'}, ...
    'Interpreter', 'latex', 'Box', 'off', 'FontSize', 10, 'Location', 'southeast');

% --- 3 (tiles 3-4): step response, depression ---
% Handles are named for the legends: three steady-state segments sit between the
% drive and the traces, so a legend given only labels attaches them wrongly.
ax3 = nexttile(tl, 3, [1 2]); hold(ax3, 'on');
h3r = plot(ax3, t, r_t, '-', 'LineWidth', 1.5, 'Color', identity_color);
h3s = gobjects(1, numel(ss_b));
for k = 1:numel(ss_b)
    h3s(k) = plot(ax3, seg_start([2*k 2*k+1]), ss_b([k k]), '-', ...
        'LineWidth', 1.5, 'Color', single_color);
end
h3o = plot(ax3, t, b_t(1, :), '--', 'LineWidth', 1, 'Color', single_color);
h3p = plot(ax3, t, b_prod_t, 'LineWidth', lw, 'Color', prod_color);
box(ax3, 'off'); set(ax3, 'FontSize', tick_fs);
xlim(ax3, [0 t(end)]); ylim(ax3, [0 1.05]);
xlabel(ax3, 'time (s)', 'FontSize', label_fs);
ylabel(ax3, 'depression  $\prod_k b_k$', 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax3, 'Step response: depression', 'FontWeight', 'normal', 'FontSize', title_fs);
legend(ax3, [h3r h3s(1) h3o h3p], ...
    {'rate  $r(t)$', 'steady state', 'single $b_k$', 'depression  $\prod_k b_k$'}, ...
    'Interpreter', 'latex', 'Box', 'off', 'FontSize', 10, 'Location', 'southeast');

% --- 4 (tiles 5-6): step response, synaptic output ---
ax4 = nexttile(tl, 5, [1 2]); hold(ax4, 'on');
h4r = plot(ax4, t, r_t, '-', 'LineWidth', 1.5, 'Color', identity_color);
h4s = gobjects(1, numel(ss_theta));
for k = 1:numel(ss_theta)
    h4s(k) = plot(ax4, seg_start([2*k 2*k+1]), ss_theta([k k]), '-', ...
        'LineWidth', 1.5, 'Color', single_color);
end
h4o = plot(ax4, t, b_t(1, :) .* r_t, '--', 'LineWidth', 1, 'Color', single_color);
h4p = plot(ax4, t, theta_t, 'LineWidth', lw, 'Color', prod_color);
box(ax4, 'off'); set(ax4, 'FontSize', tick_fs);
xlim(ax4, [0 t(end)]); ylim(ax4, [0 max([theta_t, r_t]) * 1.12]);
xlabel(ax4, 'time (s)', 'FontSize', label_fs);
ylabel(ax4, 'synaptic output  $\prod_k b_k \cdot r$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax4, 'Step response: synaptic output', 'FontWeight', 'normal', 'FontSize', title_fs);
legend(ax4, [h4r h4s(1) h4o h4p], ...
    {'rate  $r(t)$', 'steady state', 'single $b_k \cdot r$', ...
     'synaptic output  $\prod_k b_k \cdot r$'}, ...
    'Interpreter', 'latex', 'Box', 'off', 'FontSize', 10, 'Location', 'northeast');

% One time axis for both step panels: same protocol, only useful compared.
linkaxes([ax3, ax4], 'x');

%% ---- Save -----------------------------------------------------------------
if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_STD_steady_state';
out = struct('figs', fig, 'files', {{}}, 'source', ['preset: ' cfg.preset_name], ...
    'rho', rho, 'r_peak', r_peak, 'theta_peak', theta_peak, 'has_peak', has_peak, ...
    'R', theta_peak / theta(end), 'recovery_frac', recovery_frac);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

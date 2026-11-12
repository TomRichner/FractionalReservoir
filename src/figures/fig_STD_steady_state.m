function out = fig_STD_steady_state(cfg)
% FIG_STD_STEADY_STATE Short-term depression: steady state and step response.
%
%   out = FIG_STD_STEADY_STATE()
%   out = FIG_STD_STEADY_STATE('preset_name', p)
%   out = FIG_STD_STEADY_STATE('preset_name', p, 'conditions', names)
%
% TWO PATHS, one file name (Fig_STD_steady_state):
%
% * THE CONDITION FIGURE (condition_figure, below) -- used whenever
%   cfg.conditions names a selection, or the preset's most-adapted route is
%   not a plain product (a chain, a geometric-mean combine, facilitation).
%   One curve per selected condition that carries depression on the route,
%   titled and coloured by condition NAME (manuscript_style). Every number is
%   the ROUTE MODEL's: the route is compiled by
%   SRNNCellTypePairs.compile_route from the condition's synapse_config, the
%   steady state is the model's `steady_state` handle (for a chain route
%   that is (C4) of route_pool, the one place it is coded), the depression
%   factor is its `depression` handle, the synaptic output factor its
%   `output`, and the step response integrates its `rhs` (classical RK4 at
%   cfg.fs, rate constant within each step). A chain has no closed-form
%   step response (its equations are bilinear), so every model is
%   integrated the same way here. Same four panels as below; for a chained
%   route the reserve b_2 is drawn dashed in the condition's colour beside
%   its releasable pool b_1 (the depression factor). The synaptic-output
%   panel marks an interior peak (the product's turnover) with a dot, and a
%   monotone curve's high-rate ceiling (r F(r) at r -> infinity, read off
%   the model's steady state) with a dotted line. What it makes visible:
%   1TS and the chain share the releasable-pool equation and the reserve
%   starts full, so their step responses coincide at onset and separate as
%   the reserve is drawn down; the chain's sustained output is lower and
%   rises to a ceiling, where the product's peaks and falls. Conditions with
%   no depression on the route (no_adaptation) are skipped and listed in
%   out.skipped; the undepressed reference y = r is drawn in grey.
%
% * THE PRODUCT FIGURE (product_figure, below) -- no selection AND a plain
%   product route: the figure as it was before the chain existed, kept
%   number for number. Its description follows.
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
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
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
    % Condition names to draw, in order; empty = the product figure for a plain
    % product preset, every depressing condition otherwise (select_conditions).
    cfg.conditions                = {}
    % Rates at which the condition figure marks and reports the steady state
    % (out.steady_state); Plan A section 1 quotes these four.
    cfg.probe_rates (1,:) double  = [0.1 0.25 0.5 1]
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
% THE TWO PATHS (see the header). The product figure below is kept exactly as
% it was for the case it was written for: no selection, and a most-adapted
% route that is a plain product (no model field other than 'product', no
% combine other than 'product', no facilitation). Everything else -- a
% selection, a chain route, a geometric-mean or facilitated route -- is drawn
% by condition_figure from the route models themselves.
entry = sc.(cfg.route_pre).(cfg.route_post);
if is_empty_request(cfg.conditions) && is_plain_product(entry)
    out = product_figure(cfg, entry.std, out_dir, st);
else
    out = condition_figure(cfg, conditions, cond_names, out_dir, st);
end
end

%% ========================================================================
function out = product_figure(cfg, route, out_dir, st)
% The product-of-independent-pools figure, unchanged since before the chain
% existed (its numbers and its file names are what the stages regression
% tier compares). Only reached for a plain product route; see the dispatch.
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
% Created invisible and shown at the end when asked (it used to open a
% window first and hide it after drawing).
fig = figure('Color', 'white', 'Visible', 'off', ...
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
if cfg.visible; set(fig, 'Visible', 'on'); end

fig_tag = 'Fig_STD_steady_state';
out = struct('figs', fig, 'files', {{}}, 'source', ['preset: ' cfg.preset_name], ...
    'rho', rho, 'r_peak', r_peak, 'theta_peak', theta_peak, 'has_peak', has_peak, ...
    'R', theta_peak / theta(end), 'recovery_frac', recovery_frac);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

%% ========================================================================
function out = condition_figure(cfg, conditions, cond_names, out_dir, st)
% One curve per selected condition with depression on the route, every
% number from the route model (see the header).
names = select_conditions(cond_names, cfg.conditions);
label = sprintf('%s->%s', cfg.route_pre, cfg.route_post);
[preset, ~] = srnn_param_preset(cfg.preset_name);
zero_floor = isfield(preset, 'std_zero_floor') && preset.std_zero_floor;

r  = linspace(0, 1, 4000)';         % rate, over the full range of the nonlinearity
rf = linspace(0, 1, 4e5)';          % fine grid for the peak search
r_inf = 1e9;                        % "infinity" for the high-rate ceiling

% Step protocol, as the product figure's: settle, then on/off per step rate.
seg_rate = [0, reshape([cfg.step_rates; zeros(1, numel(cfg.step_rates))], 1, [])];
seg_dur  = [cfg.settle_s, repmat([cfg.on_s, cfg.off_s], 1, numel(cfg.step_rates))];
seg_n    = round(seg_dur * cfg.fs);
dt       = 1 / cfg.fs;
r_t      = repelem(seg_rate, seg_n);
t        = (1:sum(seg_n)) * dt;
seg_end  = cumsum(seg_n);           % sample index at the end of each segment
seg_start = [0, cumsum(seg_dur)];

C = struct('name', {}, 'title', {}, 'color', {}, 'model', {}, 'P', {}, 'params', {}, ...
    'has_reserve', {}, 'D_ss', {}, 'b2_ss', {}, 'theta_ss', {}, 'has_peak', {}, ...
    'r_peak', {}, 'theta_peak', {}, 'ceiling', {}, 'probe', {}, 's_t', {}, ...
    'D_t', {}, 'theta_t', {}, 'step_ss', {});
skipped = {};
for k = 1:numel(names)
    cnd = conditions{strcmp(cond_names, names{k})};
    entry = route_entry(cnd, cfg.route_pre, cfg.route_post);
    if isempty(entry)
        skipped{end+1} = names{k}; %#ok<AGROW>
        continue
    end
    stf = [];
    if isfield(entry, 'stf'); stf = entry.stf; end
    route_cfg = struct('std', entry.std, 'stf', stf);
    opts = struct('npre', 1, 'label', label, 'zero_floor', zero_floor, ...
        'err_id_prefix', 'fig_STD_steady_state');
    R1 = SRNNCellTypePairs.compile_route(route_cfg, opts);          % one neuron: the step response
    opts.npre = numel(r);
    Rr = SRNNCellTypePairs.compile_route(route_cfg, opts);          % one "neuron" per rate
    opts.npre = numel(rf);
    Rf = SRNNCellTypePairs.compile_route(route_cfg, opts);
    model = R1.model;  P = R1.P;

    c = struct();
    c.name   = names{k};
    title_k  = titles_by_name(names(k));
    c.title  = title_k{1};
    c.color  = condition_color(st, names{k}, k);
    c.model  = P.model;
    c.P      = P;
    c.params = param_string(P);
    % Only the vesicle-pool model has a reserve: its second b column.
    c.has_reserve = strcmp(P.model, 'pool');

    % Steady state, from the model's own closed form.
    ss = model.steady_state(r, Rr.P);
    c.D_ss = model.depression(ss.s(:, 1:P.n_b), Rr.P);
    c.b2_ss = [];
    if c.has_reserve; c.b2_ss = ss.s(:, 2); end
    c.theta_ss = r .* ss.F;
    ssf = model.steady_state(rf, Rf.P);
    thf = rf .* ssf.F;
    [c.theta_peak, ip] = max(thf);
    c.r_peak   = rf(ip);
    c.has_peak = ip < numel(rf);
    ss_inf = model.steady_state(r_inf, P);
    c.ceiling = r_inf * ss_inf.F;           % 0 where theta falls to zero (the product)
    if c.has_peak || c.ceiling < 1e-6; c.ceiling = NaN; end
    ssp = model.steady_state(cfg.probe_rates(:), compile_n(route_cfg, opts, numel(cfg.probe_rates)));
    c.probe = struct('r', cfg.probe_rates, ...
        'D', reshape(model.depression(ssp.s(:, 1:P.n_b), compile_n(route_cfg, opts, numel(cfg.probe_rates))), 1, []), ...
        'theta', reshape(cfg.probe_rates(:) .* ssp.F, 1, []), ...
        'b2', []);
    if c.has_reserve; c.probe.b2 = reshape(ssp.s(:, 2), 1, []); end

    % Step response: integrate the route's own equations.
    c.s_t = integrate_route(model, P, r_t, dt);
    s3 = reshape(c.s_t.', 1, P.ns, []);                            % 1 x ns x nt
    c.D_t = reshape(model.depression(s3(:, 1:P.n_b, :), P), 1, []);
    c.theta_t = reshape(model.output(s3, P), 1, []) .* r_t;
    sss = model.steady_state(cfg.step_rates(:), compile_n(route_cfg, opts, numel(cfg.step_rates)));
    c.step_ss = struct('D', reshape(model.depression(sss.s(:, 1:P.n_b), ...
        compile_n(route_cfg, opts, numel(cfg.step_rates))), 1, []), ...
        'theta', reshape(cfg.step_rates(:) .* sss.F, 1, []));
    C(end+1) = c; %#ok<AGROW>
end
if isempty(C)
    error('fig_STD_steady_state:NoDepression', ...
        'None of the conditions {%s} of preset ''%s'' carries depression on route %s.', ...
        strjoin(names, ', '), cfg.preset_name, label);
end
if ~isempty(skipped)
    vprintf(cfg.verbose, 'verbose', '[fig_STD_steady_state] no depression on %s, not drawn: %s\n', ...
        label, strjoin(skipped, ', '));
end

% Recovery between steps, measured on the integrated traces: the smallest
% fraction of the depression recovered by the end of an off period.
recovery = inf;
for k = 1:numel(C)
    for j = 1:numel(cfg.step_rates)
        d_on  = C(k).D_t(seg_end(2*j));
        d_off = C(k).D_t(seg_end(2*j + 1));
        recovery = min(recovery, (d_off - d_on) / (1 - d_on));
    end
end

%% ---- Figure ---------------------------------------------------------------
ref_color = [0.62 0.62 0.62];       % undepressed reference / the rate itself
tick_fs = st.tick_fs; label_fs = st.label_fs; title_fs = 15; lw = st.line_lw;
fig_size = [900, 950];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', 'Visible', 'off', ...
    'Position', [scr(1:2) + max((scr(3:4) - fig_size)/2, 0), fig_size]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
head = cell(1, numel(C));
for k = 1:numel(C)
    head{k} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%s:  %s', C(k).color, ...
        tex_safe(C(k).title), C(k).params);
end
title(tl, head, 'FontSize', 12, 'Interpreter', 'tex');
reserve_note = any([C.has_reserve]);

% --- 1: steady-state depression factor ---
ax1 = nexttile(tl); hold(ax1, 'on');
h1 = gobjects(1, numel(C));
for k = numel(C):-1:1                % back to front: the first condition on top
    if C(k).has_reserve
        plot(ax1, r, C(k).b2_ss, '--', 'LineWidth', 1, 'Color', C(k).color);
    end
    h1(k) = plot(ax1, r, C(k).D_ss, 'LineWidth', lw, 'Color', C(k).color);
    plot(ax1, C(k).probe.r, C(k).probe.D, 'o', 'MarkerSize', 4, ...
        'MarkerFaceColor', C(k).color, 'MarkerEdgeColor', 'none');
end
box(ax1, 'off'); set(ax1, 'FontSize', tick_fs, 'XTick', [0 1], 'YTick', [0 1]);
xlim(ax1, [0 1]); ylim(ax1, [0 1]);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'depression factor', 'FontSize', label_fs);
title(ax1, 'Steady-state depression', 'FontWeight', 'normal', 'FontSize', title_fs);
leg_h = h1; leg_t = tex_safe_all({C.title});
if reserve_note
    leg_h(end+1) = plot(ax1, NaN, NaN, '--', 'LineWidth', 1, 'Color', [0.3 0.3 0.3]);
    leg_t{end+1} = 'reserve b_2 (chain)';
end
legend(ax1, leg_h, leg_t, 'Box', 'off', 'FontSize', 9, 'Location', 'northeast');

% --- 2: steady-state synaptic output ---
ymax2 = max([arrayfun(@(c) max(c.theta_ss), C), [C.ceiling]]) * 1.25;
ax2 = nexttile(tl); hold(ax2, 'on');
h_id = plot(ax2, [0 1], [0 1], '-', 'LineWidth', 1, 'Color', ref_color);
for k = numel(C):-1:1
    plot(ax2, r, C(k).theta_ss, 'LineWidth', lw, 'Color', C(k).color);
    if C(k).has_peak
        plot(ax2, C(k).r_peak, C(k).theta_peak, 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', C(k).color, 'MarkerEdgeColor', 'none');
        text(ax2, C(k).r_peak, C(k).theta_peak, sprintf('  peak (%.2f, %.3f)', ...
            C(k).r_peak, C(k).theta_peak), 'FontSize', 10, 'Color', C(k).color, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    elseif ~isnan(C(k).ceiling)
        plot(ax2, [0 1], C(k).ceiling * [1 1], ':', 'LineWidth', 1.2, 'Color', C(k).color);
        text(ax2, 1, C(k).ceiling, sprintf('ceiling %.3f ', C(k).ceiling), 'FontSize', 10, ...
            'Color', C(k).color, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs, 'XTick', [0 1], 'YTick', [0 round(ymax2, 3)]);
xlim(ax2, [0 1]); ylim(ax2, [0 ymax2]);
xlabel(ax2, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax2, 'synaptic output  \theta = r F', 'FontSize', label_fs);
title(ax2, 'Steady-state synaptic output', 'FontWeight', 'normal', 'FontSize', title_fs);
legend(ax2, h_id, {'undepressed  y = r'}, 'Box', 'off', 'FontSize', 9, 'Location', 'northeast');

% --- 3 (tiles 3-4): step response, depression ---
ax3 = nexttile(tl, 3, [1 2]); hold(ax3, 'on');
h3r = plot(ax3, t, r_t, '-', 'LineWidth', 1.5, 'Color', ref_color);
for k = numel(C):-1:1
    for j = 1:numel(cfg.step_rates)
        plot(ax3, seg_start([2*j 2*j+1]), C(k).step_ss.D([j j]), ':', ...
            'LineWidth', 1.2, 'Color', C(k).color);
    end
    if C(k).has_reserve
        plot(ax3, t, C(k).s_t(:, 2), '--', 'LineWidth', 1, 'Color', C(k).color);
    end
    plot(ax3, t, C(k).D_t, 'LineWidth', lw, 'Color', C(k).color);
end
h3s = plot(ax3, NaN, NaN, ':', 'LineWidth', 1.2, 'Color', [0.3 0.3 0.3]);
box(ax3, 'off'); set(ax3, 'FontSize', tick_fs);
xlim(ax3, [0 t(end)]); ylim(ax3, [0 1.05]);
xlabel(ax3, 'time (s)', 'FontSize', label_fs);
ylabel(ax3, 'depression factor', 'FontSize', label_fs);
title(ax3, 'Step response: depression', 'FontWeight', 'normal', 'FontSize', title_fs);
legend(ax3, [h3r h3s], {'rate  r(t)', 'steady state'}, 'Box', 'off', 'FontSize', 9, ...
    'Location', 'southeast');

% --- 4 (tiles 5-6): step response, synaptic output ---
ax4 = nexttile(tl, 5, [1 2]); hold(ax4, 'on');
plot(ax4, t, r_t, '-', 'LineWidth', 1.5, 'Color', ref_color);
for k = numel(C):-1:1
    for j = 1:numel(cfg.step_rates)
        plot(ax4, seg_start([2*j 2*j+1]), C(k).step_ss.theta([j j]), ':', ...
            'LineWidth', 1.2, 'Color', C(k).color);
    end
    plot(ax4, t, C(k).theta_t, 'LineWidth', lw, 'Color', C(k).color);
end
box(ax4, 'off'); set(ax4, 'FontSize', tick_fs);
xlim(ax4, [0 t(end)]); ylim(ax4, [0 max([C.theta_t, r_t]) * 1.12]);
xlabel(ax4, 'time (s)', 'FontSize', label_fs);
ylabel(ax4, 'synaptic output  \theta = r F', 'FontSize', label_fs);
title(ax4, 'Step response: synaptic output', 'FontWeight', 'normal', 'FontSize', title_fs);
linkaxes([ax3, ax4], 'x');

%% ---- Save -----------------------------------------------------------------
if cfg.visible; set(fig, 'Visible', 'on'); end
fig_tag = 'Fig_STD_steady_state';
summary = rmfield(C, {'P', 's_t', 'D_t', 'theta_t', 'D_ss', 'b2_ss', 'theta_ss'});
out = struct('figs', fig, 'files', {{}}, 'source', ['preset: ' cfg.preset_name], ...
    'conditions', {{C.name}}, 'skipped', {skipped}, 'steady_state', summary, ...
    'recovery_frac', recovery, 't', t, 'r_t', r_t, ...
    'traces', struct('name', {C.name}, 's', {C.s_t}, 'D', {C.D_t}, 'theta', {C.theta_t}));
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

%% ------------------------------------------------------------------------
% integrate_route (classical RK4 on the route's rhs, the rate held at r_t(i)
% over step i; S nt x ns, the state at the END of each step) is shared with
% fig_std_1TS_vs_MTS: src/figures/helpers/integrate_route.m.

function P = compile_n(route_cfg, opts, n)
% The route compiled for n neurons (n rates evaluated at once).
opts.npre = n;
R = SRNNCellTypePairs.compile_route(route_cfg, opts);
P = R.P;
end

function entry = route_entry(cnd, pre, post)
% The condition's synapse_config entry for the route, or [] when the route
% carries no depression in this condition.
entry = [];
if ~isfield(cnd, 'synapse_config') || ~isstruct(cnd.synapse_config); return; end
sc = cnd.synapse_config;
if ~isfield(sc, pre) || ~isfield(sc.(pre), post); return; end
e = sc.(pre).(post);
if isfield(e, 'std') && ~isempty(e.std) && ~(isstruct(e.std) && isempty(fieldnames(e.std)))
    entry = e;
end
end

function s = param_string(P)
% The route's configured values, from the compiled route's nominal (TeX).
n = P.nominal;
parts = {};
add = {'tau_rec', '\tau_{rec}', ' s'; 'tau_rel', '\tau_{rel}', ' s'; 'n2', 'n_2', ''; ...
       'tau_dec', '\tau_{dec}', ' s'; 'tau_fac', '\tau_{fac}', ' s'; 'G', 'G', ''};
for k = 1:size(add, 1)
    if isfield(n, add{k, 1}) && ~isempty(n.(add{k, 1}))
        v = n.(add{k, 1});
        if isscalar(v)
            vs = sprintf('%g', v);
        else
            vs = ['[' strjoin(compose('%g', v), ' ') ']'];
        end
        parts{end+1} = sprintf('%s = %s%s', add{k, 2}, vs, add{k, 3}); %#ok<AGROW>
    end
end
if strcmp(P.model, 'pool')
    model_word = 'chained pools';
elseif P.n_b == 1
    model_word = 'one pool';
else
    model_word = sprintf('product of %d pools', P.n_b);
end
s = sprintf('%s,  %s', model_word, strjoin(parts, ',  '));
end

function col = condition_color(st, name, k)
if isKey(st.condition_color, name)
    col = st.condition_color(name);
else
    L = lines(max(k, 7));
    col = L(k, :);
end
end

function s = tex_safe(s)
% Condition titles are TeX already (srnn_condition_titles writes \tau); only
% an underscore from a raw-name fallback needs escaping.
s = strrep(s, '_', '\_');
end

function c = tex_safe_all(c)
c = cellfun(@tex_safe, c, 'UniformOutput', false);
end

function tf = is_empty_request(r)
% {}, '', "" and [] mean "no selection", as in select_conditions.
tf = isempty(r) || ((ischar(r) || isstring(r)) && all(strlength(r) == 0));
end

function tf = is_plain_product(entry)
% A route the product figure draws correctly: product depression (no model
% field, or 'product'), the plain product combine, and no facilitation.
tf = isfield(entry, 'std') && isstruct(entry.std) && ...
    (~isfield(entry, 'stf') || isempty(entry.stf)) && ...
    (~isfield(entry.std, 'model') || strcmp(char(string(entry.std.model)), 'product')) && ...
    (~isfield(entry.std, 'combine') || isempty(entry.std.combine) || ...
     strcmp(char(string(entry.std.combine)), 'product'));
end

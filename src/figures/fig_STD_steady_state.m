function out = fig_STD_steady_state(cfg)
% FIG_STD_STEADY_STATE Multi-timescale STD steady state, and what it costs.
%
%   out = FIG_STD_STEADY_STATE()
%   out = FIG_STD_STEADY_STATE('preset_name', p)
%
% Conceptual, ANALYTIC (no simulation). Setting db/dt = 0 in
%     db_k/dt = (1 - b_k)/tau_rec_k - b_k*r/tau_rel_k
% gives b_k(r) = 1/(1 + r*tau_rec_k/tau_rel_k); the synapse multiplies the
% timescales, so it sees prod_k b_k(r).
%
%   Left  : prod(b) against r, with one component b_k dashed for reference.
%   Right : prod(b)*r -- what the recurrent sum actually receives -- again with
%           one component dashed, and a faint identity line y = r marking the
%           UNDEPRESSED synapse. The vertical gap to a curve is what depression
%           costs at that rate.
%
% TWO VERSIONS are written, differing ONLY in the right panel's y range. The
% full-scale one keeps equal pixel scaling on both axes, so the identity line is
% a literal 45 degrees and a vertical distance can be compared against a
% horizontal one by eye -- honest, but both curves sit in the bottom tenth of
% the panel. The zoom drops the 1:1 aspect deliberately (it cannot survive a y
% range that much shorter than x without squashing the panel to a sliver, so no
% cross-axis angle there means anything) and buys the turnover, which is the
% point of the panel and is nearly invisible at full scale.
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
% See also: fig_SFA_steady_state, srnn_param_preset, srnn_adaptation_conditions

arguments
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.route_pre   (1,:) char    = 'E'
    cfg.route_post  (1,:) char    = 'E'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

% Pull the depression timescales out of the preset's own STD routes. They live
% on the conditions, not on the model_defaults struct: synapse_config can only
% reach the model through a condition, so that is where a preset puts them.
[~, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
sc = conditions{strcmp(cond_names, 'sfa_and_std')}.synapse_config;
route = sc.E.E.std;                 % uniform across all four routes in this preset
tau_rec = route.tau_rec(:)';
tau_rel = route.tau_rel(:)';
ratio = tau_rec ./ tau_rel;         % the only combination that sets the steady state

r = linspace(0, 1, 400);            % rate, over the full range of the nonlinearity

% b_k(r) for each timescale, then their product. rows = timescale, cols = r.
b_each = 1 ./ (1 + ratio(:) * r);
b_prod = prod(b_each, 1);

tick_fs  = st.tick_fs;
label_fs = st.label_fs;
title_fs = 16;    % panel titles
lw       = st.line_lw;

prod_color     = [0.85 0.325 0.098];   % warm, matching the E colour used elsewhere
single_color   = [0.5 0.5 0.5];
identity_color = [0.55 0.80 0.55];     % faint green: y = r, the undepressed synapse

%% ---- Derived quantities ---------------------------------------------------
output = b_prod .* r;
output_single = b_each(1, :) .* r;
[peak_val, peak_idx] = max(output);
r_peak = r(peak_idx);

% The zoomed version's ceiling: just above the largest value a SINGLE b_k ever
% delivers, which on this preset is its value at r = 1 (b_k*r rises
% monotonically to the asymptote tau_rel/tau_rec, so the right edge is the max).
% Rounded up to the next hundredth so the tick is a round number rather than
% 0.1111..., and so the curve does not touch the top of the axes.
zoom_ymax = ceil(max(output_single) * 100) / 100;

%% ---- Figures --------------------------------------------------------------
% TWO VERSIONS of the same figure, differing only in the right panel's y range:
%
%   'square' : y over [0, 1] with daspect [1 1 1]. The identity line y = r is a
%              literal 45 degrees and a vertical distance can be compared
%              against a horizontal one by eye. Honest but small: both curves
%              sit in the bottom tenth of the panel.
%   'zoom'   : y over [0, zoom_ymax], filling the panel with the curves
%              themselves. The 1:1 aspect ratio is DELIBERATELY DROPPED here --
%              it cannot survive a y range 9x shorter than x's without squashing
%              the panel to a sliver, so the identity line is no longer 45
%              degrees and no cross-axis angle in this version means anything.
%              What it buys is the turnover, which is the point of the panel and
%              is nearly invisible at full scale.
%
% Size copied from Fig_FI_curve so the two figures match, but the POSITION is
% computed rather than copied: that script hardcodes x = 4429, a second-monitor
% coordinate from the machine it was written on, which lands the window
% off-screen on a single 1920-wide display and looks like the figure never
% opened. Centre it on whatever screen is actually attached.
variants = {'square', 'zoom'};
fig_size = [623, 322];
scr = get(groot, 'ScreenSize');
figs = gobjects(1, numel(variants));

for v = 1:numel(variants)
is_zoom = strcmp(variants{v}, 'zoom');

% Offset each window so the second does not land exactly on the first.
fig_pos = [scr(1:2) + max((scr(3:4) - fig_size)/2, 0) + (v - 1)*[40, -40], fig_size];
fig = figure('Color', 'white', 'Position', fig_pos);
figs(v) = fig;
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Left panel: prod(b) vs r --------------------------------------------------
% Identical in both versions: prod(b) is a gain, already bounded by [0, 1], so
% there is nothing to zoom into.
ax1 = nexttile(tl); hold(ax1, 'on');
% One component b_k, thin and dashed, so the gap to the product reads as the
% cost of the SECOND timescale rather than as depression in general. Both
% timescales share tau_rec/tau_rel here, so the components coincide -- drawing
% one is drawing both (the code plots b_1).
plot(ax1, r, b_each(1, :), '--', 'LineWidth', 1, 'Color', single_color);
plot(ax1, r, b_prod, 'LineWidth', lw, 'Color', prod_color);
hold(ax1, 'off');
box(ax1, 'off');
set(ax1, 'FontSize', tick_fs);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'depression  $\prod_k b_k$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax1, 'Steady-state depression', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax1, [0, 1]); ylim(ax1, [0, 1]);
set(ax1, 'XTick', [0, 1], 'YTick', [0, 1]);
legend(ax1, {'single $b_k$', '$\prod_k b_k$'}, 'Interpreter', 'latex', ...
    'Box', 'off', 'FontSize', 12, 'Location', 'northeast');

%% Right panel: prod(b)*r vs r -----------------------------------------------
ax2 = nexttile(tl); hold(ax2, 'on');
% Identity first, so it sits UNDER the data. y = r is the undepressed synapse
% (b == 1): the vertical gap down to each curve is exactly what depression costs
% at that rate. In the zoomed version it leaves the top of the panel almost
% immediately -- that steep exit IS the depression, seen edge-on.
plot(ax2, [0, 1], [0, 1], '-', 'LineWidth', 1, 'Color', identity_color);
% Same pairing as the left panel: one timescale dashed, the product solid. The
% two curves differ in KIND here, not just in scale -- a single b_k rises
% monotonically toward its asymptote 1/(tau_rec/tau_rel), so more rate always
% delivers more, whereas the product turns over.
plot(ax2, r, output_single, '--', 'LineWidth', 1, 'Color', single_color);
plot(ax2, r, output, 'LineWidth', lw, 'Color', prod_color);
plot(ax2, r_peak, peak_val, 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', prod_color, 'MarkerEdgeColor', 'none');
hold(ax2, 'off');
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs);
xlabel(ax2, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax2, 'synaptic output  $\prod_k b_k \cdot r$', 'Interpreter', 'latex', ...
    'FontSize', label_fs);
title(ax2, 'Delivered output', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax2, [0, 1]);
if is_zoom
    % No daspect: see the header comment above. The y tick is zoom_ymax rather
    % than 1, which is the one number a reader must notice to avoid reading this
    % panel as if it were the square one.
    ylim(ax2, [0, zoom_ymax]);
    set(ax2, 'XTick', [0, 1], 'YTick', [0, zoom_ymax]);
else
    % Equal pixel scaling on both axes, which forces y to share x's [0, 1]
    % range: the curves then sit low in a square panel, which is the honest
    % picture -- at every rate the depressed synapse delivers a small fraction
    % of what an undepressed one would.
    ylim(ax2, [0, 1]);
    daspect(ax2, [1 1 1]);
    set(ax2, 'XTick', [0, 1], 'YTick', [0, 1]);
end

% Mark where the output peaks -- past this rate the synapse delivers less.
% Above-right of the marker: there is headroom there in the square version
% (the title is far off at y = 1) and in the zoomed one (the peak sits below
% the single-b_k curve it is being contrasted with).
text(ax2, r_peak, peak_val, sprintf('  r = %.2f', r_peak), ...
    'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

end

%% ---- Save (stable filenames) ----------------------------------------------

if ~cfg.visible; set(figs, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag_base = 'Fig_STD_steady_state';
fig_tags = {fig_tag_base, [fig_tag_base '_zoom']};
out = struct('figs', figs, 'files', {{}}, 'source', ['preset: ' cfg.preset_name]);
if cfg.save
    for v = 1:numel(figs)
        save_figure_stable(out_dir, fig_tags{v}, figs(v));
        out.files = [out.files, existing_outputs(out_dir, fig_tags{v})];
    end

end
end

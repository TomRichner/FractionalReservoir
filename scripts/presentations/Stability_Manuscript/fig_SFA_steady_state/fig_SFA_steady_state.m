function out = fig_SFA_steady_state(cfg)
% FIG_SFA_STEADY_STATE SFA steady state: one timescale against three.
%
%   out = FIG_SFA_STEADY_STATE()
%   out = FIG_SFA_STEADY_STATE('preset_name', p)
%
% Conceptual, ANALYTIC (no simulation). Setting da/dt = 0 in
%     da_k/dt = (-a_k + r)/tau_k
% gives a_k = r for EVERY timescale, so c*sum(a) = c*n_a*r.
%
%   Left  : c*sum(a) against r, for the preset's n_a timescales and for one.
%   Right : the transient response to a step r: 0 -> 1.
%
% THE RESULT: THE TWO LINES IN THE LEFT PANEL COINCIDE. Because a_k = r
% independently of tau_k, the steady-state SFA current depends only on the
% PRODUCT c*n_a -- and splitting c as a BUDGET over the timescales holds that
% product fixed. This is why the project uses c = budget/n_a: adding timescales
% does not silently move the operating point.
%
% WHERE THE COUNT DOES SHOW UP: THE TRANSIENT. Several timescales give a
% multi-exponential approach -- fast partial adaptation then a slow tail --
% against a single exponential. Same destination, different route.
%
% PARAMETERS COME FROM THE PRESET, not from literals. The BUDGET is recovered as
% c(1)*n_a(1): the preset stores the ALREADY-DIVIDED value (c = 0.5/3), not the
% budget, so it has to be multiplied back by the SFA count that the conditions
% declare. tau_a is NOT in the preset at all -- it is the class default,
% auto-filled per n_a -- so it is read off a BUILT model rather than off the
% preset struct.
%
% THE SINGLE-TIMESCALE CASE uses tau(1) and c = budget/1. The model's own rule,
% logspace(log10(0.25), log10(10), n_a), is deliberately NOT used at n_a = 1:
% logspace(a, b, 1) returns 10^b, the SLOW end, so the comparison would be
% against the slowest timescale alone and would be a statement about the slow
% tail rather than about the count.
%
% CONTRAST WITH STD (see fig_STD_steady_state): SFA enters the dynamics as a
% SUM, so a budget-split c makes the count invisible at steady state. STD enters
% as a PRODUCT, so n_b = 2 squares the depression however the taus are chosen --
% there is no budget-split that would make dual STD match single STD.
%
% See also: fig_STD_steady_state, srnn_param_preset, manuscript_style

arguments
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

%% Timescales and the c budget, read off the preset
[tau_3, c_budget, n_a] = sfa_from_preset(cfg.preset_name);
tau_1 = tau_3(1);          % the FAST component, not logspace(...,1); see header

c_3 = c_budget / numel(tau_3);
c_1 = c_budget / numel(tau_1);

fprintf('[fig_SFA_steady_state] preset=%s  n_a=%d  tau_a=%s  budget=%.4g\n', ...
    cfg.preset_name, n_a, mat2str(tau_3, 4), c_budget);

r = linspace(0, 1, 400);       % rate, over the full range of the nonlinearity

% Steady state: a_k = r for every k, so sum(a) = n_a * r.
sfa_3 = c_3 * numel(tau_3) * r;
sfa_1 = c_1 * numel(tau_1) * r;

r_step = 1;                     % rate step for the transient panel: the full range
                               % of the nonlinearity, so the y-axes of the two
                               % panels share a scale and the transient lands on
                               % the left panel's endpoint.
t = linspace(0, 20, 800);

tick_fs  = st.tick_fs;
label_fs = st.label_fs;
title_fs = 16;    % panel titles
lw       = st.line_lw;

three_color = [0.85 0.325 0.098];   % warm, matching the E colour used elsewhere
one_color   = [0.5 0.5 0.5];

%% ---- Figure ---------------------------------------------------------------
% Position computed rather than hardcoded: Fig_FI_curve's [4429 ...] is a
% second-monitor coordinate and lands off-screen on a single display.
fig_size = [623, 322];
scr = get(groot, 'ScreenSize');
fig = figure('Color', 'white', ...
    'Position', [scr(1:2) + max((scr(3:4) - fig_size)/2, 0), fig_size]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Left panel: steady-state SFA current vs rate ------------------------------
ax1 = nexttile(tl); hold(ax1, 'on');
% Drawn thick-then-dashed BECAUSE they coincide exactly. A single line would
% look like one case was forgotten; the dashed overlay is the evidence that both
% were computed and landed on top of each other.
plot(ax1, r, sfa_3, 'LineWidth', lw + 1, 'Color', three_color);
plot(ax1, r, sfa_1, '--', 'LineWidth', 1.2, 'Color', one_color);
hold(ax1, 'off');
box(ax1, 'off');
set(ax1, 'FontSize', tick_fs);
xlabel(ax1, 'firing rate  r', 'FontSize', label_fs);
ylabel(ax1, 'SFA current  c\Sigmaa', 'FontSize', label_fs);
title(ax1, 'Steady-state SFA', 'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax1, [0, 1]); ylim(ax1, [0, 1]);
set(ax1, 'XTick', [0, 1], 'YTick', [0, 1]);
legend(ax1, ...
    {sprintf('n_a = 3,  c = %.3g', c_3), sprintf('n_a = 1,  c = %.3g', c_1)}, ...
    'Box', 'off', 'FontSize', 12, 'Location', 'northeast');

%% Right panel: the transient, where the count DOES matter ------------------
% Step r: 0 -> r_step at t = 0, from a = 0. Then a_k(t) = r_step*(1 - e^{-t/tau_k}).
sfa_3_t = c_3 * r_step * sum(1 - exp(-t ./ tau_3(:)), 1);
sfa_1_t = c_1 * r_step * sum(1 - exp(-t ./ tau_1(:)), 1);

ax2 = nexttile(tl); hold(ax2, 'on');
plot(ax2, t, sfa_3_t, 'LineWidth', lw, 'Color', three_color);
plot(ax2, t, sfa_1_t, '--', 'LineWidth', 1.2, 'Color', one_color);
yline(ax2, c_3 * numel(tau_3) * r_step, ':', 'Color', [0.6 0.6 0.6]);
hold(ax2, 'off');
box(ax2, 'off');
set(ax2, 'FontSize', tick_fs);
xlabel(ax2, 'time (s)', 'FontSize', label_fs);
ylabel(ax2, 'SFA current  c\Sigmaa', 'FontSize', label_fs);
title(ax2, sprintf('Transient (step to r = %.2g)', r_step), ...
    'FontWeight', 'normal', 'FontSize', title_fs);
xlim(ax2, [0, t(end)]);
% Same [0, 1] range as the left panel, so the two are directly comparable and
% the transient visibly lands on the left panel's endpoint.
ylim(ax2, [0, 1]);
set(ax2, 'XTick', [0, t(end)], 'YTick', [0, 1]);

%% ---- Save (stable filenames) ----------------------------------------------

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_SFA_steady_state';
out = struct('figs', fig, 'files', {{}}, 'source', ['preset: ' cfg.preset_name]);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    capture_git_provenance(out_dir, project_root);

    result_note = sprintf([ ...
        'THE TWO LINES IN THE LEFT PANEL ARE IDENTICAL: both are %.4g*r exactly. ' ...
        'Because a_k = r independently of tau_k, the steady-state SFA current ' ...
        'depends only on the PRODUCT c*n_a, and splitting c as a budget holds ' ...
        'that product fixed. The dashed line is drawn over the solid one to show ' ...
        'both were computed; they coincide to machine precision. This is why the ' ...
        'project splits c over its timescales -- adding timescales does not ' ...
        'silently move the operating point.'], c_budget);

    transient_note = sprintf([ ...
        'WHERE THE COUNT SHOWS UP: THE TRANSIENT. tau_a (n_a = %d) = %s against ' ...
        'a single tau_a = %.4g, which is the FAST component of that set. The ' ...
        'model rule logspace(log10(0.25), log10(10), n_a) is deliberately not ' ...
        'used at n_a = 1: logspace(a, b, 1) returns 10^b, the slow end, so the ' ...
        'comparison would be against the slowest timescale alone. Several ' ...
        'timescales give a multi-exponential approach -- fast partial adaptation ' ...
        'then a slow tail -- against a single exponential. Both settle to %.4f.'], ...
        numel(tau_3), mat2str(tau_3, 4), tau_1, c_budget);

    contrast_note = [ ...
        'CONTRAST WITH STD (fig_STD_steady_state). SFA enters the dynamics as a ' ...
        'SUM, so a budget-split c makes the count invisible at steady state. STD ' ...
        'enters as a PRODUCT, so two timescales square the depression however ' ...
        'the taus are chosen -- there is no budget-split of tau that would make ' ...
        'dual STD match single STD.'];

    provenance_note = [ ...
        'PARAMETERS COME FROM THE PRESET. The c BUDGET is recovered as ' ...
        'c(1)*n_a(1), because the preset stores the already-divided value, not ' ...
        'the budget. tau_a is not in the preset at all -- it is the class ' ...
        'default, auto-filled per n_a -- so it is read off a BUILT model.'];

    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_SFA_steady_state', ...
        'title',  'Stability_Manuscript figure: SFA steady state, one timescale vs several', ...
        'script', 'fig_SFA_steady_state.m', ...
        'what',   ['Conceptual, ANALYTIC figure (no simulation). Setting ' ...
                   'da_k/dt = 0 in da_k/dt = (-a_k + r)/tau_k gives a_k = r for ' ...
                   'every timescale, so c*sum(a) = c*n_a*r. Left: c*sum(a) ' ...
                   'against r, for the preset''s timescale count and for one. ' ...
                   'Right: the transient response to a step r: 0 -> 1.'], ...
        'source', struct('preset', cfg.preset_name), ...
        'settings', struct('tau_a', tau_3, 'n_a', numel(tau_3), ...
                           'c_budget', c_budget, 'c_per_timescale', c_3, ...
                           'tau_single', tau_1, 'c_single', c_1), ...
        'figures', {out.files}, ...
        'sections', struct( ...
            'heading', {'the result', 'the transient', 'contrast with std', 'provenance'}, ...
            'body',    {result_note, transient_note, contrast_note, provenance_note})));
end
end

%% ------------------------------------------------------------------------
function [tau_a, c_budget, n_a] = sfa_from_preset(preset_name)
% The SFA timescales and the c BUDGET a preset actually runs at.
%
% Two things have to be recovered rather than read:
%
%   * The BUDGET. A preset stores c already divided by the timescale count
%     (c = 0.5/3), so the budget is c(1)*n_a(1). n_a comes from the sfa_and_std
%     CONDITION, since adaptation counts are condition-owned and never live in
%     a preset.
%   * tau_a. It is usually absent from the preset, being the class default
%     auto-filled per n_a inside build(). Reading it off a BUILT model is the
%     only way to get the values that would actually be integrated.
model = build_from_preset(preset_name, 'sfa_and_std');

if isprop(model, 'tau_a')            % SRNNCellTypePairs: 1 x C cell
    tau_a = model.tau_a{1};
    c_vec = model.c;
    n_a_v = model.n_a;
    n_a   = n_a_v(1);
    c_budget = c_vec(1) * n_a;
else                                  % SRNNModel2
    tau_a = model.tau_a_E;
    n_a   = model.n_a_E;
    c_budget = model.c_E * n_a;
end

assert(n_a >= 1, 'The sfa_and_std condition of ''%s'' switches on no SFA.', preset_name);
end

function out = fig_sine_paired_stim(cfg)
% FIG_SINE_PAIRED_STIM How SFA and STD make the evoked response lead a sine-modulated pulse train.
%
%   out = FIG_SINE_PAIRED_STIM('run_dir', d)                  % d holds stim_evoked/
%   out = FIG_SINE_PAIRED_STIM('run_dir', d, 'out_dir', o)
%
% The mechanism figure for the Minnesota Neuromodulation poster (TR,
% 2026-10-20), from the simulations of run_stim_evoked
% (sine_paired_stim_config): the paper's network driven by the human SEEG
% protocols at a fixed set of E neurons, read out at the UNDRIVEN neurons.
%
%   (A) The multiple-timescale network (sfa3_std2) under the 2, 5 and 10 s
%       sine trains, one column per period. Top: the stimulus envelope with
%       its pulses (left axis, % of peak) and the cycle-averaged response size
%       per pulse with the fitted sinusoid (right axis, the human figure's
%       layout); dashed verticals at the envelope peak (grey) and the fit's
%       peak (blue), the lead printed. Under it the cycle-averaged population
%       means over one cycle: the undriven neurons' dendritic potential x
%       (the CCEP-like signal, E and I), and for the DRIVEN E neurons, where
%       the pulse acts, their rate r, their synaptic output r prod_m b_m, the
%       SFA feedback (c/K) sum_k a_k and the depression product prod_m b_m.
%       Rows 1-3 print the lead of their per-pulse response size; rows 4-5
%       mark (triangle) and print the lag behind the envelope peak of the SFA
%       feedback's peak and of the depression product's trough (fundamental),
%       i.e. when each gain reduction is strongest. The quantities are those
%       of docs/EquationsParametersDocs/Equations_stability_paper.md.
%   (B) Lead against modulation period for the three conditions (95 %
%       bootstrap intervals over cycles, clipped at lead_ylim) with the human
%       medians (9.0, 6.6, 3.9 deg at 2, 5, 10 s) for reference.
%   (C) Paired pulse in the multiple-timescale network: the undriven mean x
%       after the single pulse and after the 100 ms pair, with the LTI
%       prediction (single response plus its 100 ms-shifted copy) dashed;
%       under it the driven E neurons' SFA feedback and depression over the
%       same window, each scaled to its largest change; right, the ratio
%       second / single response for 2, 10 and 100 ms per condition with the
%       LTI ratio as a black tick on each bar (the human definition:
%       stim_evoked_measures).
%
% Measures: stim_evoked_measures (epochs -20..190 ms, baseline -20..-2 ms,
% RMS 10..60 ms of the cycle average, least-squares sinusoid, lead positive
% when the response peaks first).
%
% LAYOUT IS BY HAND, in inches, not tiledlayout. tiledlayout reserves the
% largest label inset of the grid in EVERY tile, which squeezed the
% one-quantity trace rows to 0.5 in (and to zero at one grid row each). The
% canvas is 17.2 x 10.6 in: a figure larger than the screen is clamped to it
% (this PC's is 17.8 x 11.2 in), visible or not. Drawn invisible, shown
% after saving if asked. Fonts 15 pt (ticks 13): at about 17 in on the
% poster that is the style guide's 14 pt.
%
% out.figs / out.files / out.source as every figure; plus out.leads (table:
% condition, period, lead, CI, r2 of the primary readout), out.ratios
% (table), out.lead_all (every readout) and out.measures.
%
% See also: run_stim_evoked, stim_evoked_measures, sine_paired_stim_config,
%           fig_example_timeseries

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.mat_file    (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''      % unused; accepted for a uniform call
    cfg.show_condition (1,:) char = 'sfa3_std2'
    cfg.discard_s   (1,1) double  = 10
    cfg.human_periods (1,:) double = [2 5 10]
    cfg.human_leads   (1,:) double = [9.0 6.6 3.9]
    cfg.font_size   (1,1) double  = 15
    cfg.lead_ylim   (1,2) double  = [-60 60]   % intervals beyond it are clipped (no adaptation)
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.verbose     (1,:) char    = 'minimal'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st = manuscript_style();
fsz = cfg.font_size;

if isempty(cfg.mat_file)
    mat_file = resolve_data_file('', cfg.run_dir, {fullfile(cfg.run_dir, 'stim_evoked')}, ...
        'stim_evoked_data.mat', 'run_stim_evoked(sine_paired_stim_config())');
else
    mat_file = cfg.mat_file;
end
S = load(mat_file, 'results', 'settings');
M = stim_evoked_measures(S.results, 'discard_s', cfg.discard_s);

conds = S.settings.cond_names;
is_sine = strcmp({M.kind}, 'sine');
periods = unique([M(is_sine).mod_period]);
nP = numel(periods);
tcol = SRNNCellTypePairs.type_colors(2);
show_col = st.condition_color(cfg.show_condition);
grey = [0.2 0.2 0.2];
pick = @(cond, proto) find(strcmp({M.condition}, cond) & strcmp({M.protocol}, proto), 1);

%% Canvas and the hand layout (inches from the lower left)
W = 17.2; H = 10.6;
fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'inches', 'Position', [0.2 0.2 W H]);
mkax = @(l, b, w, h) axes(fig, 'Units', 'inches', 'Position', [l b w h]);
colA_x = 1.25 + (0:nP - 1) * 3.6;   colA_w = 2.75;
top_b = 7.3;  top_h = 2.45;
tr_h = 0.98;  tr_b = 0.65 + (4:-1:0) * (tr_h + 0.33);   % rows 1..5, top to bottom
rx = 12.95;                                              % right block

% ---- (A) the multiple-timescale network under the sine trains ----------------
E_drv = tcol(1, :) * 0.8;
% Greek letters as Unicode characters, never TeX commands: a single \Pi,
% \Sigma, \Delta or \it anywhere in the figure makes the vector SVG export
% (R2026a) outline EVERY text object as paths, so the poster SVG would have no
% editable text. Subscripts (_m) and \color are safe.
PI = char(928); SIG = char(931);
trace_rows = { ...
    'x',   {'I_nonstim', 'E_nonstim'}, {'Undriven', 'mean x'},             'x_nonstim'; ...
    'r',   {'E_stim'},                 {'Driven E', 'rate r'},             'r_stim'; ...
    'syn', {'E_stim'},                 {'Driven E', ['r' PI '_m b_m']},    'syn_stim'; ...
    'sfa', {'E_stim'},                 {'SFA', ['(c/K) ' SIG '_k a_k']},   ''; ...
    'std', {'E_stim'},                 {'STD', [PI '_m b_m']},             ''};
group_color = struct('E_stim', E_drv, 'E_nonstim', tcol(1, :), 'I_nonstim', tcol(2, :));
ax_tr = gobjects(size(trace_rows, 1), nP);
ax_top = gobjects(1, nP);
for p = 1:nP
    Tm = periods(p);
    k = pick(cfg.show_condition, sprintf('sine_%gs', Tm));
    m = M(k); R = S.results(k); f = m.fits(1);

    % Stimulus (left) and response size (right)
    ax = mkax(colA_x(p), top_b, colA_w, top_h); ax_top(p) = ax;
    yyaxis(ax, 'left'); hold(ax, 'on');
    pt = m.t_in_cycle; pa = 100 * R.info.pulse_amps(1:numel(pt)) / R.info.amp;
    for j = 1:numel(pt)
        plot(ax, [pt(j) pt(j)], [60 pa(j)], '-', 'Color', [0.62 0.62 0.62], 'LineWidth', 0.9 + 0.5 * (Tm <= 2));
    end
    plot(ax, m.envelope_t, 100 * m.envelope / R.info.amp, '-', 'Color', grey, 'LineWidth', 2.5);
    ylim(ax, [60 160]); set(ax, 'YTick', [70 100], 'YColor', grey);
    if p == 1, ylabel(ax, {'Stimulus', '(% of peak)'}, 'FontSize', fsz); end
    xline(ax, Tm / 2, '--', 'Color', grey, 'LineWidth', 1.5);
    yyaxis(ax, 'right'); hold(ax, 'on');
    tt = linspace(0, Tm, 400); ph = 2 * pi * tt / Tm;
    fit_y = f.a * cos(ph) + f.b * sin(ph) + f.offset;
    sz = m.sizes{1};
    plot(ax, tt, fit_y, '-', 'Color', show_col, 'LineWidth', 3);
    t_peak = Tm * f.peak_phase / (2 * pi);
    xline(ax, t_peak, '--', 'Color', show_col, 'LineWidth', 1.5);
    lo = min([sz fit_y]); hi = max([sz fit_y]); span = hi - lo;
    ylim(ax, [lo - 1.5 * span, hi + 0.75 * span]);
    set(ax, 'YColor', show_col, 'YTick', nice_ticks(lo, hi));
    % The per-pulse sizes as filled PATCH discs, not line markers: on the
    % right side of a yyaxis axes the 'o' markers vanished from the 600-dpi
    % PNG export (present at 75 dpi). Radius in inches, set after the limits.
    discs(ax, pt, sz, 0.045 + 0.02 * (Tm <= 2), colA_w, top_h, show_col);
    if p == nP, ylabel(ax, {'Response size', '(RMS 10-60 ms)'}, 'FontSize', fsz); end
    text(ax, 0.98, 0.99, sprintf('lead %.0f%s', f.lead_deg, char(176)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', show_col, 'FontSize', fsz + 2);
    text(ax, 0.98, 0.84, sprintf('r^2 = %.2f', f.r2), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', [0.3 0.3 0.3], 'FontSize', fsz - 3);
    title(ax, sprintf('%g s period  (%d cycles)', Tm, m.n_cycles_used), 'FontWeight', 'normal', 'FontSize', fsz + 1);
    xlim(ax, [0 Tm]); set(ax, 'XTickLabel', []);
    style_axes(ax, fsz);

    % Cycle-averaged traces
    for row = 1:size(trace_rows, 1)
        ax = mkax(colA_x(p), tr_b(row), colA_w, tr_h); ax_tr(row, p) = ax; hold(ax, 'on');
        var = trace_rows{row, 1}; groups = trace_rows{row, 2};
        for g = 1:numel(groups)
            plot(ax, m.cycle_t, m.cycle_traces.(var).(groups{g}), '-', ...
                'Color', group_color.(groups{g}), 'LineWidth', 1.4);
        end
        xline(ax, Tm / 2, '--', 'Color', grey, 'LineWidth', 1.2);
        ro = trace_rows{row, 4};
        if ~isempty(ro)
            fi = m.fits(strcmp(m.readout_names, ro));
            label = sprintf('lead %.0f%s', fi.lead_deg, char(176));
        else
            % When the gain reduction is strongest: peak of the SFA feedback,
            % trough of the depression product (fundamental), driven E.
            lag = m.trace_lag_deg.(var).E_stim;
            t_ext = mod(Tm / 2 + Tm * lag / 360, Tm);
            y = m.cycle_traces.(var).E_stim;
            if strcmp(var, 'sfa')
                plot(ax, t_ext, max(y), 'v', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
                label = sprintf('peak lags %.0f%s', lag, char(176));
            else
                plot(ax, t_ext, min(y), '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
                label = sprintf('trough lags %.0f%s', lag, char(176));
            end
        end
        text(ax, 0.99, 1.02, label, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', fsz - 2);
        xlim(ax, [0 Tm]);
        if p == 1, ylabel(ax, trace_rows{row, 3}, 'FontSize', fsz - 1); end
        if row < size(trace_rows, 1)
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'Time in modulation cycle (s)', 'FontSize', fsz);
        end
        style_axes(ax, fsz);
    end
end
for row = 1:size(trace_rows, 1)
    linkaxes(ax_tr(row, :), 'y');
    yl = ylim(ax_tr(row, 1));
    set(ax_tr(row, :), 'YTick', nice_ticks(yl(1), yl(2)));
end
% E and I of the undriven row, named in their colours rather than a legend box.
text(ax_tr(1, 1), 0.01, 1.02, sprintf('\\color[rgb]{%g %g %g}E   \\color[rgb]{%g %g %g}I', ...
    tcol(1, :), tcol(2, :)), 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontSize', fsz - 1);
% The heading is written in PIECES: a single text object wider than about
% 6.8 in (4096 px at 600 dpi) silently vanished from the 600-dpi PNG export
% (it was there at 60 dpi), as an annotation textbox and as axes text alike.
pieces = {'(A)   ', [cond_label(cfg.show_condition) ' network, '], 'sine-modulated 5 Hz pulse trains'};
xh = (0.1 - colA_x(1)) / colA_w;
for i = 1:numel(pieces)
    th = text(ax_top(1), xh, (H - 0.3 - top_b) / top_h, pieces{i}, 'Units', 'normalized', ...
        'FontSize', fsz + 2, 'Color', 'k', 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', 'Clipping', 'off');
    xh = xh + th.Extent(3);
end

% ---- (B) lead against period ------------------------------------------------
axB = mkax(rx, 5.75, W - rx - 0.25, 3.0); hold(axB, 'on');
hb = gobjects(1, numel(conds) + 1); lead_rows = {};
jit = linspace(-0.035, 0.035, numel(conds));
for c = 1:numel(conds)
    L = nan(1, nP); CI = nan(2, nP); r2 = nan(1, nP);
    for p = 1:nP
        k = pick(conds{c}, sprintf('sine_%gs', periods(p)));
        L(p) = M(k).fits(1).lead_deg; CI(:, p) = M(k).fits(1).lead_ci(:); r2(p) = M(k).fits(1).r2;
        lead_rows(end + 1, :) = {conds{c}, periods(p), L(p), CI(1, p), CI(2, p), r2(p), ...
            M(k).fits(1).rel_mod, M(k).n_cycles_used}; %#ok<AGROW>
    end
    col = st.condition_color(conds{c});
    xp = periods .* 10 .^ jit(c);
    errorbar(axB, xp, L, L - CI(1, :), CI(2, :) - L, '-', 'Color', col, 'LineWidth', 2, 'CapSize', 0);
    hb(c) = plot(axB, xp, L, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 9);
end
hb(end) = plot(axB, cfg.human_periods, cfg.human_leads, 'd--', 'Color', [0.45 0.45 0.45], ...
    'MarkerFaceColor', [0.45 0.45 0.45], 'MarkerSize', 10, 'LineWidth', 1.5);
yline(axB, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
set(axB, 'XScale', 'log', 'XTick', periods, 'XLim', [periods(1) / 1.3, periods(end) * 1.3], ...
    'YLim', cfg.lead_ylim);
xlabel(axB, 'Modulation period (s)', 'FontSize', fsz);
ylabel(axB, ['Lead (' char(176) ')'], 'FontSize', fsz);
style_axes(axB, fsz);
leg_names = cellfun(@cond_label, conds, 'UniformOutput', false);
% The non-adapting network fluctuates strongly and its fits are poor; the legend
% says so where the reader meets its long (clipped) intervals.
for c = 1:numel(conds)
    rc = [lead_rows{strcmp(lead_rows(:, 1), conds{c}), 6}];
    if min(rc) < 0.5, leg_names{c} = sprintf('%s (r^2 %.2f-%.2f)', leg_names{c}, min(rc), max(rc)); end
end
lg = legend(axB, hb, [leg_names, {'Human SEEG (median)'}], 'FontSize', fsz - 3, 'Box', 'off');
lg.Units = 'inches'; lg.Position(1:2) = [rx + 0.05, 8.95];
annotation(fig, 'textbox', [(rx - 0.85) / W, (H - 0.5) / H, 0.6 / W, 0.42 / H], 'String', '(B)', ...
    'EdgeColor', 'none', 'FontSize', fsz + 2, 'VerticalAlignment', 'middle', 'Margin', 0);

% ---- (C) paired pulse --------------------------------------------------------
kp = pick(cfg.show_condition, 'paired_pulse');
mp = M(kp); P = mp.paired; et = 1000 * mp.epoch_t;
j100 = find(abs(P.isis - 0.1) < 1e-9, 1);
cw = 1.85;
axC = mkax(rx, 2.35, cw, 1.85); hold(axC, 'on');
s1 = mp.epochs{1}(1, :); pr = mp.epochs{1}(j100, :); lti = P.lti{1}(j100, :);
h1 = plot(axC, et, s1, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
h2 = plot(axC, et, pr, '-', 'Color', show_col, 'LineWidth', 2.5);
h3 = plot(axC, et, lti, '--', 'Color', show_col * 0.55, 'LineWidth', 2);
top = max([s1 pr lti]);
ylim(axC, [-0.22 * top, 1.95 * top]);      % room for the legend above the traces
plot(axC, [0 100], [-0.12 -0.12] * top, '^', 'MarkerSize', 7, 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none');                               % the pulses, under the baseline
xlim(axC, [-40 350]); set(axC, 'XTickLabel', [], 'YTick', nice_ticks(0, top));
ylabel(axC, {'Undriven', 'mean x'}, 'FontSize', fsz - 1);
style_axes(axC, fsz);
legend(axC, [h1 h2 h3], {'single pulse', '100 ms pair', 'LTI prediction'}, ...
    'Location', 'northeast', 'FontSize', fsz - 4, 'Box', 'off');
title(axC, sprintf('Paired pulse, %s', lower(cond_label(cfg.show_condition))), ...
    'FontWeight', 'normal', 'FontSize', fsz - 1, 'Color', show_col, 'HorizontalAlignment', 'left', ...
    'Units', 'normalized', 'Position', [-0.35 1.03 0]);

axC2 = mkax(rx, 0.65, cw, 1.4); hold(axC2, 'on');
i_sfa = find(strcmp(mp.readout_names, 'sfa_stim'), 1);
i_std = find(strcmp(mp.readout_names, 'std_stim'), 1);
d_sfa = mp.epochs{i_sfa}(j100, :); d_std = -mp.epochs{i_std}(j100, :);
std_col = [0.42 0.16 0.60];
hs1 = plot(axC2, et, d_sfa / max(abs(d_sfa)), '-', 'Color', E_drv, 'LineWidth', 2.2);
hs2 = plot(axC2, et, d_std / max(abs(d_std)), '-', 'Color', std_col, 'LineWidth', 2.2);
plot(axC2, [0 100], [-0.12 -0.12], '^', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
xlim(axC2, [-40 350]); ylim(axC2, [-0.25 1.75]); set(axC2, 'YTick', [0 1], 'XTick', [0 100 200 300]);
ylabel(axC2, {'Driven E', '(scaled)'}, 'FontSize', fsz - 1);
xlabel(axC2, 'Time after first pulse (ms)', 'FontSize', fsz);
legend(axC2, [hs1 hs2], {sprintf('SFA (+%.2f)', max(d_sfa)), sprintf('depression (%.2f)', -max(d_std))}, ...
    'Location', 'northwest', 'FontSize', fsz - 4, 'Box', 'off');
style_axes(axC2, fsz);
annotation(fig, 'textbox', [(rx - 0.85) / W, 4.55 / H, 0.6 / W, 0.42 / H], 'String', '(C)', ...
    'EdgeColor', 'none', 'FontSize', fsz + 2, 'VerticalAlignment', 'middle', 'Margin', 0);

axC3 = mkax(rx + cw + 0.85, 0.65, W - (rx + cw + 0.85) - 0.2, 3.55); hold(axC3, 'on');
isi_ix = find(P.isis > 0);
nc = numel(conds); bw = 0.8 / nc;
ratio_rows = {};
for c = 1:nc
    kc = pick(conds{c}, 'paired_pulse'); Pc = M(kc).paired;
    col = st.condition_color(conds{c});
    for jj = 1:numel(isi_ix)
        j = isi_ix(jj); xc = jj + (c - (nc + 1) / 2) * bw;
        bar(axC3, xc, Pc.ratio(1, j), bw * 0.9, 'FaceColor', col, 'EdgeColor', 'none');
        plot(axC3, xc + [-0.5 0.5] * bw, Pc.ratio_lti(1, j) * [1 1], 'k-', 'LineWidth', 2.5);
        ratio_rows(end + 1, :) = {conds{c}, 1000 * P.isis(j), Pc.ratio(1, j), Pc.ratio_lti(1, j), ...
            Pc.ratio(1, j) / Pc.ratio_lti(1, j), Pc.ratio_isolated(1, j), Pc.n_events(j)}; %#ok<AGROW>
    end
end
set(axC3, 'XTick', 1:numel(isi_ix), 'XTickLabel', arrayfun(@(d) sprintf('%g', 1000 * d), P.isis(isi_ix), ...
    'UniformOutput', false), 'XLim', [0.5 numel(isi_ix) + 0.5], 'YLim', [0 2.6], 'YTick', 0:0.5:2.5);
xlabel(axC3, 'Interval (ms)', 'FontSize', fsz);
ylabel(axC3, 'Second / single response', 'FontSize', fsz);
hk = plot(axC3, NaN, NaN, 'k-', 'LineWidth', 2.5);
legend(axC3, hk, {'LTI'}, 'Location', 'north', 'FontSize', fsz - 3, 'Box', 'off');
style_axes(axC3, fsz);

%% Tables and provenance
leads = cell2table(lead_rows, 'VariableNames', {'condition', 'period_s', 'lead_deg', ...
    'ci_lo', 'ci_hi', 'r2', 'rel_mod', 'n_cycles'});
ratios = cell2table(ratio_rows, 'VariableNames', {'condition', 'isi_ms', 'ratio', ...
    'ratio_lti', 'ratio_over_lti', 'ratio_isolated', 'n_events'});
lead_all = {};
for k = find(is_sine)
    for i = 1:numel(M(k).fits)
        lead_all(end + 1, :) = {M(k).condition, M(k).mod_period, M(k).readout_names{i}, ...
            M(k).fits(i).lead_deg, M(k).fits(i).r2}; %#ok<AGROW>
    end
end
lead_all = cell2table(lead_all, 'VariableNames', {'condition', 'period_s', 'readout', 'lead_deg', 'r2'});
c = S.settings.cfg;
source = sprintf(['%s; run_stim_evoked, preset %s, fs %g Hz, seeds %s, sigma_u_noise %s; ', ...
    'E-only drive (%g of E), amp %g, width %g ms, mod_depth %g, settle %g s, ', ...
    'first %g s of stimulation discarded; readout mean x of undriven neurons; git %s'], ...
    mat_file, c.preset_name, c.fs, mat2str(c.rng_seeds), num2str_or_preset(c.sigma_u_noise), ...
    c.stim.neuron_frac, c.stim.amp, 1000 * c.stim.width_s, c.stim.mod_depth, c.settle_s, ...
    cfg.discard_s, S.settings.git.commit_short);
out = struct('figs', fig, 'files', {{}}, 'source', source, 'leads', leads, ...
    'ratios', ratios, 'lead_all', lead_all, 'measures', M);
vprintf(cfg.verbose, 'minimal', '  [sine_paired_stim] leads (deg):\n');
for i = 1:height(leads)
    vprintf(cfg.verbose, 'minimal', '    %-14s %4.0f s  %6.1f  [%6.1f %6.1f]  r2 %.2f\n', leads.condition{i}, ...
        leads.period_s(i), leads.lead_deg(i), leads.ci_lo(i), leads.ci_hi(i), leads.r2(i));
end
if cfg.save
    save_figure_stable(out_dir, 'fig_sine_paired_stim', fig);
    out.files = existing_outputs(out_dir, 'fig_sine_paired_stim');
end
if cfg.visible, set(fig, 'Visible', 'on'); end
end

function style_axes(ax, fsz)
set(ax, 'FontSize', fsz - 2, 'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'off');
end

function s = num2str_or_preset(v)
if isempty(v), s = 'preset (0.025)'; else, s = num2str(v); end
end

function s = cond_label(name)
% Full condition names in the style guide's sentence case.
labels = containers.Map({'no_adaptation', 'sfa1_std1', 'sfa3_std2'}, ...
    {'No adaptation', 'Single-timescale adaptation', 'Multiple-timescale adaptation'});
if isKey(labels, name)
    s = labels(name);
else
    s = srnn_condition_titles();
    s = s(name);
end
end

function t = nice_ticks(lo, hi)
% Two or three round tick values inside [lo, hi] (1-2-5 steps).
span = hi - lo;
if ~(span > 0), t = lo; return; end
steps = 10 .^ floor(log10(span)) * [0.1 0.2 0.25 0.5 1 2 5];
for s = steps(end:-1:1)
    t = (ceil(lo / s) * s):s:hi;
    if numel(t) >= 2, break; end
end
if numel(t) > 3, t = t(1:2:end); end
t = round(t, 10);
end

function discs(ax, x, y, r_in, w_in, h_in, col)
% Filled circles of radius r_in inches at (x, y), drawn as patches so every
% export path renders them; the axes' current limits and size give the scale.
xl = xlim(ax); yl = ylim(ax);
rx = r_in * diff(xl) / w_in; ry = r_in * diff(yl) / h_in;
th = linspace(0, 2 * pi, 40)';
for i = 1:numel(x)
    patch(ax, x(i) + rx * cos(th), y(i) + ry * sin(th), col, 'EdgeColor', 'w', 'LineWidth', 0.75);
end
end

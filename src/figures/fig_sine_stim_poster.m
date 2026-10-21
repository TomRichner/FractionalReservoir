function out = fig_sine_stim_poster(cfg)
% FIG_SINE_STIM_POSTER Column-width poster version of fig_sine_paired_stim: the 2 s mechanism and the lead against period.
%
%   out = FIG_SINE_STIM_POSTER('run_dir', d)                 % d holds stim_evoked/
%   out = FIG_SINE_STIM_POSTER('run_dir', d, 'out_dir', o)
%
% For the Minnesota Neuromodulation poster (TR via the lead, 2026-10-20),
% drawn from the same saved simulations as fig_sine_paired_stim
% (run_stim_evoked, sine_paired_stim_config); nothing is simulated here.
%
%   Left: the multiple-timescale network (sfa3_std2) at the 2 s period,
%   four rows over one modulation cycle sharing the time axis:
%     1. the stimulus envelope with its pulses (left axis, % of peak) and
%        the per-pulse response size of the undriven neurons' mean x with
%        its fitted sinusoid (right axis, arbitrary units, unlabelled), the
%        lead printed;
%     2. the driven E neurons' synaptic output r prod_m b_m, its per-pulse
%        lead printed;
%     3. their SFA feedback (c/K) sum_k a_k, the lag of its peak;
%     4. their depression product prod_m b_m, the lag of its trough;
%   a dashed line at the envelope peak runs through all four.
%   Right: lead against modulation period for the single- and
%   multiple-timescale networks (95 % bootstrap intervals over cycles; the
%   no-adaptation network only if show_no_adaptation, faint, as its fits are
%   poor: r^2 0.03-0.65), and the
%   human medians as grey diamonds; legend below the axes.
%
% Measures as in fig_sine_paired_stim (stim_evoked_measures); the quantities
% are those of docs/EquationsParametersDocs/Equations_stability_paper.md.
%
% SIZE AND TYPE. Canvas 10.0 x 9.2 in, every text 20 pt or more (titles
% 23 pt), Arial. 10.0 rather than the requested 10.6 in wide because
% save_figure_stable exports PNGs at a fixed 600 dpi and the PNG must stay
% under 6000 px on its long side for previewing. Greek letters are Unicode
% characters: a TeX \Pi or \it anywhere makes the R2026a vector SVG outline
% every text object (see fig_sine_paired_stim). Layout by hand in inches,
% for the reasons given there.
%
% out.figs / out.files / out.source as every figure; plus out.leads.
%
% See also: fig_sine_paired_stim, stim_evoked_measures, run_stim_evoked

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.mat_file    (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''      % unused; accepted for a uniform call
    cfg.show_condition (1,:) char = 'sfa3_std2'
    cfg.period      (1,1) double  = 2
    cfg.discard_s   (1,1) double  = 10
    cfg.human_periods (1,:) double = [2 5 10]
    cfg.human_leads   (1,:) double = [10.0 7.7 6.2]   % wider cohort medians (22 subjects)
    cfg.human_label (1,:) char    = 'Human SEEG (22 subjects)'
    cfg.show_no_adaptation (1,1) logical = false   % true: drawn faint, no intervals
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.verbose     (1,:) char    = 'minimal'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st = manuscript_style();
FN = 'Arial'; fs_txt = 20; fs_ttl = 23;

if isempty(cfg.mat_file)
    mat_file = resolve_data_file('', cfg.run_dir, {fullfile(cfg.run_dir, 'stim_evoked')}, ...
        'stim_evoked_data.mat', 'run_stim_evoked(sine_paired_stim_config())');
else
    mat_file = cfg.mat_file;
end
S = load(mat_file, 'results', 'settings');
M = stim_evoked_measures(S.results, 'discard_s', cfg.discard_s);
pick = @(cond, proto) find(strcmp({M.condition}, cond) & strcmp({M.protocol}, proto), 1);
is_sine = strcmp({M.kind}, 'sine');
periods = unique([M(is_sine).mod_period]);

tcol = SRNNCellTypePairs.type_colors(2);
E_drv = tcol(1, :) * 0.8;
show_col = st.condition_color(cfg.show_condition);
grey = [0.2 0.2 0.2];
PI = char(928); SIG = char(931); DEG = char(176);

%% Canvas (inches from the lower left)
W = 10.0; H = 9.2;
fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'inches', 'Position', [0.2 0.2 W H]);
mkax = @(l, b, w, h) axes(fig, 'Units', 'inches', 'Position', [l b w h], 'FontName', FN);
lx = 1.6; lw = 2.65;                       % left column
row_b = [5.55, 3.85, 2.45, 1.05]; row_h = [2.25, 1.05, 1.05, 1.05];

k = pick(cfg.show_condition, sprintf('sine_%gs', cfg.period));
m = M(k); R = S.results(k); f = m.fits(1); Tm = cfg.period;

% ---- row 1: stimulus and response size --------------------------------------
ax1 = mkax(lx, row_b(1), lw, row_h(1));
yyaxis(ax1, 'left'); hold(ax1, 'on');
pt = m.t_in_cycle; pa = 100 * R.info.pulse_amps(1:numel(pt)) / R.info.amp;
for j = 1:numel(pt)
    plot(ax1, [pt(j) pt(j)], [60 pa(j)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6);
end
plot(ax1, m.envelope_t, 100 * m.envelope / R.info.amp, '-', 'Color', grey, 'LineWidth', 3);
ylim(ax1, [60 190]); set(ax1, 'YTick', [70 100], 'YColor', grey);
ylabel(ax1, 'Stimulus (%)', 'FontSize', fs_txt, 'FontName', FN);
yyaxis(ax1, 'right'); hold(ax1, 'on');
tt = linspace(0, Tm, 400); ph = 2 * pi * tt / Tm;
fit_y = f.a * cos(ph) + f.b * sin(ph) + f.offset;
sz = m.sizes{1};
plot(ax1, tt, fit_y, '-', 'Color', show_col, 'LineWidth', 3.5);
t_peak = Tm * f.peak_phase / (2 * pi);
xline(ax1, t_peak, '--', 'Color', show_col, 'LineWidth', 2);
lo = min([sz fit_y]); hi = max([sz fit_y]); span = hi - lo;
ylim(ax1, [lo - 2.6 * span, hi + 1.5 * span]);   % the sizes above the envelope, labels above them
set(ax1, 'YTick', []); ax1.YAxis(2).Visible = 'off';   % size scale off: arbitrary units
text(ax1, 0.98, 0.85, 'response size', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'Color', show_col, 'FontSize', fs_txt, 'FontName', FN);
xlim(ax1, [0 Tm]);
discs(ax1, pt, sz, 0.075, lw, row_h(1), show_col);      % patches: see fig_sine_paired_stim
yyaxis(ax1, 'left');
text(ax1, 0.98, 0.99, sprintf('lead %.0f%s', f.lead_deg, DEG), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', show_col, ...
    'FontSize', fs_txt + 2, 'FontName', FN);

% ---- rows 2-4: driven E synaptic output, SFA, STD -----------------------------
rows = { ...
    'syn', {'Synaptic', ['output r' PI 'b']}, 'syn_stim'; ...
    'sfa', {'SFA', ['(c/K)' SIG 'a']},        ''; ...
    'std', {'STD', [PI 'b']},                 ''};
ax_r = gobjects(1, 3);
for i = 1:3
    ax = mkax(lx, row_b(i + 1), lw, row_h(i + 1)); ax_r(i) = ax; hold(ax, 'on');
    var = rows{i, 1};
    y = m.cycle_traces.(var).E_stim;
    plot(ax, m.cycle_t, y, '-', 'Color', E_drv, 'LineWidth', 2);
    yr = [min(y), max(y)]; pad = diff(yr);
    if isempty(rows{i, 3})
        lag = m.trace_lag_deg.(var).E_stim;
        t_ext = mod(Tm / 2 + Tm * lag / 360, Tm);
        if strcmp(var, 'sfa')
            plot(ax, t_ext, max(y), 'v', 'MarkerSize', 13, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
            label = sprintf('peak lags %.0f%s', lag, DEG);
        else
            plot(ax, t_ext, min(y), '^', 'MarkerSize', 13, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
            label = sprintf('trough lags %.0f%s', lag, DEG);
        end
    else
        fi = m.fits(strcmp(m.readout_names, rows{i, 3}));
        label = sprintf('lead %.0f%s', fi.lead_deg, DEG);
    end
    ylim(ax, [yr(1) - 0.08 * pad, yr(2) + 0.75 * pad]);     % headroom for the label
    set(ax, 'YTick', nice_ticks(yr(1), yr(2)));
    text(ax, 0.99, 0.99, label, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', fs_txt, 'FontName', FN);
    ylabel(ax, rows{i, 2}, 'FontSize', fs_txt, 'FontName', FN);
end
all_left = [ax1, ax_r];
for ax = all_left
    xlim(ax, [0 Tm]);
    xline(ax, Tm / 2, '--', 'Color', grey, 'LineWidth', 1.8);   % the stimulus peak, every row
    style(ax, fs_txt, FN);
    if ax ~= ax_r(end), set(ax, 'XTickLabel', []); end
end
set(ax_r(end), 'XTick', 0:1:Tm);
xlabel(ax_r(end), 'Time in cycle (s)', 'FontSize', fs_txt, 'FontName', FN);
text(ax1, 0.5, 1.06, sprintf('Model, %g s period', Tm), 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs_ttl, 'FontName', FN);
% Who the traces belong to, under the heading of rows 2-4.
text(ax_r(1), 0.01, 0.99, 'driven E', 'Units', 'normalized', 'Color', E_drv, ...
    'VerticalAlignment', 'top', 'FontSize', fs_txt, 'FontName', FN);

% ---- right: lead against period ------------------------------------------------
axB = mkax(6.25, 3.4, 3.55, 4.75); hold(axB, 'on');
series = {'sfa1_std1', 'Single timescale'; 'sfa3_std2', 'Multiple timescales'};
if cfg.show_no_adaptation, series = [{'no_adaptation', 'No adaptation'}; series]; end
hb = gobjects(1, 0); names = {}; lead_rows = {};
jit = linspace(-0.03, 0.03, size(series, 1));
for c = 1:size(series, 1)
    L = nan(1, numel(periods)); CI = nan(2, numel(periods));
    for p = 1:numel(periods)
        kk = pick(series{c, 1}, sprintf('sine_%gs', periods(p)));
        L(p) = M(kk).fits(1).lead_deg; CI(:, p) = M(kk).fits(1).lead_ci(:);
        lead_rows(end + 1, :) = {series{c, 1}, periods(p), L(p), CI(1, p), CI(2, p), M(kk).fits(1).r2}; %#ok<AGROW>
    end
    col = st.condition_color(series{c, 1});
    xp = periods .* 10 .^ jit(c);
    if strcmp(series{c, 1}, 'no_adaptation')
        % Faint and without intervals: its fits are poor (r^2 0.03-0.65).
        plot(axB, xp, L, 'o-', 'Color', [0.75 0.75 0.75], 'MarkerFaceColor', [0.75 0.75 0.75], ...
            'MarkerSize', 7, 'LineWidth', 1.5);
        continue
    end
    errorbar(axB, xp, L, L - CI(1, :), CI(2, :) - L, '-', 'Color', col, 'LineWidth', 3, 'CapSize', 0);
    hb(end + 1) = plot(axB, xp, L, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 12); %#ok<AGROW>
    names{end + 1} = series{c, 2}; %#ok<AGROW>
end
hb(end + 1) = plot(axB, cfg.human_periods, cfg.human_leads, 'd--', 'Color', [0.45 0.45 0.45], ...
    'MarkerFaceColor', [0.45 0.45 0.45], 'MarkerSize', 13, 'LineWidth', 2);
names{end + 1} = cfg.human_label;
yline(axB, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
set(axB, 'XScale', 'log', 'XTick', periods, 'XLim', [periods(1) / 1.3, periods(end) * 1.3], ...
    'YLim', [0 45], 'YTick', 0:20:40);
xlabel(axB, 'Modulation period (s)', 'FontSize', fs_txt, 'FontName', FN);
ylabel(axB, ['Lead (' DEG ')'], 'FontSize', fs_txt, 'FontName', FN);
title(axB, 'Lead vs period', 'FontWeight', 'normal', 'FontSize', fs_ttl, 'FontName', FN);
style(axB, fs_txt, FN);
lg = legend(axB, hb, names, 'FontSize', fs_txt, 'FontName', FN, 'Box', 'off');
lg.Units = 'inches';
lg.Position(1:2) = [5.45, 1.0];
if cfg.show_no_adaptation
    text(axB, 0.98, 0.02, 'grey: no adaptation', 'Units', 'normalized', 'Color', [0.55 0.55 0.55], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs_txt, 'FontName', FN);
end

%% Save
leads = cell2table(lead_rows, 'VariableNames', {'condition', 'period_s', 'lead_deg', 'ci_lo', 'ci_hi', 'r2'});
c = S.settings.cfg;
source = sprintf(['%s; run_stim_evoked, preset %s, fs %g Hz, seeds %s; E-only drive, amp %g; ', ...
    'readout mean x of undriven neurons; human medians %s deg (%s); git %s'], ...
    mat_file, c.preset_name, c.fs, mat2str(c.rng_seeds), c.stim.amp, mat2str(cfg.human_leads), ...
    cfg.human_label, S.settings.git.commit_short);
out = struct('figs', fig, 'files', {{}}, 'source', source, 'leads', leads);
if cfg.save
    save_figure_stable(out_dir, 'fig_sine_stim_poster', fig);
    out.files = existing_outputs(out_dir, 'fig_sine_stim_poster');
end
vprintf(cfg.verbose, 'minimal', '  [sine_stim_poster] %d files in %s\n', numel(out.files), out_dir);
if cfg.visible, set(fig, 'Visible', 'on'); end
end

function style(ax, fs, fn)
set(ax, 'FontSize', fs, 'FontName', fn, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off', ...
    'TickLength', [0.02 0.02]);
end

function discs(ax, x, y, r_in, w_in, h_in, col)
% Filled circles of radius r_in inches, as patches (line markers on a
% yyaxis right side vanished from the 600-dpi PNG; fig_sine_paired_stim).
xl = xlim(ax); yl = ylim(ax);
rx = r_in * diff(xl) / w_in; ry = r_in * diff(yl) / h_in;
th = linspace(0, 2 * pi, 40)';
for i = 1:numel(x)
    patch(ax, x(i) + rx * cos(th), y(i) + ry * sin(th), col, 'EdgeColor', 'w', 'LineWidth', 1);
end
end

function t = nice_ticks(lo, hi)
% Two or three round tick values inside [lo, hi] (1-2-5 steps); the same
% rule as fig_sine_paired_stim's.
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

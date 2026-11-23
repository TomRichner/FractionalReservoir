function files = amp_outlier_figs(varargin)
%AMP_OUTLIER_FIGS Figures and the table of medians for amp_outlier_compute's results.
%
%   files = amp_outlier_figs()
%   files = amp_outlier_figs('res_file', f, 'fig_dir', d)
%
% Fig_AmpOutlier_budget   per condition: Re lambda_O (the rightmost
%     eigenvalue of the mean part P, units of W) against the feed-forward
%     weight w_FF of P, one point per calibration state; filled = measured
%     gains, open = uniform gain (every gain at the rms), a thin line joins
%     the two of each state; dotted arcs lambda_O^2 + w_FF^2 = const through
%     each network's median uniform-gain point.
% Fig_AmpOutlier_margin   by network (E-to-E dose on the horizontal axis):
%     (A) non-normal margin omega - alpha of J_xx, (B) peak of
%     ||expm(J t)||_2 over [0, 2] s, (C) ||expm(J tau_d)||_2; medians and
%     the range over the 12 states, measured solid, uniform gain dashed.
% Each PNG (300 dpi) has a .md sidecar; amp_outlier_table.md holds the
% table of medians [min, max] for every quantity.
%
% See also: amp_outlier_compute, manuscript_style

p = inputParser;
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
p.addParameter('res_file', fullfile(root, 'data', 'explorations', 'amplification_outlier', 'amp_outlier_results.mat'));
p.addParameter('fig_dir', fullfile(root, 'figs', 'narrative_amplification_outlier'));
p.parse(varargin{:});
o = p.Results;
if ~isfolder(o.fig_dir); mkdir(o.fig_dir); end
L = load(o.res_file); T = L.T; tg = L.t_grid;
T.reO = real(T.lamO); T.reO_u = real(T.lamO_u);
T.l2 = real(T.lam2); T.l2_u = real(T.lam2_u);
T.Gc = cellfun(@(g, a) max(g .* exp(-a * tg)), T.Gt, num2cell(T.alpha));
T.Gc_u = cellfun(@(g, a) max(g .* exp(-a * tg)), T.Gt_u, num2cell(T.alpha_u));

st = manuscript_style();
nets = {'reference', 'EEx1p1', 'EEx1p25', 'EEx1p5', 'EEx1p75'};
dose = [1 1.1 1.25 1.5 1.75];
conds = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
mk = {'o', 's', 'd', '^', 'v'};
fs = st.tick_fs;
files = {};

%% Figure 1: the budget
f1 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 1300 430]);
tl = tiledlayout(f1, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
xl = [min([T.reO; T.reO_u]) - 0.3, max([T.reO; T.reO_u]) + 0.3];
yl = [0, max([T.wFF; T.wFF_u]) * 1.08];
hl = gobjects(1, numel(nets));
for c = 1:numel(conds)
    ax = nexttile(tl); hold(ax, 'on');
    cc = st.condition_color(conds{c});
    for q = 1:numel(nets)
        S = T(strcmp(T.network, nets{q}) & strcmp(T.condition, conds{c}), :);
        col = shade(cc, q, numel(nets));
        r0 = hypot(median(S.reO_u), median(S.wFF_u));
        th = linspace(0, pi, 200);
        plot(ax, r0 * cos(th), r0 * sin(th), ':', 'Color', col, 'LineWidth', 1.0);
        for k = 1:height(S)
            plot(ax, [S.reO_u(k) S.reO(k)], [S.wFF_u(k) S.wFF(k)], '-', 'Color', [col 0.35], 'LineWidth', 0.6);
        end
        plot(ax, S.reO_u, S.wFF_u, mk{q}, 'MarkerEdgeColor', col, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.0);
        hl(q) = plot(ax, S.reO, S.wFF, mk{q}, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'MarkerSize', 6);
    end
    xline(ax, 1, '-', 'Color', st.zeroline_color, 'LineWidth', 0.5);
    set(ax, 'XLim', xl, 'YLim', yl, 'FontSize', fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    title(ax, st.condition_title(conds{c}), 'Color', cc, 'FontSize', fs, 'FontWeight', 'normal');
    xlabel(ax, 'Re \lambda_O', 'FontSize', fs);
    if c == 1; ylabel(ax, 'w_{FF}', 'FontSize', fs); end
    if c == numel(conds)
        lg = legend(ax, hl, {'reference', 'EE \times 1.1', 'EE \times 1.25', 'EE \times 1.5', 'EE \times 1.75'}, ...
            'Location', 'northeast', 'FontSize', fs - 2, 'Box', 'off');
        lg.ItemTokenSize = [10 10];
    end
    text(ax, -0.12, 1.02, sprintf('(%s)', char('A' + c - 1)), 'Units', 'normalized', 'FontSize', fs);
end
files{end + 1} = save_png(f1, fullfile(o.fig_dir, 'Fig_AmpOutlier_budget.png'));
close(f1);

%% Figure 2: margin and gain by dose
f2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 1700 450]);
tl = tiledlayout(f2, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ys = {{'margin', 'margin_u'}, {'Gpeak', 'Gpeak_u'}, {'G_tau', 'G_tau_u'}, {'Gc', 'Gc_u'}};
ylab = {'\omega - \alpha  (1/s)', 'max_{t \leq 2 s} ||e^{Jt}||', '||e^{J\tau_d}||', 'max_{t \leq 2 s} e^{-\alpha t}||e^{Jt}||'};
islog = [false true true true];
hc = gobjects(1, numel(conds));
for pnl = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    for c = 1:numel(conds)
        cc = st.condition_color(conds{c});
        off = (c - 2) * 0.012;
        for w = 1:2
            v = ys{pnl}{w};
            md = nan(1, 5); lo = md; hi = md;
            for q = 1:5
                x = T.(v)(strcmp(T.network, nets{q}) & strcmp(T.condition, conds{c}));
                md(q) = median(x); lo(q) = min(x); hi(q) = max(x);
            end
            xs = dose + off + (w - 1) * 0.006;
            ls = '-'; if w == 2; ls = '--'; end
            for q = 1:5
                plot(ax, [xs(q) xs(q)], [lo(q) hi(q)], ls, 'Color', [cc 0.6], 'LineWidth', 1.0);
            end
            h = plot(ax, xs, md, ls, 'Color', cc, 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', ifelse(w == 1, cc, [1 1 1]), 'MarkerEdgeColor', cc);
            if w == 1; hc(c) = h; end
        end
    end
    if islog(pnl); set(ax, 'YScale', 'log'); end
    set(ax, 'XLim', [0.95 1.8], 'XTick', dose, 'XTickLabel', {'1', '1.1', '1.25', '1.5', '1.75'}, ...
        'FontSize', fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    xlabel(ax, 'E-to-E mean (\times reference)', 'FontSize', fs);
    ylabel(ax, ylab{pnl}, 'FontSize', fs);
    text(ax, -0.17, 1.02, sprintf('(%s)', char('A' + pnl - 1)), 'Units', 'normalized', 'FontSize', fs);
    if pnl == 1
        hm = plot(ax, NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        hu = plot(ax, NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        lg = legend(ax, [hc hm hu], [cellfun(@(k) st.condition_title(k), conds, 'UniformOutput', false), ...
            {'measured gains', 'uniform gain'}], 'Location', 'northeast', 'FontSize', fs - 3, 'Box', 'off');
        lg.ItemTokenSize = [18 10];
    end
end
files{end + 1} = save_png(f2, fullfile(o.fig_dir, 'Fig_AmpOutlier_margin.png'));
close(f2);

%% Tables
vars = {'d_rms', 'dE', 'dI', 'frac_on_E', 'R_eff', ...
    'alpha', 'omega', 'margin', 'Gpeak', 'G_tau', 'G_0p3', 'Gc', ...
    'alpha_u', 'omega_u', 'margin_u', 'Gpeak_u', 'G_tau_u', 'G_0p3_u', 'Gc_u', ...
    'reO', 'l2', 'PF', 'wFF', 'reO_u', 'l2_u', 'PF_u', 'wFF_u', 'henrici', 'henrici_u'};
md_tab = med_table(T, nets, conds, vars);
n_pos = sprintf('States with alpha > 0 (peak gain over [0, 2] s is then modal growth): measured %d of %d, uniform gain %d of %d.', ...
    sum(T.alpha > 0), height(T), sum(T.alpha_u > 0), height(T));
legend_txt = [ ...
    '- d_rms: rms of every route gain (the uniform control''s gain); dE, dI: mean route gain of E, I neurons; frac_on_E: fraction of E with nonzero slope; R_eff: bulk radius (units of W)\n' ...
    '- alpha, omega (1/s): spectral and numerical abscissa of J = (-I + M)/tau_d; margin = omega - alpha\n' ...
    '- Gpeak: max over t in [0, 2] s (step 0.05 s) of ||expm(J t)||_2; G_tau at t = tau_d = 0.1 s; G_0p3 at 0.3 s; Gc: max_t e^{-alpha t} ||expm(J t)||_2 (gain beyond the modal envelope)\n' ...
    '- suffix _u: the uniform-gain control M_u = W rms(d)\n' ...
    '- reO, l2: real parts of the two nonzero eigenvalues of the mean part P (rightmost, other), units of W; PF = ||P||_F; wFF = sqrt(PF^2 - sum |lambda_k|^2)\n' ...
    '- henrici: sqrt(||M||_F^2 - sum |eig M|^2)\n'];
fid = fopen(fullfile(o.fig_dir, 'amp_outlier_table.md'), 'w');
fprintf(fid, '# Amplification and the effective outlier: medians [min, max] over 12 calibration states\n\n');
fprintf(fid, 'Source: %s (amp_outlier_compute, from data/narrative_essential_spectral/outlier_spectra; one seed, n = 500).\n\n', o.res_file);
fprintf(fid, [legend_txt '\n']);
fprintf(fid, '%s\n\n%s\n', n_pos, md_tab);
fclose(fid);
files{end + 1} = fullfile(o.fig_dir, 'amp_outlier_table.md');

side = {
    'Fig_AmpOutlier_budget', ...
    ['Per condition (A no adaptation, B single-timescale, C multiple-timescale): Re lambda_O (rightmost eigenvalue of the mean part P of M, units of W; the green line is the stability boundary 1) against the feed-forward weight w_FF of P. ' ...
     'One point per calibration state; filled = measured gains, open = uniform gain (rms), thin lines join each state''s pair. Shades/markers: reference, EE x1.1, x1.25, x1.5, x1.75 (light to dark). ' ...
     'Dotted arcs: Re lambda_O^2 + w_FF^2 = const through each network''s median uniform-gain point (exact budget only for rank one; P here has rank two, its second eigenvalue is in the table as l2).'], ...
     {'reO', 'wFF', 'PF', 'l2', 'reO_u', 'wFF_u', 'PF_u', 'l2_u', 'd_rms', 'dE', 'dI'};
    'Fig_AmpOutlier_margin', ...
    ['By network (E-to-E dose), conditions in their colors; measured gains solid with filled markers, uniform gain dashed with open markers; markers at the median over 12 states, bars the range. ' ...
     '(A) non-normal margin omega - alpha of J_xx (1/s); (B) peak ||expm(J t)||_2 over [0, 2] s (log; modal growth where alpha > 0); (C) ||expm(J tau_d)||_2 (log); (D) max over [0, 2] s of e^{-alpha t} ||expm(J t)||_2, the gain beyond the rightmost mode''s own growth (log).'], ...
     {'margin', 'margin_u', 'Gpeak', 'Gpeak_u', 'G_tau', 'G_tau_u', 'Gc', 'Gc_u', 'alpha', 'alpha_u'}};
for s = 1:size(side, 1)
    fid = fopen(fullfile(o.fig_dir, [side{s, 1} '.md']), 'w');
    fprintf(fid, '# %s\n\n%s\n\nSource: %s. Code: scripts/explorations/amplification_outlier/amp_outlier_figs.m. Full table: amp_outlier_table.md.\n\n', ...
        side{s, 1}, side{s, 2}, o.res_file);
    fprintf(fid, [legend_txt '\n']);
    fprintf(fid, '%s\n\n%s\n', n_pos, med_table(T, nets, conds, side{s, 3}));
    fclose(fid);
    files{end + 1} = fullfile(o.fig_dir, [side{s, 1} '.md']); %#ok<AGROW>
end
end

function t = med_table(T, nets, conds, vars)
t = sprintf('| network | condition | %s |\n|---|---|%s\n', strjoin(vars, ' | '), repmat('---|', 1, numel(vars)));
for q = 1:numel(nets)
    for c = 1:numel(conds)
        S = T(strcmp(T.network, nets{q}) & strcmp(T.condition, conds{c}), :);
        cells = cellfun(@(v) fmt(S.(v)), vars, 'UniformOutput', false);
        t = [t sprintf('| %s | %s | %s |\n', nets{q}, conds{c}, strjoin(cells, ' | '))]; %#ok<AGROW>
    end
end
end

function s = fmt(x)
s = sprintf('%s [%s, %s]', g3(median(x)), g3(min(x)), g3(max(x)));
end

function s = g3(v)
if abs(v) >= 1e4; s = sprintf('%.2e', v); else; s = sprintf('%.3g', v); end
end

function col = shade(cc, q, nq)
w = 0.65 * (1 - (q - 1) / (nq - 1));     % reference lightest, EE x1.75 darkest
col = cc * (1 - w) + w * [1 1 1];
col = col * (1 - 0.25 * (q - 1) / (nq - 1));
end

function f = save_png(fig, f)
exportgraphics(fig, f, 'Resolution', 300);
end

function v = ifelse(c, a, b)
if c; v = a; else; v = b; end
end

function out = fig_transient_gain(cfg)
% FIG_TRANSIENT_GAIN Transient amplification with adaptation frozen vs active.
%
%   out = FIG_TRANSIENT_GAIN('run_dir', d)
%
% Row 1, one panel per adaptation regime: the worst-case x-in / x-out
% transient gain G(t) = ||P_x Phi(t) P_x'||_2, median with an interquartile
% band over the REGULAR samples of run_transient_gain, for the three
% propagators -- the dendritic block frozen at the state (dotted; the
% conventional rate-network picture), the full Jacobian frozen at the state
% (dashed; adaptation's feedback without drift) and the tangent flow along
% the trajectory (solid; what a perturbation actually does). Log y, a G = 1
% line. The gap between dotted and dashed is adaptation's dynamic negative
% feedback; between dashed and solid, nonstationarity.
%
% Row 2: for the ACTIVE propagator, the four direction readings taken from
% the same n x n block -- worst case, noise average ||.||_F/sqrt(n), the E/I
% difference mode (balanced amplification), the E/I sum mode and the leading
% Lyapunov direction at the sample -- as medians.
%
% A markdown table beside the figure (and on the console) gives, per regime
% and propagator, G_max median [min, max], t_peak, and the peak optimal
% direction's E fraction, participation ratio and |cos| with the leading
% Lyapunov direction.
%
% See also: run_transient_gain, fig_transient_gain_excursions,
%           fig_transient_amplification, SRNNCellTypePairs.transient_gain

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'transient_gain')}, ...
    'transient_gain_data.mat', ...
    'Run run_transient_gain first');
D = load(data_file);
R = D.results; n_cond = numel(R);
variants = D.settings.variants;
v_style = containers.Map({'frozen_x', 'frozen_full', 'active'}, {':', '--', '-'});
v_label = containers.Map({'frozen_x', 'frozen_full', 'active'}, ...
    {'J_{xx} frozen', 'J frozen', 'active'});
readings = {'G_worst', 'G_noise', 'G_ei_diff', 'G_ei_sum', 'G_lyap'};
r_label  = {'worst case', 'noise avg', 'E/I diff', 'E/I sum', 'Lyapunov dir'};
r_grey   = [0 0.15 0.35 0.5 0.65];
r_style  = {'-', '-', '--', '-.', ':'};
i_act = find(strcmp(variants, 'active'));

fig = figure('Color', 'w', 'Position', [80 80 380 * n_cond, 600]);
tl = tiledlayout(fig, 2, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
rows = {};
ax_top = gobjects(1, n_cond); ax_bot = gobjects(1, n_cond);
for i = 1:n_cond
    col = st.condition_color(R(i).name);
    smp = R(i).samples(strcmp({R(i).samples.kind}, 'regular'));
    t = smp(1).t;

    ax = nexttile(tl, [1 1]); hold(ax, 'on');
    h = gobjects(1, numel(variants));
    for v = 1:numel(variants)
        Gm = cell2mat(arrayfun(@(s) s.G_worst(v, :), smp(:), 'UniformOutput', false));
        med = median(Gm, 1); q1 = prctile(Gm, 25, 1); q3 = prctile(Gm, 75, 1);
        fill(ax, [t, fliplr(t)], [q1, fliplr(q3)], col, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        h(v) = plot(ax, t, med, v_style(variants{v}), 'Color', col, 'LineWidth', 1.8);
        rows{end + 1} = table_row(R(i).title, variants{v}, smp, v); %#ok<AGROW>
    end
    yline(ax, 1, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax, 'off');
    set(ax, 'YScale', 'log', 'FontSize', st.tick_fs); box(ax, 'off');
    title(ax, R(i).title, 'FontWeight', 'normal', 'FontSize', st.title_fs);
    xlabel(ax, 'time after perturbation (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'worst-case gain G(t)', 'FontSize', st.label_fs); end
    legend(ax, h, cellfun(@(v) v_label(v), variants, 'UniformOutput', false), ...
        'Location', 'northeast', 'FontSize', 9, 'Box', 'off');
    text(ax, 0.03, 0.03, sprintf('%d states, %d seed(s)', numel(smp), D.settings.n_seeds), ...
        'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'bottom');
    ax_top(i) = ax;

    ax2 = nexttile(tl, [1 1]); hold(ax2, 'on');
    h2 = gobjects(1, numel(readings));
    for k = 1:numel(readings)
        Gm = cell2mat(arrayfun(@(s) s.(readings{k})(i_act, :), smp(:), 'UniformOutput', false));
        if k == 1; c = col; else; c = r_grey(k) * [1 1 1]; end
        if k == 1; lw = 1.8; else; lw = 1.2; end
        h2(k) = plot(ax2, t, median(Gm, 1), r_style{k}, 'Color', c, 'LineWidth', lw);
    end
    yline(ax2, 1, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax2, 'off');
    set(ax2, 'YScale', 'log', 'FontSize', st.tick_fs); box(ax2, 'off');
    xlabel(ax2, 'time after perturbation (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax2, 'active gain by direction', 'FontSize', st.label_fs); end
    legend(ax2, h2, r_label, 'Location', 'northeast', 'FontSize', 9, 'Box', 'off');
    ax_bot(i) = ax2;
end
linkaxes(ax_top, 'y'); linkaxes(ax_bot, 'y');
title(tl, sprintf('Transient gain of a dendritic perturbation, n = %d, T = %g s, %d seed(s), horizon %g s', ...
    R(1).n, D.settings.T, D.settings.n_seeds, D.settings.horizon_s), ...
    'FontWeight', 'normal', 'FontSize', 11);
hdr = '| Condition | Propagator | G_max | t_peak (s) | frac_E(v_opt) | participation | cos(v_opt, v_lyap) |';
vprintf(cfg.verbose, 'verbose', '%s\n|---|---|---|---|---|---|---|\n', hdr);
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});
% The frozen operating point against the trajectory's own rate.
hdr2 = '| Condition | lambda_1 (trajectory) | alpha(J_xx) at state | omega(J_xx) at state | alpha(J) at state |';
rows2 = arrayfun(@(r) op_row(r), R, 'UniformOutput', false);
vprintf(cfg.verbose, 'verbose', '\n%s\n|---|---|---|---|---|\n', hdr2);
vprintf(cfg.verbose, 'verbose', '%s\n', rows2{:});

if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_Transient_Gain';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n|---|---|---|---|---|---|---|\n', hdr);
        fprintf(fid, '%s\n', rows{:});
        fprintf(fid, '\n%s\n|---|---|---|---|---|\n', hdr2);
        fprintf(fid, '%s\n', rows2{:});
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function s = table_row(title, variant, smp, v)
g  = arrayfun(@(x) x.G_max(v), smp);
tp = arrayfun(@(x) x.t_peak(v), smp);
fe = arrayfun(@(x) x.frac_E_opt(v), smp);
pr = arrayfun(@(x) x.participation_opt(v), smp);
al = arrayfun(@(x) x.align_opt_lyap(v), smp);
s = sprintf('| %s | %s | %s | %s | %s | %s | %s |', title, strrep(variant, '_', ' '), ...
    mmm(g, '%.2f'), mmm(tp, '%.2f'), mmm(fe, '%.2f'), mmm(pr, '%.0f'), mmm(al, '%.2f'));
end

function s = op_row(r)
% The frozen operating point vs the trajectory: lambda_1 of the top-K run
% against the spectral abscissa of J_xx and of the full J at the sampled
% states (medians over every sample of the condition; NaN-tolerant).
smp = r.samples;
lam = [r.trials.LLE];
if isfield(smp, 'alpha_full')
    ax = [smp.alpha_xx]; ox = [smp.omega_xx]; af = [smp.alpha_full];
else
    ax = NaN; ox = NaN; af = NaN;
end
s = sprintf('| %s | %s | %s | %s | %s |', r.title, mmm(lam, '%+.3f'), mmm(ax(~isnan(ax)), '%+.2f'), ...
    mmm(ox(~isnan(ox)), '%.1f'), mmm(af(~isnan(af)), '%+.2f'));
end

function s = mmm(v, fmt)
if isempty(v)
    s = '-';
elseif isscalar(v) || max(v) == min(v)
    s = sprintf(fmt, median(v));
else
    s = sprintf([fmt ' [' fmt ', ' fmt ']'], median(v), min(v), max(v));
end
end

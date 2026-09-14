function out = fig_transient_gain_excursions(cfg)
% FIG_TRANSIENT_GAIN_EXCURSIONS Transient gain at excursion onsets vs quiet states.
%
%   out = FIG_TRANSIENT_GAIN_EXCURSIONS('run_dir', d)
%
% The intermittency question: is the network's own transient divergence
% non-normal amplification triggered by its fluctuations? run_transient_gain
% samples states at the ONSETS of local-Lyapunov-rate excursions and at the
% midpoints of QUIET stretches. Row 1, one panel per regime: the ACTIVE
% worst-case gain G(t), median with interquartile band, onset states solid,
% quiet states dashed. Row 2: the |cos| between the optimal input direction
% at the peak and the leading Lyapunov direction at that state, onset vs
% quiet, as box plots. A regime with fewer than two states of a class says so
% in the panel instead of drawing it (the no-adaptation regime rarely has a
% 1-s quiet stretch; the multiple-timescale regime rarely has an onset).
%
% The markdown table beside the figure: per regime, onsets and quiets found
% and used, G_max and alignment medians per class, and a Wilcoxon rank-sum p
% on G_max (onset vs quiet) when both classes have at least three states.
%
% See also: run_transient_gain, fig_transient_gain,
%           SRNNCellTypePairs.excursion_samples

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
i_act = find(strcmp(D.settings.variants, 'active'));
classes = {'onset', 'quiet'}; c_style = {'-', '--'};

fig = figure('Color', 'w', 'Position', [80 80 380 * n_cond, 600]);
tl = tiledlayout(fig, 2, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
rows = cell(1, n_cond);
ax_top = gobjects(1, n_cond);
for i = 1:n_cond
    col = st.condition_color(R(i).name);
    smp = R(i).samples;
    S = struct();
    for c = 1:2
        S.(classes{c}) = smp(strcmp({smp.kind}, classes{c}));
    end
    n_on = numel(S.onset); n_qu = numel(S.quiet);

    ax = nexttile(tl, [1 1]); hold(ax, 'on');
    h = gobjects(1, 2); drawn = false(1, 2);
    for c = 1:2
        cs = S.(classes{c});
        if numel(cs) < 2; continue; end
        t = cs(1).t;
        Gm = cell2mat(arrayfun(@(s) s.G_worst(i_act, :), cs(:), 'UniformOutput', false));
        fill(ax, [t, fliplr(t)], [prctile(Gm, 25, 1), fliplr(prctile(Gm, 75, 1))], col, ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none');
        h(c) = plot(ax, t, median(Gm, 1), c_style{c}, 'Color', col, 'LineWidth', 1.8);
        drawn(c) = true;
    end
    yline(ax, 1, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax, 'off');
    set(ax, 'YScale', 'log', 'FontSize', st.tick_fs); box(ax, 'off');
    title(ax, R(i).title, 'FontWeight', 'normal', 'FontSize', st.title_fs);
    xlabel(ax, 'time after perturbation (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'active worst-case gain G(t)', 'FontSize', st.label_fs); end
    lab = {sprintf('onset (n = %d)', n_on), sprintf('quiet (n = %d)', n_qu)};
    if any(drawn)
        legend(ax, h(drawn), lab(drawn), 'Location', 'northeast', 'FontSize', 9, 'Box', 'off');
    end
    if ~all(drawn)
        miss = classes(~drawn);
        text(ax, 0.03, 0.03, sprintf('fewer than 2 %s states', strjoin(miss, ' / ')), ...
            'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'bottom');
    end
    ax_top(i) = ax;

    ax2 = nexttile(tl, [1 1]); hold(ax2, 'on');
    al = struct('onset', arrayfun(@(s) s.align_opt_lyap(i_act), S.onset), ...
        'quiet', arrayfun(@(s) s.align_opt_lyap(i_act), S.quiet));
    if any(drawn)
        x = []; y = [];
        for c = 1:2
            if ~drawn(c); continue; end
            y = [y; al.(classes{c})(:)]; x = [x; c * ones(numel(al.(classes{c})), 1)]; %#ok<AGROW>
        end
        boxchart(ax2, x, y, 'BoxFaceColor', col, 'MarkerColor', col, 'BoxWidth', 0.5);
        scatter(ax2, x + 0.08 * randn(size(x)), y, 14, col, 'filled', 'MarkerFaceAlpha', 0.5);
    else
        text(ax2, 0.5, 0.5, 'no excursion samples', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    hold(ax2, 'off');
    set(ax2, 'FontSize', st.tick_fs, 'YLim', [0 1], 'XTick', 1:2, 'XTickLabel', classes, 'XLim', [0.4 2.6]); box(ax2, 'off');
    if i == 1; ylabel(ax2, '|cos(v_{opt}, v_{Lyap})| at peak', 'FontSize', st.label_fs); end

    rows{i} = table_row(R(i), S, al, i_act);
end
linkaxes(ax_top(isgraphics(ax_top)), 'y');
title(tl, sprintf('Active transient gain at excursion onsets vs quiet states, n = %d, T = %g s, %d seed(s)', ...
    R(1).n, D.settings.T, D.settings.n_seeds), 'FontWeight', 'normal', 'FontSize', 11);
hdr = '| Condition | onsets found / used | quiets found / used | G_max onset | G_max quiet | p (rank-sum) | align onset | align quiet |';
vprintf(cfg.verbose, 'verbose', '%s\n|---|---|---|---|---|---|---|---|\n', hdr);
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});

if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_Transient_Gain_Excursions';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n|---|---|---|---|---|---|---|---|\n', hdr);
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function s = table_row(Ri, S, al, i_act)
g_on = arrayfun(@(x) x.G_max(i_act), S.onset);
g_qu = arrayfun(@(x) x.G_max(i_act), S.quiet);
p = 'n/a';
if numel(g_on) >= 3 && numel(g_qu) >= 3 && exist('ranksum', 'file') == 2
    p = sprintf('%.3g', ranksum(g_on, g_qu));
end
s = sprintf('| %s | %d / %d | %d / %d | %s | %s | %s | %s | %s |', Ri.title, ...
    Ri.n_onset_found, numel(S.onset), Ri.n_quiet_found, numel(S.quiet), ...
    mmm(g_on, '%.2f'), mmm(g_qu, '%.2f'), p, mmm(al.onset, '%.2f'), mmm(al.quiet, '%.2f'));
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

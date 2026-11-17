function out = fig_mc_noise_sweeps(cfg)
% FIG_MC_NOISE_SWEEPS Memory capacity (and lambda_1) of the preset network against each kind of noise: one panel per 1-D noise sweep (block C, FD-029).
%
%   out = FIG_MC_NOISE_SWEEPS('run_dir', d)        % the grids noise1d_wiener and noise1d_release
%   out = FIG_MC_NOISE_SWEEPS('run_dir', d, 'grids', {'a', 'b'})
%   out = FIG_MC_NOISE_SWEEPS('data', {D1, D2}, 'first_level', 1, 'save', false)
%
% Tom (2026-11-16): "a figure at just the preset but with each type of
% noise varied ... two 1D sweeps, which would be much faster and still tell
% important info about the effect of noise on MC". Each sweep is a grid of
% the block C stage (run_mc_gain_sweep) with ONE level on its first axis
% (the preset's gain) and a noise amplitude on its second; this figure
% draws one column per sweep:
%   row 1   memory capacity against the noise level
%   row 2   lambda_1 with the task input on, the same way (cfg.show_lle,
%           default true), with a thin gray line at zero
% One line per condition in the condition colors (mc_condition_colors,
% titles by titles_by_name): the median over networks (seeds), a band from
% the smallest to the largest network (with two to five networks a range
% says more than quartiles), and every network as a small point. Failed
% jobs are left out and counted in the caption.
% Two options change that for a figure that wants it (the defaults are the
% above): cfg.band = 'iqr' draws the band from the 25th to the 75th
% percentile over networks instead, prctile(v, [25 75]) (base MATLAB;
% it interpolates linearly between the sorted values, the k-th of
% n sitting at the 100 (k - 0.5) / n percentile, so with 5 networks the
% band runs from 3/4 of the way from the smallest to the second smallest
% up to 1/4 of the way from the fourth to the largest), and
% cfg.show_points = false leaves the per-network points out. The caption
% says which band is drawn; the table keeps the smallest and largest.
%
% THE NOISE AXIS is drawn as the maps draw it (mc_axis_spacing): a zero
% level followed by a log ladder ('zero_log') at its log ladder, the zero one
% log step below the lowest positive level and labeled 0, the axis broken
% there (said in the caption); a log ladder in log; anything else linear.
% The ticks are the sweep's levels.
%
% WHAT EACH PANEL HOLDS FIXED is read from the result's settings (its
% model_args), so the caption states it rather than assuming it: a
% spontaneous release rate above zero (release_rate_all or release_rate,
% the "leak" r_s) is named in the panel's title, and the caption says the
% zero of that panel is the leak alone, not the noise-free network.
%
% ONE LEVEL OF THE FIRST AXIS. A 1-D sweep has one. A result with several
% (a 2-D grid, for example the noise grids of FD-027) is drawn at the level
% nearest cfg.first_level, which must then be given, and the caption says
% which level was used; this is how a preview is drawn from a 2-D grid.
%
% THE TABLE <fig_tag>_table.md beside the figure (when saved): per sweep,
% condition and noise level, the number of networks, the median, smallest
% and largest memory capacity, and the median lambda_1.
%
% The results are taken through resolve_mc_gain_sweep (cfg.data, a cell of
% results or tables, or cfg.run_dir with cfg.grids naming the grids).
% cfg.note adds a line to the caption (for example "preview from ...").
%
% See also: fig_mc_gain_bias_maps, mc_gain_sweep_grid, mu7chain_gainbias_config,
%           mc_axis_spacing, mc_sweep_quantity

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.run_dir     (1,:) char    = ''
    cfg.grids                     = {'noise1d_wiener', 'noise1d_release'}
    cfg.data                      = {}
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                = {}
    cfg.first_level double        = []   % with several first-axis levels: the one drawn (nearest)
    cfg.show_lle    (1,1) logical = true
    cfg.band        (1,:) char {mustBeMember(cfg.band, {'range', 'iqr'})} = 'range'   % the band over networks
    cfg.show_points (1,1) logical = true   % every network as a small point
    cfg.note        (1,:) char    = ''
    cfg.fig_tag     (1,:) char    = 'Fig_MC_noise_sweeps'
end

setup_paths();
LEAK_NAMES = {'release_rate_all', 'release_rate'};
BAND_ALPHA = 0.18;
POINT_SIZE = 12;

%% The sweeps
if isempty(cfg.data)
    grids = cellstr(cfg.grids);
    sources = cellfun(@(g) resolve_mc_gain_sweep('run_dir', cfg.run_dir, 'grid', g), grids, ...
        'UniformOutput', false);
else
    data = cfg.data;
    if ~iscell(data); data = {data}; end
    sources = cellfun(@(d) resolve_mc_gain_sweep('data', d), data, 'UniformOutput', false);
end
np = numel(sources);
S = cell(1, np);
for p = 1:np
    S{p} = one_sweep(sources{p}, cfg.first_level, LEAK_NAMES);
end
conds = select_conditions(S{1}.D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();
cond_titles = titles_by_name(conds);
cond_colors = mc_condition_colors(conds);

%% Draw
rows = {'mc_total'};
row_labels = {'Memory capacity'};
if cfg.show_lle
    rows{end+1} = 'lle_input';
    row_labels{end+1} = '\lambda_1 (1/s), task input on';
end
nr = numel(rows);
fig = figure('Color', 'w', 'Position', [60 60 430 * np + 160, 300 * nr + 170]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, nr, np, 'TileSpacing', 'compact', 'Padding', 'compact');
h = gobjects(1, nc);
md = {'| sweep | condition | level | networks | MC median | MC min | MC max | lambda_1 median (1/s) |', ...
      '|---|---|---:|---:|---:|---:|---:|---:|'};
n_failed = 0;
for r = 1:nr
    for p = 1:np
        P = S{p};
        ax = nexttile(tl, (r - 1) * np + p); hold(ax, 'on');
        if strcmp(rows{r}, 'lle_input')
            plot(ax, P.x([1 end]), [0 0], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
        for c = 1:nc
            T = P.T(strcmp(P.T.condition, conds{c}), :);
            if r == 1; n_failed = n_failed + nnz(T.failed); end
            T = T(~T.failed, :);
            [med, lo, hi, blo, bhi, n] = deal(nan(size(P.levels)));   % lo, hi: the range; blo, bhi: the band
            for k = 1:numel(P.levels)
                v = T.(rows{r})(T.(P.column) == P.levels(k));
                v = v(isfinite(v));
                n(k) = numel(v);
                if isempty(v); continue; end
                med(k) = median(v); lo(k) = min(v); hi(k) = max(v);
                if strcmp(cfg.band, 'iqr')
                    q = prctile(v, [25 75]); blo(k) = q(1); bhi(k) = q(2);
                else
                    blo(k) = lo(k); bhi(k) = hi(k);
                end
                if cfg.show_points
                    scatter(ax, repmat(P.x(k), size(v)), v, POINT_SIZE, cond_colors(c, :), 'filled', ...
                        'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
                end
            end
            ok = isfinite(med);
            if nnz(ok) > 1
                fill(ax, [P.x(ok), fliplr(P.x(ok))], [blo(ok), fliplr(bhi(ok))], cond_colors(c, :), ...
                    'FaceAlpha', BAND_ALPHA, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
            h(c) = plot(ax, P.x(ok), med(ok), '-o', 'Color', cond_colors(c, :), 'LineWidth', 2, ...
                'MarkerSize', 4, 'MarkerFaceColor', cond_colors(c, :));
            if r == 1
                lle = arrayfun(@(v) median(T.lle_input(T.(P.column) == v & isfinite(T.lle_input))), P.levels);
                for k = 1:numel(P.levels)
                    md{end+1} = sprintf('| %s | %s | %.4g | %d | %.3f | %.3f | %.3f | %.3f |', ...
                        P.name, conds{c}, P.levels(k), n(k), med(k), lo(k), hi(k), lle(k)); %#ok<AGROW>
                end
            end
        end
        set(ax, 'XTick', P.x, 'XTickLabel', P.ticklabels, 'XTickLabelRotation', 45, ...
            'TickDir', 'out', 'Box', 'off', 'FontSize', st.tick_fs - 4);
        xlim(ax, P.x([1 end]) + [-0.5 0.5] * P.step);
        if r == 1
            title(ax, P.title, 'FontSize', st.title_fs - 2, 'FontWeight', 'normal');
            ylim(ax, [0, max([ax.YLim(2), eps])]);
        end
        if r == nr; xlabel(ax, P.xlabel, 'FontSize', st.label_fs - 3); end
        if p == 1; ylabel(ax, row_labels{r}, 'FontSize', st.label_fs - 3); end
    end
end
lg = legend(h(isgraphics(h)), cond_titles(isgraphics(h)), 'Orientation', 'horizontal', ...
    'FontSize', st.tick_fs - 4, 'Box', 'off');
lg.Layout.Tile = 'south';

%% Caption
n_networks = max(cellfun(@(P) numel(unique(P.T.seed)), S));
band_text = 'smallest to largest network';
if strcmp(cfg.band, 'iqr'); band_text = 'interquartile range over networks (prctile, 25 to 75)'; end
points_text = '';
if cfg.show_points; points_text = '; points: each network'; end
caption = {sprintf('Preset network at %s. Line: median over %d networks; band: %s%s.', ...
    S{1}.at_text, n_networks, band_text, points_text)};
broken = cellfun(@(P) P.zero_log, S);
if any(broken)
    caption{end+1} = 'Noise axes: 0, then log-spaced levels; the axis is broken between 0 and the lowest positive level.';
end
for p = find(cellfun(@(P) ~isempty(P.leak), S))
    caption{end+1} = sprintf(['%s: every level includes the leak r_s = %g (spontaneous release), so its 0 ' ...
        'is the leak alone, not the noise-free network.'], capitalized(S{p}.Q.symbol), S{p}.leak); %#ok<AGROW>
end
if n_failed > 0
    caption{end+1} = sprintf('%d failed job(s) left out.', n_failed);
end
if ~isempty(cfg.note); caption{end+1} = cfg.note; end
title(tl, caption, 'FontSize', st.tick_fs - 5);

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', {cellfun(@(P) P.D.source, S, 'UniformOutput', false)}, ...
    'table', {md'}, 'n_networks', n_networks);
vprintf(cfg.verbose, 'verbose', '%s\n', md{:});
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);   % deletes <fig_tag>*, so the table is written after
    table_name = [cfg.fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    fprintf(fid, '# %s\n\n%s\n\n', cfg.fig_tag, strjoin(caption, ' '));
    fprintf(fid, '%s\n', md{:});
    fclose(fid);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
    if ~any(endsWith(out.files, table_name)); out.files{end+1} = table_name; end
end
end

%% ------------------------------------------------------------------------
function P = one_sweep(D, first_level, leak_names)
% One panel's data: the rows at one level of the first axis, the noise
% axis's positions and labels, and what the run held fixed.
P = struct('D', D);
first_col = D.first_axis.column;
levels1 = D.first_levels;
if isscalar(levels1)
    chosen = levels1;
elseif isempty(first_level)
    error('fig_mc_noise_sweeps:SeveralFirstLevels', ...
        ['The result has %d levels of %s; a noise sweep has one. Give first_level to draw ' ...
         'the one nearest it.'], numel(levels1), first_col);
else
    [~, k] = min(abs(levels1 - first_level));
    chosen = levels1(k);
end
P.T = D.T(D.T.(first_col) == chosen, :);
Q1 = mc_sweep_quantity(first_col, D);
P.at_text = Q1.at(chosen);
if ~isscalar(levels1)
    % Not Q1.at, which rounds to three digits: the level used must read as
    % itself (1.003 is not 1).
    P.at_text = sprintf('%s %.4g, its level nearest %g', Q1.symbol, chosen, first_level);
end
P.column = D.second_axis.column;
P.Q = mc_sweep_quantity(P.column, D);
P.name = D.second_axis.property;
P.levels = D.second_levels(:)';
[spacing, to_position] = mc_axis_spacing(P.levels);
switch spacing
    case 'zero_log'
    case 'log';    to_position = @log10;
    otherwise;     to_position = @(v) v;
end
P.zero_log = strcmp(spacing, 'zero_log');
P.x = to_position(P.levels);
P.step = max([median(diff(P.x)), eps]);
P.ticklabels = arrayfun(@(v) sprintf('%.3g', v), P.levels, 'UniformOutput', false);
P.xlabel = capitalized(P.Q.label);
% The leak, from the run's own settings.
P.leak = [];
if isfield(D.settings, 'model_args')
    args = D.settings.model_args;
    for name = leak_names
        k = find(strcmp(args(1:2:end), name{1}), 1, 'last');
        if ~isempty(k) && isnumeric(args{2 * k}) && any(args{2 * k} > 0)
            P.leak = max(args{2 * k});
        end
    end
end
P.title = capitalized(P.Q.symbol);
if ~isempty(P.leak)
    P.title = sprintf('%s, at the leak r_s = %g', P.title, P.leak);
end
end

function s = capitalized(s)
if ~isempty(s); s(1) = upper(s(1)); end
end

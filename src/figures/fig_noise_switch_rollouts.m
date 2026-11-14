function out = fig_noise_switch_rollouts(cfg)
% FIG_NOISE_SWITCH_ROLLOUTS Copies of one network that share the noise, then stop sharing it.
%
%   out = FIG_NOISE_SWITCH_ROLLOUTS('run_dir', d)
%   out = FIG_NOISE_SWITCH_ROLLOUTS('data_file', f, 'out_dir', o)
%
% Draws what run_noise_switch_rollouts saved: copies of one network that see
% the identical Wiener path up to a switch time and independent noise seeds
% after it. It reads one run directory (<run_dir>/noise_switch_rollouts/)
% or the file given, and simulates nothing.
%
% ONE TIME-SERIES FIGURE PER CONDITION, tag
% Fig_noise_switch_timeseries_<condition>_copies. Rows: x, r, the synaptic
% output theta, the adaptation term (c/K) sum_k a_k and the depression
% factor (b_1 on a chain route), whichever the condition has; a condition
% without SFA or depression has fewer rows. Columns: the cell types, each
% with the neurons the stage stored (the first 25 of each type, as
% fig_example_timeseries chooses them). The TWO COPIES are overlaid in two
% color families, copy 1 blue to green and copy 2 yellow through orange to
% red, each neuron keeping one shade of its copy's family so a trace can be
% followed across the switch. Copy 2 is drawn on top, so before the switch,
% where the copies are bit-identical, only its colors show; the title says
% so. A vertical line marks the switch. The time axis covers the last
% cfg.t_before seconds of the shared part and the whole rollout.
%
% THE SUMMARY FIGURE, tag Fig_noise_switch_divergence. Top: per state
% block, the RMS difference between copies divided by the block's own
% temporal standard deviation (rms_diff_rel) against time since the
% switch, one line per condition in its usual color, both axes
% logarithmic (the early rise spans two decades of time; a diffusive rise
% is a straight line of slope 1/2 there). The dotted line at sqrt(2) is
% where two independent samples of a stationary fluctuation sit. Bottom:
% per block and condition the fitted early growth rate, with each
% condition's lambda_1 as a dashed line in its color, the time to half
% saturation, and the saturation level; an open marker means the block had
% not saturated within the rollout. Every number is also written to
% Fig_noise_switch_divergence_table.md.
%
% The tag of a time-series figure ends in _copies because
% save_figure_stable deletes <tag>* before saving and one condition name
% can be a prefix of another (sfa1_std1 and sfa1_std1_stf1); with the
% suffix no tag is a prefix of another.
%
% The two copies' colors are defined here (copy_colors), not in
% src/plotting: they mean "copy", which nothing else in the repo draws.
%
% cfg.conditions: names drawn, in order; empty draws all the run holds
% (select_conditions). Titles and colors come from the condition NAME.
%
% See also: run_noise_switch_rollouts, fig_example_timeseries, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''      % unused; the preset is recorded in the .mat
    cfg.conditions                = {}      % names to draw, in order; empty = all the run holds
    cfg.t_before    (1,1) double {mustBePositive} = 3   % seconds of the shared part shown
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'noise_switch_rollouts')}, ...
    'noise_switch_rollouts_data.mat', 'Run run_noise_switch_rollouts first');
D = load(data_file);
[~, idx] = select_conditions({D.results.name}, cfg.conditions);
R = D.results(idx);
titles = titles_by_name({R.name}, {R.title});
vis = 'on'; if ~cfg.visible, vis = 'off'; end

out = struct('figs', gobjects(1, 0), 'files', {{}}, 'source', data_file, 'tags', {{}});

%% One time-series figure per condition
for i = 1:numel(R)
    fig = timeseries_figure(R(i), titles{i}, D.settings, cfg.t_before, st, vis);
    out.figs(end + 1) = fig;
    out.tags{end + 1} = sprintf('Fig_noise_switch_timeseries_%s_copies', R(i).name);
end

%% The summary figure and its table
[fig, table_lines] = divergence_figure(R, titles, D.settings, st, vis);
out.figs(end + 1) = fig;
out.tags{end + 1} = 'Fig_noise_switch_divergence';
vprintf(cfg.verbose, 'verbose', '%s\n', table_lines{:});

%% Save
if cfg.save
    for k = 1:numel(out.figs)
        save_figure_stable(out_dir, out.tags{k}, out.figs(k));
        out.files = [out.files, existing_outputs(out_dir, out.tags{k})];
    end
    table_name = 'Fig_noise_switch_divergence_table.md';
    fid = fopen(fullfile(out_dir, table_name), 'w');
    fprintf(fid, '%s\n', table_lines{:});
    fclose(fid);
    out.files{end + 1} = table_name;
end
end

% =========================================================================
function fig = timeseries_figure(R, ttl, settings, t_before, st, vis)
% Blocks down, cell types across, two copies overlaid per neuron.
E = R.example;
blocks = R.blocks;
names  = E.type_names;
C      = numel(names);
n_rows = numel(blocks);
t_sw   = settings.t_shared;
t      = E.t;
show   = t >= t_sw - t_before - 1e-9;

fig = figure('Color', 'w', 'Visible', vis, 'Position', [60 40 560 * C + 120, 190 * n_rows + 150]);
tl  = tiledlayout(fig, n_rows, C, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_all = gobjects(n_rows, C);
for row = 1:n_rows
    blk = blocks{row};
    S   = E.series.(blk);                        % neurons x time x 2 copies
    for q = 1:C
        ax = nexttile(tl); ax_all(row, q) = ax; hold(ax, 'on');
        rows_q = find(E.type.(blk) == q);
        if isempty(rows_q)
            text(ax, 0.5, 0.5, sprintf('%s neurons have none', names{q}), 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
            set(ax, 'YTick', []);
        end
        col = copy_colors(numel(rows_q));
        for copy = 1:2                           % copy 2 last, so it lies on top
            for k = 1:numel(rows_q)
                plot(ax, t(show), S(rows_q(k), show, copy), '-', ...
                    'Color', col{copy}(k, :), 'LineWidth', 0.6);
            end
        end
        xline(ax, t_sw, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1);
        hold(ax, 'off');
        xlim(ax, [t_sw - t_before, t(end)]);
        set(ax, 'FontSize', st.tick_fs - 3); box(ax, 'off');
        if q == 1, ylabel(ax, R.block_labels.(blk), 'FontSize', st.label_fs - 2); end
        if row < n_rows, set(ax, 'XTickLabel', []); else, xlabel(ax, 'time (s)', 'FontSize', st.label_fs - 2); end
        if row == 1
            title(ax, sprintf('%s neurons (%d of %d shown)', names{q}, numel(rows_q), R.n_per_type(q)), ...
                'FontWeight', 'normal', 'FontSize', st.title_fs - 2);
            % The switch, labeled on both sides at the top of the panel
            % (on white, since the traces fill it).
            x_sw = t_before / (t(end) - t_sw + t_before);
            text(ax, x_sw, 0.99, ' own noise seeds \rightarrow', 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.25 0.25 0.25], ...
                'BackgroundColor', 'w', 'Margin', 1);
            text(ax, x_sw, 0.99, '\leftarrow shared noise ', 'Units', 'normalized', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8, ...
                'Color', [0.25 0.25 0.25], 'BackgroundColor', 'w', 'Margin', 1);
        end
    end
    linkaxes(ax_all(row, :), 'y');
end
% Legend: one proxy line per copy, mid-family.
ax = ax_all(1, C); hold(ax, 'on');
mid = copy_colors(3);
h1 = plot(ax, NaN, NaN, '-', 'Color', mid{1}(2, :), 'LineWidth', 2);
h2 = plot(ax, NaN, NaN, '-', 'Color', mid{2}(2, :), 'LineWidth', 2);
hold(ax, 'off');
lg = legend(ax, [h1 h2], {'copy 1', 'copy 2 (on top)'}, 'Orientation', 'horizontal', ...
    'FontSize', 9, 'Box', 'off');
lg.Layout.Tile = 'south';

title(tl, {sprintf('%s:  two copies of one network, \\lambda_1 = %+.2f s^{-1}', ttl, R.lambda_1), ...
    sprintf(['identical noise until t = %g s, then each copy on its own noise seed; ' ...
    'before the switch the copies coincide bit for bit, so only copy 2 shows'], t_sw), ...
    sprintf('%s   (n = %d, \\sigma_u = %g, seeds %s)', strrep(settings.preset_name, '_', '\_'), ...
    R.n, settings.sigma_u_noise, mat2str(settings.rng_seeds))}, ...
    'FontWeight', 'normal', 'FontSize', 9);
end

% =========================================================================
function col = copy_colors(n)
% Per-neuron colors of the two copies: col{1} runs blue to green, col{2}
% yellow through orange to red, n rows each. Neuron k has the same position
% along both ramps, dealt out by the golden ratio so that neighboring
% neurons land far apart and a trace can be told from the next.
ramp1 = [0.05 0.20 0.70; 0.00 0.55 0.75; 0.15 0.70 0.30];
ramp2 = [0.95 0.80 0.10; 0.95 0.50 0.10; 0.75 0.10 0.10];
pos = mod((0:n - 1)' * 0.6180339887, 1);
anchors = linspace(0, 1, 3);
col = {interp1(anchors, ramp1, pos), interp1(anchors, ramp2, pos)};
if n == 0, col = {zeros(0, 3), zeros(0, 3)}; end
end

% =========================================================================
function [fig, lines_md] = divergence_figure(R, titles, settings, st, vis)
% Top row: rms_diff_rel against time since the switch, one panel per block.
% Bottom row: growth rate (with lambda_1), half time and saturation.
all_blocks = {'x', 'r', 'theta', 'a', 'b'};
have = cellfun(@(b) any(arrayfun(@(r) ismember(b, r.blocks), R)), all_blocks);
blocks = all_blocks(have);
block_ticks = strrep(blocks, 'theta', '\theta');
n_blk  = numel(blocks);
n_cond = numel(R);
colors = zeros(n_cond, 3); short = cell(1, n_cond);
for i = 1:n_cond
    colors(i, :) = [0.5 0.5 0.5]; short{i} = titles{i};
    if isKey(st.condition_color, R(i).name), colors(i, :) = st.condition_color(R(i).name); end
    if isKey(st.condition_short, R(i).name), short{i} = st.condition_short(R(i).name); end
end

n_col = max(n_blk, 3);
fig = figure('Color', 'w', 'Visible', vis, 'Position', [60 60 300 * n_col + 80, 720]);
tl  = tiledlayout(fig, 2, 6 * n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- top: the curves ------------------------------------------------------
ax_top = gobjects(1, n_blk);
h_cond = gobjects(1, n_cond);
for k = 1:n_blk
    ax = nexttile(tl, 6 * (k - 1) + 1, [1 6]); ax_top(k) = ax; hold(ax, 'on');
    % The block's label, or a general one where the conditions drawn label
    % it differently (b is b_1 on a chain route, a product elsewhere).
    label = '';
    for i = 1:n_cond
        if ~ismember(blocks{k}, R(i).blocks), continue; end
        S = R(i).stats.(blocks{k});
        after = R(i).t > settings.t_shared + 1e-9;
        h_cond(i) = plot(ax, R(i).t(after) - settings.t_shared, S.rms_diff_rel(after), '-', ...
            'Color', colors(i, :), 'LineWidth', 1.6);
        this = R(i).block_labels.(blocks{k});
        if isempty(label), label = this; end
        if ~strcmp(label, this), label = 'depression factor b'; end
    end
    yline(ax, sqrt(2), ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    hold(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', st.tick_fs - 3); box(ax, 'off');
    xlim(ax, [1 / settings.store_fs, settings.t_rollout]);
    title(ax, label, 'FontWeight', 'normal', 'FontSize', st.title_fs - 2);
    xlabel(ax, 'time since the switch (s)', 'FontSize', st.label_fs - 3);
    if k == 1
        ylabel(ax, {'RMS difference between copies', '/ temporal std of the block'}, ...
            'FontSize', st.label_fs - 3);
    end
end
linkaxes(ax_top, 'y');
text(ax_top(1), 1 / settings.store_fs, sqrt(2), ' independent: \surd2', 'VerticalAlignment', 'bottom', ...
    'FontSize', 8, 'Color', [0.3 0.3 0.3]);
ok = isgraphics(h_cond);
entries = arrayfun(@(i) sprintf('%s  (\\lambda_1 = %+.2f s^{-1})', short{i}, R(i).lambda_1), ...
    1:n_cond, 'UniformOutput', false);
lg = legend(ax_top(1), h_cond(ok), entries(ok), 'Orientation', 'horizontal', 'FontSize', 9, 'Box', 'off');
lg.Layout.Tile = 'south';

% --- bottom: the fitted numbers, per block and condition --------------------
quantity = {'growth_rate', 't_half', 'saturation'};
q_label  = {'early growth rate (s^{-1})', 'time to half saturation (s)', 'saturation level'};
span = 2 * n_col;                                % three panels across the row
offsets = linspace(-0.25, 0.25, n_cond);
for p = 1:3
    ax = nexttile(tl, 6 * n_col + (p - 1) * span + 1, [1 span]); hold(ax, 'on');
    for i = 1:n_cond
        for k = 1:n_blk
            if ~ismember(blocks{k}, R(i).blocks), continue; end
            F = R(i).stats.(blocks{k}).fit;
            face = colors(i, :);
            if ~F.saturated, face = 'none'; end
            plot(ax, k + offsets(i), F.(quantity{p}), 'o', 'Color', colors(i, :), ...
                'MarkerFaceColor', face, 'MarkerSize', 7, 'LineWidth', 1.2);
        end
        if p == 1
            yline(ax, R(i).lambda_1, '--', 'Color', colors(i, :), 'LineWidth', 1);
        end
    end
    if p == 3, yline(ax, sqrt(2), ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1); end
    hold(ax, 'off');
    set(ax, 'XTick', 1:n_blk, 'XTickLabel', block_ticks, 'FontSize', st.tick_fs - 3); box(ax, 'off');
    xlim(ax, [0.5, n_blk + 0.5]);
    if p == 2, set(ax, 'YScale', 'log'); end
    ylabel(ax, q_label{p}, 'FontSize', st.label_fs - 3);
    if p == 1
        title(ax, 'dashed: \lambda_1 of each condition', 'FontWeight', 'normal', 'FontSize', 9);
    elseif p == 2
        title(ax, 'open marker: not saturated within the rollout', 'FontWeight', 'normal', 'FontSize', 9);
    else
        title(ax, 'dotted: \surd2', 'FontWeight', 'normal', 'FontSize', 9);
    end
end

title(tl, {sprintf(['How fast copies of one network drift apart once their noise is independent ' ...
    '(%d copies, shared %g s, rollout %g s)'], settings.n_copies, settings.t_shared, settings.t_rollout), ...
    sprintf('%s   (n = %d, \\sigma_u = %g, %s at %g Hz)', strrep(settings.preset_name, '_', '\_'), ...
    settings.n, settings.sigma_u_noise, settings.ode_solver, settings.fs)}, ...
    'FontWeight', 'normal', 'FontSize', 10);

lines_md = table_lines(R, titles, blocks, settings);
end

% =========================================================================
function L = table_lines(R, titles, blocks, settings)
% The numbers behind the summary figure, as a markdown table.
L = {sprintf('# Noise-switch rollouts: %s', settings.preset_name), '', ...
    sprintf(['%d copies; identical noise for %g s, then independent seeds for %g s; %s at %g Hz; ' ...
    'statistics at %g Hz. The RMS difference between copies is divided by the block''s temporal ' ...
    'standard deviation over [%g, %g] s; two independent samples of a stationary fluctuation ' ...
    'give 1.414. Growth rate: slope of the log of that ratio against time, from the first sample ' ...
    'after the switch to half saturation; rise exponent: the same samples in log-log ' ...
    '(0.5 is diffusion). Saturation: mean over the last fifth of the rollout.'], ...
    settings.n_copies, settings.t_shared, settings.t_rollout, settings.ode_solver, settings.fs, ...
    settings.store_fs, settings.reference_window(1), settings.reference_window(2)), '', ...
    '| condition | lambda_1 (1/s) | block | temporal std | growth rate (1/s) | rise exponent | t_half (s) | saturation | saturated | correlation at end |', ...
    '|---|---|---|---|---|---|---|---|---|---|'};
yes_no = {'no', 'yes'};
for i = 1:numel(R)
    for k = 1:numel(blocks)
        if ~ismember(blocks{k}, R(i).blocks), continue; end
        S = R(i).stats.(blocks{k}); F = S.fit;
        L{end + 1} = sprintf('| %s (%s) | %+.3f | %s | %.4g | %.3g | %.2f | %.3g | %.3g | %s | %.2f |', ...
            titles{i}, R(i).name, R(i).lambda_1, blocks{k}, S.scale, F.growth_rate, F.rise_exponent, ...
            F.t_half, F.saturation, yes_no{F.saturated + 1}, S.corr(end)); %#ok<AGROW>
    end
end
end

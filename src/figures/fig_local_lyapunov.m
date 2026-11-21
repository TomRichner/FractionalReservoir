function out = fig_local_lyapunov(cfg)
% FIG_LOCAL_LYAPUNOV Local Lyapunov exponents and expansive transients under a stimulus staircase.
%
%   out = FIG_LOCAL_LYAPUNOV('run_dir', d)
%
% One column per condition of run_local_lyapunov's output, six rows: the
% external input (every neuron); the K local rates with Benettin's over the
% top; the K accumulating exponents with Benettin's finite-time curve; the
% share of the accumulation window each local rate is positive; how many of
% the K are positive at each moment; and the local positive QR-rate sum
% sum_k max(local_k, 0)/log(2) in bits/s. K is the number of directions the
% run tracked (settings.K), whatever it is.
%
% THE WINDOW (Plan 2, package C, 2026-11-21). Every share, median and mean,
% and the x range of the accumulating row, is taken over the stage's own
% accumulation window, settings.lya_T_interval, the window of the exponents
% the table prints. Until then the panels and percentages used [T/2, T]
% (20 to 40 s) beside exponents accumulated over 0 to 40 s. A run saved
% before the stage recorded its window falls back to [T/2, T], and the
% table says so. A dotted vertical line marks the start of the window only
% where it opens after the first displayed time.
%
% NO TEXT ON THE PANELS (figure plan, 2026-11-21): the per-panel titles
% with percentages, the key sentence and the window means, which collided
% with each other and with the y labels, are in Fig_Local_Lyapunov_table.md
% with the run's settings. One key for the line styles sits under the
% panels. Condition titles in sentence case and in the condition's color.
% The local sum describes transient growth in tracked directions, not KS
% entropy. The legacy h_KS_bits field is retained for saved-data
% compatibility.
%
% PNG and .fig only, NO SVG: the figure holds ~15 million points at K = 100
% and the vector export took minutes (2026-09-14).
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions (FD-006), condition names drawn
% as columns in the order given; empty draws every condition the .mat holds,
% as before. Resolved by select_conditions; titles from the NAME, not the
% title the stage saved.
%
% See also: run_local_lyapunov, fig_local_vs_finite_lle, fig_lyapunov_spectrum

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; recorded in the .mat
    cfg.conditions                = {}    % names to draw, in order; empty = all the run holds (select_conditions)
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
data_file = resolve_data_file(cfg.data_file, cfg.run_dir, {fullfile(cfg.run_dir, 'local_lyapunov')}, ...
    'local_lyapunov_data.mat', 'Run run_local_lyapunov first');
D = load(data_file);
% The selected conditions in the selected order (FD-006), titled by NAME.
[names, idx] = select_conditions({D.results.name}, cfg.conditions);
R = D.results(idx); n_cond = numel(R); T = D.settings.T; K = D.settings.K;
titles = grouped_condition_label(names, 'sentence');
% The stage's accumulation window; [T/2, T] only for a run that did not record it.
window_recorded = isfield(D.settings, 'lya_T_interval') && numel(D.settings.lya_T_interval) == 2;
if window_recorded; win = D.settings.lya_T_interval(:)'; else; win = [T/2, T]; end

% 1000 px tall: a figure taller than the screen (1073 px on the workstation,
% 2026-11-21) is shrunk by MATLAB at creation or at export.
fig = figure('Color', 'w', 'Position', [40 40 560 * n_cond + 120, 1000]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 6, n_cond, 'TileSpacing', 'compact', 'Padding', 'loose');
tl.TileIndexing = 'columnmajor';
rows = cell(1, n_cond);
tcol = SRNNCellTypePairs.type_colors(numel(R(1).cell_type_names));
ylabels = {{'External', 'input'}, {'Local rate', '(s^{-1})'}, {'Accumulating', '\lambda_k (s^{-1})'}, ...
    {'Time with local', 'rate > 0 (%)'}, {'Local rates > 0', sprintf('(of %d)', K)}, {'Expansion', '(bits/s)'}};
for i = 1:n_cond
    tk = R(i).topk; bn = R(i).ben;
    col = st.condition_color(R(i).name);
    t  = tk.t_lya(:); L = tk.local_LE_spectrum_t; F = tk.finite_LE_spectrum_t;
    assert(size(L, 2) == K, 'fig_local_lyapunov:K', ...
        'The run''s settings say K = %d but %s holds %d local rates.', K, R(i).name, size(L, 2));
    inwin = t >= win(1) & t <= win(2);
    bwin = bn.t_lya >= win(1) & bn.t_lya <= win(2);
    greys = 0.85 - 0.6 * (0:K-1)' / max(1, K - 1);
    axs = gobjects(1, 6);

    ax = nexttile(tl); axs(1) = ax; hold(ax, 'on');
    ti = R(i).type_indices;
    for q = numel(ti):-1:1
        plot(ax, R(i).t_ex, R(i).u_ex(ti{q}, :)', '-', 'Color', tcol(q, :), 'LineWidth', 0.4);
    end
    xlim(ax, [0 T]);
    title(ax, titles{i}, 'FontWeight', 'normal', 'FontSize', st.title_fs, 'Color', col);

    ax = nexttile(tl); axs(2) = ax; hold(ax, 'on');
    for k = K:-1:2; plot(ax, t, L(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.5); end
    h_one = plot(ax, t, L(:, 1), '-', 'Color', col, 'LineWidth', 1.2);
    h_ben = plot(ax, bn.t_lya, bn.local_lya, '-', 'Color', 'k', 'LineWidth', 0.8);
    plot(ax, [0 T], [0 0], '--', 'Color', st.zeroline_color, 'LineWidth', 1);
    mark_window_start(ax); xlim(ax, [0 T]);

    ax = nexttile(tl); axs(3) = ax; hold(ax, 'on');
    for k = K:-1:2; plot(ax, t, F(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.6); end
    plot(ax, t, F(:, 1), '-', 'Color', col, 'LineWidth', 1.6);
    plot(ax, bn.t_lya, bn.finite_lya, '--', 'Color', 'k', 'LineWidth', 1.2);
    plot(ax, win, [0 0], '--', 'Color', st.zeroline_color, 'LineWidth', 1);
    xlim(ax, win);
    % The first tenth of the window is the estimate settling: it is drawn, but
    % does not set the y range, which would hide the spread at the end.
    settled = t >= win(1) + 0.1 * diff(win) & t <= win(2);
    v = [reshape(F(settled, :), [], 1); reshape(bn.finite_lya(bn.t_lya >= win(1) + 0.1 * diff(win) & bn.t_lya <= win(2)), [], 1); 0];
    v = v(isfinite(v));
    ylim(ax, [min(v), max(v)] + 0.06 * max(max(v) - min(v), eps) * [-1 1]);

    share = mean(L(inwin, :) > 0, 1); share_ben = mean(bn.local_lya(bwin) > 0);
    ax = nexttile(tl); axs(4) = ax; hold(ax, 'on');
    bar(ax, 1:K, 100 * share, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    xlim(ax, [0.5 K + 0.5]); ylim(ax, [0 100]);
    xlabel(ax, 'Exponent index k', 'FontSize', st.tick_fs);

    n_pos_t = sum(L > 0, 2); any_pos = mean(n_pos_t(inwin) > 0);
    ax = nexttile(tl); axs(5) = ax; hold(ax, 'on');
    stairs(ax, t, n_pos_t, '-', 'Color', col, 'LineWidth', 0.8);
    mark_window_start(ax); xlim(ax, [0 T]); ylim(ax, [0 K]);

    h_loc = sum(max(L, 0), 2) / log(2); h_mean = mean(h_loc(inwin));
    ax = nexttile(tl); axs(6) = ax; hold(ax, 'on');
    plot(ax, t, h_loc, '-', 'Color', col, 'LineWidth', 0.9);
    mark_window_start(ax); xlim(ax, [0 T]);
    xlabel(ax, 'Time (s)', 'FontSize', st.tick_fs);

    set(axs, 'FontSize', st.tick_fs, 'Box', 'off', 'TickDir', 'out', 'LineWidth', st.axis_lw);
    if i == 1
        for q = 1:6; ylabel(axs(q), ylabels{q}, 'FontSize', st.tick_fs); end
        h_rest = plot(axs(2), NaN, NaN, '-', 'Color', 0.55 * [1 1 1], 'LineWidth', 0.8);
        lg = legend(axs(2), [h_one h_rest h_ben], {'k = 1 (in the condition''s color)', sprintf('k = 2 to %d (darker: larger k)', K), ...
            'Benettin (dashed where accumulating)'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', st.tick_fs);
        lg.Layout.Tile = 'south';
    end

    rows{i} = sprintf('| %s | %+.4f | %+.4f | %+.4f | %.0f%% | %.0f%% | %d | %.0f%% | %d | %.2f | %.2f |', titles{i}, ...
        tk.LE_spectrum(1), tk.LE_spectrum(end), bn.LLE, 100 * share(1), 100 * share_ben, ...
        round(median(n_pos_t(inwin))), 100 * any_pos, nnz(tk.LE_spectrum > 0), tk.h_KS_bits, h_mean);
end
if window_recorded
    win_note = 'the stage''s accumulation window (settings.lya_T_interval)';
else
    win_note = 'the second half of the run, because this run did not record its accumulation window';
end
hdr = {'# Local Lyapunov exponents under the switching input', '', ...
    sprintf(['Preset `%s`, seeds %s, T = %g s, K = %d tracked directions (top-K QR) and Benettin''s leading exponent. ' ...
    'Window of every exponent, share, median and mean below, and of the accumulating row: [%g, %g] s, %s.'], ...
    D.settings.preset_name, mat2str(D.settings.seeds), T, K, win(1), win(2), win_note), '', ...
    ['Rows of the figure: external input (every neuron, by cell type); the K local rates (color: k = 1; gray, darker ' ...
    'for larger k: the others; black: Benettin); the accumulating exponents (dashed black: Benettin); share of the ' ...
    'window each local rate is positive; number of local rates positive at each moment; sum of the positive local ' ...
    'rates over log 2 (bits/s; transient growth in the tracked directions, a lower bound when more than K expand; ' ...
    'not KS entropy). Green dashed: zero. Dotted vertical: start of the window, where it opens after t = 0.'], '', ...
    sprintf(['| Condition | top-K lambda_1 (1/s) | top-K lambda_%d (1/s) | Benettin lambda_1 (1/s) | leading local rate > 0 | ' ...
    'Benettin local > 0 | median number of positive local rates | at least one positive | positive exponents | ' ...
    'positive accumulated-rate sum (bits/s) | mean local positive-rate sum (bits/s) |'], K), ...
    '|---|---|---|---|---|---|---|---|---|---|---|'};
table_text = strjoin([hdr, rows], newline);
vprintf(cfg.verbose, 'verbose', '%s\n', table_text);

fig_tag = 'Fig_Local_Lyapunov';
out = struct('figs', fig, 'files', {{}}, 'source', data_file, 'window', win, 'K', K, 'table_text', table_text);
if cfg.save
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n', table_text);
        fclose(fid);
    end
    for ext = {'.png', '.fig'}
        f = fullfile(out_dir, [fig_tag ext{1}]); if isfile(f); delete(f); end
    end
    save_some_figs_to_folder_2(out_dir, fig_tag, fig.Number, {'png', 'fig'});
    num = num2str(fig.Number);
    movefile(fullfile(out_dir, [fig_tag '_figure_' num '.png']), fullfile(out_dir, [fig_tag '.png']), 'f');
    movefile(fullfile(out_dir, [fig_tag '_f_' num '.fig']),      fullfile(out_dir, [fig_tag '.fig']), 'f');
    out.files = existing_outputs(out_dir, fig_tag);
end
    function mark_window_start(a)
        if win(1) > 0
            xline(a, win(1), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        end
    end
end

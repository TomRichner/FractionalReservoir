function out = fig_local_lyapunov(cfg)
% FIG_LOCAL_LYAPUNOV Local Lyapunov exponents and local KS entropy under a stimulus staircase.
%
%   out = FIG_LOCAL_LYAPUNOV('run_dir', d)
%
% One column per condition of run_local_lyapunov's output, six rows: the
% external input (every neuron); the K local rates with Benettin's over the
% top; the K accumulating exponents with Benettin's finite-time curve; the
% share of the window each local exponent is positive; how many of the K are
% positive at each moment; and the local KS entropy rate sum_k max(local_k, 0)
% in bit/s. Table: final lambda_1 / lambda_K / Benettin, positive shares,
% number of positive exponents, accumulated and window-mean h_KS.
%
% PNG and .fig only, NO SVG: the figure holds ~15 million points at K = 100
% and the vector export took minutes (2026-09-14).
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
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
data_file = resolve_data_file(cfg.data_file, cfg.run_dir, {fullfile(cfg.run_dir, 'local_lyapunov')}, ...
    'local_lyapunov_data.mat', 'Run run_local_lyapunov first');
D = load(data_file);
R = D.results; n_cond = numel(R); T = D.settings.T; K = D.settings.K;

fig = figure('Color', 'w', 'Position', [40 40 620 * n_cond, 1500]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 6, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
rows = cell(1, n_cond);
tcol = SRNNCellTypePairs.type_colors(numel(R(1).cell_type_names));
for i = 1:n_cond
    tk = R(i).topk; bn = R(i).ben;
    col = st.condition_color(R(i).name);
    t  = tk.t_lya(:); L = tk.local_LE_spectrum_t; F = tk.finite_LE_spectrum_t;
    inwin = t >= T/2;
    greys = 0.85 - 0.6 * (0:K-1)' / max(1, K - 1);

    ax = nexttile(tl); hold(ax, 'on');
    ti = R(i).type_indices;
    for q = numel(ti):-1:1
        plot(ax, R(i).t_ex, R(i).u_ex(ti{q}, :)', '-', 'Color', tcol(q, :), 'LineWidth', 0.4);
    end
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    title(ax, R(i).title, 'FontWeight', 'normal', 'FontSize', st.title_fs);
    if i == 1; ylabel(ax, 'external input u', 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs);

    ax = nexttile(tl); hold(ax, 'on');
    for k = K:-1:2; plot(ax, t, L(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.5); end
    plot(ax, t, L(:, 1), '-', 'Color', col, 'LineWidth', 1.2);
    plot(ax, bn.t_lya, bn.local_lya, '-', 'Color', 'k', 'LineWidth', 0.8);
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]); xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    if i == 1; ylabel(ax, sprintf('local rate, top %d (s^{-1})', K), 'FontSize', st.label_fs); end
    text(ax, 0.02, 0.97, 'grey: top-K local rates (darker = larger k); colour: k = 1; black: Benettin', ...
        'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top');
    set(ax, 'FontSize', st.tick_fs);

    ax = nexttile(tl); hold(ax, 'on');
    for k = K:-1:2; plot(ax, t, F(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.6); end
    plot(ax, t, F(:, 1), '-', 'Color', col, 'LineWidth', 1.6);
    plot(ax, bn.t_lya, bn.finite_lya, '--', 'Color', 'k', 'LineWidth', 1.2);
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [T/2 T]);
    if i == 1; ylabel(ax, '\lambda_k(t), accumulating (s^{-1})', 'FontSize', st.label_fs); end
    text(ax, 0.98, 0.97, sprintf('\\lambda_1 = %+.3f, \\lambda_{%d} = %+.3f\nBenettin \\lambda_1 = %+.3f', ...
        tk.LE_spectrum(1), K, tk.LE_spectrum(end), bn.LLE), 'Units', 'normalized', ...
        'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    set(ax, 'FontSize', st.tick_fs);

    share = mean(L(inwin, :) > 0, 1); share_ben = mean(bn.local_lya(bn.t_lya >= T/2) > 0);
    ax = nexttile(tl); hold(ax, 'on');
    bar(ax, 1:K, 100 * share, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0.5 K + 0.5]); ylim(ax, [0 100]);
    xlabel(ax, 'exponent index k', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'time with local rate > 0 (%)', 'FontSize', st.label_fs); end
    title(ax, sprintf('leading local rate > 0: %.0f%% of the window', 100 * share(1)), 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'FontSize', st.tick_fs);

    n_pos_t = sum(L > 0, 2); any_pos = mean(n_pos_t(inwin) > 0);
    ax = nexttile(tl); hold(ax, 'on');
    stairs(ax, t, n_pos_t, '-', 'Color', col, 'LineWidth', 0.8);
    xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]); ylim(ax, [0 K]);
    if i == 1; ylabel(ax, sprintf('local rates > 0 (of %d)', K), 'FontSize', st.label_fs); end
    title(ax, sprintf('median %d positive; at least one positive %.0f%% of the window', ...
        round(median(n_pos_t(inwin))), 100 * any_pos), 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'FontSize', st.tick_fs);

    h_loc = sum(max(L, 0), 2) / log(2); h_mean = mean(h_loc(inwin));
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, t, h_loc, '-', 'Color', col, 'LineWidth', 0.9);
    yline(ax, h_mean, '--', sprintf('window mean %.1f', h_mean), 'Color', 'k', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    xlabel(ax, 'time (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'local h_{KS} (bit/s)', 'FontSize', st.label_fs); end
    title(ax, sprintf('local KS entropy rate; accumulated h_{KS} = %.2f bit/s', tk.h_KS_bits), 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'FontSize', st.tick_fs);

    rows{i} = sprintf('| %s | %+.4f | %+.4f | %+.4f | %.0f%% | %.0f%% | %d | %.2f | %.2f |', R(i).title, ...
        tk.LE_spectrum(1), tk.LE_spectrum(end), bn.LLE, 100 * share(1), 100 * share_ben, ...
        nnz(tk.LE_spectrum > 0), tk.h_KS_bits, h_mean);
end
title(tl, sprintf('%s, seeds %s, T = %g s, window [%g %g] s, top-%d QR + Benettin', ...
    strrep(D.settings.preset_name, '_', '\_'), mat2str(D.settings.seeds), T, T/2, T, K), ...
    'FontWeight', 'normal', 'FontSize', 10);
hdr = sprintf('| Condition | top-K lambda_1 | top-K lambda_%d | Benettin lambda_1 | leading local rate > 0 | Benettin local > 0 | n positive exponents | h_KS (bit/s) | window-mean local h_KS (bit/s) |', K);
vprintf(cfg.verbose, 'verbose', '%s\n|---|---|---|---|---|---|---|---|---|\n%s\n', hdr, strjoin(rows, newline));

fig_tag = 'Fig_Local_Lyapunov';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '%s\n|---|---|---|---|---|---|---|---|---|\n%s\n', hdr, strjoin(rows, newline));
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
end

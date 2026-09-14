% explore_local_lyapunov_exponents.m - how often are the local Lyapunov
% exponents of the paper's base model positive?
%
% One network (rng_seeds [1 2]) of the sfaEI_fast preset
% (celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25, n = 500,
% noise on, SRA1) under each of its three adaptation conditions, at the
% lyapunov_spectrum stage's MEDIUM settings: T = 40 s, fs 400, the exponents
% accumulated over the last 20 s after a 10-s alignment, lya_dt 0.05 s.
%
% Each condition is run TWICE on the same trajectory seed: once with the
% top-K QR method at K = 30 (no automatic retry, so it is exactly 30) and once
% with Benettin reshooting (K = 1, the manuscript's estimator), so the two can
% be drawn on top of each other. Noise increments are regenerable from the
% seed, so both runs see the same Brownian path.
%
% Figure, one column per condition:
%   row 1  the LOCAL rate of every one of the 30 exponents against time (grey,
%          darker = larger index; the leading one in the condition colour) with
%          Benettin's local rate over the top (black).
%   row 2  the ACCUMULATING finite-time exponents lambda_k(t) for the same 30,
%          with Benettin's finite-time lambda_1(t) over the top (black dashed),
%          zero line, final lambda_1 / lambda_30 and Benettin printed.
%   row 3  the share of the accumulation window in which each local exponent
%          is positive, k = 1..30 (bars), with Benettin's share as a line, and
%          the share of time the LEADING local rate is positive in the title.
%
% The point: "transient expansion" is a statement about the local rates; the
% finite-time exponents say what survives. In the multiple-timescale regime
% the leading local rate was positive ~3% of the time at medium in the sweeps
% (fig_local_vs_finite_lle); this shows the whole top of the spectrum.
%
% Output: figs/explorations/local_lyapunov_exponents/ (png, svg, fig) and a
% markdown table of the final exponents and positive shares. ~10-15 min.
% Assumes setup_paths has run.
%
% See also: fig_local_vs_finite_lle, fig_lyapunov_spectrum, lyapunov_topk,
%           SRNNCellTypePairs.lya_summary

setup_paths();

P        = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
seeds    = [1 2];
T        = 40;
K        = 30;
common   = {'rng_seeds', seeds, 'fs', 400, 'T_range', [0 T], ...
            'lya_T_interval', [T/2 T], 'lya_warmup', T/4, 'verbose', 'minimal'};
out_dir  = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', 'local_lyapunov_exponents');
if ~isfolder(out_dir); mkdir(out_dir); end

[~, ~, conditions] = srnn_param_preset(P);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
n_cond = numel(cond_names);
st = manuscript_style();

R = struct('name', cond_names, 'topk', [], 'ben', [], 'seconds', []);
for i = 1:n_cond
    t0 = tic;
    m = build_from_preset(P, cond_names{i}, common{:}, ...
        'lya_method', 'topk', 'lya_K', K, 'lya_K_auto', false, 'lya_dt', 0.05);
    m.run();
    R(i).topk = m.lya_results;
    m = build_from_preset(P, cond_names{i}, common{:}, 'lya_method', 'benettin');
    m.run();
    R(i).ben = m.lya_results;
    R(i).seconds = toc(t0);
    fprintf('%-14s top-%d lambda_1 %+.4f, lambda_%d %+.4f | Benettin %+.4f | %.0f s\n', ...
        cond_names{i}, K, R(i).topk.LE_spectrum(1), K, R(i).topk.LE_spectrum(end), ...
        R(i).ben.LLE, R(i).seconds);
end

%% Figure
fig = figure('Color', 'w', 'Position', [60 60 420 * n_cond, 950]);
tl = tiledlayout(fig, 3, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
rows = cell(1, n_cond);
for i = 1:n_cond
    tk = R(i).topk; bn = R(i).ben;
    col = st.condition_color(cond_names{i});
    t  = tk.t_lya(:);
    L  = tk.local_LE_spectrum_t;        % nt x K local rates
    F  = tk.finite_LE_spectrum_t;       % nt x K accumulating exponents (NaN before the window)
    inwin = t >= T/2;
    greys = 0.85 - 0.6 * (0:K-1)' / max(1, K - 1);   % index 1 darkest of the greys, drawn under

    % row 1: local rates
    ax = nexttile(tl); hold(ax, 'on');
    for k = K:-1:2
        plot(ax, t, L(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.5);
    end
    plot(ax, t, L(:, 1), '-', 'Color', col, 'LineWidth', 1.2);
    plot(ax, bn.t_lya, bn.local_lya, '-', 'Color', 'k', 'LineWidth', 0.8);
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
    xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    title(ax, st.condition_title(cond_names{i}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
    if i == 1; ylabel(ax, sprintf('local rate, top %d (s^{-1})', K), 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs);
    text(ax, 0.02, 0.97, 'grey: top-K local rates (darker = larger k); colour: k = 1; black: Benettin', ...
        'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top');

    % row 2: accumulating exponents
    ax = nexttile(tl); hold(ax, 'on');
    for k = K:-1:2
        plot(ax, t, F(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.6);
    end
    plot(ax, t, F(:, 1), '-', 'Color', col, 'LineWidth', 1.6);
    plot(ax, bn.t_lya, bn.finite_lya, '--', 'Color', 'k', 'LineWidth', 1.2);
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [T/2 T]);
    if i == 1; ylabel(ax, '\lambda_k(t), accumulating (s^{-1})', 'FontSize', st.label_fs); end
    xlabel(ax, 'time (s)', 'FontSize', st.label_fs);
    set(ax, 'FontSize', st.tick_fs);
    text(ax, 0.98, 0.97, sprintf('\\lambda_1 = %+.3f, \\lambda_{%d} = %+.3f\nBenettin \\lambda_1 = %+.3f', ...
        tk.LE_spectrum(1), K, tk.LE_spectrum(end), bn.LLE), 'Units', 'normalized', ...
        'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    % row 3: share of the window with a positive local rate, per exponent
    share = mean(L(inwin, :) > 0, 1);
    bwin  = bn.t_lya >= T/2;
    share_ben = mean(bn.local_lya(bwin) > 0);
    ax = nexttile(tl); hold(ax, 'on');
    bar(ax, 1:K, 100 * share, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    yline(ax, 100 * share_ben, '--', sprintf('Benettin %.0f%%', 100 * share_ben), 'Color', 'k', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0.5 K + 0.5]); ylim(ax, [0 100]);
    xlabel(ax, 'exponent index k', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'time with local rate > 0 (%)', 'FontSize', st.label_fs); end
    title(ax, sprintf('leading local rate > 0: %.0f%% of the window', 100 * share(1)), ...
        'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'FontSize', st.tick_fs);

    rows{i} = sprintf('| %s | %+.4f | %+.4f | %+.4f | %.0f%% | %.0f%% | %d |', ...
        st.condition_title(cond_names{i}), tk.LE_spectrum(1), tk.LE_spectrum(end), bn.LLE, ...
        100 * share(1), 100 * share_ben, nnz(tk.LE_spectrum > 0));
end
title(tl, sprintf('%s, seeds %s, T = %g s, window [%g %g] s, top-%d QR + Benettin', ...
    strrep(P, '_', '\_'), mat2str(seeds), T, T/2, T, K), 'FontWeight', 'normal', 'FontSize', 10);

save_figure_stable(out_dir, 'local_lyapunov_exponents', fig);
fid = fopen(fullfile(out_dir, 'local_lyapunov_exponents_table.md'), 'w');
fprintf(fid, '# Local Lyapunov exponents, %s, seeds %s, T = %g s, window [%g %g] s\n\n', P, mat2str(seeds), T, T/2, T);
fprintf(fid, '| Condition | top-K lambda_1 | top-K lambda_%d | Benettin lambda_1 | leading local rate > 0 | Benettin local > 0 | n positive exponents |\n|---|---|---|---|---|---|---|\n', K);
fprintf(fid, '%s\n', rows{:});
fclose(fid);
save(fullfile(out_dir, 'local_lyapunov_exponents_data.mat'), 'R', 'P', 'seeds', 'T', 'K', '-v7.3');
fprintf('saved to %s\n', out_dir);

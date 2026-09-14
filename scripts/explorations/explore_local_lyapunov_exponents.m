function out_dir = explore_local_lyapunov_exponents(preset_name, K, opts)
% EXPLORE_LOCAL_LYAPUNOV_EXPONENTS How often are the local Lyapunov exponents positive?
%
%   explore_local_lyapunov_exponents()                 % the paper preset (sfaEI_fast), K = 30
%   explore_local_lyapunov_exponents(preset_name, K)   % e.g. the tauSpread0p25 preset, top-80
%   explore_local_lyapunov_exponents(preset_name, K, 'T', 60, 'input_config', ic, 'tag', 'steps6')
%
% One network (rng_seeds [1 2]) of the preset (default
% celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25, n = 500,
% noise on, SRA1) under each of its three adaptation conditions, at the
% lyapunov_spectrum stage's MEDIUM settings scaled to T: fs 400, the exponents
% accumulated over the last T/2 after a T/4 alignment, lya_dt 0.05 s.
%
% Each condition is run TWICE on the same trajectory seed: once with the
% top-K QR method (no automatic retry, so it is exactly K) and once with
% Benettin reshooting (K = 1, the manuscript's estimator), so the two can be
% drawn on top of each other. Noise increments are regenerable from the seed,
% so both runs see the same Brownian path.
%
% Options:
%   'T'             simulation length in s (default 40); window [T/2 T], warmup T/4
%   'input_config'  a full input_config struct to REPLACE the preset's (default
%                   [] = the preset's own three-step pattern). Use it for a
%                   staircase of random steps: n_steps = T/10 with
%                   no_stim_pattern all false gives a new random step every
%                   10 s that never returns to zero.
%   'tag'           suffix for the output folder (default '')
%
% Figure, one column per condition, six rows:
%   row 1  the external input u(t), every neuron, E warm / I cool.
%   row 2  the LOCAL rate of every one of the K exponents against time (grey,
%          darker = larger index; the leading one in the condition colour) with
%          Benettin's local rate over the top (black).
%   row 3  the ACCUMULATING finite-time exponents lambda_k(t) for the same K,
%          with Benettin's finite-time lambda_1(t) over the top (black dashed),
%          zero line, final lambda_1 / lambda_K and Benettin printed.
%   row 4  the share of the accumulation window in which each local exponent
%          is positive, k = 1..K (bars), with Benettin's share as a line, and
%          the share of time the LEADING local rate is positive in the title.
%   row 5  how many of the K local rates are positive at each moment (a step
%          plot), with the share of the window at which at least one is.
%   row 6  the LOCAL Kolmogorov-Sinai entropy rate, sum_k max(local_k, 0) in
%          bit/s at each moment, with its window mean and the accumulated h_KS
%          of the top-K run in the title. Expect it to spike at every step of
%          the input.
% Plus, per condition, the class's own time-series summary from the Benettin
% run (timeseries_<condition>.png): u, x, r, synaptic output, SFA, STD and the
% local Lyapunov exponent, every neuron drawn.
%
% The point: "transient expansion" is a statement about the local rates; the
% finite-time exponents say what survives. In the multiple-timescale regime
% the leading local rate was positive ~3% of the time at medium in the sweeps
% (fig_local_vs_finite_lle); this shows the whole top of the spectrum.
%
% Output: figs/explorations/local_lyapunov_exponents/<preset>_K<K><tag>/ (png,
% svg, fig) and a markdown table of the final exponents and positive shares.
% ~4 min at K = 30, T = 40; ~20 min at K = 100, T = 60.
%
% See also: fig_local_vs_finite_lle, fig_lyapunov_spectrum, lyapunov_topk,
%           SRNNCellTypePairs.lya_summary, SRNNCellTypePairs.generate_external_input

arguments
    preset_name (1,:) char = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25'
    K           (1,1) double = 30
    opts.T            (1,1) double = 40
    opts.input_config              = []
    opts.tag          (1,:) char   = ''
end
setup_paths();

P        = preset_name;
seeds    = [1 2];
T        = opts.T;
common   = {'rng_seeds', seeds, 'fs', 400, 'T_range', [0 T], ...
            'lya_T_interval', [T/2 T], 'lya_warmup', T/4, 'verbose', 'minimal'};
if ~isempty(opts.input_config)
    common = [common, {'input_config', opts.input_config}];
end
out_dir  = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', ...
    'local_lyapunov_exponents', sprintf('%s_K%d%s', P, K, opts.tag));
if ~isfolder(out_dir); mkdir(out_dir); end

[~, ~, conditions] = srnn_param_preset(P);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
n_cond = numel(cond_names);
st = manuscript_style();

R = struct('name', cond_names, 'topk', [], 'ben', [], 'u_ex', [], 't_ex', [], ...
    'type_indices', [], 'cell_type_names', [], 'seconds', []);
for i = 1:n_cond
    t0 = tic;
    m = build_from_preset(P, cond_names{i}, common{:}, ...
        'lya_method', 'topk', 'lya_K', K, 'lya_K_auto', false, 'lya_dt', 0.05);
    m.run();
    R(i).topk = m.lya_results;
    R(i).u_ex = m.u_ex; R(i).t_ex = m.t_ex;
    R(i).type_indices = m.type_indices; R(i).cell_type_names = m.cell_type_names;
    m = build_from_preset(P, cond_names{i}, common{:}, 'lya_method', 'benettin');
    m.run();
    R(i).ben = m.lya_results;
    % The class's own summary figure from the Benettin run: u, x, r, synaptic
    % output, SFA, STD and the local Lyapunov exponent, every neuron drawn.
    [fh, ~] = m.plot();
    set(fh, 'Visible', 'off');
    save_figure_stable(out_dir, sprintf('timeseries_%s', cond_names{i}), fh);
    close(fh);
    R(i).seconds = toc(t0);
    fprintf('%-14s top-%d lambda_1 %+.4f, lambda_%d %+.4f | Benettin %+.4f | h_KS %.2f bit/s | %.0f s\n', ...
        cond_names{i}, K, R(i).topk.LE_spectrum(1), K, R(i).topk.LE_spectrum(end), ...
        R(i).ben.LLE, R(i).topk.h_KS_bits, R(i).seconds);
end

%% Figure
n_rows = 6;
fig = figure('Color', 'w', 'Position', [40 40 620 * n_cond, 1500]);
tl = tiledlayout(fig, n_rows, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
rows = cell(1, n_cond);
tcol = SRNNCellTypePairs.type_colors(numel(R(1).cell_type_names));
for i = 1:n_cond
    tk = R(i).topk; bn = R(i).ben;
    col = st.condition_color(cond_names{i});
    t  = tk.t_lya(:);
    L  = tk.local_LE_spectrum_t;        % nt x K local rates
    F  = tk.finite_LE_spectrum_t;       % nt x K accumulating exponents (NaN before the window)
    inwin = t >= T/2;
    greys = 0.85 - 0.6 * (0:K-1)' / max(1, K - 1);   % index 1 darkest of the greys, drawn under

    % row 1: the external input, every neuron, type colours, E on top
    ax = nexttile(tl); hold(ax, 'on');
    ti = R(i).type_indices; names = R(i).cell_type_names;
    for q = numel(names):-1:1
        plot(ax, R(i).t_ex, R(i).u_ex(ti{q}, :)', '-', 'Color', tcol(q, :), 'LineWidth', 0.4);
    end
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    title(ax, st.condition_title(cond_names{i}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
    if i == 1; ylabel(ax, 'external input u', 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs);

    % row 2: local rates
    ax = nexttile(tl); hold(ax, 'on');
    for k = K:-1:2
        plot(ax, t, L(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.5);
    end
    plot(ax, t, L(:, 1), '-', 'Color', col, 'LineWidth', 1.2);
    plot(ax, bn.t_lya, bn.local_lya, '-', 'Color', 'k', 'LineWidth', 0.8);
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
    xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    if i == 1; ylabel(ax, sprintf('local rate, top %d (s^{-1})', K), 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs);
    text(ax, 0.02, 0.97, 'grey: top-K local rates (darker = larger k); colour: k = 1; black: Benettin', ...
        'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top');

    % row 3: accumulating exponents
    ax = nexttile(tl); hold(ax, 'on');
    for k = K:-1:2
        plot(ax, t, F(:, k), '-', 'Color', greys(k) * [1 1 1], 'LineWidth', 0.6);
    end
    plot(ax, t, F(:, 1), '-', 'Color', col, 'LineWidth', 1.6);
    plot(ax, bn.t_lya, bn.finite_lya, '--', 'Color', 'k', 'LineWidth', 1.2);
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [T/2 T]);
    if i == 1; ylabel(ax, '\lambda_k(t), accumulating (s^{-1})', 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs);
    text(ax, 0.98, 0.97, sprintf('\\lambda_1 = %+.3f, \\lambda_{%d} = %+.3f\nBenettin \\lambda_1 = %+.3f', ...
        tk.LE_spectrum(1), K, tk.LE_spectrum(end), bn.LLE), 'Units', 'normalized', ...
        'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    % row 4: share of the window with a positive local rate, per exponent
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

    % row 5: how many local exponents are positive at each moment
    n_pos_t = sum(L > 0, 2);
    any_pos = mean(n_pos_t(inwin) > 0);
    ax = nexttile(tl); hold(ax, 'on');
    stairs(ax, t, n_pos_t, '-', 'Color', col, 'LineWidth', 0.8);
    xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]); ylim(ax, [0 K]);
    if i == 1; ylabel(ax, sprintf('local rates > 0 (of %d)', K), 'FontSize', st.label_fs); end
    title(ax, sprintf('median %d positive; at least one positive %.0f%% of the window', ...
        round(median(n_pos_t(inwin))), 100 * any_pos), 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'FontSize', st.tick_fs);

    % row 6: local KS entropy rate, sum of the positive local rates, bit/s
    h_loc = sum(max(L, 0), 2) / log(2);
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, t, h_loc, '-', 'Color', col, 'LineWidth', 0.9);
    yline(ax, mean(h_loc(inwin)), '--', sprintf('window mean %.1f', mean(h_loc(inwin))), 'Color', 'k', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    xline(ax, T/2, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off'); xlim(ax, [0 T]);
    xlabel(ax, 'time (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'local h_{KS} (bit/s)', 'FontSize', st.label_fs); end
    title(ax, sprintf('local KS entropy rate; accumulated h_{KS} = %.2f bit/s', tk.h_KS_bits), ...
        'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'FontSize', st.tick_fs);

    rows{i} = sprintf('| %s | %+.4f | %+.4f | %+.4f | %.0f%% | %.0f%% | %d | %.2f |', ...
        st.condition_title(cond_names{i}), tk.LE_spectrum(1), tk.LE_spectrum(end), bn.LLE, ...
        100 * share(1), 100 * share_ben, nnz(tk.LE_spectrum > 0), tk.h_KS_bits);
end
title(tl, sprintf('%s, seeds %s, T = %g s, window [%g %g] s, top-%d QR + Benettin%s', ...
    strrep(P, '_', '\_'), mat2str(seeds), T, T/2, T, K, ...
    tern(isempty(opts.input_config), '', ', input_config overridden')), ...
    'FontWeight', 'normal', 'FontSize', 10);

save_figure_stable(out_dir, 'local_lyapunov_exponents', fig);
fid = fopen(fullfile(out_dir, 'local_lyapunov_exponents_table.md'), 'w');
fprintf(fid, '# Local Lyapunov exponents, %s, seeds %s, T = %g s, window [%g %g] s, K = %d%s\n\n', ...
    P, mat2str(seeds), T, T/2, T, K, tern(isempty(opts.input_config), '', ' (input_config overridden)'));
fprintf(fid, '| Condition | top-K lambda_1 | top-K lambda_%d | Benettin lambda_1 | leading local rate > 0 | Benettin local > 0 | n positive exponents | h_KS (bit/s) |\n|---|---|---|---|---|---|---|---|\n', K);
fprintf(fid, '%s\n', rows{:});
fclose(fid);
input_config = opts.input_config; %#ok<NASGU>  saved with the data
save(fullfile(out_dir, 'local_lyapunov_exponents_data.mat'), 'R', 'P', 'seeds', 'T', 'K', 'input_config', '-v7.3');
fprintf('saved to %s\n', out_dir);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

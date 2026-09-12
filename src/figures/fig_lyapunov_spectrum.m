function out = fig_lyapunov_spectrum(cfg)
% FIG_LYAPUNOV_SPECTRUM The top-K Lyapunov spectrum per adaptation regime.
%
%   out = FIG_LYAPUNOV_SPECTRUM('run_dir', d)
%
% One panel per condition: the sorted exponents lambda_i against their index
% i, one line per network seed in the condition colour (preset noise on,
% solid; noise off, dashed), a zero line, and the x axis cut at the largest
% K used. Below the panels a text table gives, per condition and noise
% variant, the median [min, max] over seeds of n_+ (positive exponents),
% h_KS (bit/s) and D_KY (or "> K" when unresolved), plus the convergence
% drift of lambda_1. Data from run_lyapunov_spectrum.
%
% What to read off it: whether the regimes differ in the NUMBER of unstable
% directions and in dimension, not only in the sign of lambda_1; whether
% noise lowers h_KS and D_KY (Engelken et al. 2023); and, in the stable
% regime, the flat band of exponents near -1/tau_a,max that every neuron's
% slow adaptation contributes.
%
% See also: run_lyapunov_spectrum, lyapunov_topk, manuscript_style

arguments
    cfg.data_file   (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.x_log       (1,1) logical = false % log x axis (index) for wide spectra
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'lyapunov_spectrum')}, ...
    'lyapunov_spectrum_data.mat', ...
    'Run run_lyapunov_spectrum first');
D = load(data_file);
R = D.results; n_cond = numel(R);
variants = D.settings.variants;
styles = containers.Map({'noise_on', 'noise_off'}, {'-', '--'});

fig = figure('Color', 'w', 'Position', [80 80 380 * n_cond, 560]);
tl = tiledlayout(fig, 2, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
K_max_used = 0;
rows = {};
for i = 1:n_cond
    ax = nexttile(tl, [1 1]); hold(ax, 'on');
    col = st.condition_color(R(i).name);
    for v = 1:numel(variants)
        var = variants{v};
        if ~isfield(R(i), var) || isempty(R(i).(var)); continue; end
        runs = R(i).(var);
        for s = 1:numel(runs)
            lam = runs(s).LE_spectrum;
            plot(ax, 1:numel(lam), lam, styles(var), 'Color', [col, 0.75], 'LineWidth', 1.4);
            K_max_used = max(K_max_used, numel(lam));
        end
        rows{end + 1} = summary_row(R(i).title, var, runs); %#ok<AGROW>
    end
    yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax, 'off');
    title(ax, R(i).title, 'FontWeight', 'normal', 'FontSize', st.title_fs);
    xlabel(ax, 'index i', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, '\lambda_i (1/s)', 'FontSize', st.label_fs); end
    set(ax, 'FontSize', st.tick_fs); box(ax, 'off');
    if cfg.x_log; set(ax, 'XScale', 'log'); end
    ax_top(i) = ax; %#ok<AGROW>

    % Below: finite-time convergence of the first three exponents (noise on
    % if present), so the reader can judge whether the window sufficed.
    ax2 = nexttile(tl, [1 1]); hold(ax2, 'on');
    var = variants{1};
    if isfield(R(i), var) && ~isempty(R(i).(var))
        runs = R(i).(var);
        for s = 1:numel(runs)
            F = runs(s).finite_LE_spectrum_t; t = runs(s).t_lya;
            for j = 1:min(3, size(F, 2))
                plot(ax2, t, F(:, j), '-', 'Color', [col, 0.35 + 0.3 * (j == 1)], 'LineWidth', 1 + (j == 1));
            end
        end
    end
    yline(ax2, 0, ':', 'Color', [0.3 0.3 0.3]);
    hold(ax2, 'off');
    xlabel(ax2, 'time (s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax2, 'finite-time \lambda_{1..3}', 'FontSize', st.label_fs); end
    set(ax2, 'FontSize', st.tick_fs); box(ax2, 'off');
end
for i = 1:n_cond
    xlim(ax_top(i), [1, max(2, K_max_used)]);
end
linkaxes(ax_top, 'y');
if numel(variants) > 1
    legend(ax_top(1), {'noise on', 'noise off'}, 'Location', 'southwest', 'FontSize', 10);
end
title(tl, {sprintf('Top-K Lyapunov spectrum, %s, n = %d, T = %g s, %d seed(s)', ...
    strrep(D.settings.preset_name, '_', '\_'), D.settings.n, D.settings.T, D.settings.n_seeds), ...
    strjoin(rows, '   |   ')}, 'FontWeight', 'normal', 'FontSize', 10);

if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_Lyapunov_Spectrum';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    % The table as text too: the title line is cramped for many conditions.
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '| Condition | Noise | n+ | h_KS (bit/s) | D_KY | K used | lambda_1 | drift |\n|---|---|---|---|---|---|---|---|\n');
        for i = 1:n_cond
            for v = 1:numel(variants)
                var = variants{v};
                if ~isfield(R(i), var) || isempty(R(i).(var)); continue; end
                fprintf(fid, '%s\n', table_row(R(i).title, var, R(i).(var)));
            end
        end
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function s = summary_row(title, var, runs)
s = sprintf('%s (%s): n_+ %s, h_{KS} %s bit/s, D_{KY} %s', title, strrep(var, '_', ' '), ...
    mmm([runs.n_positive], '%d'), mmm([runs.h_KS_bits], '%.2f'), dky_mmm(runs));
end

function s = table_row(title, var, runs)
s = sprintf('| %s | %s | %s | %s | %s | %s | %s | %s |', title, strrep(var, '_', ' '), ...
    mmm([runs.n_positive], '%d'), mmm([runs.h_KS_bits], '%.2f'), dky_mmm(runs), ...
    mmm([runs.K_used], '%d'), mmm([runs.LLE], '%+.4f'), mmm(abs([runs.lambda_1_drift]), '%.3g'));
end

function s = mmm(v, fmt)
% median [min, max] over seeds, collapsed when they coincide.
if numel(v) == 1 || max(v) == min(v)
    s = sprintf(fmt, median(v));
else
    s = sprintf([fmt ' [' fmt ', ' fmt ']'], median(v), min(v), max(v));
end
end

function s = dky_mmm(runs)
res = [runs.D_KY_resolved] == 1;
if all(res)
    s = mmm([runs.D_KY], '%.1f');
elseif ~any(res)
    s = sprintf('> %d', min([runs.K_used]));
else
    s = sprintf('%s (%d of %d unresolved)', mmm([runs(res).D_KY], '%.1f'), nnz(~res), numel(runs));
end
end

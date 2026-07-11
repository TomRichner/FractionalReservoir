function visualize_campagnola_matrices(save_figs)
%VISUALIZE_CAMPAGNOLA_MATRICES  Labeled heatmaps of the Campagnola 2022 cell-type matrices.
%
%   VISUALIZE_CAMPAGNOLA_MATRICES() loads campagnola_matrices.mat (co-located; falls
%   back to the CSVs via load_campagnola_matrices) and draws a set of heatmap figures,
%   one thematic figure per group:
%       1. Connectivity        (raw + distance-adjusted probability, coverage)
%       2. Strength & kinetics (PSP amplitudes, latency, rise/decay)
%       3. Short-term plasticity (induction, recovery, variability)
%       4. Stochastic-release model (facilitation/depression, Pr, sites, mini amp)
%       5. Spike-frequency adaptation (per-type bar chart)
%
%   Every matrix is pre (rows) x post (cols), ordered {Pyr,Pvalb,Sst,Vip}. Signed
%   metrics (PSP amplitude, STP) use a diverging blue-white-red map centered at 0;
%   magnitudes use parula. Under-sampled elements (n<2) are gray "n/a".
%
%   VISUALIZE_CAMPAGNOLA_MATRICES(true) also saves PNGs into a figures/ subdirectory.
%
%   See src/connectivity/load_campagnola_matrices.m and campagnola/PROVENANCE.md.

    if nargin < 1 || isempty(save_figs), save_figs = false; end

    this_dir = fileparts(mfilename('fullpath'));
    mat_file = fullfile(this_dir, 'campagnola_matrices.mat');
    if isfile(mat_file)
        S = load(mat_file);
        C = S.campagnola;
        C.types = reshape(cellstr(C.types), 1, []);
    else
        addpath(fileparts(this_dir));            % src/connectivity, for the CSV loader
        C = load_campagnola_matrices();
    end

    % Pretty type labels: pyr -> Pyr, etc.
    lab = cellfun(@(t)[upper(t(1)) t(2:end)], C.types, 'UniformOutput', false);

    figs = gobjects(0);

    % --- Figure 1: connectivity ------------------------------------------
    panels1 = { ...
        'conn_prob',      'Connection prob (raw)',        'seq', 1,   '%.2f'; ...
        'conn_prob_adj',  'Connection prob (adjusted)',   'seq', 1,   '%.2f'; ...
        'n_connected',    '# connections found',          'seq', 1,   '%d'};
    figs(end+1) = draw_figure(C, lab, panels1, 'Campagnola 2022 — mouse connectivity', [1 3]);

    % --- Figure 2: strength & kinetics -----------------------------------
    panels2 = { ...
        'psp_amplitude',   'Resting PSP amp (mV)',   'div', 1e3, '%.2f'; ...
        'pulse_amp_90pct', '90th-pct PSP amp (mV)',  'div', 1e3, '%.2f'; ...
        'latency',         'Latency (ms)',           'seq', 1e3, '%.2f'; ...
        'psc_rise_time',   'PSC rise time (ms)',     'seq', 1e3, '%.2f'; ...
        'psc_decay_tau',   'PSC decay tau (ms)',     'seq', 1e3, '%.1f'};
    figs(end+1) = draw_figure(C, lab, panels2, 'Campagnola 2022 — synaptic strength & kinetics', [2 3]);

    % --- Figure 3: short-term plasticity ---------------------------------
    panels3 = { ...
        'stp_induction_50hz',  'STP induction (+facil / -depress)', 'div', 1, '%.2f'; ...
        'stp_recovery_250ms',  'STP recovery @250ms',               'div', 1, '%.2f'; ...
        'variability_resting', 'Resting variability (aCV)',         'seq', 1, '%.2f'};
    figs(end+1) = draw_figure(C, lab, panels3, 'Campagnola 2022 — short-term plasticity', [1 3]);

    % --- Figure 4: stochastic-release model ------------------------------
    panels4 = { ...
        'ml_facilitation_amount', 'Facilitation amount',   'seq', 1,   '%.2f'; ...
        'ml_facilitation_tau',    'Facilitation tau (s)',  'seq', 1,   '%.2f'; ...
        'ml_depression_amount',   'Depression amount',     'div', 1,   '%.2f'; ...
        'ml_depression_tau',      'Depression tau (s)',    'seq', 1,   '%.2f'; ...
        'ml_release_prob',        'Base release prob',     'seq', 1,   '%.2f'; ...
        'ml_n_release_sites',     '# release sites',       'seq', 1,   '%.1f'; ...
        'ml_mini_amplitude',      'Mini amplitude (mV)',   'div', 1e3, '%.2f'};
    figs(end+1) = draw_figure(C, lab, panels4, 'Campagnola 2022 — stochastic-release model (STF/STD source)', [2 4]);

    % --- Figure 5: SFA bar chart -----------------------------------------
    f5 = figure('Name', 'SFA', 'Color', 'w', 'Position', [200 200 560 420]);
    ax = axes(f5);
    b = bar(ax, categorical(lab, lab), C.sfa_adaptation_index, 0.6, 'FaceColor', [0.30 0.45 0.70]);
    ylabel(ax, 'adaptation index (median)');
    title(ax, 'Spike-frequency adaptation per type');
    grid(ax, 'on'); ax.YGrid = 'on'; ax.XGrid = 'off';
    % annotate with sample counts
    xt = b.XEndPoints; yt = b.YEndPoints;
    for i = 1:numel(xt)
        text(ax, xt(i), yt(i), sprintf('  n=%d', round(C.sfa_adaptation_index_n(i))), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    ylim(ax, [0, max(C.sfa_adaptation_index) * 1.25]);
    figs(end+1) = f5;

    % --- optional save ---------------------------------------------------
    if save_figs
        fig_dir = fullfile(this_dir, 'figures');
        if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
        names = {'connectivity', 'strength_kinetics', 'stp', 'release_model', 'sfa'};
        for k = 1:numel(figs)
            exportgraphics(figs(k), fullfile(fig_dir, ['campagnola_' names{k} '.png']), 'Resolution', 150);
        end
        fprintf('Saved %d figures to %s\n', numel(figs), fig_dir);
    end
end

% ======================================================================= %
function f = draw_figure(C, lab, panels, ttl, layout)
    f = figure('Name', ttl, 'Color', 'w', ...
        'Position', [100 100 360*layout(2), 340*layout(1)]);
    tl = tiledlayout(f, layout(1), layout(2), 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:size(panels, 1)
        field = panels{p, 1};
        if ~isfield(C, field)
            warning('visualize_campagnola:missingField', 'skipping missing field %s', field);
            continue;
        end
        ax = nexttile(tl);
        heatmap_panel(ax, C.(field), lab, panels{p, 2}, panels{p, 3}, panels{p, 4}, panels{p, 5});
    end
    sgtitle(f, ttl, 'FontWeight', 'bold', 'Interpreter', 'none');
end

% ----------------------------------------------------------------------- %
function heatmap_panel(ax, M, lab, ttl, ctype, scale, fmt)
    n = size(M, 1);
    Ms = M * scale;
    valid = ~isnan(Ms);

    imagesc(ax, Ms, 'AlphaData', valid);
    ax.Color = [0.9 0.9 0.9];                 % background shows through NaN cells
    axis(ax, 'equal'); axis(ax, 'tight');
    ax.XTick = 1:n; ax.YTick = 1:n;
    ax.XTickLabel = lab; ax.YTickLabel = lab;
    ax.XAxisLocation = 'top';
    ax.TickLength = [0 0];
    xlabel(ax, 'postsynaptic'); ylabel(ax, 'presynaptic');
    title(ax, ttl, 'Interpreter', 'none', 'FontSize', 10);

    if strcmp(ctype, 'div')
        cmap = diverging_bwr(256);
        a = max(abs(Ms(valid)));
        if isempty(a) || a == 0, a = 1; end
        lim = [-a, a];
    else
        cmap = parula(256);
        lo = min(Ms(valid)); hi = max(Ms(valid));
        if isempty(lo), lo = 0; hi = 1; end
        if lo == hi, hi = lo + eps(lo) + 1; end
        lim = [lo, hi];
    end
    colormap(ax, cmap);
    clim(ax, lim);
    cb = colorbar(ax); cb.FontSize = 8;

    % annotate cells with values (luminance-aware text color)
    for i = 1:n
        for j = 1:n
            v = Ms(i, j);
            if isnan(v)
                text(ax, j, i, 'n/a', 'HorizontalAlignment', 'center', ...
                    'Color', [0.45 0.45 0.45], 'FontSize', 8);
            else
                tc = text_color(v, lim, cmap);
                text(ax, j, i, sprintf(fmt, v), 'HorizontalAlignment', 'center', ...
                    'Color', tc, 'FontSize', 8, 'FontWeight', 'bold');
            end
        end
    end
end

% ----------------------------------------------------------------------- %
function tc = text_color(v, lim, cmap)
    % pick black/white text by the luminance of the cell's mapped color
    t = (v - lim(1)) / (lim(2) - lim(1));
    t = min(max(t, 0), 1);
    idx = round(t * (size(cmap, 1) - 1)) + 1;
    rgb = cmap(idx, :);
    lum = 0.299 * rgb(1) + 0.587 * rgb(2) + 0.114 * rgb(3);
    if lum < 0.5, tc = [1 1 1]; else, tc = [0 0 0]; end
end

% ----------------------------------------------------------------------- %
function c = diverging_bwr(n)
    % colorblind-aware blue -> white -> red diverging colormap
    h = floor(n / 2);
    blue = [0.23 0.30 0.75];
    red  = [0.79 0.10 0.15];
    b2w = [linspace(blue(1), 1, h)', linspace(blue(2), 1, h)', linspace(blue(3), 1, h)'];
    w2r = [linspace(1, red(1), n - h)', linspace(1, red(2), n - h)', linspace(1, red(3), n - h)'];
    c = [b2w; w2r];
end

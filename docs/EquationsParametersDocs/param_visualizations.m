function param_visualizations(save_figs)
%PARAM_VISUALIZATIONS  Heatmaps of the Campagnola 2022 -> SRNN parameter mapping.
%
%   PARAM_VISUALIZATIONS() draws (and by default saves) the figures used by
%   docs/EquationsParametersDocs/SFA_STD_STF_parameter_mapping.md:
%
%     1. Campagnola 2022 source data actually used (connectivity, PSP amplitude,
%        release-model STF/STD metrics) as pre x post heatmaps.
%     2. The derived SRNNModelCellTypes synaptic parameters (block weight,
%        release prob p0, facilitation rate kappa and tau_f, STD tau_rec/tau_rel).
%     3. Spike-frequency adaptation per cell type (adaptation index and the
%        derived SFA strength c).
%
%   Rows = PREsynaptic type, cols = POSTsynaptic type, order {Pyr,Pvalb,Sst,Vip}.
%   Signed quantities (PSP amplitude, depression amount, block weight) use a
%   diverging blue-white-red map centered at 0; magnitudes use parula.
%
%   The RAW source matrices come from load_campagnola_matrices(); the DERIVED
%   model parameters are read/recomputed from a built SRNNModelCellTypes using
%   the exact expressions in its load_parameter_tables()/get_params() (constants
%   c_gain, tau_b_rel_ref are read from the object so nothing drifts).
%
%   PARAM_VISUALIZATIONS(false) draws without saving.

    if nargin < 1 || isempty(save_figs), save_figs = true; end

    this_dir     = fileparts(mfilename('fullpath'));         % docs/EquationsParametersDocs
    project_root = fileparts(fileparts(this_dir));
    addpath(genpath(fullfile(project_root, 'src')));         % loader + model classes

    % --- raw Campagnola source (rows=pre, cols=post) ---------------------
    C   = load_campagnola_matrices();
    lab = cellfun(@(t)[upper(t(1)) t(2:end)], C.types, 'UniformOutput', false);

    % --- derived model parameters (authoritative: from a built model) ----
    m = SRNNModelCellTypes('n', 8);
    m.build();                                               % fills the K x K tables
    tau_b_rec = m.dep_tau;                                              % = ml_depression_tau
    tau_b_rel = max(m.tau_b_rel_ref * (1 - min(max(m.dep_amount,0),0.95)), 0.05);
    p0        = min(max(m.rel_prob, 0.05), 0.95);                       % clamped release prob
    tau_f     = m.fac_tau;                                              % = ml_facilitation_tau
    kappa     = max(m.kappa, 0);                                        % = ml_facilitation_amount
    c_type    = m.c_gain * m.adapt_index(:);                            % SFA strength per type
    w_mean    = m.conn_prob .* m.psp_amp;                              % mean block weight seed (V)

    figs = gobjects(0);

    % --- Figure 1: Campagnola inputs used --------------------------------
    p1 = { ...
        C.conn_prob_adj,          'Connection prob (adjusted)', 'seq', 1,   '%.2f'; ...
        C.psp_amplitude,          'PSP amplitude (mV)',         'div', 1e3, '%.2f'; ...
        C.ml_release_prob,        'Release prob (raw)',         'seq', 1,   '%.2f'; ...
        C.sfa_adaptation_index,   'SFA index (per type)',       'col', 1,   '%.3f'; ...
        C.ml_facilitation_amount, 'Facilitation amount',        'seq', 1,   '%.2f'; ...
        C.ml_facilitation_tau,    'Facilitation tau (s)',       'seq', 1,   '%.3f'; ...
        C.ml_depression_amount,   'Depression amount (signed)', 'div', 1,   '%.2f'; ...
        C.ml_depression_tau,      'Depression tau (s)',         'seq', 1,   '%.3f'};
    figs(end+1) = draw_matrices(p1, lab, ...
        'Campagnola 2022 source data used by the model', [2 4]); %#ok<AGROW>

    % --- Figure 2: derived model synaptic parameters ---------------------
    p2 = { ...
        w_mean,    'Mean block weight w (mV)',   'div', 1e3, '%.2f'; ...
        p0,        'Release prob p_0 (model)',   'seq', 1,   '%.2f'; ...
        kappa,     'Facilitation rate \kappa',   'seq', 1,   '%.2f'; ...
        tau_f,     'Facilitation \tau_f (s)',    'seq', 1,   '%.3f'; ...
        tau_b_rec, 'STD recovery \tau_{rec} (s)','seq', 1,   '%.3f'; ...
        tau_b_rel, 'STD release \tau_{rel} (s)', 'seq', 1,   '%.3f'};
    figs(end+1) = draw_matrices(p2, lab, ...
        'Derived SRNNModelCellTypes synaptic parameters (per pre \times post)', [2 3]); %#ok<AGROW>

    % --- Figure 3: SFA per type ------------------------------------------
    f3 = figure('Name', 'SFA per type', 'Color', 'w', 'Position', [200 200 820 380]);
    tl = tiledlayout(f3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(tl);
    bar(ax, categorical(lab, lab), C.sfa_adaptation_index, 0.6, 'FaceColor', [0.30 0.45 0.70]);
    ylabel(ax, 'adaptation index (median)'); title(ax, 'Campagnola SFA index');
    grid(ax, 'on'); ax.XGrid = 'off';
    for i = 1:numel(lab)
        text(ax, i, C.sfa_adaptation_index(i), sprintf('  n=%d', round(C.sfa_adaptation_index_n(i))), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    ax2 = nexttile(tl);
    bar(ax2, categorical(lab, lab), c_type, 0.6, 'FaceColor', [0.70 0.45 0.30]);
    ylabel(ax2, 'SFA strength c = c_{gain}\cdot index'); title(ax2, sprintf('Model SFA strength (c_{gain}=%.2g)', m.c_gain));
    grid(ax2, 'on'); ax2.XGrid = 'off';
    sgtitle(f3, 'Spike-frequency adaptation per cell type', 'FontWeight', 'bold');
    figs(end+1) = f3; %#ok<AGROW>

    % --- save ------------------------------------------------------------
    if save_figs
        fig_dir = fullfile(this_dir, 'figures');
        if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
        names = {'param_map_campagnola_inputs', 'param_map_model_std_stf', 'param_map_sfa'};
        for k = 1:numel(figs)
            exportgraphics(figs(k), fullfile(fig_dir, [names{k} '.png']), 'Resolution', 150);
        end
        fprintf('Saved %d figures to %s\n', numel(figs), fig_dir);
    end
end

% ======================================================================= %
function f = draw_matrices(panels, lab, ttl, layout)
% DRAW_MATRICES  Tiled heatmaps. panels rows: {M, title, ctype, scale, fmt},
% ctype in {'seq','div','col'} ('col' = per-type 4x1 shown as a column).
    f = figure('Name', ttl, 'Color', 'w', 'Position', [80 80 340*layout(2), 340*layout(1)]);
    tl = tiledlayout(f, layout(1), layout(2), 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:size(panels, 1)
        ax = nexttile(tl);
        heatmap_panel(ax, panels{p,1}, lab, panels{p,2}, panels{p,3}, panels{p,4}, panels{p,5});
    end
    sgtitle(f, ttl, 'FontWeight', 'bold');
end

% ----------------------------------------------------------------------- %
function heatmap_panel(ax, M, lab, ttl, ctype, scale, fmt)
    M = M(:, :);                                   % 4x4, or 4x1 for 'col'
    is_col = strcmp(ctype, 'col') || size(M, 2) == 1;
    Ms = M * scale;
    valid = ~isnan(Ms);
    n = size(Ms, 1); ncol = size(Ms, 2);

    imagesc(ax, Ms, 'AlphaData', valid);
    ax.Color = [0.9 0.9 0.9];
    axis(ax, 'equal'); axis(ax, 'tight');
    ax.YTick = 1:n; ax.YTickLabel = lab;
    if is_col
        ax.XTick = 1; ax.XTickLabel = {''};
        ylabel(ax, 'cell type');
    else
        ax.XTick = 1:ncol; ax.XTickLabel = lab; ax.XAxisLocation = 'top';
        xlabel(ax, 'postsynaptic'); ylabel(ax, 'presynaptic');
    end
    ax.TickLength = [0 0];
    title(ax, ttl, 'FontSize', 10);

    if strcmp(ctype, 'div')
        cmap = diverging_bwr(256);
        a = max(abs(Ms(valid))); if isempty(a) || a == 0, a = 1; end
        lim = [-a, a];
    else
        cmap = parula(256);
        lo = min(Ms(valid)); hi = max(Ms(valid));
        if isempty(lo), lo = 0; hi = 1; end
        if lo == hi, hi = lo + eps(lo) + 1; end
        lim = [lo, hi];
    end
    colormap(ax, cmap); clim(ax, lim);
    cb = colorbar(ax); cb.FontSize = 8;

    for i = 1:n
        for j = 1:ncol
            v = Ms(i, j);
            if isnan(v)
                text(ax, j, i, 'n/a', 'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.45], 'FontSize', 8);
            else
                text(ax, j, i, sprintf(fmt, v), 'HorizontalAlignment', 'center', ...
                    'Color', text_color(v, lim, cmap), 'FontSize', 8, 'FontWeight', 'bold');
            end
        end
    end
end

% ----------------------------------------------------------------------- %
function tc = text_color(v, lim, cmap)
    t = min(max((v - lim(1)) / (lim(2) - lim(1)), 0), 1);
    rgb = cmap(round(t * (size(cmap,1) - 1)) + 1, :);
    if 0.299*rgb(1) + 0.587*rgb(2) + 0.114*rgb(3) < 0.5, tc = [1 1 1]; else, tc = [0 0 0]; end
end

% ----------------------------------------------------------------------- %
function c = diverging_bwr(n)
    h = floor(n/2); blue = [0.23 0.30 0.75]; red = [0.79 0.10 0.15];
    b2w = [linspace(blue(1),1,h)', linspace(blue(2),1,h)', linspace(blue(3),1,h)'];
    w2r = [linspace(1,red(1),n-h)', linspace(1,red(2),n-h)', linspace(1,red(3),n-h)'];
    c = [b2w; w2r];
end

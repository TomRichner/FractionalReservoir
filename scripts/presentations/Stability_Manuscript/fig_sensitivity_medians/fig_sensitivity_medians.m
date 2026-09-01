function out = fig_sensitivity_medians(cfg)
% FIG_SENSITIVITY_MEDIANS Sensitivity medians, conditions overlaid, 2 x 3.
%
%   out = FIG_SENSITIVITY_MEDIANS()
%   out = FIG_SENSITIVITY_MEDIANS('run_dir', d)
%
% The compact counterpart of fig_sensitivity_analysis_allStd. That figure tiles
% the four adaptation conditions as COLUMNS and shows the full rep distribution
% per panel (28 panels over two sheets per metric). Here the conditions are
% OVERLAID in one axes per swept parameter, as median lines only, so all six
% parameters fit one 2 x 3 sheet and the conditions can be compared directly
% rather than across columns. No simulation is re-run.
%
% THREE CONSEQUENCES of that choice:
%   1. NO DISTRIBUTIONS -- only the median across reps is drawn. The percentile
%      machinery is written generally (see pcts / band_pcts) so adding an IQR
%      band later is a two-line change.
%   2. SIX PARAMETERS, not seven. level_of_chaos ("Synaptic Gain") is dropped so
%      the rest fit 2 x 3; it is the least surprising of the seven, since it
%      simply scales W.
%   3. NO REPLOT DETOUR. replot_sensitivity / plot_sensitivity /
%      assemble_sensitivity_figure exist to build the imagesc sheets; harvesting
%      median lines back out of their saved .fig files would be fragile. The
%      medians are computed straight off the saved PSA objects, from the same
%      ParamSpaceAnalysis2.collect_level_values that plot_sensitivity uses.
%
% See also: fig_sensitivity_analysis_allStd, resolve_run_dir,
%           preset_default_values, manuscript_style

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% (no 'close all force' -- destroyed sibling figures in a batch; see header)

%% -------------------- What goes in which panel --------------------
% Panel order IS the tile order (row-major over 2 x 3), and is what the panel
% letters and the styling lookups are indexed by. Names are the raw property
% names, matching psa.grid_params.
%
% level_of_chaos is deliberately absent: seven parameters do not tile 2 x 3, and
% the synaptic-gain sweep is the least surprising of the seven (it simply scales
% W, so every condition rises monotonically through the edge of chaos).
panel_params = {'f_E', 'n', 'mu_EE_relative', ...
                'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'};
n_rows = 2;
n_cols = 3;

% Per-parameter x-axis styling. Empty ticks means "leave MATLAB's alone".
%
% f_E: the sweep varies the fraction excitatory, but because mu_EI and mu_IE are
% held fixed it is the E:I NEURON COUNT that the axis is really reporting. This
% run has n = 500, so f_E = 0.2/0.5/0.8 is 100:400 / 250:250 / 400:100 neurons --
% spelled out in counts rather than as the reduced ratios 1:4 / 1:1 / 4:1, which
% would hide the network size the counts are drawn from.
% mu_*: RMT block means, indexed (post, pre) -- mu_EE is E<-E. Rendered in tex.
% The four mu axes are re-expressed as PERCENT OF THE PRESET DEFAULT below, so
% their raw ticks are left empty here.
row_style = containers.Map( ...
    {'f_E', 'n', 'mu_EE_relative', 'mu_EI_relative', ...
     'mu_IE_relative', 'mu_II_relative'}, ...
    {struct('xlabel', 'E:I neuron ratio', 'xticks', [0.2 0.5 0.8], ...
                                          'xticklabels', {{'100:400','250:250','400:100'}}), ...
     struct('xlabel', 'Network Size',     'xticks', [100 500 1000], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{EE}',         'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{EI}',         'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{IE}',         'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{II}',         'xticks', [], 'xticklabels', {{}})});

% Axes shown as percent departure from the preset default, (value/default - 1)*100.
% Absolute mu_tilde values mean little on their own; what the sweep varies is the
% departure from the network the preset defines. See apply_percent_axis for the
% XDir reversal that keeps the two negative (inhibitory) blocks reading
% left-to-right like the other two.
pct_params = {'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'};

% Preset defaults for every panel, not just the percent ones: they also place the
% default marker, and on the f_E and Network Size panels there is no 0% tick to
% carry that information, which is exactly where the marker earns its keep.
default_value = preset_default_values(run_dir, panel_params);

%% -------------------- Conditions: colour, order, labels --------------------
% Okabe-Ito colorblind-safe palette, copied from
% ../fig_memory_capacity_example/Fig_memory_capacity_example.m so the two
% presentation figures name the same condition with the same colour. Reddish
% purple (instead of bluish green) keeps all four hues well separated.
%
% Matched BY NAME, not by position in psa.conditions, so a run whose conditions
% were declared in a different order cannot silently recolour the figure. The
% display names match the condition_titles map in
% ParamSpaceAnalysis2.plot_sensitivity, so the legend here and the column titles
% of the allStd sheets read identically.
% Driven by the PRESET's conditions, in its order, rather than a hardcoded four.
% A 7-regime preset would otherwise have plotted only the four this list knew
% about, with no indication that three were missing.
[~, ~, preset_conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, preset_conditions, 'UniformOutput', false);
cond_spec  = struct( ...
    'name',    cond_names, ...
    'display', cellfun(@(n) st.condition_title(n), cond_names, 'UniformOutput', false), ...
    'color',   cellfun(@(n) st.condition_color(n), cond_names, 'UniformOutput', false));

%% -------------------- Percentiles --------------------
% pcts(median_col) is the curve that gets drawn. Everything downstream indexes
% the percentile dimension, so plotting the 25th/75th later means adding them to
% pcts and setting band_pcts to their column indices -- the shaded-band branch in
% the plotting loop is already written, just inactive while band_pcts is empty.
pcts       = 50;
median_col = find(pcts == 50, 1);
assert(~isempty(median_col), 'pcts must include 50 (the median is the plotted curve).');
band_pcts  = [];   % e.g. [find(pcts==25) find(pcts==75)] to shade the IQR

%% -------------------- Styling --------------------
tick_fs   = st.tick_fs;
label_fs  = st.label_fs;
legend_fs = 14;
letter_fs = 18;    % panel letters
line_lw   = 2.5;   % opaque: four overlaid curves, so the allStd alpha would muddy them
band_alpha = 0.15; % only used when band_pcts is non-empty

zeroline_lw = 2;                        % green dashed lambda_1 = 0 line
zeroline_color = [0 0.7 0];

% Short reddish-gray tick rising off the x-axis at the preset's default -- the
% "you are here" mark for the network the sweep departs from. Deliberately
% low-contrast: it is a reference, not data.
default_mark_color = [0.80 0.60 0.60];
default_mark_lw    = 2;
default_mark_frac  = 0.05;

%% -------------------- Load the sweeps --------------------
% One 1D_sensitivity_* subfolder per swept parameter. Which parameter a folder
% swept is read off the PSA itself (the single non-reps grid param), exactly as
% replot_sensitivity does -- never guessed from the folder name.
sens_listing = dir(fullfile(run_dir, '1D_sensitivity_*'));
sens_listing = sens_listing([sens_listing.isdir]);
if isempty(sens_listing)
    error('fig_sensitivity_medians:NotFound', ...
        'No 1D_sensitivity_* subfolder found in:\n  %s', run_dir);
end

psa_of_param = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(sens_listing)
    src_dir = fullfile(sens_listing(k).folder, sens_listing(k).name);
    if ~isfile(fullfile(src_dir, 'psa_object.mat'))
        warning('Fig_sensitivity_medians:NoPsaFile', ...
            'Skipping (no psa_object.mat): %s', src_dir);
        continue;
    end
    psa_k = ParamSpaceAnalysis2.from_dir(src_dir);
    swept = setdiff(psa_k.grid_params, {'reps'}, 'stable');
    if isempty(swept)
        warning('Fig_sensitivity_medians:NoSweptParam', ...
            'Skipping (no swept param found): %s', src_dir);
        continue;
    end
    if ~ismember(swept{1}, panel_params)
        continue;   % e.g. level_of_chaos, which this figure omits
    end
    psa_of_param(swept{1}) = psa_k;
    fprintf('Loaded sweep %-16s from %s\n', swept{1}, sens_listing(k).name);
end

missing = panel_params(~cellfun(@(p) isKey(psa_of_param, p), panel_params));
if ~isempty(missing)
    error('Fig_sensitivity_medians:MissingSweep', ...
        ['No 1D_sensitivity_* sweep found for: %s\n' ...
         'Source run: %s'], strjoin(missing, ', '), run_dir);
end

%% -------------------- Collect the percentile curves --------------------
% curves.(param).x                     : 1 x n_levels swept-parameter values
% curves.(param).y.(cond)              : n_levels x numel(pcts)
% Computed once and reused by both metrics' figures.
metric_names = {'LLE', 'mean_rate'};
curves = struct();
for mi = 1:numel(metric_names)
    metric = metric_names{mi};
    for pi = 1:numel(panel_params)
        param = panel_params{pi};
        psa_p = psa_of_param(param);
        param_dim = find(strcmp(psa_p.grid_params, param));
        x_values  = psa_p.param_vectors{param_dim};
        n_levels  = numel(x_values);

        curves.(metric).(param).x = x_values;
        for ci = 1:numel(cond_spec)
            cond = cond_spec(ci).name;
            y = NaN(n_levels, numel(pcts));
            for level_idx = 1:n_levels
                % Same accessor plot_sensitivity uses: successful, non-NaN reps
                % for this (level, condition, metric).
                vals = ParamSpaceAnalysis2.collect_level_values( ...
                    psa_p, param, level_idx, cond, metric);
                if ~isempty(vals)
                    y(level_idx, :) = prctile(vals, pcts);
                end
            end
            curves.(metric).(param).y.(cond) = y;
        end
    end
end

%% -------------------- One figure per metric --------------------
% Same shape as the allStd script's metric_specs. Only the ylabel, the y window,
% the zero line and the output name differ:
%   LLE       -> lambda_1, window kept at [-1.75 1.75] to match the allStd
%                sheets. This run's LLEs span roughly p1 = -10 to p99 = +3.7, so
%                some medians (No Adaptation especially) run off the bottom and
%                are CLIPPED rather than rescaling every panel around them.
%   mean_rate -> [0 1] by construction; zero line dropped (it would sit on the
%                bottom axis and carry no meaning for a rate).
metric_specs = struct( ...
    'name',      {'LLE',                         'mean_rate'}, ...
    'ylabel',    {'\lambda_1',                   'Mean Firing Rate'}, ...
    'fig_tag',   {'Fig_Sensitivity_LLE_medians', 'Fig_sensitivity_mean_rate_medians'}, ...
    'ylim',      {[-1.75, 1.75],                 [0, 1]}, ...
    'yticks',    {[],                            [0, 1]}, ...
    'zero_line', {true,                          false});

made_tags = {};
for mi = 1:numel(metric_specs)
    spec = metric_specs(mi);

    fh = figure('Name', sprintf('%s Sensitivity medians', spec.name), ...
        'Position', [50, 50, 1300, 680], 'Color', 'w');
    % 'loose' horizontal spacing, not 'compact': each panel has its own x-axis
    % with labels running to both ends (e.g. "+100%", "1000"), and at compact
    % spacing the last label of one tile collides with the next tile's y-axis.
    % Padding 'loose' as well as spacing: with 'compact' the last x-tick label of
    % the right-hand column ("+100%") ends flush with the figure's right edge and
    % gets shaved off by the canvas boundary on export.
    tl = tiledlayout(fh, n_rows, n_cols, 'TileSpacing', 'loose', 'Padding', 'loose');

    ax_cell   = cell(1, numel(panel_params));
    leg_lines = gobjects(1, numel(cond_spec));

    for pi = 1:numel(panel_params)
        param = panel_params{pi};
        ax = nexttile(tl);
        ax_cell{pi} = ax;
        hold(ax, 'on');

        cx = curves.(spec.name).(param);

        % Zero line first, under the data. HandleVisibility off keeps it and the
        % default marker out of the legend.
        if spec.zero_line
            yline(ax, 0, '--', 'Color', zeroline_color, 'LineWidth', zeroline_lw, ...
                'Alpha', 0.5, 'HandleVisibility', 'off');
        end

        % Optional percentile band (inactive while band_pcts is empty).
        if numel(band_pcts) == 2
            for ci = 1:numel(cond_spec)
                y = cx.y.(cond_spec(ci).name);
                ok = all(isfinite(y(:, band_pcts)), 2);
                fill(ax, [cx.x(ok), fliplr(cx.x(ok))], ...
                     [y(ok, band_pcts(1))', fliplr(y(ok, band_pcts(2))')], ...
                     cond_spec(ci).color, 'FaceAlpha', band_alpha, ...
                     'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end

        % Median curves. Drawn in condition order, so the black baseline goes
        % down FIRST and the three coloured curves sit on top of it rather than
        % being hidden under it.
        for ci = 1:numel(cond_spec)
            y = cx.y.(cond_spec(ci).name)(:, median_col);
            h = plot(ax, cx.x, y, '-', 'Color', cond_spec(ci).color, ...
                'LineWidth', line_lw, 'DisplayName', cond_spec(ci).display);
            if pi == 1
                leg_lines(ci) = h;
            end
        end
        hold(ax, 'off');

        % --- Axes limits + ticks ------------------------------------------
        xlim(ax, [min(cx.x), max(cx.x)]);
        ylim(ax, spec.ylim);
        if ~isempty(spec.yticks)
            set(ax, 'YTick', spec.yticks);
        end
        set(ax, 'FontSize', tick_fs);
        box(ax, 'off');

        % ylabel on the left column only.
        if mod(pi - 1, n_cols) == 0
            ylabel(ax, spec.ylabel, 'Interpreter', 'tex', 'FontSize', label_fs);
        end

        % --- x-axis: percent scale where the default gives one -------------
        rs = row_style(param);
        has_default = isKey(default_value, param);
        % A zero default has no percent scale (x/0), so such an axis stays in raw
        % units even if it is on the percent list.
        use_pct = has_default && ismember(param, pct_params) && default_value(param) ~= 0;
        if use_pct
            apply_percent_axis(ax, default_value(param), rs.xlabel, label_fs);
        else
            xlabel(ax, rs.xlabel, 'Interpreter', 'tex', 'FontSize', label_fs);
            if ~isempty(rs.xticks)
                set(ax, 'XTick', rs.xticks);
                if ~isempty(rs.xticklabels)
                    set(ax, 'XTickLabel', rs.xticklabels);
                end
            end
        end

        % Default marker LAST: it re-pins xlim/ylim, so everything that sets
        % them must already have run.
        if has_default
            mark_default_value(ax, default_value(param), ...
                default_mark_color, default_mark_lw, default_mark_frac);
        end
    end

    % One legend for the whole sheet, as a horizontal strip ABOVE the tiles.
    % Inside a panel it would sit on top of the curves -- every panel's upper
    % half carries data somewhere across the four conditions -- so it goes in the
    % layout's own north tile, where it covers nothing.
    lg = legend(leg_lines, {cond_spec.display}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', legend_fs);
    lg.Layout.Tile = 'north';

    % Panel letters. Passing the ORDERED axes cell (not the figure) is what
    % guarantees (a)..(f) follow panel_params rather than a position sort.
    letters = panel_letters(numel(ax_cell));   % (a).. and past (z) correctly
    % HShift/VShift are larger than the allStd figure's -0.03/-0.04: the topmost
    % y-tick label sits at the axes' top-left corner (on mean_rate it is the "1"
    % of the [0 1] ticks, right at the top), so a letter placed just outside that
    % corner lands on top of it. These offsets clear it on both metrics.
    AddLetters2Plots(ax_cell, letters, ...
        'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.075, 'VShift', -0.05);

    save_figure_stable(out_dir, spec.fig_tag, fh);
    made_tags(end+1) = {spec.fig_tag}; %#ok<SAGROW>
end

% Log the git state alongside the figures so this presentation output can be

%% --- Record ------------------------------------------------------------------
out = struct('figs', gobjects(0), 'files', {{}}, 'source', run_dir);
if cfg.save
    for k = 1:numel(made_tags)
        out.files = [out.files, existing_outputs(out_dir, made_tags{k})];
    end






end
end








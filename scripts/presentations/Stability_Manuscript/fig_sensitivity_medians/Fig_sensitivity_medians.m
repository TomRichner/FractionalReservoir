close all
clc

% Fig_sensitivity_medians.m
% Stability_Manuscript presentation figure: the 1D sensitivity sweeps COLLAPSED
% ACROSS CONDITIONS. One subplot per swept parameter, four median lines per
% subplot (one per adaptation condition), 2 x 3 on a single sheet per metric.
%
% Sibling of ../fig_sensitivity_analysis_allStd/Fig_sensitivity_analysis_allStd.m,
% which is left untouched. That script lays the same run out as rows = swept
% parameter, columns = condition, each panel an imagesc of the full rep
% distribution -- 28 panels over two sheets per metric. Here the conditions are
% overlaid instead of tiled, which is what makes the comparison between them
% direct and fits everything on one sheet. Three consequences:
%
%   1. NO DISTRIBUTIONS. Only the median across reps is drawn. The percentile
%      machinery is written generally (see pcts / band_pcts below) so adding an
%      IQR band later is a two-line change, but nothing but the median is
%      plotted now.
%   2. SIX PARAMETERS, not seven. level_of_chaos ("Synaptic Gain") is dropped so
%      the remaining six fit 2 x 3.
%   3. NO REPLOT DETOUR. replot_sensitivity / plot_sensitivity /
%      assemble_sensitivity_figure exist to build the imagesc sheets; harvesting
%      median lines back out of their saved .fig files would be fragile. The
%      medians are computed straight off the saved PSA objects here, from the
%      same ParamSpaceAnalysis2.collect_level_values that plot_sensitivity
%      medians internally.
%
% NOTE ON SHARED HELPERS. preset_default_values / apply_percent_axis /
% mark_default_value / save_figure_stable now live as standalone files in
% scripts/run_all_analyses/replot/. They were lifted verbatim from
% Fig_sensitivity_analysis_allStd.m, which still carries its own local copies so
% its committed figures cannot shift; deleting those copies in favour of these is
% a safe follow-up once that figure is next regenerated.
%
% See also: Fig_sensitivity_analysis_allStd, Fig_memory_capacity_example

this_dir = fileparts(mfilename('fullpath'));
setup_paths();
project_root = fileparts(which('setup_paths'));

% Source run (a run_all_<dt> folder with 1D_sensitivity_* subdirs). Same run the
% allStd figures are built from -- swap this one line to re-point both.
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_aug_14_26_17_25');
out_dir   = this_dir;   % write the final figures next to this script

close all force;

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
default_value = preset_default_values(data_root, panel_params);

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
cond_spec = struct( ...
    'name',    {'no_adaptation', 'sfa_only',  'std_only',  'sfa_and_std'}, ...
    'display', {'No Adaptation', 'SFA Only',  'STD Only',  'SFA + STD'}, ...
    'color',   {[0.000 0.000 0.000], ...   % black          #000000
                [0.902 0.624 0.000], ...   % orange         #E69F00
                [0.337 0.706 0.914], ...   % sky blue       #56B4E9
                [0.800 0.475 0.655]});     % reddish purple #CC79A7

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
tick_fs   = 14;    % tick numbers -- matches the allStd + MC figures
label_fs  = 15.4;  % axis labels
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
sens_listing = dir(fullfile(data_root, '1D_sensitivity_*'));
sens_listing = sens_listing([sens_listing.isdir]);
if isempty(sens_listing)
    error('Fig_sensitivity_medians:NotFound', ...
        'No 1D_sensitivity_* subfolder found in:\n  %s', data_root);
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
         'Source run: %s'], strjoin(missing, ', '), data_root);
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
    panel_letters = arrayfun(@(ch) sprintf('(%c)', ch), ...
        char('a' + (0:numel(ax_cell)-1)), 'UniformOutput', false);
    % HShift/VShift are larger than the allStd figure's -0.03/-0.04: the topmost
    % y-tick label sits at the axes' top-left corner (on mean_rate it is the "1"
    % of the [0 1] ticks, right at the top), so a letter placed just outside that
    % corner lands on top of it. These offsets clear it on both metrics.
    AddLetters2Plots(ax_cell, panel_letters, ...
        'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.075, 'VShift', -0.05);

    save_figure_stable(out_dir, spec.fig_tag, fh);
    made_tags(end+1) = {spec.fig_tag}; %#ok<SAGROW>
end

% Log the git state alongside the figures so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(out_dir, 'README_fig_sensitivity_medians.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: sensitivity MEDIANS, collapsed across conditions\n');
fprintf(fid, '============================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'WHAT THIS IS\n');
fprintf(fid, '  The compact counterpart of ../fig_sensitivity_analysis_allStd. That figure\n');
fprintf(fid, '  tiles the four adaptation conditions as COLUMNS and shows the full rep\n');
fprintf(fid, '  distribution per panel as an imagesc histogram (28 panels over two sheets\n');
fprintf(fid, '  per metric). Here the conditions are OVERLAID in one axes per swept\n');
fprintf(fid, '  parameter, as median lines only, so all six parameters fit one 2 x 3 sheet\n');
fprintf(fid, '  and the conditions can be compared directly rather than across columns.\n');
fprintf(fid, '  No simulation is re-run.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  1D_sensitivity subfolders used:\n');
for pi = 1:numel(panel_params)
    fprintf(fid, '    %-16s\n', panel_params{pi});
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for k = 1:numel(made_tags)
    fprintf(fid, '  %s.png / .svg / .fig\n', made_tags{k});
end
fprintf(fid, '\n');

fprintf(fid, 'PANELS (row-major, 2 x 3)\n');
for pi = 1:numel(panel_params)
    rs = row_style(panel_params{pi});
    fprintf(fid, '  (%c) %-16s  x-axis: %s\n', char('a' + pi - 1), panel_params{pi}, rs.xlabel);
end
fprintf(fid, '\n');
fprintf(fid, '  level_of_chaos ("Synaptic Gain") is deliberately DROPPED: the allStd run\n');
fprintf(fid, '  sweeps seven parameters, which do not tile 2 x 3, and the gain sweep is the\n');
fprintf(fid, '  least surprising of the seven (it simply scales W).\n\n');

fprintf(fid, 'E:I NEURON RATIO. The f_E sweep varies the fraction excitatory with mu_EI and\n');
fprintf(fid, '  mu_IE held fixed, so what the axis really reports is the E:I neuron COUNT.\n');
fprintf(fid, '  This run has n = 500, so the ticks read 100:400 / 250:250 / 400:100 rather\n');
fprintf(fid, '  than the reduced ratios 1:4 / 1:1 / 4:1, which would hide the network size\n');
fprintf(fid, '  the counts come from.\n\n');

fprintf(fid, 'PERCENT AXES. The four mu axes are shown as percent departure from the preset\n');
fprintf(fid, '  default, (value/default - 1)*100. mu_EI and mu_II have NEGATIVE defaults, so\n');
fprintf(fid, '  "+100%%" means twice as inhibitory; those two panels have their x-direction\n');
fprintf(fid, '  REVERSED so that on all four rightward means "stronger synapse of this type".\n');
fprintf(fid, '  Only the ruler changes -- the plotted data is untouched.\n\n');

fprintf(fid, 'DEFAULT MARKER. Every panel carries a short reddish-gray tick rising from the\n');
fprintf(fid, '  x-axis (%g of the y-range, %g pt) at the preset default for that parameter.\n', ...
    default_mark_frac, default_mark_lw);
fprintf(fid, '  Resolved by preset_default_values() from the run''s own run_manifest.mat, not\n');
fprintf(fid, '  hardcoded. For this run:\n');
for pi = 1:numel(panel_params)
    if isKey(default_value, panel_params{pi})
        fprintf(fid, '    %-16s default %g\n', panel_params{pi}, default_value(panel_params{pi}));
    end
end
fprintf(fid, '\n');

fprintf(fid, 'CONDITION COLOURS (Okabe-Ito, colorblind-safe; same as fig_memory_capacity_example)\n');
for ci = 1:numel(cond_spec)
    c = cond_spec(ci).color;
    fprintf(fid, '  %-14s %-14s [%.3f %.3f %.3f]\n', ...
        cond_spec(ci).name, cond_spec(ci).display, c(1), c(2), c(3));
end
fprintf(fid, '  Conditions are matched to colours BY NAME, so a run declaring them in a\n');
fprintf(fid, '  different order cannot silently recolour the figure.\n\n');

fprintf(fid, 'STATISTIC. Median across reps at each swept-parameter level, per condition\n');
fprintf(fid, '  (ParamSpaceAnalysis2.collect_level_values -> prctile), i.e. exactly the blue\n');
fprintf(fid, '  median line of the allStd sheets. Failed / NaN reps are excluded first. The\n');
fprintf(fid, '  script computes prctile(vals, pcts) with pcts = %s and plots only the median;\n', mat2str(pcts));
fprintf(fid, '  adding e.g. an IQR band means extending pcts and setting band_pcts, whose\n');
fprintf(fid, '  shading branch is already written.\n\n');

fprintf(fid, 'PER-METRIC DIFFERENCES\n');
fprintf(fid, '  %s: ylabel \\lambda_1; y window [%.2f, %.2f], kept identical to the\n', ...
    metric_specs(1).fig_tag, metric_specs(1).ylim(1), metric_specs(1).ylim(2));
fprintf(fid, '    allStd sheets; green dashed zero line kept (the sign of lambda_1 marks the\n');
fprintf(fid, '    edge of chaos). NOTE: this preset''s LLEs span roughly p1 = -10.0 to\n');
fprintf(fid, '    p99 = +3.7, so several medians -- No Adaptation above all -- leave the\n');
fprintf(fid, '    window and are CLIPPED rather than rescaling every panel around them.\n');
fprintf(fid, '  %s: ylabel "Mean Firing Rate"; y window [0, 1] (nothing can\n', metric_specs(2).fig_tag);
fprintf(fid, '    fall outside it); y ticks at 0 and 1 only; zero line removed.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);

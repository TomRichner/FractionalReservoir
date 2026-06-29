function combined_fig = assemble_sensitivity_figure(replot_dir, metric)
% ASSEMBLE_SENSITIVITY_FIGURE Combine 1D sensitivity figs into one stacked plot.
%
% Loads every saved sensitivity figure under <replot_dir>/figures/ whose
% figure name starts with "<metric> Sensitivity - " (one per swept param),
% copies their per-condition axes into a single new figure with one row
% per swept param, and saves the result.
%
%   combined_fig = assemble_sensitivity_figure(replot_dir)
%   combined_fig = assemble_sensitivity_figure(replot_dir, metric)
%
% Inputs:
%   replot_dir : path to a replot_sensitivity_<dt> folder (containing figures/)
%   metric     : metric prefix used in figure Name (default 'LLE')

    if nargin < 2 || isempty(metric)
        metric = 'LLE';
    end

    fig_dir = fullfile(replot_dir, 'figures');
    listing = dir(fullfile(fig_dir, '*.fig'));
    if isempty(listing)
        warning('assemble_sensitivity_figure:NoFigs', ...
            'No .fig files in %s', fig_dir);
        combined_fig = [];
        return;
    end

    name_prefix = sprintf('%s Sensitivity - ', metric);

    % Open each .fig invisibly, keep only those matching the metric
    matched_figs = gobjects(0);
    matched_params = {};
    for k = 1:length(listing)
        fpath = fullfile(listing(k).folder, listing(k).name);
        f = openfig(fpath, 'invisible');
        fname = get(f, 'Name');
        if startsWith(fname, name_prefix)
            param_name = extractAfter(fname, name_prefix);
            matched_figs(end+1) = f; %#ok<AGROW>
            matched_params{end+1} = param_name; %#ok<AGROW>
        else
            close(f);
        end
    end

    if isempty(matched_figs)
        warning('assemble_sensitivity_figure:NoMatch', ...
            'No figures found with name starting "%s" in %s', name_prefix, fig_dir);
        combined_fig = [];
        return;
    end

    % Sort param rows alphabetically for determinism
    [matched_params, sort_idx] = sort(matched_params);
    matched_figs = matched_figs(sort_idx);

    % Collect per-fig axes ordered left-to-right
    nRows = length(matched_figs);
    src_axes = cell(nRows, 1);
    for r = 1:nRows
        ax = findobj(matched_figs(r), 'Type', 'axes');
        positions = cell2mat(get(ax, 'Position'));
        [~, order] = sort(positions(:, 1));
        src_axes{r} = ax(order);
    end
    nCols = length(src_axes{1});

    % Build combined figure
    combined_fig = figure('Name', sprintf('%s Sensitivity - combined', metric), ...
        'Position', [50, 50, 350 * nCols, 300 * nRows]);
    for r = 1:nRows
        for c = 1:nCols
            placeholder = subplot(nRows, nCols, (r - 1) * nCols + c, 'Parent', combined_fig);
            target_pos = get(placeholder, 'Position');
            delete(placeholder);
            copied_ax = copyobj(src_axes{r}(c), combined_fig);
            set(copied_ax, 'Position', target_pos);
        end
    end

    % Close the source (invisible) figures
    close(matched_figs);

    % Save combined figure
    out_base = fullfile(fig_dir, sprintf('sensitivity_%s_combined', metric));
    saveas(combined_fig, [out_base '.fig']);
    exportgraphics(combined_fig, [out_base '.png'], 'Resolution', 300);
    fprintf('Combined %s sensitivity figure saved to:\n  %s.{fig,png}\n', ...
        metric, out_base);

    % Close so the figure isn't picked up by later save_some_figs_to_folder_2 calls
    close(combined_fig);
    if nargout < 1
        clear combined_fig;
    end
end

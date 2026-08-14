function combined_fig = assemble_sensitivity_figure(replot_dir, metric, params, tag)
% ASSEMBLE_SENSITIVITY_FIGURE Combine 1D sensitivity figs into one stacked plot.
%
% Loads every saved sensitivity figure under <replot_dir>/figures/ whose
% figure name starts with "<metric> Sensitivity - " (one per swept param),
% copies their per-condition axes into a single new figure with one row
% per swept param, and saves the result.
%
%   combined_fig = assemble_sensitivity_figure(replot_dir)
%   combined_fig = assemble_sensitivity_figure(replot_dir, metric)
%   combined_fig = assemble_sensitivity_figure(replot_dir, metric, params, tag)
%
% Inputs:
%   replot_dir : path to a replot_sensitivity_<dt> folder (containing figures/)
%   metric     : metric prefix used in figure Name (default 'LLE')
%   params     : cellstr of swept-parameter names to include, in the ROW ORDER
%                given. Default {} means every parameter found, sorted
%                alphabetically. Use this to put a related group on one sheet --
%                the four mu blocks, say, which are only comparable side by side.
%   tag        : suffix for the output filename, so a filtered figure does not
%                overwrite the all-parameters one. Default '' -> "..._combined".
%
% A name in PARAMS with no matching figure is skipped with a warning rather
% than erroring: which parameters were swept depends on the model class, so
% asking for the mu blocks on an SRNNModel2 run is a no-op, not a failure.

    if nargin < 2 || isempty(metric)
        metric = 'LLE';
    end
    if nargin < 3
        params = {};
    end
    if nargin < 4 || isempty(tag)
        tag = '';
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

    % Row order. Only sort_idx is used from here on -- the source axes already
    % carry their own labels, so the names are needed for ordering, not drawing.
    if isempty(params)
        % Sort param rows alphabetically for determinism
        [~, sort_idx] = sort(matched_params);
    else
        % Keep only the requested params, in the order requested -- for a group
        % like the mu blocks the natural reading order is not the alphabetical
        % one. Missing names are reported, not fatal.
        sort_idx = [];
        for k = 1:numel(params)
            hit = find(strcmp(matched_params, params{k}), 1);
            if isempty(hit)
                warning('assemble_sensitivity_figure:MissingParam', ...
                    'No "%s" figure for parameter ''%s''; skipping that row.', ...
                    metric, params{k});
            else
                sort_idx(end+1) = hit; %#ok<AGROW>
            end
        end
        if isempty(sort_idx)
            close(matched_figs);
            warning('assemble_sensitivity_figure:NoRequestedParams', ...
                'None of the requested parameters were swept; no figure made.');
            combined_fig = [];
            return;
        end
        % Close the ones not wanted, before reindexing.
        close(matched_figs(setdiff(1:numel(matched_figs), sort_idx)));
    end
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
    if isempty(tag)
        combined_name = sprintf('%s Sensitivity - combined', metric);
    else
        combined_name = sprintf('%s Sensitivity - combined %s', metric, tag);
    end
    combined_fig = figure('Name', combined_name, ...
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

    % Save combined figure. The tag keeps a filtered figure from overwriting the
    % all-parameters one written by the earlier call.
    if isempty(tag)
        out_base = fullfile(fig_dir, sprintf('sensitivity_%s_combined', metric));
    else
        out_base = fullfile(fig_dir, sprintf('sensitivity_%s_combined_%s', metric, tag));
    end
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

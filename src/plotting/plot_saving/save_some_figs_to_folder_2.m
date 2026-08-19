function [] = save_some_figs_to_folder_2(save_folder, save_name, fig_vec, fig_type)

if not(exist(save_folder,'dir'))
    mkdir(save_folder)
else
    warning('folder already exists')
end

if isempty(fig_vec)
    figHandles = findobj('Type', 'figure')
    for i_f = 1:length(figHandles)
        fig_vec(i_f) = figHandles(i_f).Number;
    end
end

if isempty(fig_type)
    fig_type = {'fig','png'};
end
    
h = get(0,'children');
h = flipud(h);

% Initialize PDF report if requested
if any(strcmpi(fig_type, 'pdf'))
    pdf_file = fullfile(save_folder, [save_name '_report.pdf']);
    if exist(pdf_file, 'file')
        delete(pdf_file);
    end
end

% Each figure is saved inside a try/catch so that ONE failed export cannot
% abort the caller. This matters because run_all_analyses saves figures at the
% end of each stage: an export that throws there discards a completed
% multi-minute sweep and everything downstream of it, to lose a picture that can
% be regenerated from the saved .mat at any time.
%
% The failure this guards against is real and not understood. Roughly one run in
% four, a LIVE figure reaches a state where every rasterizing path --
% exportgraphics at any resolution, print -dpng, saveas png -- throws
% MATLAB:class:InvalidHandle, while the figure handle and all of its descendants
% report valid, saveas to .fig still works, a freshly built figure exports fine,
% and rebuilding the same figure from the same saved data exports fine. So it is
% neither a stale handle in this function nor a broken renderer nor anything
% about the data. Until it is diagnosed, warn and carry on.
failed = {};
for i=fig_vec
    try
        set(i,'PaperPositionMode','auto')
        if any(strcmpi(fig_type,'fig'))
            saveas(i, fullfile(save_folder, [save_name '_f_' num2str(i)]), 'fig');
        end
        if any(strcmpi(fig_type,'png'))
            exportgraphics(figure(i), fullfile(save_folder, [save_name '_figure_' num2str(i) '.png']), 'Resolution', 600)
        end
        if any(strcmpi(fig_type,'svg'))
            set(gcf, 'Renderer', 'painters');
            exportgraphics(figure(i), fullfile(save_folder, [save_name '_figure_' num2str(i) '.svg']), 'BackgroundColor', 'none', 'ContentType', 'vector');
        end
        if any(strcmpi(fig_type, 'pdf'))
            % Append to PDF report
            exportgraphics(figure(i), pdf_file, 'Append', true, 'ContentType', 'vector');
        end
    catch save_err
        failed{end+1} = sprintf('%d (%s)', i, save_err.identifier); %#ok<AGROW>
        warning('save_some_figs_to_folder_2:FigureSaveFailed', ...
            'Could not save figure %d as ''%s'': %s\nContinuing with the rest.', ...
            i, save_name, save_err.message);
    end
end
if ~isempty(failed)
    warning('save_some_figs_to_folder_2:SomeFiguresNotSaved', ...
        '%d of %d figures did not save for ''%s'': %s', ...
        numel(failed), numel(fig_vec), save_name, strjoin(failed, ', '));
end
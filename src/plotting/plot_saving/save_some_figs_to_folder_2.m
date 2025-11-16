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
    fig_type = {'fig','svg','png'};
end
    
h = get(0,'children');
h = flipud(h);

for i=fig_vec
    set(i,'PaperPositionMode','auto')  
    if any(strcmpi(fig_type,'fig'))
        saveas(i, fullfile(save_folder, [save_name '_f_' num2str(i)]), 'fig');
    end
    if any(strcmpi(fig_type,'png'))
        exportgraphics(figure(i), fullfile(save_folder, [save_name '_figure_' num2str(i) '.png']), 'Resolution', 300)
    end
    if any(strcmpi(fig_type,'svg'))
        set(gcf, 'Renderer', 'painters');
        saveas(i, fullfile(save_folder, [save_name '_figure_' num2str(i)]), 'svg');
    end
end
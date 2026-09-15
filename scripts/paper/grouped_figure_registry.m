function F = grouped_figure_registry(preset_name, source_fig_root, pytorch_file, human_psd_file)
% GROUPED_FIGURE_REGISTRY Append these entries after their source figures.
arguments
    preset_name char
    source_fig_root char
    pytorch_file char = ''
    human_psd_file char = ''
end
F=cell(1,8);
for k=1:8
    name=sprintf('fig_main%d_grouped',k);
    F{k}=struct('name',name,'fn',str2func(name),'in_paper',true,'args', ...
        {{'preset_name',preset_name,'source_fig_root',source_fig_root, ...
        'pytorch_file',pytorch_file,'human_psd_file',human_psd_file}});
end
end

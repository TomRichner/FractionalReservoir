% Compare fig_eig_heatmap's two density_scale options on one eig_heatmap_data.mat:
%   'log'    -> log10(1 + D)              (the default, the paper's figure)
%   'loglog' -> log10(1 + log10(1 + D))   (compresses the peaks further)
% Nothing is saved; one figure window pops up per option. Edit mat_file to
% point at another run's eig_heatmap folder.

setup_paths();
project_root = fileparts(which('setup_paths'));

mat_file = fullfile(project_root, 'data', 'single_multi_TS_independent_med', ...
    'eig_heatmap', 'eig_heatmap_data.mat');

scales = {'log', 'loglog'};
for k = 1:numel(scales)
    out = fig_eig_heatmap('data_file', mat_file, 'density_scale', scales{k}, ...
        'save', false);
    set(out.figs, 'Name', sprintf('density_scale = %s', scales{k}), 'NumberTitle', 'off');
end

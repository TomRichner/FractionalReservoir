function F = eigdens_figure_data(cfg, caller)
% EIGDENS_FIGURE_DATA What the eigenvalue-density figures share: the stage's file, one network, the conditions, their colors and names.
%
%   F = EIGDENS_FIGURE_DATA(cfg, caller)
%
%   cfg     the figure's cfg: data_file, run_dir, network, conditions
%   caller  the figure's name, for the error identifiers
%
% Reads <run_dir>/eigenvalue_density/eigenvalue_density_data.mat (the
% opt-in stage eigenvalue_density, run_eigenvalue_density_stage) through
% resolve_data_file: ONE run directory or an explicit file, no fallback.
% Keeps the runs of cfg.network (an error, <caller>:NoNetwork, listing the
% networks the file holds, when it has none) and the conditions
% select_conditions gives for cfg.conditions among the conditions that
% network holds, in the order of the stage's index.
%
% F: file, D (the whole file), network, net (its entry of D.networks),
% cond_names, by_cond (cell of run struct arrays, one per condition),
% colors (cell of RGB, manuscript_style's condition_color, so no
% adaptation follows manuscript_style; FD-044), titles
% (grouped_condition_label 'sentence': 'No adaptation', 'Single-timescale
% adaptation', ...), ticks (grouped_condition_label 'tick'), suffix (''
% for the reference network, '_<network>' otherwise, for the file tag),
% st (manuscript_style).
%
% See also: fig_eigdens_radial, fig_eigdens_gains, fig_eigdens_measures,
%           fig_eigdens_plane, run_eigenvalue_density_stage

F = struct();
F.file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'eigenvalue_density')}, 'eigenvalue_density_data.mat', ...
    'Run the opt-in stage eigenvalue_density (run_eigenvalue_density_stage) first');
F.D = load(F.file);
nets = {F.D.networks.name};
if ~ismember(cfg.network, nets)
    error([caller ':NoNetwork'], 'No runs of network ''%s'' in %s; it holds: %s.', ...
        cfg.network, F.file, strjoin(nets, ', '));
end
F.network = cfg.network;
F.net = F.D.networks(strcmp(nets, cfg.network));
runs = F.D.runs(strcmp({F.D.runs.network}, cfg.network));
available = unique({runs.condition}, 'stable');
F.cond_names = select_conditions(available, cfg.conditions);
F.by_cond = cellfun(@(c) runs(strcmp({runs.condition}, c)), F.cond_names, 'UniformOutput', false);
F.st = manuscript_style();
F.colors = cellfun(@(c) F.st.condition_color(c), F.cond_names, 'UniformOutput', false);
F.titles = grouped_condition_label(F.cond_names, 'sentence');
F.ticks = grouped_condition_label(F.cond_names, 'tick');
F.suffix = '';
if ~strcmp(cfg.network, 'reference'); F.suffix = ['_' cfg.network]; end
end

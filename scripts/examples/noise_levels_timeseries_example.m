% noise_levels_timeseries_example.m - the mu7chain network's four adaptation
% conditions at three levels of the Wiener noise sigma_u (0, 0.025, 0.05),
% one network and one noise seed in every panel.
%
% Draws fig_noise_levels_timeseries twice, for the four-condition set and the
% three-condition set (chained MTS in the product's place, as mu7chain_config
% orders them). About 40 s each; no pool needed.
%
% Output: figs/noise_levels/wiener_4cond/ and figs/noise_levels/wiener_3cond/
% under the project root, tag Fig_noise_levels_timeseries.
%
% The release-noise (gamma) version, once the model has release_rate and
% release_sigma (FD-027), is the same function with another level:
%   fig_noise_levels_timeseries('level_name', 'release_sigma', ...
%       'levels', [0 0.007 0.014], 'level_label', '\sigma_\gamma', ...
%       'model_args', {'sigma_u_noise', 0, 'release_rate', 0.01}, ...
%       'conditions', sets{k}, 'out_dir', fullfile(root, ['gamma_' set_names{k}]))
%
% See also: fig_noise_levels_timeseries, mu7chain_config

setup_paths();
project_root = fileparts(which('setup_paths'));
root = fullfile(project_root, 'figs', 'noise_levels');

set_names = {'4cond', '3cond'};
sets = {{'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}, ...
        {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}};
for k = 1:numel(sets)
    out = fig_noise_levels_timeseries('level_name', 'sigma_u_noise', ...
        'levels', [0 0.025 0.05], 'level_label', '\sigma_u', ...
        'conditions', sets{k}, 'out_dir', fullfile(root, ['wiener_' set_names{k}]), ...
        'visible', false);
    close(out.figs);
    fprintf('%s: std(x), rows = levels, columns = %s\n', set_names{k}, strjoin(out.cond_names, ', '));
    disp(out.x_std);
end

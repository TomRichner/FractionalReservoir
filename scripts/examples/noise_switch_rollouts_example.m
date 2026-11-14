% noise_switch_rollouts_example.m - copies of one network that share the
% noise for 10 s and then continue on independent noise seeds for 5 s.
%
% Runs run_noise_switch_rollouts and fig_noise_switch_rollouts with their
% defaults: the mu7chain network (n = 500, sigma_u_noise = 0.025, no
% stimulus), all four adaptation conditions, 8 copies each. About a minute
% on 4 pool workers; without a pool the copies run serially in the client.
% Give the session a pool first if one is wanted: pool_set(N).
%
% Output: data/noise_switch_rollouts/noise_switch_rollouts_data.mat and
% figs/noise_switch_rollouts/ (one time-series figure per condition, the
% summary figure and its table), both under the project root.
%
% See also: run_noise_switch_rollouts, fig_noise_switch_rollouts

setup_paths();
project_root = fileparts(which('setup_paths'));

data_dir = run_noise_switch_rollouts();

fig_dir = fullfile(project_root, 'figs', 'noise_switch_rollouts');
out = fig_noise_switch_rollouts('run_dir', fileparts(data_dir), 'out_dir', fig_dir, 'visible', false);
close(out.figs);

fprintf('data:    %s\n', out.source);
fprintf('figures: %s\n', fig_dir);
fprintf('  %s\n', out.files{:});

function [run_dir, fig_results] = ied_stim_run(run_mode, varargin)
% IED_STIM_RUN The IED stimulation block run and its figures, Fig 8 included (FD-031 items 4, 5).
%
%   ied_stim_run                  % 'essential': 2 networks, 1 + 1 blocks of 60 s, about 2 min on 3 workers
%   ied_stim_run('fast')          % 3 networks, 2 + 2 blocks of 120 s
%   ied_stim_run('production')    % the paper's run: 12 networks, 4 + 4 five-minute blocks (Tom's to start)
%   ied_stim_run(mode, 'clinical_dir', d, ...)   % name-value overrides of ied_stim_config
%   [run_dir, fig_results] = ied_stim_run(...)
%
% setup_paths is called on the first line. Every setting comes from
% ied_stim_config(run_mode). Simulation (ied_stim_blocks) needs the
% session's pool allocation first (pool_set(N)); it uses pool_ensure. Then
% ied_stim_figures draws every registered figure into cfg.fig_root:
%   fig_ied_stim_blocks/             model traces, normalised rate, PSD (Fig 8 D-F source)
%   fig_ied_stim_blocks_rate/        the rate panel alone
%   fig_ied_stim_blocks_supp_state/  supplement: rate, SFA feedback, depression 1 - b_1
%   Fig_Stim_Combined/               Fig 8 (clinical recordings: never git add -f)
%
% The simulation resumes (finished blocks are skipped) and refuses a run
% directory written with other settings. To redraw a finished run without
% simulating, call ied_stim_figures(ied_stim_config(run_mode)).
arguments
    run_mode (1,:) char = 'essential'
end
arguments (Repeating)
    varargin
end
setup_paths();
cfg = ied_stim_config(run_mode, varargin{:});
assert(~strcmp(run_mode, 'stim_blocks_20261117') || isfile(fullfile(cfg.run_dir, 'summary.mat')), ...
    'ied_stim_run:Legacy', ['stim_blocks_20261117 is the finished 400-Hz run; it is redrawn, ' ...
    'not resimulated. Use ied_stim_figures(ied_stim_config(''stim_blocks_20261117'')).']);
t = tic;
sim_args = namedargs2cell(cfg.sim);
if ~strcmp(run_mode, 'stim_blocks_20261117')
    ied_stim_blocks('data_dir', cfg.run_dir, sim_args{:});
end
fprintf('[ied_stim_run] simulation step done in %.1f min\n', toc(t)/60);
fig_results = ied_stim_figures(cfg);
run_dir = cfg.run_dir;

fprintf('\n========================================================\n');
fprintf('IED STIMULATION BLOCK RUN COMPLETE (%s)\n', run_mode);
fprintf('  fs      : %g Hz\n', cfg.sim.fs);
fprintf('  run_dir : %s\n', cfg.run_dir);
fprintf('  figures : %d of %d succeeded  (%s)\n', sum([fig_results.ok]), numel(fig_results), cfg.fig_root);
fprintf('  clinical: %s\n', cfg.clinical_dir);
fprintf('  total   : %.1f min\n', toc(t)/60);
fprintf('========================================================\n');
if nargout == 0, clear run_dir fig_results; end
end

% SINE_PAIRED_STIM_RUN Simulate the sine / paired-pulse stimulation experiment and draw its figure.
%
%   Open this file and press Run, or launch it through the MATLAB MCP server.
%   Every setting comes from sine_paired_stim_config().
%
%   run_dir   data/sine_paired_stim      fig_root  figs/sine_paired_stim
%
% RERUNNING THE SIMULATIONS REQUIRES DELETING data/sine_paired_stim FIRST.
% To redraw the figure only, call fig_sine_paired_stim as in the last lines.
% About 20 minutes of compute at 2 kHz (12 runs, single process).
setup_paths();
cfg = sine_paired_stim_config();
project_root = fileparts(which('setup_paths'));
run_dir = fullfile(project_root, cfg.run_dir);
if ~isfile(fullfile(run_dir, 'stim_evoked', 'stim_evoked_data.mat'))
    run_stim_evoked(cfg);
end
out = fig_sine_paired_stim('run_dir', run_dir, ...
    'out_dir', fullfile(project_root, cfg.fig_root), 'discard_s', cfg.discard_s, ...
    'verbose', cfg.verbose);
vprintf(cfg.verbose, 'near-none', 'sine_paired_stim: %d files in %s\n', numel(out.files), fullfile(project_root, cfg.fig_root));

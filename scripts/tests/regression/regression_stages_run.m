function [out_dir, seconds] = regression_stages_run(out_dir)
% REGRESSION_STAGES_RUN The middle regression tier: three analysis stages and their figures, in a few minutes.
%
%   out_dir = regression_stages_run()          % a new folder under tempdir
%   [out_dir, seconds] = regression_stages_run(out_dir)
%
% THREE TIERS (TR, 2026-11-10):
%   ultrarapid  regression_rapid_run    ~30 s    the three classes, directly
%   stages      THIS                    ~3 min   stage functions and figures
%   full        regression_fast_run     ~50 min  the whole bundle
%
% The ultrarapid tier never calls a stage function or a figure, and the full
% tier is too long for routine use. This one fills the gap: it calls three
% of the manuscript bundle's stage functions directly, on the 100-neuron
% regression network at run_mode 'fast', laid out as they would be inside a
% run directory, and then builds the figures that read them. So it exercises
% the stage wrappers, the files they save, and the figure code on top.
%
%   run_lyapunov_spectrum         -> OUT_DIR/lyapunov_spectrum/
%   run_transient_gain            -> OUT_DIR/transient_gain/
%   run_memory_capacity_example   -> OUT_DIR/mc_example/     (noise-free twin)
%   fig_lyapunov_spectrum, fig_transient_gain,
%   fig_memory_capacity_example, fig_STD_steady_state   -> OUT_DIR/figs/<name>/
%
% These three stages were chosen because they are the cheapest that still
% reach every part of the model class the stages use: the top-K spectrum
% with noise on and off, the tangent propagators with adaptation frozen and
% active, and the reservoir with its synaptic readout. fig_STD_steady_state
% needs no stage; it is here because it draws the depression model itself.
%
% regression_digest reduces OUT_DIR and regression_compare compares, as for
% the other tiers. The figure FILES are compared by name only, not by pixel.
%
% DO NOT RETUNE the settings below; see regression_rapid_run.
%
% See also: regression_stages_check, regression_rapid_run, regression_fast_run
arguments
    out_dir (1,:) char = ''
end
t_all = tic;
if isempty(out_dir)
    out_dir = fullfile(tempdir, ['regression_stages_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
end
if isfolder(out_dir) && numel(dir(out_dir)) > 2
    error('regression_stages_run:NotEmpty', 'Output folder is not empty: %s', out_dir);
end
if ~isfolder(out_dir); mkdir(out_dir); end

P    = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
P_mc = 'regression_pairs_noStim_noise0_dualStd_3cond';
quiet = 'near-none';

%% Stages
run_lyapunov_spectrum('preset_name', P, 'run_mode', 'fast', 'verbose', quiet, ...
    'out_dir', fullfile(out_dir, 'lyapunov_spectrum'));
run_transient_gain('preset_name', P, 'run_mode', 'fast', 'verbose', quiet, ...
    'out_dir', fullfile(out_dir, 'transient_gain'));
run_memory_capacity_example('preset_name', P_mc, 'run_mode', 'fast', 'verbose', quiet, ...
    'ode_solver', 'sra1', 'output_dir', fullfile(out_dir, 'mc_example'));

%% Figures that read them (and one that draws the depression model)
fig_dir = fullfile(out_dir, 'figs');
figs = { ...
    'fig_lyapunov_spectrum',       {'run_dir', out_dir}
    'fig_transient_gain',          {'run_dir', out_dir}
    'fig_memory_capacity_example', {'run_dir', out_dir}
    'fig_STD_steady_state',        {'preset_name', P}};
for k = 1:size(figs, 1)
    out = feval(figs{k, 1}, figs{k, 2}{:}, 'out_dir', fullfile(fig_dir, figs{k, 1}), ...
        'save', true, 'visible', false, 'verbose', quiet);
    if isfield(out, 'figs')
        close(out.figs(isgraphics(out.figs)));
    end
end

seconds = toc(t_all);
fprintf('regression_stages_run: %.1f s -> %s\n', seconds, out_dir);
end

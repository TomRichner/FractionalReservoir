function out_dir = run_param_space_analysis(ctx)
% RUN_PARAM_SPACE_ANALYSIS Multi-dimensional grid sweep over n, f and gain.
%
%   out_dir = RUN_PARAM_SPACE_ANALYSIS(ctx)
%   out_dir = RUN_PARAM_SPACE_ANALYSIS()           % standalone, class defaults
%
% Joint grid over network size, fraction excitatory and synaptic gain, run under
% all four adaptation conditions on the SAME W per grid point so the regimes are
% directly comparable. No reps axis: total simulations are n_levels^3 x 4.
%
% ctx comes from resolve_run_context('param_space', ...).
%
% RENAMED from run_param_space_analysis2.m. The trailing 2 was inherited from
% ParamSpaceAnalysis2 and meant nothing here -- there was never a
% run_param_space_analysis.m for it to disambiguate from.
%
% WAS A SCRIPT reading the master_* base-workspace protocol; see
% resolve_run_context.
%
% See also: resolve_run_context, ParamSpaceAnalysis2, run_all_analyses

arguments
    ctx struct = resolve_run_context('param_space')
end

setup_paths();

psa = ParamSpaceAnalysis2( ...
    'n_levels', ctx.n_levels, ...   % set by run_mode (fast=3, production=5)
    'batch_size', 25, ...           % configs per batch (for checkpointing)
    'note', '', ...                 % EMPTY, not 'param_space'. The folder name is
    ...                             % <folder_prefix>_<note>_nLevs_N_<dt>, and
    ...                             % folder_prefix already defaults to
    ...                             % 'param_space' -- so a note here doubles it
    ...                             % (param_space_param_space_nLevs_3_...).
    ...                             % ParamSpaceAnalysis2 drops the separator when
    ...                             % note is empty, giving param_space_nLevs_3_<dt>.
    ...                             % It used to be 'test_refactor', which is why
    ...                             % every such folder on disk is called
    ...                             % param_space_test_refactor_nLevs_*.
    'verbose', ctx.verbose);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = ctx.model_defaults;
if ~isempty(ctx.output_dir)
    psa.output_dir = ctx.output_dir;
end

%% Grid axes
% Same swept variables as run_sensitivity_analysis (n, f, level_of_chaos), but
% here as a joint multi-dimensional grid. n and level_of_chaos use the same
% ranges as the 1-D sweeps; only f is narrower here, 0.25-0.75 against 0.2-0.8,
% since a joint grid has far fewer levels per axis to spend.
psa.add_grid_parameter('n',              [100, 1000]);   % network size
psa.add_grid_parameter(ctx.f_param,      [0.25, 0.75]);  % fraction excitatory
psa.add_grid_parameter('level_of_chaos', [0.25, 2.5]);   % W scaling (edge of chaos)
% Widened from [0.5, 1.5] to match the 1-D sweep, which was rebased on measured
% edge-of-chaos crossings (see run_sensitivity_analysis for the numbers). The
% cost is real and worth knowing: this grid spends the SAME number of levels
% over a 2.25-wide gain axis instead of a 1.0-wide one, so it samples the region
% around gain = 1 roughly half as finely and puts more of the grid deep in the
% chaotic regime. Accepted so the two analyses keep describing the same span.

%% Conditions
% From the preset rather than PSA's built-in defaults, which are spelled in
% SRNNModel2's n_a_E / n_b_E vocabulary and would be passed verbatim to whatever
% model_class is set -- silently wrong for any other class.
psa.set_conditions(ctx.conditions);

%% Run
% Generates all combinations, randomizes execution order (so early stopping is
% representative), runs batched parfor with checkpoints, consolidates into
% per-condition MAT files. Results are saved during execution.
psa.run();

copyfile([mfilename('fullpath') '.m'], psa.output_dir);

% run() writes psa_object.mat itself -- once before batching so an interrupted
% run stays recoverable, and again on completion.
fprintf('PSA object saved to: %s\n', fullfile(psa.output_dir, 'psa_object.mat'));

%% Plot
% Colour by the fraction-excitatory axis. This has to be named explicitly for
% SRNNCellTypePairs: plot's default color_by is 'f', which there is a 1 x C row
% and blows up the scalar assignment the histogram colouring does.
psa.plot('metric', 'LLE',       'color_by', ctx.f_param);
psa.plot('metric', 'mean_rate', 'color_by', ctx.f_param);

if ctx.save_figs
    fig_dir = fullfile(psa.output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'param_space', [], {'fig', 'png'});
    fprintf('Figures saved to %s\n', fig_dir);
end

out_dir = psa.output_dir;

%% Summary
fprintf('\n=== Parameter Space Analysis Summary ===\n');
fprintf('Output directory: %s\n', out_dir);
fprintf('Grid parameters: %s\n', strjoin(psa.grid_params, ', '));
fprintf('Levels per parameter: %d\n', psa.n_levels);
fprintf('Total combinations: %d^%d = %d\n', ...
    psa.n_levels, numel(psa.grid_params), psa.n_levels^numel(psa.grid_params));
fprintf('Conditions: %s\n', ...
    strjoin(cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false), ', '));
end

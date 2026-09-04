function ctx = resolve_run_context(analysis, opts)
% RESOLVE_RUN_CONTEXT Settings one analysis needs, resolved from named arguments.
%
%   ctx = RESOLVE_RUN_CONTEXT('sensitivity')
%   ctx = RESOLVE_RUN_CONTEXT('param_space', 'preset_name', 'celltype_pairs', ...
%                             'run_mode', 'fast', 'output_dir', d, 'save_figs', true)
%
% Replaces two things at once.
%
% FIRST, the ~60 lines of preamble that were copy-pasted across
% run_sensitivity_analysis, run_tau_sensitivity_analysis and
% run_param_space_analysis: save_figs resolution, the run_mode default, the
% preset lookup, the model class, the conditions, and the
% merge_struct(preset, cfg.model) layering. Three copies meant three chances to
% drift, and they had already drifted in their defaults (`production` in two of
% them, `medium` in the orchestrator).
%
% SECOND, and the reason this exists at all, the `master_*` BASE-WORKSPACE
% PROTOCOL. The sub-scripts used to read master_output_dir, master_save_figs,
% master_model_overrides, master_model_class, master_conditions, run_mode and
% save_figs out of whatever workspace happened to be calling them, with
% exist(...,'var'). That had three costs:
%
%   * You could not tell what a sub-script needed without grepping for exist().
%   * A variable left behind by one run silently applied to the next --
%     run_overnight_queue.m existed ONLY to scrub them between pipelines, which
%     its own header admitted was the reason it was a script at all.
%   * The sub-scripts had to skip their own clear/clc when master_output_dir was
%     set, so that "am I being orchestrated?" leaked into their cleanup logic.
%
% All three vanish once the settings are arguments: a function has its own
% scope, so there is nothing to leak and nothing to scrub.
%
% INPUTS
%   analysis        which table row of analysis_run_config to read. 'dc_lle' is
%                   accepted but has no row -- it carries its own run_mode
%                   switch, because the DC staircase fixes T_range by
%                   construction (dc_levels x hold_dur) and the fs/T_range
%                   tuning the PSA analyses share does not apply. ctx.cfg is
%                   left empty in that case.
%   preset_name     srnn_param_preset key. Supplies the model overrides, the
%                   MODEL CLASS they are written for, and the ADAPTATION
%                   CONDITIONS -- a preset with its own depression routes can
%                   only express them through a condition.
%   run_mode        analysis_run_config mode; see run_mode_names() for the
%                   canonical list ('fast' | 'medium' | 'medium2' | 'production').
%   output_dir      shared run directory. '' means "let ParamSpaceAnalysis2
%                   create its own dated folder", which is what a standalone
%                   run wants.
%   save_figs       whether the caller wants extra figure formats written.
%
% OUTPUT ctx, with the fields every caller needs:
%   .analysis .run_mode .preset_name .verbose .save_figs .output_dir
%   .preset_defaults  the preset struct as srnn_param_preset returned it
%   .model_class      'SRNNModel2' | 'SRNNCellTypePairs'
%   .conditions       the four adaptation regimes, spelled for that class
%   .cfg              analysis_run_config output (empty for 'dc_lle')
%   .model_defaults   merge_struct(preset_defaults, cfg.model) -- preset FIRST,
%                     so run_mode keeps final say over its own knobs
%   .n_levels .n_reps sweep size (n_reps absent for 'param_space')
%   .integer_params   {'n','indegree'} -- NOT the class default, which lists
%                     SRNNModel2's adaptation counts and is wrong for Pairs
%   .f_param          'f' on SRNNModel2, 'f_E' on SRNNCellTypePairs. The
%                     fraction-excitatory axis is a scalar property on one and a
%                     scalar ALIAS onto a 1 x C row on the other, so the axis
%                     name differs even though the quantity does not.
%
% See also: analysis_run_config, srnn_param_preset, validate_preset_conditions,
%           merge_struct, run_all_analyses

arguments
    analysis (1,:) char {mustBeMember(analysis, ...
        {'sensitivity', 'tau_sensitivity', 'param_space', 'dc_lle'})}
    opts.preset_name (1,:) char = 'default'
    opts.run_mode    (1,:) char = 'production'
    opts.output_dir  (1,:) char = ''
    opts.save_figs   (1,1) logical = false
    opts.verbose     (1,1) logical = true
end

ctx = struct();
ctx.analysis    = analysis;
ctx.run_mode    = opts.run_mode;
ctx.preset_name = opts.preset_name;
ctx.output_dir  = opts.output_dir;
ctx.save_figs   = opts.save_figs;
ctx.verbose     = opts.verbose;

% The preset carries three things, not one: the overrides, the class they are
% written for, and the conditions. Taking all three from the same call is what
% stops a Pairs preset being swept with SRNNModel2-shaped conditions.
[ctx.preset_defaults, ctx.model_class, ctx.conditions] = ...
    srnn_param_preset(opts.preset_name);

% Resolved BEFORE analysis_run_config, which needs it: a preset carrying
% sigma_u_noise > 0 selects that mode's STOCHASTIC integrator in place of its
% deterministic one.
if strcmp(analysis, 'dc_lle')
    ctx.cfg            = struct([]);
    ctx.model_defaults = ctx.preset_defaults;
else
    ctx.cfg = analysis_run_config(analysis, opts.run_mode, ctx.preset_defaults);
    % Preset first, run_mode timings second, so run_mode keeps final say over
    % ode_solver / fs / T_range / lya_T_interval, and so a whole-struct preset
    % assignment cannot clobber them.
    ctx.model_defaults = merge_struct(ctx.preset_defaults, ctx.cfg.model);

    ctx.n_levels = ctx.cfg.n_levels;
    if isfield(ctx.cfg, 'n_reps')
        ctx.n_reps = ctx.cfg.n_reps;
    end
end

% ParamSpaceAnalysis2's class default for integer_params lists SRNNModel2's
% adaptation counts, which are not SRNNCellTypePairs properties at all. Only n
% and indegree are integer-valued axes common to both classes.
ctx.integer_params = {'n', 'indegree'};

if strcmp(ctx.model_class, 'SRNNCellTypePairs')
    ctx.f_param = 'f_E';
else
    ctx.f_param = 'f';
end

if ctx.verbose
    report(ctx);
end
end

function report(ctx)
% One line per run, naming everything that decides what the numbers mean.
if isempty(ctx.cfg)
    fprintf('[%s] run_mode=%s, preset=%s (%s)\n', ...
        ctx.analysis, ctx.run_mode, ctx.preset_name, ctx.model_class);
    return
end

if isfield(ctx, 'n_reps')
    reps_str = sprintf(', n_reps=%d', ctx.n_reps);
else
    reps_str = '';       % param_space has no reps axis
end
fprintf('[%s] run_mode=%s, preset=%s (%s), n_levels=%d%s\n', ...
    ctx.analysis, ctx.run_mode, ctx.preset_name, ctx.model_class, ...
    ctx.n_levels, reps_str);
fprintf('[%s] ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    ctx.analysis, ctx.cfg.model.ode_solver, ctx.cfg.model.fs, ...
    ctx.cfg.model.T_range(1), ctx.cfg.model.T_range(2));
if ctx.cfg.is_stochastic
    % Worth its own line: the integrator was chosen by the PRESET's noise, not
    % by the run mode, which is the one place the two knobs are not orthogonal.
    fprintf('[%s] stochastic: sigma_u_noise=%g, integrator=%s\n', ...
        ctx.analysis, ctx.preset_defaults.sigma_u_noise, ctx.cfg.model.ode_solver);
end
end

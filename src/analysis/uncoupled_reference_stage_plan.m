function plan = uncoupled_reference_stage_plan(cfg)
% UNCOUPLED_REFERENCE_STAGE_PLAN The runs a bundle's uncoupled_reference stage makes, validated before any compute.
%
%   plan = UNCOUPLED_REFERENCE_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'uncoupled_reference'
% (paper_opt_in_stages; Plan 2, package C, 2026-11-21). For each condition,
% the leading growth rate of the bundle's own network with every recurrent
% weight scaled to near zero: what a neuron's adaptive variables do by
% themselves, under the same input and noise. It is the dotted line of
% Figure 4, panel C (fig_leading_direction): the joint sample's contracting
% configurations sit at this rate, which is why the panel can say that what
% decays slowly there is single-neuron adaptation and not a network mode.
% Config in, plan out; nothing is computed and no file is touched.
% run_uncoupled_reference_stage carries the plan out.
%
% WHY A STAGE AND NOT THE GAIN SWEEP. The 1-D sweep over overall gain
% (level_of_chaos) stops at 0.25 of the preset's weights. In the paper run
% of 2026-11-15 the adapting conditions had nearly settled there (median
% lambda_1 -0.49 1/s at 0.25 and -0.46 at 0.40 with one timescale; -0.124
% and -0.122 with the chain), but no adaptation had not: -1.94 1/s at 0.25
% against -1/tau_d = -10 1/s uncoupled. One more low level on the sweep
% would have changed the sweep's axis, and with it the sensitivity figures
% and their statistics. This stage adds a few runs and changes nothing
% else.
%
% WHAT IT RUNS. One ParamSpaceAnalysis2 with the joint sample's own model
% settings (resolve_run_context('param_space', ...): the same network size,
% input, noise, run length, integrator, rate and Lyapunov window as the
% points of the panel), a single axis level_of_chaos at the listed gains,
% and a repetition axis (each repetition a different network seed). Nothing
% else is varied.
%
% THE CONFIG
%   cfg.preset_name, cfg.run_mode   the bundle's network and size. Required.
%   cfg.run_conditions   the conditions simulated (bundle_run_conditions);
%                        empty = every condition the preset states
%   cfg.fs, cfg.ode_solver   the bundle's rate and integrator, as the joint
%                        sample takes them
%   cfg.uncoupled_reference   optional struct of overrides:
%     gains    the values of level_of_chaos, default [0.01 0.05]. The FIRST
%              is the reference the figure draws; the others are there so
%              the table can show that the rate has stopped changing.
%              Each in (0, 0.25): at or above the sweep's lowest level it
%              would not be a reference.
%     n_reps   repetitions per gain (network seeds); default by run mode
%     model    a struct merged over the model settings, for tests at tiny
%              size (n, T_range, lya_T_interval, ...)
%   A field the stage does not know is refused.
%
% REPETITIONS by run mode: essential 2, fast 4, fastmedium, medium,
% medium2 8, production 12. A run of the joint sample took 6 to 13 s at
% fastmedium (2026-11-15), so the stage is minutes at any size.
%
% THE PLAN
%   plan.folder    where it writes, relative to the run directory
%                  ('uncoupled_reference')
%   plan.gains, plan.n_reps, plan.conditions, plan.preset_name
%   plan.ctx       the joint sample's run context with the overrides merged
%   plan.verbose
%
% See also: run_uncoupled_reference_stage, fig_leading_direction,
%           resolve_run_context, run_param_space_analysis, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

%            mode           repetitions
REPS = {    'essential',    2; ...
            'fast',         4; ...
            'fastmedium',   8; ...
            'medium',       8; ...
            'medium2',      8; ...
            'production',  12};
GAINS_DEFAULT = [0.01 0.05];
SWEEP_LOWEST = 0.25;       % the lowest level of the 1-D sweep and the joint sample
ALLOWED = {'gains', 'n_reps', 'model'};

assert(isequal(sort(REPS(:, 1)), sort(run_mode_names())'), 'uncoupled_reference_stage_plan:TableIncomplete', ...
    'The repetition table must name every run mode of run_mode_names.');
if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, REPS(:, 1))
    bad_config('The uncoupled_reference stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(REPS(:, 1)', ', '));
end
preset_name = field_or_empty(cfg, 'preset_name');
if isempty(preset_name)
    bad_config('The uncoupled_reference stage runs on cfg.preset_name, which this bundle does not set.');
end

opts = struct();
if isfield(cfg, 'uncoupled_reference') && ~isempty(cfg.uncoupled_reference)
    opts = cfg.uncoupled_reference;
end
if ~(isstruct(opts) && isscalar(opts))
    bad_config('cfg.uncoupled_reference must be a scalar struct.');
end
unknown = setdiff(fieldnames(opts), ALLOWED);
if ~isempty(unknown)
    bad_config('cfg.uncoupled_reference has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(ALLOWED, ', '));
end

row = REPS(strcmp(REPS(:, 1), cfg.run_mode), :);
gains  = override(opts, 'gains', GAINS_DEFAULT);
n_reps = override(opts, 'n_reps', row{2});
model  = override(opts, 'model', struct());
if ~(isnumeric(gains) && isvector(gains) && all(isfinite(gains)) && all(gains > 0) && all(gains < SWEEP_LOWEST) && ...
        numel(unique(gains)) == numel(gains))
    bad_config(['cfg.uncoupled_reference.gains must be distinct values of level_of_chaos in (0, %g): ' ...
        'the sweep''s lowest level is %g, and the reference has to lie below it.'], SWEEP_LOWEST, SWEEP_LOWEST);
end
if ~(isscalar(n_reps) && isnumeric(n_reps) && isfinite(n_reps) && n_reps == round(n_reps) && n_reps >= 1)
    bad_config('cfg.uncoupled_reference.n_reps must be a whole number of at least 1.');
end
if ~(isstruct(model) && isscalar(model))
    bad_config('cfg.uncoupled_reference.model must be a scalar struct of model settings.');
end
if isfield(model, 'level_of_chaos')
    bad_config('cfg.uncoupled_reference.model must not set level_of_chaos: that is the stage''s axis (gains).');
end

%% The joint sample's context: the same model settings as the panel's points
verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
ctx_args = {'preset_name', preset_name, 'run_mode', cfg.run_mode, 'verbose', 'near-none', ...
    'conditions', bundle_run_conditions(cfg)};
fs = field_or_empty(cfg, 'fs');
if ~isempty(fs) && fs > 0; ctx_args = [ctx_args, {'fs', fs}]; end
solver = field_or_empty(cfg, 'ode_solver');
if ~isempty(solver); ctx_args = [ctx_args, {'ode_solver', solver}]; end
ctx = resolve_run_context('param_space', ctx_args{:});
if ~strcmp(ctx.model_class, 'SRNNCellTypePairs')
    bad_config('The uncoupled_reference stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        preset_name, ctx.model_class);
end
fn = fieldnames(model);
for k = 1:numel(fn)
    ctx.model_defaults.(fn{k}) = model.(fn{k});
end

plan = struct();
plan.folder      = 'uncoupled_reference';
plan.preset_name = preset_name;
plan.gains       = reshape(gains, 1, []);
plan.n_reps      = n_reps;
plan.conditions  = cellfun(@(c) c.name, ctx.conditions, 'UniformOutput', false);
plan.ctx         = ctx;
plan.verbose     = verbose;
end

%% ------------------------------------------------------------------------
function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('uncoupled_reference_stage_plan:BadConfig', varargin{:});
end

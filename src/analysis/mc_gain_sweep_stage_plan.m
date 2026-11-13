function plan = mc_gain_sweep_stage_plan(cfg)
% MC_GAIN_SWEEP_STAGE_PLAN The runs a bundle's mc_gain_sweep stage makes, validated before any compute.
%
%   plan = MC_GAIN_SWEEP_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'mc_gain_sweep' (block C,
% FD-010; paper_opt_in_stages). Config in, plan out; nothing is computed
% and no file is touched. run_mc_gain_sweep_stage carries the plan out.
%
% THE CONFIG
%   cfg.mc_preset       the network: a noise-free memory-capacity preset
%   cfg.mc_ode_solver   its integrator; absent or empty = run_mc_gain_sweep's
%   cfg.run_mode        sets the grid's size (run_mc_gain_sweep's table)
%   cfg.mc_gain_sweep   optional struct. Every field is optional, and one
%                       that is absent or empty leaves run_mc_gain_sweep's
%                       own value for the mode:
%     gain_levels, bias_levels, seeds, conditions
%     model_args, protocol, lya_K, lya_T_sec, lya_acc_sec, segment_len
%                       run_mc_gain_sweep's options of those names (the
%                       second line is for tests)
%     variants          struct array with fields tag and preset: further
%                       networks, each swept along GAIN ONLY, at its own
%                       bias (mc_gain_sweep_preset_point), with the main
%                       run's gains, conditions and options, into
%                       <run_dir>/mc_gain_sweep_<tag>
%                       (mc_gain_sweep_variant_folder)
%     variant_seeds     seeds of the variant runs; absent = the main run's
%   A field the stage does not know is refused, so a typo cannot be ignored.
%
% THE PLAN
%   plan.runs   struct array, the main run first, then the variants:
%     label     'main', or the variant's tag
%     folder    where it writes, relative to the run directory ('' for the
%               main run, whose results land in <run_dir>/mc_gain_sweep)
%     args      name-value arguments for run_mc_gain_sweep, all but
%               output_dir
%
% See also: run_mc_gain_sweep_stage, run_mc_gain_sweep, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

GRID_OPTIONS  = {'gain_levels', 'bias_levels', 'seeds', 'conditions'};
OTHER_OPTIONS = {'model_args', 'protocol', 'lya_K', 'lya_T_sec', 'lya_acc_sec', 'segment_len'};
STAGE_FIELDS  = {'variants', 'variant_seeds'};

if ~isfield(cfg, 'mc_preset') || isempty(cfg.mc_preset)
    error('mc_gain_sweep_stage_plan:BadConfig', ...
        'The mc_gain_sweep stage runs on cfg.mc_preset, which this bundle does not set.');
end
opts = struct();
if isfield(cfg, 'mc_gain_sweep') && ~isempty(cfg.mc_gain_sweep); opts = cfg.mc_gain_sweep; end
if ~(isstruct(opts) && isscalar(opts))
    error('mc_gain_sweep_stage_plan:BadConfig', 'cfg.mc_gain_sweep must be a scalar struct.');
end
allowed = [GRID_OPTIONS, OTHER_OPTIONS, STAGE_FIELDS];
unknown = setdiff(fieldnames(opts), allowed);
if ~isempty(unknown)
    error('mc_gain_sweep_stage_plan:BadConfig', ...
        'cfg.mc_gain_sweep has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(allowed, ', '));
end

% What every run of the stage shares.
shared = {'run_mode', cfg.run_mode, 'verbose', cfg.verbose};
if isfield(cfg, 'mc_ode_solver') && ~isempty(cfg.mc_ode_solver)
    shared = [shared, {'ode_solver', cfg.mc_ode_solver}];
end
shared = [shared, given(opts, [setdiff(GRID_OPTIONS, {'bias_levels', 'seeds'}, 'stable'), OTHER_OPTIONS])];

runs = struct('label', 'main', 'folder', '', 'args', ...
    {[{'preset_name', cfg.mc_preset}, shared, given(opts, {'bias_levels', 'seeds'})]});

variants = struct('tag', {}, 'preset', {});
if isfield(opts, 'variants') && ~isempty(opts.variants); variants = opts.variants; end
if ~(isstruct(variants) && all(isfield(variants, {'tag', 'preset'})))
    error('mc_gain_sweep_stage_plan:BadConfig', ...
        'cfg.mc_gain_sweep.variants must be a struct array with fields tag and preset.');
end
variant_seeds = given(opts, {'seeds'});
if isfield(opts, 'variant_seeds') && ~isempty(opts.variant_seeds)
    variant_seeds = {'seeds', opts.variant_seeds};
end
model_args = {};
if isfield(opts, 'model_args') && ~isempty(opts.model_args); model_args = opts.model_args; end
for k = 1:numel(variants)
    v = variants(k);
    [~, own_bias] = mc_gain_sweep_preset_point(v.preset, model_args);   % errors on an unknown preset
    runs(end+1) = struct('label', v.tag, 'folder', mc_gain_sweep_variant_folder(v.tag), 'args', ...
        {[{'preset_name', v.preset}, shared, variant_seeds, {'bias_levels', own_bias}]}); %#ok<AGROW>
end
labels = {runs.label};
if numel(unique(labels)) < numel(labels)
    error('mc_gain_sweep_stage_plan:BadConfig', ...
        'cfg.mc_gain_sweep.variants must have distinct tags, none of them ''main''; got %s.', ...
        strjoin(labels(2:end), ', '));
end
plan = struct('runs', runs);
end

%% ------------------------------------------------------------------------
function args = given(opts, names)
% Name-value pairs of the named fields that opts sets (present and non-empty).
args = {};
for k = 1:numel(names)
    if isfield(opts, names{k}) && ~isempty(opts.(names{k}))
        args = [args, {names{k}, opts.(names{k})}]; %#ok<AGROW>
    end
end
end

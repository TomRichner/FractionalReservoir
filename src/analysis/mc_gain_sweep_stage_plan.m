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
%   cfg.run_mode        sets every grid's size (mc_gain_sweep_grid's table:
%                       levels per axis and networks)
%   cfg.mc_gain_sweep   optional struct. Every field is optional; one that
%                       is absent or empty leaves the default.
%     grids             struct array, ONE ELEMENT PER GRID the stage sweeps,
%                       each one scalar model property against another
%                       (mc_sweep_axis; mc_gain_sweep_grid has the defaults):
%       tag             names the grid: its folder in the run directory
%                       (mc_gain_sweep_folder) and how a figure asks for it
%                       (resolve_mc_gain_sweep's 'grid'). Required.
%       first_axis      the model property of its first axis (default
%                       'level_of_chaos', the gain)
%       second_axis     the model property of its second axis (default
%                       'S_c', the bias)
%       first_range, first_spacing   the first axis: [lo hi] and 'log' or
%                       'linear'; for the gain also gain_range and
%                       gain_spacing (default [0.1 4], 'log'). Another
%                       property must give a range or its levels, and is
%                       spaced linearly unless told.
%       second_range, second_spacing   the second axis: [lo hi] and 'linear'
%                       (the default), 'log' or 'zero_log' (the bias has the
%                       default range [-0.6 0.6]; another axis must give a
%                       range or its levels). Either axis takes any of the
%                       three spacings; mc_gain_sweep_grid defines them.
%       first_levels, second_levels   explicit levels, in place of the range
%                       and the mode's count (gain_levels and bias_levels on
%                       the default axes)
%       model_args      name-value overrides of the preset for THIS grid's
%                       runs (its variants' too), applied after the stage's
%                       model_args: the arguments the grid holds fixed (a
%                       noise grid's leak, the other noise off)
%       readout_signal  the memory readout of this grid's runs, in place of
%                       the stage's ('rate', 'synaptic', 'release_filtered';
%                       run_mc_gain_sweep)
%       variants        struct array with fields tag and preset: further
%                       networks, each swept along THE FIRST AXIS ONLY, on
%                       this grid's levels of it, at its own value of the
%                       grid's second axis (its bias, by default;
%                       mc_gain_sweep_preset_point)
%       variant_seeds   seeds of the variant runs; absent = the grid's
%     run_grids         cell of grid tags: the grids that RUN, in the order
%                       of grids (FD-029). Absent or empty = every grid
%                       listed. A grid not named stays defined (a bundle's
%                       figure entries can follow the same list, and one
%                       tag added brings the grid back) and is validated
%                       like any other, but its runs, its variants' too,
%                       are left out of the plan and of its estimate. A
%                       tag that names no grid is refused; it needs the
%                       list form (grids).
%     seeds             network seeds of every grid (default: the mode's)
%     conditions, model_args, readout_signal, protocol, lya_K, lya_T_sec,
%     lya_acc_sec, segment_len   run_mc_gain_sweep's options of those names,
%                       the same for every run (the last five are for
%                       tests); a grid adds its own model_args after these
%                       and may name its own readout_signal
%   ONE UNNAMED GRID, the form from before a stage could hold several, is
%   still read: with no grids field, the grid's own fields (first_axis,
%   second_axis, gain_range, ..., variants, variant_seeds) are fields of
%   cfg.mc_gain_sweep itself (its model_args and readout_signal are then
%   the stage's).
%   That grid has the empty tag: it writes to <run_dir>/mc_gain_sweep and
%   its variants to <run_dir>/mc_gain_sweep_<variant>, where the results of
%   such a bundle have always been.
%   A field the stage does not know is refused, so a typo cannot be ignored.
%
% THE PLAN
%   plan.runs   struct array: each grid, followed by its variants
%     label     the grid's tag, or <grid>_<variant>; for the unnamed grid
%               'main' and the variant's tag
%     grid, variant   the two tags ('' where there is none)
%     folder    where it writes, relative to the run directory
%               (mc_gain_sweep_folder)
%     args      name-value arguments for run_mc_gain_sweep, all but
%               output_dir; both axes, their levels and the seeds are
%               explicit
%     n_jobs    its jobs (first-axis levels x second-axis levels x seeds)
%   plan.n_jobs, plan.worker_minutes   the stage's total, and the estimate
%               of its compute (mc_gain_sweep_grid's worker-minutes a job,
%               at the release-noise cost for a grid whose axis is a
%               release-noise amplitude, release_sigma or release_sigma_all,
%               or whose model_args set one above zero; a preset that
%               carries release noise of its own is not looked at)
%
% See also: run_mc_gain_sweep_stage, run_mc_gain_sweep_set, mc_gain_sweep_grid,
%           mc_gain_sweep_folder, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

RUN_OPTIONS  = {'conditions', 'model_args', 'readout_signal', 'protocol', 'lya_K', 'lya_T_sec', ...
                'lya_acc_sec', 'segment_len'};
GRID_SHAPE   = {'first_axis', 'second_axis', 'first_range', 'first_spacing', 'first_levels', ...
                'gain_range', 'gain_spacing', 'gain_levels', 'second_range', 'second_spacing', ...
                'second_levels', 'bias_levels'};
% A grid's own run options: added to (model_args) or replacing (readout_signal)
% the stage's, for that grid and its variants only.
GRID_RUN     = {'model_args', 'readout_signal'};
GRID_FIELDS  = [{'tag'}, GRID_SHAPE, GRID_RUN, {'variants', 'variant_seeds'}];
% The options that go by a default axis's own name: refused on another axis.
DEFAULT_NAMED = struct('first', {{'gain_range', 'gain_spacing', 'gain_levels'}}, ...
                       'second', {{'bias_levels'}});
STAGE_FIELDS = {'grids', 'run_grids', 'seeds'};

if ~isfield(cfg, 'mc_preset') || isempty(cfg.mc_preset)
    bad_config('The mc_gain_sweep stage runs on cfg.mc_preset, which this bundle does not set.');
end
opts = struct();
if isfield(cfg, 'mc_gain_sweep') && ~isempty(cfg.mc_gain_sweep); opts = cfg.mc_gain_sweep; end
if ~(isstruct(opts) && isscalar(opts))
    bad_config('cfg.mc_gain_sweep must be a scalar struct.');
end
% The unnamed grid's fields sit beside the stage's; its run options ARE the
% stage's (model_args, readout_signal), so they are not listed twice.
single_grid_fields = setdiff(GRID_FIELDS, [{'tag'}, GRID_RUN], 'stable');
allowed = [STAGE_FIELDS, RUN_OPTIONS, single_grid_fields];
unknown = setdiff(fieldnames(opts), allowed);
if ~isempty(unknown)
    bad_config('cfg.mc_gain_sweep has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(allowed, ', '));
end

%% The grids: the list, or the one unnamed grid written at the top level
if isfield(opts, 'grids') && ~isempty(opts.grids)
    grids = opts.grids;
    if ~(isstruct(grids) && isfield(grids, 'tag'))
        bad_config('cfg.mc_gain_sweep.grids must be a struct array with a field tag.');
    end
    unknown = setdiff(fieldnames(grids), GRID_FIELDS);
    if ~isempty(unknown)
        bad_config('cfg.mc_gain_sweep.grids has unknown field(s) %s. Allowed: %s.', ...
            strjoin(unknown', ', '), strjoin(GRID_FIELDS, ', '));
    end
    misplaced = names_given(opts, single_grid_fields);
    if ~isempty(misplaced)
        bad_config(['cfg.mc_gain_sweep lists grids, so %s belong(s) in a grid, ' ...
            'not beside the list.'], strjoin(misplaced, ', '));
    end
    if any(cellfun(@isempty, {grids.tag}))
        bad_config('Every grid of cfg.mc_gain_sweep.grids needs a tag.');
    end
    run_grids = field_or_empty(opts, 'run_grids');
    if ~isempty(run_grids)
        if ischar(run_grids) || isstring(run_grids); run_grids = cellstr(run_grids); end
        if ~iscellstr(run_grids)
            bad_config('cfg.mc_gain_sweep.run_grids must be a cell of grid tags.');
        end
        unknown = setdiff(run_grids, {grids.tag});
        if ~isempty(unknown)
            bad_config('cfg.mc_gain_sweep.run_grids names %s, which no grid has as its tag (the grids: %s).', ...
                strjoin(unknown, ', '), strjoin({grids.tag}, ', '));
        end
    else
        run_grids = {grids.tag};
    end
else
    if ~isempty(field_or_empty(opts, 'run_grids'))
        bad_config('cfg.mc_gain_sweep.run_grids selects among listed grids; this config lists none (no grids field).');
    end
    grids = cell2struct([{''}; cellfun(@(f) field_or_empty(opts, f), single_grid_fields(:), ...
        'UniformOutput', false)], [{'tag'}; single_grid_fields(:)], 1);
    run_grids = {''};
end
% EVERY grid is planned, so a grid that does not run is still validated (a
% typo in it fails now, not on the day it is brought back); only the runs
% of the grids in run_grids are kept, and only they count in the estimate.

%% What every run of the stage shares
shared = {'run_mode', cfg.run_mode, 'verbose', cfg.verbose};
if isfield(cfg, 'mc_ode_solver') && ~isempty(cfg.mc_ode_solver)
    shared = [shared, {'ode_solver', cfg.mc_ode_solver}];
end
shared = [shared, given(opts, setdiff(RUN_OPTIONS, GRID_RUN, 'stable'))];
stage_model_args = field_or_empty(opts, 'model_args');
if isempty(stage_model_args); stage_model_args = {}; end
stage_readout = field_or_empty(opts, 'readout_signal');
seeds = given(opts, {'seeds'});

%% One run per grid, then one per variant network of that grid
runs = struct('label', {}, 'grid', {}, 'variant', {}, 'folder', {}, 'args', {}, 'n_jobs', {});
run_minutes = zeros(1, 0);   % the estimate of each run, in the order of runs
for g = 1:numel(grids)
    grid = grids(g);
    first_axis = mc_sweep_axis('first', field_or_empty(grid, 'first_axis'));
    second_axis = mc_sweep_axis('second', field_or_empty(grid, 'second_axis'));
    % An option under a default axis's own name on another axis is a config
    % mistake; mc_gain_sweep_grid refuses the rest (both names of one option).
    for A = [first_axis, second_axis]
        own = names_given(grid, DEFAULT_NAMED.(A.position));
        if ~A.is_default && ~isempty(own)
            bad_config('Grid ''%s'': %s belong(s) to the %s axis when it is %s; its %s axis is %s. Use the %s_* options.', ...
                grid.tag, strjoin(own, ', '), A.position, mc_sweep_axis(A.position).property, ...
                A.position, A.property, A.generic_option);
        end
    end
    % The grid's own run options: its model_args after the stage's (the run
    % applies them in order, so the grid's win), its readout in place of the
    % stage's. Neither is passed when nothing sets it, so a grid that sets
    % neither plans the arguments it always has.
    grid_model_args = field_or_empty(grid, 'model_args');
    if isempty(grid_model_args); grid_model_args = {}; end
    if ~(iscell(grid_model_args) && mod(numel(grid_model_args), 2) == 0)
        bad_config('Grid ''%s'': model_args must be a cell of name-value pairs.', grid.tag);
    end
    model_args = [stage_model_args, grid_model_args];
    % The estimate only: a grid with the release noise on costs more a job.
    release_noise = has_release_noise({first_axis.property, second_axis.property}, model_args);
    shape = given(grid, GRID_SHAPE);
    G = mc_gain_sweep_grid(cfg.run_mode, shape{:}, seeds{:}, 'release_noise', release_noise);
    readout = field_or_empty(grid, 'readout_signal');
    if isempty(readout); readout = stage_readout; end
    run_options = {};
    if ~isempty(model_args); run_options = [run_options, {'model_args', model_args}]; end %#ok<AGROW>
    if ~isempty(readout); run_options = [run_options, {'readout_signal', readout}]; end %#ok<AGROW>
    % The run takes each axis's levels under the axis's own option name
    % (gain_levels and bias_levels for the default axes), as it always has.
    on_axis = [shared, run_options, {'first_axis', first_axis.property, 'second_axis', second_axis.property, ...
        first_axis.levels_option, G.first_levels}];
    runs(end+1) = new_run(grid.tag, '', G, [{'preset_name', cfg.mc_preset}, on_axis, ...
        {second_axis.levels_option, G.second_levels, 'seeds', G.seeds}]); %#ok<AGROW>
    run_minutes(end+1) = G.worker_minutes; %#ok<AGROW>

    variants = field_or_empty(grid, 'variants');
    if isempty(variants); variants = struct('tag', {}, 'preset', {}); end
    if ~(isstruct(variants) && all(isfield(variants, {'tag', 'preset'})))
        bad_config('The variants of a grid must be a struct array with fields tag and preset.');
    end
    variant_seeds = field_or_empty(grid, 'variant_seeds');
    if isempty(variant_seeds); variant_seeds = G.seeds; end
    for k = 1:numel(variants)
        v = variants(k);
        if isempty(v.tag)
            bad_config('Grid ''%s'': every variant needs a tag.', grid.tag);
        end
        [~, own_value] = mc_gain_sweep_preset_point(v.preset, model_args, second_axis.property, ...
            first_axis.property);   % errors on an unknown preset
        V = mc_gain_sweep_grid(cfg.run_mode, 'first_axis', first_axis.property, ...
            'second_axis', second_axis.property, ...
            'first_levels', G.first_levels, 'second_levels', own_value, 'seeds', variant_seeds, ...
            'release_noise', release_noise);
        runs(end+1) = new_run(grid.tag, v.tag, V, [{'preset_name', v.preset}, on_axis, ...
            {second_axis.levels_option, V.second_levels, 'seeds', V.seeds}]); %#ok<AGROW>
        run_minutes(end+1) = V.worker_minutes; %#ok<AGROW>
    end
end
if numel(unique({runs.folder})) < numel(runs)
    bad_config(['Two runs of the stage would write to one folder (%s): grid tags must be ' ...
        'distinct, and so must the variant tags of a grid.'], strjoin({runs.folder}, ', '));
end
keep = ismember({runs.grid}, run_grids);
runs = runs(keep);
plan = struct('runs', runs, 'n_jobs', sum([runs.n_jobs]), 'worker_minutes', sum(run_minutes(keep)));
end

%% ------------------------------------------------------------------------
function run = new_run(grid_tag, variant_tag, G, args)
% One run of the plan. The unnamed grid's own run is labelled 'main'.
tags = {grid_tag, variant_tag};
label = strjoin(tags(~cellfun(@isempty, tags)), '_');
if isempty(label); label = 'main'; end
run = struct('label', label, 'grid', grid_tag, 'variant', variant_tag, ...
    'folder', mc_gain_sweep_folder(grid_tag, variant_tag), 'args', {args}, 'n_jobs', G.n_jobs);
end

function tf = has_release_noise(axis_properties, model_args)
% Whether a grid's runs have the release noise on, for the estimate: an axis
% is a release-noise amplitude, or the fixed arguments set one above zero
% (the last value given for a name is the one the run applies).
RELEASE_SIGMA = {'release_sigma', 'release_sigma_all'};
tf = any(ismember(axis_properties, RELEASE_SIGMA));
names = model_args(1:2:end);
for p = RELEASE_SIGMA
    k = find(strcmp(names, p{1}), 1, 'last');
    if ~isempty(k)
        v = model_args{2 * k};
        tf = tf || (isnumeric(v) && any(v(:) > 0));
    end
end
end

function args = given(s, names)
% Name-value pairs of the named fields that s sets (present and non-empty).
names = names_given(s, names);
args = cell(1, 2 * numel(names));
args(1:2:end) = names;
args(2:2:end) = cellfun(@(f) s.(f), names, 'UniformOutput', false);
end

function names = names_given(s, names)
% Those of the names that are fields of s with a non-empty value.
names = names(cellfun(@(f) isfield(s, f) && ~isempty(s.(f)), names));
end

function v = field_or_empty(s, name)
v = [];
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('mc_gain_sweep_stage_plan:BadConfig', varargin{:});
end

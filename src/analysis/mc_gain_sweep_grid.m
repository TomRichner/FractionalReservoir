function G = mc_gain_sweep_grid(run_mode, spec)
% MC_GAIN_SWEEP_GRID The levels and seeds of one block C sweep of memory capacity at a run mode: the stage's one size table.
%
%   G = MC_GAIN_SWEEP_GRID('fast')                       % gain x bias, the default grid
%   G = MC_GAIN_SWEEP_GRID('fast', 'second_axis', 'self_coupling', ...
%           'gain_spacing', 'linear', 'second_range', [-4 4])
%   G = MC_GAIN_SWEEP_GRID('fast', 'first_axis', 'self_coupling', 'first_range', [-4 4], ...
%           'second_axis', 'sigma_u_noise', 'second_range', [0 0.05])
%   G = MC_GAIN_SWEEP_GRID(run_mode, 'gain_levels', g, 'second_levels', b, 'seeds', s)
%
% A grid is two axes, each a scalar model property with a RANGE and a
% SPACING (mc_sweep_axis names them); the run mode says how many levels
% each axis gets and how many networks (seeds) run at every point.
% run_mc_gain_sweep takes its default grid from here, and a bundle's stage
% (mc_gain_sweep_stage_plan) resolves every grid it lists through here, so
% the sizes are written once.
%
% THE SIZE TABLE (Tom, 2026-11-13): levels per axis x networks, the same for
% every grid of a mode.
%   mode         levels    networks   one grid   the mu7chain bundle
%   essential     3 x  3       1         9 jobs     24 jobs, about 20 min
%   fast          7 x  7       3       147 jobs    336 jobs, about 4 h 25
%   medium(2)    13 x 13      10      1690 jobs   3640 jobs, about 48 h
%   production   17 x 17      10      2890 jobs   6120 jobs, about 81 h
% "The mu7chain bundle" is mu7chain_gainbias_config: two grids (gain x bias,
% gain x self-coupling) and two networks swept along gain only, so
% (2 n^2 + 2 n) x networks jobs for n levels. Times are for 12 pool workers
% at WORKER_MINUTES_PER_JOB = 9.5 worker-minutes per job of four conditions
% (measured on the first block C runs: 9 to 12), with every job of the stage
% in one parfor (run_mc_gain_sweep_set), so the time is jobs x 9.5 / 12.
% Agents run 'essential' and 'fast'; 'medium' and larger are Tom's to start.
% 'production' is the lead's choice of a size above medium (17 levels keeps
% the count odd, so zero stays a level of a symmetric axis).
%
% THE AXES (each named by a model property; mc_sweep_axis)
%   first axis    first_axis, default 'level_of_chaos' (the gain).
%                 first_range, and first_spacing 'log' or 'linear'. For the
%                 gain the range defaults to [0.1 4] and the spacing to 'log',
%                 and the options also go by gain_range, gain_spacing and
%                 gain_levels; another property has no default range (give
%                 one, or its levels) and is spaced linearly unless told.
%   second axis   second_axis, default 'S_c' (the bias). second_range,
%                 linear spacing. For the bias the range defaults to
%                 [-0.6 0.6], and its levels also go by bias_levels; another
%                 property has no default range.
% An option under an axis's default name (gain_*, bias_levels) is refused
% on another axis, and the two names of one option are not both given.
% WHY THESE RANGES (FD-010, decisions 14, 18 and 21): above bias 0.6 the
% network is silent; the no-adaptation peak lies at low bias, so the axis
% runs down to -0.6; and at low bias the adapting conditions are still
% stable, with their best memory capacity, at gain 2, so the gain axis runs
% to 4 to bracket their collapse.
% LEVELS ARE ROUNDED TO 10 DECIMALS, so a level is the number it reads as:
% 0.2 and not linspace's 0.19999999999999998, and the middle of a symmetric
% axis with an odd count is exactly 0. Job files are named by the exact
% value, and a figure asks for "the level 0.2".
% Seeds are mc_protocol's seed_net_base + 1, 2, ... (3001, 3002, ...).
%
% EXPLICIT LEVELS OR SEEDS (first_levels / gain_levels, second_levels /
% bias_levels, seeds) replace the table's for that axis, as given (not
% rounded).
%
% G: first_axis, second_axis (mc_sweep_axis descriptors), first_levels,
% second_levels, seeds (rows), n_jobs (their product), and worker_minutes,
% the estimate n_jobs x WORKER_MINUTES_PER_JOB.
%
% Errors carry the stage's identifiers: run_mc_gain_sweep:badMode (an
% unknown run mode), :NoFirstLevels / :NoSecondLevels (an axis with no
% default range given neither a range nor levels), :TwoLevelOptions /
% :TwoAxisOptions (both names of one
% option), and :<Default><Option>OnOtherAxis (GainLevelsOnOtherAxis,
% BiasLevelsOnOtherAxis, ...: an option under the default axis's name on
% another axis).
%
% See also: run_mc_gain_sweep, mc_gain_sweep_stage_plan, mc_sweep_axis,
%           mc_protocol, run_mode_names

arguments
    run_mode (1,:) char
    spec.first_axis    (1,:) char   = 'level_of_chaos'
    spec.second_axis   (1,:) char   = 'S_c'
    spec.first_range   double {mustBeEmptyOrRange} = []
    spec.gain_range    double {mustBeEmptyOrRange} = []
    spec.first_spacing (1,:) char   = ''
    spec.gain_spacing  (1,:) char   = ''
    spec.first_levels  (1,:) double = []
    spec.gain_levels   (1,:) double = []
    spec.second_range  double {mustBeEmptyOrRange} = []
    spec.second_levels (1,:) double = []
    spec.bias_levels   (1,:) double = []
    spec.seeds         (1,:) double = []
end

WORKER_MINUTES_PER_JOB = 9.5;      % one job of four conditions; see the header
DEFAULT_RANGE   = struct('first', [0.1 4], 'second', [-0.6 0.6]);   % of the default axes
DEFAULT_SPACING = struct('first', 'log', 'second', 'linear');       % of the default axes
LEVEL_DECIMALS  = 10;

[n_levels, n_seeds] = size_for_mode(run_mode);
axes = {mc_sweep_axis('first', spec.first_axis), mc_sweep_axis('second', spec.second_axis)};

G = struct('first_axis', axes{1}, 'second_axis', axes{2});
for a = 1:2
    A = axes{a};
    levels = axis_option(spec, A, 'levels');
    range = axis_option(spec, A, 'range');
    spacing = axis_option(spec, A, 'spacing');
    if isempty(spacing)
        spacing = 'linear';
        if A.is_default; spacing = DEFAULT_SPACING.(A.position); end
    end
    if isempty(levels)
        if isempty(range) && A.is_default; range = DEFAULT_RANGE.(A.position); end
        if isempty(range)
            error(['run_mc_gain_sweep:No' capitalized(A.position) 'Levels'], ...
                ['The %s axis ''%s'' needs its levels (%s_levels) or a range (%s_range): ' ...
                 'the size table has a default range for the gain and the bias only.'], ...
                A.position, A.property, A.generic_option, A.generic_option);
        end
        levels = round(spaced(range, n_levels, spacing, A), LEVEL_DECIMALS);
    end
    G.([A.position '_levels']) = levels;
end
G.seeds = spec.seeds;
if isempty(G.seeds)
    P = mc_protocol();
    G.seeds = P.seed_net_base + (1:n_seeds);
end
G.n_jobs = numel(G.first_levels) * numel(G.second_levels) * numel(G.seeds);
G.worker_minutes = G.n_jobs * WORKER_MINUTES_PER_JOB;
end

%% ------------------------------------------------------------------------
function v = axis_option(spec, A, what)
% The value of one option of axis A ('levels', 'range' or 'spacing'), under
% its generic name (first_levels) or the default axis's name (gain_levels).
% The default name is refused on another axis, and both names together are
% refused. Empty when neither is given (or the option does not exist for
% this axis: the second axis has no spacing).
generic = [A.generic_option '_' what];
default_axis = mc_sweep_axis(A.position);
own = [default_axis.option '_' what];
v = field_value(spec, generic);
v_own = field_value(spec, own);
if isempty(v_own); return; end
if ~A.is_default
    error(['run_mc_gain_sweep:' capitalized(default_axis.option) capitalized(what) 'OnOtherAxis'], ...
        '%s is the option of the %s axis when it is %s; this %s axis is %s. Use %s.', ...
        own, A.position, default_axis.property, A.position, A.property, generic);
end
if ~isempty(v)
    if strcmp(what, 'levels'); id = 'TwoLevelOptions'; else; id = 'TwoAxisOptions'; end
    error(['run_mc_gain_sweep:' id], 'Give %s or %s (its name for the default axis), not both.', ...
        generic, own);
end
v = v_own;
end

function v = field_value(s, name)
v = [];
if isfield(s, name); v = s.(name); end
end

function levels = spaced(range, n, spacing, A)
switch spacing
    case 'log'
        if range(1) <= 0
            error('mc_gain_sweep_grid:BadRange', ...
                'A log-spaced axis needs a positive range; the %s axis %s has [%g %g].', ...
                A.position, A.property, range);
        end
        levels = logspace(log10(range(1)), log10(range(2)), n);
    case 'linear'
        levels = linspace(range(1), range(2), n);
    otherwise
        error('mc_gain_sweep_grid:BadSpacing', ...
            'The spacing of the %s axis is ''log'' or ''linear''; got ''%s''.', A.position, spacing);
end
end

function s = capitalized(s)
s(1) = upper(s(1));
end

function [n_levels, n_seeds] = size_for_mode(run_mode)
% The size table: levels per axis and networks per grid point.
switch run_mode
    case 'essential'
        n_levels = 3;   n_seeds = 1;
    case 'fast'
        n_levels = 7;   n_seeds = 3;
    case {'medium', 'medium2'}
        n_levels = 13;  n_seeds = 10;
    case 'production'
        n_levels = 17;  n_seeds = 10;
    otherwise
        error('run_mc_gain_sweep:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
            run_mode, strjoin(run_mode_names(), ', '));
end
end

function mustBeEmptyOrRange(v)
% Empty, or [lo hi] with lo < hi.
if isempty(v); return; end
if ~(numel(v) == 2 && all(isfinite(v)) && v(1) < v(2))
    error('mc_gain_sweep_grid:BadRange', 'A range is [lo hi] with lo < hi.');
end
end

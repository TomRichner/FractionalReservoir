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
%   essential     3 x  3       1         9 jobs     51 jobs, about 50 min (58.9 measured)
%   fast          7 x  7       3       147 jobs    777 jobs, about 12 h 30
%   fastmedium    9 x  9       2       162 jobs    846 jobs, about 13 h 40
%   medium(2)    13 x 13      10      1690 jobs   8710 jobs, about 141 h
%   production   17 x 17      10      2890 jobs  14790 jobs, about 239 h
% "The mu7chain bundle" is mu7chain_gainbias_config: five grids (gain x
% bias, gain x self-coupling, and since FD-027 gain x release noise,
% self-coupling x release noise, gain x Wiener noise) and two networks swept
% along gain only, so (5 n^2 + 2 n) x networks jobs for n levels, of which
% the two release-noise grids' 2 n^2 x networks. Times are for 12 pool
% workers, with every job of the stage in one parfor (run_mc_gain_sweep_set),
% at WORKER_MINUTES_PER_JOB = 10.24 worker-minutes per job of four
% conditions (the fast run of 2026-11-13: 336 jobs without release noise in
% 286.7 min on 12 workers; the first block C runs gave 9 to 12) and
% RELEASE_NOISE_COST = 1.35 times that for a job of a grid with the release
% noise on (its condition-runs took about 157 s against 116 s, FD-027 N2;
% the calibration runs of 2026-11-14 170 to 200 s against 124 to 130 s).
% G.worker_minutes counts a job at that cost when the grid is told
% release_noise = true (mc_gain_sweep_stage_plan decides it from the grid's
% axes and fixed arguments); until 2026-11-14 it counted 9.5 for every job.
% 'fastmedium' (Tom, 2026-11-14) is sized for mu7chain_AC_run to finish in
% about 19 h: 9 levels keep the count odd (zero stays a level of the
% symmetric axes), and 3 networks at 9 x 9 would be about 20.5 h alone.
% Agents run 'essential' and 'fast'; 'fastmedium' and larger are Tom's to
% start. 'production' is the lead's choice of a size above medium (17
% levels keeps the count odd, so zero stays a level of a symmetric axis).
%
% THE AXES (each named by a model property; mc_sweep_axis)
%   first axis    first_axis, default 'level_of_chaos' (the gain).
%                 first_range, and first_spacing 'log' or 'linear'. For the
%                 gain the range defaults to [0.1 4] and the spacing to 'log',
%                 and the options also go by gain_range, gain_spacing and
%                 gain_levels; another property has no default range (give
%                 one, or its levels) and is spaced linearly unless told.
%   second axis   second_axis, default 'S_c' (the bias). second_range, and
%                 second_spacing (default 'linear'). For the bias the range
%                 defaults to [-0.6 0.6], and its levels also go by
%                 bias_levels; another property has no default range.
% THE SPACINGS, the same on either axis:
%   'linear'    n levels evenly spaced from lo to hi
%   'log'       n levels evenly spaced in log from lo to hi (lo > 0)
%   'zero_log'  a ZERO LEVEL, then log-spaced levels: the n-level log ladder
%               from lo to hi with its lowest rung, lo, replaced by 0. The
%               top level is hi, the zero stands one log step below the
%               lowest positive level (where the maps draw it:
%               mc_axis_spacing), and at 3 levels the axis is zero, the
%               middle of the ladder and the top. Written for a noise
%               amplitude (FD-027, P3), whose zero is "off" and whose effect
%               grows over decades: linear levels would be mostly floor, and
%               log levels would lack the noise-free reference. The lowest
%               positive level is lo (hi/lo)^(1/(n-1)), so it moves with the
%               mode's n; lo itself is never a level.
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
% THE COST. release_noise (default false) says the grid's runs have the
% release noise on; it changes the estimate only.
%
% G: first_axis, second_axis (mc_sweep_axis descriptors), first_levels,
% second_levels, seeds (rows), n_jobs (their product), and worker_minutes,
% the estimate n_jobs x WORKER_MINUTES_PER_JOB, times RELEASE_NOISE_COST
% with release_noise.
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
    spec.second_spacing (1,:) char  = ''
    spec.second_levels (1,:) double = []
    spec.bias_levels   (1,:) double = []
    spec.seeds         (1,:) double = []
    spec.release_noise (1,1) logical = false
end

WORKER_MINUTES_PER_JOB = 10.24;    % one job of four conditions; see the header
RELEASE_NOISE_COST     = 1.35;     % a job with the release noise on, relative to one without
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
minutes_per_job = WORKER_MINUTES_PER_JOB;
if spec.release_noise; minutes_per_job = minutes_per_job * RELEASE_NOISE_COST; end
G.worker_minutes = G.n_jobs * minutes_per_job;
end

%% ------------------------------------------------------------------------
function v = axis_option(spec, A, what)
% The value of one option of axis A ('levels', 'range' or 'spacing'), under
% its generic name (first_levels) or the default axis's name (gain_levels).
% The default name is refused on another axis, and both names together are
% refused. Empty when neither is given. (The second axis's spacing and
% range have no default-axis name: there is no bias_spacing or bias_range,
% so for it only second_spacing and second_range count.)
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
                A.position, A.property, range(1), range(2));
        end
        levels = logspace(log10(range(1)), log10(range(2)), n);
    case 'linear'
        levels = linspace(range(1), range(2), n);
    case 'zero_log'
        % The n-level log ladder from lo to hi, its lowest rung replaced by 0.
        if range(1) <= 0
            error('mc_gain_sweep_grid:BadRange', ...
                ['A zero_log axis needs a positive range (the ends of its log ladder); ' ...
                 'the %s axis %s has [%g %g].'], A.position, A.property, range(1), range(2));
        end
        if n < 2
            error('mc_gain_sweep_grid:BadSpacing', 'A zero_log axis needs at least two levels.');
        end
        levels = logspace(log10(range(1)), log10(range(2)), n);
        levels(1) = 0;
    otherwise
        error('mc_gain_sweep_grid:BadSpacing', ...
            'The spacing of the %s axis is ''linear'', ''log'' or ''zero_log''; got ''%s''.', ...
            A.position, spacing);
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
    case 'fastmedium'
        n_levels = 9;   n_seeds = 2;
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

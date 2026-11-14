function G = mc_gain_sweep_grid(run_mode, spec)
% MC_GAIN_SWEEP_GRID The levels and seeds of one gain sweep of memory capacity at a run mode: the stage's one size table.
%
%   G = MC_GAIN_SWEEP_GRID('fast')                       % gain x bias, the default grid
%   G = MC_GAIN_SWEEP_GRID('fast', 'second_axis', 'self_coupling', ...
%           'gain_spacing', 'linear', 'second_range', [-4 4])
%   G = MC_GAIN_SWEEP_GRID(run_mode, 'gain_levels', g, 'second_levels', b, 'seeds', s)
%
% A grid is a RANGE and a SPACING per axis; the run mode says how many
% levels each axis gets and how many networks (seeds) run at every point.
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
% THE AXES
%   gain          gain_range (default [0.1 4]), 'log' (default) or 'linear'
%                 spacing
%   second axis   second_range, linear spacing. For the bias (second_axis
%                 'S_c', the default) the range defaults to [-0.6 0.6];
%                 another axis has no default range and must give one, or
%                 its levels.
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
% EXPLICIT LEVELS OR SEEDS (gain_levels, second_levels, seeds) replace the
% table's for that axis, as given (not rounded).
%
% G: gain_levels, second_levels, seeds (rows), n_jobs (their product), and
% worker_minutes, the estimate n_jobs x WORKER_MINUTES_PER_JOB.
%
% Errors carry the stage's identifiers: run_mc_gain_sweep:badMode (an
% unknown run mode) and run_mc_gain_sweep:NoSecondLevels (an axis other
% than the bias with neither a range nor levels).
%
% See also: run_mc_gain_sweep, mc_gain_sweep_stage_plan, mc_sweep_second_axis,
%           mc_protocol, run_mode_names

arguments
    run_mode (1,:) char
    spec.second_axis   (1,:) char   = 'S_c'
    spec.gain_range    double {mustBeEmptyOrRange} = []
    spec.gain_spacing  (1,:) char   = ''
    spec.second_range  double {mustBeEmptyOrRange} = []
    spec.gain_levels   (1,:) double = []
    spec.second_levels (1,:) double = []
    spec.seeds         (1,:) double = []
end

WORKER_MINUTES_PER_JOB = 9.5;      % one job of four conditions; see the header
DEFAULT_GAIN_RANGE     = [0.1 4];
DEFAULT_GAIN_SPACING   = 'log';
DEFAULT_BIAS_RANGE     = [-0.6 0.6];
LEVEL_DECIMALS         = 10;

[n_levels, n_seeds] = size_for_mode(run_mode);
axis2 = mc_sweep_second_axis(spec.second_axis);

gain_range = spec.gain_range;
if isempty(gain_range); gain_range = DEFAULT_GAIN_RANGE; end
gain_spacing = spec.gain_spacing;
if isempty(gain_spacing); gain_spacing = DEFAULT_GAIN_SPACING; end
second_range = spec.second_range;
if isempty(second_range) && axis2.is_default; second_range = DEFAULT_BIAS_RANGE; end

G = struct();
G.gain_levels = spec.gain_levels;
if isempty(G.gain_levels)
    switch gain_spacing
        case 'log'
            if gain_range(1) <= 0
                error('mc_gain_sweep_grid:BadGainRange', ...
                    'A log-spaced gain axis needs a positive range; got [%g %g].', gain_range);
            end
            levels = logspace(log10(gain_range(1)), log10(gain_range(2)), n_levels);
        case 'linear'
            levels = linspace(gain_range(1), gain_range(2), n_levels);
        otherwise
            error('mc_gain_sweep_grid:BadGainSpacing', ...
                'gain_spacing is ''log'' or ''linear''; got ''%s''.', gain_spacing);
    end
    G.gain_levels = round(levels, LEVEL_DECIMALS);
end
G.second_levels = spec.second_levels;
if isempty(G.second_levels)
    if isempty(second_range)
        error('run_mc_gain_sweep:NoSecondLevels', ...
            ['The second axis ''%s'' needs its levels (second_levels) or a range ' ...
             '(second_range): the size table has a default range for the bias only.'], axis2.property);
    end
    G.second_levels = round(linspace(second_range(1), second_range(2), n_levels), LEVEL_DECIMALS);
end
G.seeds = spec.seeds;
if isempty(G.seeds)
    P = mc_protocol();
    G.seeds = P.seed_net_base + (1:n_seeds);
end
G.n_jobs = numel(G.gain_levels) * numel(G.second_levels) * numel(G.seeds);
G.worker_minutes = G.n_jobs * WORKER_MINUTES_PER_JOB;
end

%% ------------------------------------------------------------------------
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

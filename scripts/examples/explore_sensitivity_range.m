function [out_dir, summary] = explore_sensitivity_range(opts)
% EXPLORE_SENSITIVITY_RANGE One 1-D sweep, to choose the range for the real one.
%
%   explore_sensitivity_range()                       % the defaults below
%   explore_sensitivity_range('param', 'n', 'range', [50 2000])
%   explore_sensitivity_range('param', 'f_E', 'range', [0.1 0.9], ...
%                             'run_mode', 'fast2')
%   [d, s] = explore_sensitivity_range('run_mode', 'medium');
%
% WHAT THIS IS FOR. run_sensitivity_analysis hardcodes the seven axes the paper
% sweeps and the range of each -- level_of_chaos over [0.5 1.5], n over
% [100 1000], and so on. Those numbers have to come from somewhere, and picking
% them by eye from a production sweep is expensive and circular. This runs ONE
% axis over whatever range you name, on the same preset, conditions, timings and
% solver the real sweep would use, and prints a per-level table so the useful
% span is visible directly.
%
% It is deliberately a THIN wrapper. Everything that decides what the numbers
% mean -- the preset, the model class, the four adaptation conditions, the
% run_mode timings, the noise-driven integrator choice -- comes from
% resolve_run_context, exactly as run_sensitivity_analysis gets it. If this
% file and that one ever disagree about how a sweep is set up, the range you
% pick here will not be the range you get there, and the whole exercise is
% worthless. Two places are load-bearing for that and are marked KEEP IN SYNC.
%
% OUTPUT
%   out_dir   the PSA output directory (also holds psa_object.mat and figures)
%   summary   struct array, one entry per level, with .value and a field per
%             condition holding [median, min, max, frac_positive, n] of the metric
%
% Written into data/param_space/explore_sensitivity_<param>_<timestamp>/ under
% its own folder_prefix, so exploration runs never sit alongside the dated
% run_all_* directories the figures resolve against.
%
% COST is n_levels x n_reps x 4 conditions simulations, printed before it starts.
% 'medium' is 11 x 15 x 4 = 660. Use 'fast' or 'fast2' to scout first.
%
% See also: run_sensitivity_analysis, resolve_run_context, analysis_run_config,
%           srnn_param_preset, ParamSpaceAnalysis2/plot_sensitivity

arguments
    opts.param       (1,:) char    = 'level_of_chaos'
    opts.range       (1,2) double  = [0.25, 2.5]
    opts.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd'
    opts.run_mode    (1,:) char    = 'medium'
    opts.metric      (1,:) char    = 'LLE'      % the metric the table reports
    opts.n_levels    double        = []         % [] -> the run_mode's
    opts.n_reps      double        = []         % [] -> the run_mode's
    opts.hist_range  (1,2) double  = [-2, 2]    % LLE heatmap y-axis
    opts.save_figs   (1,1) logical = true
    opts.output_dir  (1,:) char    = ''
    opts.verbose     (1,1) logical = true
end

setup_paths();

%% Context: preset, class, conditions, timings -- all from one place
ctx = resolve_run_context('sensitivity', ...
    'preset_name', opts.preset_name, ...
    'run_mode',    opts.run_mode, ...
    'output_dir',  opts.output_dir, ...
    'save_figs',   opts.save_figs, ...
    'verbose',     opts.verbose);

% 'f' and 'f_E' name the same quantity on the two model classes -- a scalar
% property on SRNNModel2, a scalar alias onto a 1 x C row on SRNNCellTypePairs.
% Accept either spelling and use whichever the class actually has, so the same
% call works against either preset family.
param = opts.param;
if any(strcmp(param, {'f', 'f_E'})) && ~strcmp(param, ctx.f_param)
    fprintf('[explore] ''%s'' -> ''%s'' for %s\n', param, ctx.f_param, ctx.model_class);
    param = ctx.f_param;
end

n_levels = pick(opts.n_levels, ctx.n_levels);
n_reps   = pick(opts.n_reps,   ctx.n_reps);

% reps is a GRID AXIS, and add_grid_parameter needs at least two values, so
% n_reps = 1 fails deep inside it complaining about 'param_range' -- which names
% neither reps nor the real constraint. Say it here instead. n_levels goes the
% same way, and a single level would not be a sweep anyway.
if n_reps < 2
    error('explore_sensitivity_range:TooFewReps', ...
        ['n_reps must be at least 2: reps is a grid axis and a grid axis ' ...
        'needs two or more values. Got %d.'], n_reps);
end
if n_levels < 2
    error('explore_sensitivity_range:TooFewLevels', ...
        'n_levels must be at least 2 to sweep anything. Got %d.', n_levels);
end

param_range = sort(opts.range);   % add_grid_parameter rejects a descending range
n_cond = numel(ctx.conditions);

fprintf('\n========================================\n');
fprintf('=== Explore: %s over [%.4g, %.4g] ===\n', param, param_range(1), param_range(2));
fprintf('========================================\n');
fprintf('  preset      %s (%s)\n', ctx.preset_name, ctx.model_class);
fprintf('  run_mode    %s\n', ctx.run_mode);
% resolve_run_context has already printed the run_mode's OWN n_levels/n_reps
% above; say so when these differ, or the two lines look like a contradiction.
fprintf('  levels      %d%s   reps %d%s   conditions %d\n', ...
    n_levels, overridden(opts.n_levels, ctx.n_levels), ...
    n_reps,   overridden(opts.n_reps,   ctx.n_reps), n_cond);
fprintf('  simulations %d\n', n_levels * n_reps * n_cond);
fprintf('  conditions  %s\n', strjoin(cellfun(@(c) c.name, ctx.conditions, ...
    'UniformOutput', false), ', '));

%% The sweep
psa = ParamSpaceAnalysis2( ...
    'n_levels', n_levels, ...
    'batch_size', 25, ...
    'note', param, ...
    'randomize_order', false, ...   % ordered axes matter for a sensitivity sweep
    'verbose', ctx.verbose);
psa.folder_prefix  = 'explore_sensitivity';
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = ctx.model_defaults;
if ~isempty(ctx.output_dir)
    psa.output_dir = ctx.output_dir;
end

psa.set_conditions(ctx.conditions);

% KEEP IN SYNC with run_sensitivity_analysis. On SRNNModel2 the real sweep
% upgrades the STD conditions to TWO recovery timescales; a range chosen against
% single-timescale depression would not transfer. n_b_E must come from the
% condition (PSA ignores it in model_defaults) while tau_b_E_rec flows through
% model_defaults. SRNNCellTypePairs says the same thing per route inside its
% synapse_config, so there is nothing to override -- and tau_b_E_rec is not one
% of its properties, so setting it would be a hard validate_model_defaults error.
if strcmp(ctx.model_class, 'SRNNModel2')
    for ci = 1:numel(psa.conditions)
        if psa.conditions{ci}.n_b_E > 0
            psa.conditions{ci}.n_b_E = 2;
        end
    end
    psa.model_defaults.tau_b_E_rec = [0.5, 5];
end

psa.add_grid_parameter(param, param_range);
psa.add_grid_parameter('reps', 1:n_reps);

psa.run();

copyfile([mfilename('fullpath') '.m'], psa.output_dir);
out_dir = psa.output_dir;

%% Per-level table -- the actual point of the script
summary = level_summary(psa, param, ctx.conditions, opts.metric);
print_summary(summary, param, opts.metric, ctx.conditions);

%% Figures
% KEEP IN SYNC with run_sensitivity_analysis: same two sheets, same hist_range,
% so an exploration sheet is directly comparable with a production one.
%
% Only THIS call's figures are saved. save_some_figs_to_folder_2 given an empty
% handle list saves every figure currently open, and plot_sensitivity returns no
% handles to pass instead -- so a second call in the same MATLAB session writes
% the first call's sheets into the second call's folder. Observed directly: a
% [0.25 5] run wrote five sheets, three of them left over from a [0.25 2.5] run
% minutes earlier. Diffing the open-figure list around the plotting calls is
% what isolates them, and it also leaves any figures the user already had open
% untouched -- which passing 'all' cannot do and closing 'all' would destroy.
before = findobj(0, 'Type', 'figure');
psa.plot_sensitivity('metric', 'LLE', 'hist_range', opts.hist_range);
psa.plot_sensitivity('metric', 'mean_rate');
mine = setdiff(findobj(0, 'Type', 'figure'), before);

if ctx.save_figs
    fig_dir = fullfile(out_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, sprintf('explore_%s', param), ...
        [mine.Number], {'fig', 'png'});
    fprintf('Figures saved to %s (%d sheet(s))\n', fig_dir, numel(mine));
end

fprintf('\n=== Explore complete ===\n');
fprintf('  %s\n', out_dir);
end

%% ------------------------------------------------------------------------
function v = pick(a, b)
if isempty(a); v = b; else; v = a; end
end

function s = overridden(given, from_run_mode)
% Flags a value the caller overrode, naming what the run_mode would have used.
if isempty(given) || isequal(given, from_run_mode)
    s = '';
else
    s = sprintf(' (run_mode: %d)', from_run_mode);
end
end

function s = level_summary(psa, param, conditions, metric)
% One entry per swept level, with per-condition statistics of the metric.
%
% Uses ParamSpaceAnalysis2.collect_level_values, which is public precisely so a
% caller can walk the reps axis without reimplementing the sub2ind indexing --
% it also drops failed and NaN jobs, which a naive reshape would not.
levels = psa.param_vectors{strcmp(psa.grid_params, param)};
s = struct('value', {}, 'stats', {});
for k = 1:numel(levels)
    e = struct('value', levels(k), 'stats', struct());
    for ci = 1:numel(conditions)
        name = conditions{ci}.name;
        v = ParamSpaceAnalysis2.collect_level_values(psa, param, k, name, metric);
        if isempty(v)
            e.stats.(name) = struct('median', NaN, 'min', NaN, 'max', NaN, ...
                                    'frac_positive', NaN, 'n', 0);
        else
            e.stats.(name) = struct('median', median(v), 'min', min(v), ...
                'max', max(v), 'frac_positive', mean(v > 0), 'n', numel(v));
        end
    end
    s(end+1) = e; %#ok<AGROW>
end
end

function print_summary(s, param, metric, conditions)
% Median metric per level per condition, plus the fraction of reps above zero.
%
% frac>0 is the column to read when choosing a range for LLE: it is the share of
% repetitions that came out chaotic, so the edge of chaos is where it crosses
% 0.5, and a level whose whole column reads 0 or 1 is telling you the range
% extends past anything informative.
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);

fprintf('\n=== %s by %s (median over reps) ===\n', metric, param);
fprintf('%14s', param);
for ci = 1:numel(names); fprintf(' | %20s', names{ci}); end
fprintf('\n%14s', '');
for ci = 1:numel(names); fprintf(' | %10s %9s', 'median', 'frac>0'); end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 14 + numel(names)*23));

for k = 1:numel(s)
    fprintf('%14.4g', s(k).value);
    for ci = 1:numel(names)
        st = s(k).stats.(names{ci});
        fprintf(' | %10.4f %9.2f', st.median, st.frac_positive);
    end
    fprintf('\n');
end

if strcmp(metric, 'LLE')
    fprintf('\nEdge of chaos (median %s crosses zero):\n', metric);
    for ci = 1:numel(names)
        med = arrayfun(@(e) e.stats.(names{ci}).median, s);
        val = [s.value];
        fprintf('  %-20s %s\n', names{ci}, crossing_note(val, med));
    end
end
end

function txt = crossing_note(val, med)
% Where the median first goes from <= 0 to > 0, linearly interpolated. Reported
% as an interval as well as a point, because with a handful of levels the
% interpolation is a guide to where to look, not a measurement.
ok = ~isnan(med);
val = val(ok); med = med(ok);
if numel(val) < 2;            txt = 'too few levels';                     return; end
if all(med > 0);              txt = 'positive throughout - lower the range'; return; end
if all(med <= 0);             txt = 'never positive - raise the range';   return; end
i = find(med(1:end-1) <= 0 & med(2:end) > 0, 1);
if isempty(i)
    txt = 'non-monotonic - read the table';
    return
end
x = val(i) + (val(i+1) - val(i)) * (-med(i)) / (med(i+1) - med(i));
txt = sprintf('~%.4g   (between %.4g and %.4g)', x, val(i), val(i+1));
end

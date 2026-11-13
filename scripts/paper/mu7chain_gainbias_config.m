function cfg = mu7chain_gainbias_config(run_mode, figure_set)
% MU7CHAIN_GAINBIAS_CONFIG Block C's bundle: memory capacity over gain and bias on the mu7chain network (FD-010).
%
%   cfg = mu7chain_gainbias_config()                       % 'essential', figure set '4cond'
%   cfg = mu7chain_gainbias_config('essential', '3cond')   % the paper's three conditions
%   cfg = mu7chain_gainbias_config('fast', '4cond')        % with product MTS added
%
% WHAT IT RUNS. One stage, the opt-in 'mc_gain_sweep' of
% run_all_paper_analyses (paper_opt_in_stages): run_mc_gain_sweep on the
% bundle's noise-free memory-capacity twin, over the gain x bias grid of the
% run mode, into <run_dir>/mc_gain_sweep; then the two E:I variants of the
% same twin (excitation- and inhibition-dominated mean weights, 3 : 2 and
% 2 : 3), each along gain only at its own bias, into
% <run_dir>/mc_gain_sweep_EI3to2 and _EI2to3. Nothing else: the sweeps stage
% only creates the run directory and its manifest, which record the other
% stages as skipped.
%
% HOW LONG, on 12 pool workers (run_mc_gain_sweep's SIZES):
%   essential   about 30 minutes: the main grid (12 jobs) and the two
%               variants (4 jobs each) run one after another, each taking
%               about the 10 minutes of one job. The size agents run, and
%               the default.
%   fast        about 3 h 20 (154 + 2 x 22 jobs). A run for Tom to start,
%               like 'medium' and larger.
%
% DERIVED FROM mu7chain_config, so nothing about the network is restated:
% the preset, the memory-capacity twin (cfg.mc_preset), its integrator
% (cfg.mc_ode_solver), the conditions drawn per set and the sets drawn per
% mode (cfg.figure_sets) are mu7chain_config's. What differs:
%   cfg.stages        {'mc_gain_sweep'}; cfg.sweeps and eig examples empty
%   cfg.mc_gain_sweep the stage's settings: the two variants. Add
%                     gain_levels, bias_levels, seeds or conditions to
%                     depart from the stage's table for the mode
%                     (mc_gain_sweep_stage_plan lists every field)
%   cfg.run_dir       data/mu7chain_gainbias_<run_mode>
%   cfg.fig_root      figs/mu7chain_gainbias_<run_mode>_<figure_set>
%   cfg.figures       the block C figures only: fig_mc_gain_bias_maps,
%                     fig_mc_gain_curves (the main network, then a row per
%                     variant), fig_mc_robustness, and the statistics
%                     tables (mc_gain_bias_statistics) written beside them
%
% See also: mu7chain_gainbias_run, mu7chain_config, run_mc_gain_sweep,
%           mc_gain_sweep_stage_plan, run_all_paper_analyses

arguments
    run_mode (1,:) char = 'essential'
    figure_set (1,:) char {mustBeMember(figure_set, {'3cond', '4cond'})} = '4cond'
end
if ~ismember(run_mode, run_mode_names())
    error('mu7chain_gainbias_config:BadRunMode', 'run_mode ''%s'' is not one of: %s.', ...
        run_mode, strjoin(run_mode_names(), ', '));
end

cfg = mu7chain_config(run_mode, figure_set);

%% Which stage runs: the opt-in block C stage alone
cfg.stages = {'mc_gain_sweep'};
cfg.sweeps = {};
cfg.eig_heatmap_examples = [];
% The variant networks. A tag names the run's folder and the figure's row;
% its preset is the memory-capacity twin's variant of that name.
variant_tags = {'EI3to2', 'EI2to3'};
cfg.mc_gain_sweep = struct('variants', ...
    struct('tag', variant_tags, 'preset', strcat(cfg.mc_preset, '_', variant_tags)));

%% Where things land -- fixed, and apart from mu7chain_config's
cfg.run_dir  = ['data/mu7chain_gainbias_' run_mode];
cfg.fig_root = ['figs/mu7chain_gainbias_' run_mode '_' figure_set];

%% The figure registry: block C only
% Each entry gets run_dir, out_dir, save, visible, verbose and the set's
% conditions from make_all_paper_figures, then its own arguments.
cfg.figures = { ...
    registry_entry('fig_mc_gain_bias_maps',   @fig_mc_gain_bias_maps,   true,  {}), ...
    registry_entry('fig_mc_gain_curves',      @fig_mc_gain_curves,      true,  {'variants', variant_tags}), ...
    registry_entry('fig_mc_robustness',       @fig_mc_robustness,       true,  {}), ...
    registry_entry('mc_gain_bias_statistics', @statistics_entry,        false, {})};
end

%% ------------------------------------------------------------------------
function entry = registry_entry(name, fn, in_paper, args)
entry = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end

function out = statistics_entry(opts)
% mc_gain_bias_statistics behind the figure contract: the run directory's
% main result, the set's conditions, the tables written to out_dir.
% out.files lists them; there is no figure.
arguments
    opts.run_dir     (1,:) char    = ''
    opts.data_file   (1,:) char    = ''
    opts.out_dir     (1,:) char    = ''
    opts.save        (1,1) logical = true
    opts.visible     (1,1) logical = false     % accepted for a uniform call
    opts.verbose     (1,:) char    = 'minimal'
    opts.conditions                = {}
    opts.preset_name (1,:) char    = ''        % accepted for a uniform call
    opts.reference   (1,:) char    = 'sfa3_stdc2'
end
D = resolve_mc_gain_sweep('data_file', opts.data_file, 'run_dir', opts.run_dir);
out_dir = '';
if opts.save; out_dir = opts.out_dir; end
[~, files] = mc_gain_bias_statistics('data', D, 'conditions', opts.conditions, ...
    'reference', opts.reference, 'out_dir', out_dir, 'verbose', opts.verbose);
out = struct('figs', gobjects(0), 'files', {files}, 'source', D.source);
end

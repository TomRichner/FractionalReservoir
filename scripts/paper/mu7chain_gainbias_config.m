function cfg = mu7chain_gainbias_config(run_mode, figure_set)
% MU7CHAIN_GAINBIAS_CONFIG Block C's bundle: memory capacity over gain and bias on the mu7chain network (FD-010).
%
%   cfg = mu7chain_gainbias_config()                  % run_mode 'fast', figure set '4cond'
%   cfg = mu7chain_gainbias_config('fast', '3cond')   % the paper's three conditions
%   cfg = mu7chain_gainbias_config('fast', '4cond')   % with product MTS added
%   cfg = mu7chain_gainbias_config('essential')       % a first look, one seed
%
% WHAT IT RUNS. One stage, the opt-in 'mc_gain_sweep' of
% run_all_paper_analyses (paper_stage_names, third list): run_mc_gain_sweep
% on the bundle's noise-free memory-capacity twin, over the gain x bias grid
% of the run mode (its SIZES table), into <run_dir>/mc_gain_sweep; then the
% two E:I gain-only runs (Plan C runs C3b and C3c: the _EI3to2 and _EI2to3
% variants of the same twin, at their own bias S_c = 0.2) into
% <run_dir>/mc_gain_sweep_EI3to2 and _EI2to3. Nothing else: the sweeps
% stage only creates the run directory and its manifest, which record the
% other stages as skipped.
%
% DERIVED FROM mu7chain_config, so nothing about the network is restated:
% the preset, the memory-capacity twin (cfg.mc_preset), its integrator
% (cfg.mc_ode_solver, sra1), the conditions drawn per set and the sets drawn
% per mode (cfg.figure_sets: both, except 'essential', four conditions only)
% are mu7chain_config's. What differs:
%   cfg.stages        {'mc_gain_sweep'}; cfg.sweeps and eig examples empty
%   cfg.mc_gain_sweep struct: ei_presets (the two E:I twins). Add
%                     gain_levels, bias_levels, seeds, conditions to depart
%                     from the stage's table for the mode (run_all_paper_
%                     analyses' header lists every field)
%   cfg.run_dir       data/mu7chain_gainbias_<run_mode>
%   cfg.fig_root      figs/mu7chain_gainbias_<run_mode>_<figure_set>
%   cfg.figures       the block C figures only: fig_mc_gain_bias_maps,
%                     fig_mc_gain_curves (main network, then the E:I rows
%                     when their runs are in the run directory),
%                     fig_mc_robustness, and the statistics tables
%                     (mc_gain_bias_statistics) written beside them
%
% Agents run 'essential' and 'fast' only; 'medium' and larger are Tom's.
%
% See also: mu7chain_gainbias_run, mu7chain_config, run_mc_gain_sweep,
%           run_all_paper_analyses, make_block_C_figures

arguments
    run_mode (1,:) char = 'fast'
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
cfg.mc_gain_sweep = struct('ei_presets', {{[cfg.mc_preset '_EI3to2'], [cfg.mc_preset '_EI2to3']}});

%% Where things land -- fixed, and apart from mu7chain_config's
cfg.run_dir  = ['data/mu7chain_gainbias_' run_mode];
cfg.fig_root = ['figs/mu7chain_gainbias_' run_mode '_' figure_set];

%% The figure registry: block C only
% Each entry gets run_dir, out_dir, save, visible, verbose and the set's
% conditions from make_all_paper_figures.
F = {};
F = add(F, 'fig_mc_gain_bias_maps',   @fig_mc_gain_bias_maps,   true,  {});
F = add(F, 'fig_mc_gain_curves',      @fig_mc_gain_curves,      true,  {});
F = add(F, 'fig_mc_robustness',       @fig_mc_robustness,       true,  {});
F = add(F, 'mc_gain_bias_statistics', @statistics_entry,        false, {});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end

function out = statistics_entry(opts)
% mc_gain_bias_statistics behind the figure contract: the result of
% <run_dir>/mc_gain_sweep (resolve_mc_gain_sweep: that run directory or an
% error), the set's conditions, the tables written to out_dir. out.files
% lists them; there is no figure.
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
D = resolve_mc_gain_sweep(opts.data_file, opts.run_dir);
out_dir = '';
if opts.save; out_dir = opts.out_dir; end
[~, files] = mc_gain_bias_statistics('data', D, 'conditions', opts.conditions, ...
    'reference', opts.reference, 'out_dir', out_dir, 'verbose', opts.verbose);
out = struct('figs', gobjects(0), 'files', {files}, 'source', D.source);
end

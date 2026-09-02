function cfg = paper_config(opts)
% PAPER_CONFIG The single place that names what the paper is built from.
%
%   cfg = PAPER_CONFIG()
%   cfg = PAPER_CONFIG('run_mode', 'fast')
%   cfg = PAPER_CONFIG('preset_name', 'celltype_pairs_uniform_std_n500')
%
% THIS IS THE ONE FILE TO EDIT. Change the preset here and both master scripts
% follow: run_all_paper_analyses sweeps that network, make_all_paper_figures
% resolves the matching run and rebuilds every figure from it.
%
% TWO ORTHOGONAL KNOBS, as everywhere else in this project:
%   preset_name  WHICH NETWORK   (srnn_param_preset)
%   run_mode     HOW MUCH COMPUTE (analysis_run_config)
%
% and a third that says WHERE THE OUTPUT LANDS: fig_root. See its comment below.
%
% FIGURE PRESET OVERRIDES. Five figures are DELIBERATELY a different network
% from the paper's operating point, and each names its own preset below. They
% are not oversights: two single-neuron mechanism cartoons, a Sompolinsky
% reproduction, a small bursting network and the memory-capacity reservoir are
% all making points that the 500-neuron production network cannot make. Every
% other figure inherits cfg.preset_name.
%
% A figure that needs a DIFFERENT NETWORK needs a DIFFERENT PRESET, not
% cfg.preset_name plus overrides at the call site. fig_adaptation_methods was
% handed cfg.preset_name for months and so silently built the whole 500-neuron
% recurrent network for a figure captioned "one unconnected neuron".
%
% See also: run_all_paper_analyses, make_all_paper_figures, srnn_param_preset

arguments
    opts.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    % 'medium' is the default because it is what gets run: ~3 h of compute and
    % figures that are readable. 'production' is a deliberate act -- pass it
    % explicitly, paper_config('run_mode', 'production'), for the final run.
    opts.run_mode    (1,:) char = 'medium'
    opts.run_dir     (1,:) char = ''      % '' -> resolve by preset at figure time
    % WHERE THE OUTPUT GOES. Two modes, and the empty case is not an oversight:
    %
    %   set (default)   a STABLE directory, overwritten every run. The
    %                   manuscript can cite a path that never moves, and there
    %                   is never a question of which set is current.
    %   empty           make_all_paper_figures auto-names figs/figures_<dt>/,
    %                   so a one-off cannot clobber the paper's set.
    %
    % Same convention as run_all_analyses' output_dir, where empty means "make
    % your own dated folder". A relative path resolves against the project root.
    %
    % This is the knob that lets a SECOND config file produce a separate set --
    % point it at figs/something_else and the two never collide.
    %
    % figs/ is gitignored on purpose. Figures are regenerable from run_dir plus
    % a commit, and both are recorded in the manifest written at the root; the
    % final set is force-added at submission.
    opts.fig_root    (1,:) char = 'figs/paper'
    % WHERE THE COMPUTE LANDS -- the data-side twin of fig_root, with the same
    % two modes and the same empty-means-dated convention:
    %
    %   empty (default)  run_all_analyses creates data/param_space/run_all_<dt>/,
    %                    so successive runs never overwrite each other. This is
    %                    what the paper wants: a run is a dated artefact.
    %   set              a STABLE directory, reused and overwritten. For a
    %                    smoke test or a side experiment, where a pile of dated
    %                    folders is noise rather than history.
    %
    % A relative path resolves against the project root. Note the asymmetry with
    % run_dir above, which is the READ side (which finished run the figures are
    % built from); this is the WRITE side.
    opts.output_dir  (1,:) char = ''
end

cfg = struct();
cfg.preset_name = opts.preset_name;
cfg.run_mode    = opts.run_mode;
cfg.run_dir     = opts.run_dir;
cfg.fig_root    = opts.fig_root;
cfg.output_dir  = opts.output_dir;

% Presets for the figures that are deliberately different networks.
cfg.mc_preset          = 'mc_esn';              % SRNN_ESN_reservoir; see below
cfg.bursting_preset    = 'bursting_pairs';
cfg.sompolinsky_preset = 'sompolinsky_pairs';
cfg.stf_preset         = 'single_neuron_stf';
% The single-neuron methods cartoon. It carries the PAPER'S c, SFA timescales
% and dual depression on ONE unconnected neuron, which is why it is a separate
% preset and not cfg.preset_name: passing cfg.preset_name here is exactly the
% bug that made this figure plot neuron 1 of the 500-neuron recurrent network
% while claiming to show one unconnected neuron.
cfg.single_neuron_preset = 'single_neuron_dualStd';

% The gains shared by the two halves of figure 1 panel A. They MUST match, so
% they live here rather than as a literal in each figure.
cfg.panelA_gammas = [0.9, 1.6, 2.5];

%% The figure registry
% Order is the order make_all_paper_figures runs them, and is roughly the order
% they appear in the manuscript. `in_paper` marks the ones currently used; the
% rest are kept working and regenerated, per TR, but are not in the manuscript.
%
% Every entry is  {name, handle, extra-args}  where the extra args are appended
% to the standard cfg fields (run_dir, out_dir, save, visible).
F = {};
F = add(F, 'fig_introductory_concepts',      @fig_introductory_concepts,      true, ...
        {'preset_name', cfg.sompolinsky_preset, 'gammas', cfg.panelA_gammas});
F = add(F, 'fig_energy_landscape',           @fig_energy_landscape,           false, ...
        {'gammas', cfg.panelA_gammas});
F = add(F, 'fig_example_timeseries',         @fig_example_timeseries,         true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_FI_curve',                   @fig_FI_curve,                   true, {});
F = add(F, 'fig_adaptation_methods',         @fig_adaptation_methods,         true, ...
        {'variant', 'sfa_std', 'preset_name', cfg.single_neuron_preset});
F = add(F, 'fig_adaptation_methods_stf',     @fig_adaptation_methods,         false, ...
        {'variant', 'sfa_std_stf', 'preset_name', cfg.stf_preset});
F = add(F, 'fig_SFA_steady_state',           @fig_SFA_steady_state,           false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_STD_steady_state',           @fig_STD_steady_state,           false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_stim_engages_adaptation',    @fig_stim_engages_adaptation,    true, ...
        {'preset_name', cfg.bursting_preset});
F = add(F, 'fig_sensitivity_analysis_allStd', @fig_sensitivity_analysis_allStd, true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sensitivity_medians',        @fig_sensitivity_medians,        false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_param_space_allStd',         @fig_param_space_allStd,         true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_EI_param_space',             @fig_EI_param_space,             true, ...
        {'preset_name', cfg.preset_name});
% Same sheet, coloured by the balance of the WEIGHTS instead of the balance of
% the neuron counts. Its sibling above answers "does f_E explain where a network
% lands", which stopped being the whole question once the grid started sweeping
% the four mu blocks -- an 80% excitatory network with inhibition three times as
% strong is inhibition-dominated, and f_E cannot see that. Kept immediately after
% its sibling so the pair is obvious.
F = add(F, 'fig_EI_weights_param_space',     @fig_EI_weights_param_space,     true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sfa_EOC_allStd',             @fig_sfa_EOC_allStd,             true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_memory_capacity',            @fig_memory_capacity,            true, {});
F = add(F, 'fig_memory_capacity_example',    @fig_memory_capacity_example,    true, {});
F = add(F, 'fig_eig_heatmap',                @fig_eig_heatmap,                false, {});
% The generated equation and conditions tables. An ordinary entry, not a special
% call after the loop: as a special case its failures sat outside the headline
% count, which is how the n_a refactor broke it while the run still reported
% 17/17. Lands in <fig_root>/doc_tables/ by the same name-keying rule as the
% figures.
F = add(F, 'doc_tables',                     @fig_doc_tables,                 true, ...
        {'preset_name', cfg.preset_name});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end

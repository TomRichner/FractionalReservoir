function cfg = fast4_config(opts)
% FAST4_CONFIG A fast smoke-test config: 4 conditions, into figs/ and data/fast_4.
%
%   run_dir = run_all_paper_analyses(fast4_config());
%   results = make_all_paper_figures(fast4_config('run_dir', run_dir));
%
% Same pipeline as paper_config, three knobs turned:
%
%   preset_name  the FOUR-condition dual-STD network, not the seven-condition
%                one the paper uses. Same operating point (S_c = 0.2,
%                sigma_u_noise = 0.025, dual STD); fewer adaptation regimes, so
%                every sweep does roughly half the work.
%   run_mode     'fast' -- the smoke test. Minutes, not hours: coarse grids,
%                short trajectories, rk4 or sra1 rather than ode45.
%   fig_root     figs/fast_4   } both FIXED, so a smoke test leaves one
%   output_dir   data/fast_4   } directory rather than a pile of dated ones --
%                                the opposite choice from the paper, where a run
%                                is a dated artefact.
%
% THE TWO ROOTS BEHAVE DIFFERENTLY ON A RERUN, which is not obvious and is worth
% knowing before you run this twice:
%
%   figs/fast_4   overwritten. save_figure_stable deletes <out_dir>/<tag>* before
%                 saving, so each entry's directory holds only the current set.
%   data/fast_4   ACCUMULATES. Every sub-analysis appends its own timestamped
%                 folder, so a second run adds a parallel set beside the first,
%                 while run_manifest and provenance at the top level are
%                 replaced. run_all_paper_analyses warns when it finds the
%                 directory occupied; delete it first for a clean run.
%
% WHY A SEPARATE FILE rather than paper_config('run_mode', 'fast', ...). Because
% the point is that it cannot collide with the paper's output. Both roots are
% named here, so no invocation of this config can write into figs/paper or a
% dated run directory, and no forgotten name-value pair can turn a smoke test
% into something that clobbers the real set.
%
% THE FIGURE REGISTRY IS NOT DUPLICATED. This delegates to paper_config, so the
% list of figures, the per-figure preset overrides and panelA_gammas all stay
% defined in exactly one place. Add a figure there and it appears here too.
%
% NOTE ON THE PRESET. Four conditions means the condition-keyed plotters get
% four columns instead of seven; that is the preset's business, not this file's.
% If a figure hardcodes a seven-condition layout it will show up here as a
% squashed or empty panel -- which is a reason to run this, not a reason not to.
%
% See also: paper_config, run_all_paper_analyses, make_all_paper_figures,
%           srnn_param_preset, analysis_run_config

arguments
    opts.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond'
    opts.run_mode    (1,:) char = 'fast'
    opts.run_dir     (1,:) char = ''            % '' -> resolve by preset at figure time
    opts.fig_root    (1,:) char = 'figs/fast_4'
    opts.output_dir  (1,:) char = 'data/fast_4'
end

cfg = paper_config( ...
    'preset_name', opts.preset_name, ...
    'run_mode',    opts.run_mode, ...
    'run_dir',     opts.run_dir, ...
    'fig_root',    opts.fig_root, ...
    'output_dir',  opts.output_dir);
end

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
%   run_dir      data/fast_4   } both FIXED, so a smoke test leaves one pair of
%   fig_root     figs/fast_4   } directories rather than a pile of dated ones --
%                                the opposite choice from the paper, where a run
%                                is a dated artefact.
%
% run_dir is one field used by both entry points: run_all_paper_analyses writes
% there, make_all_paper_figures reads there. So the two halves need no argument
% passed between them --
%
%     run_all_paper_analyses(fast4_config());
%     make_all_paper_figures(fast4_config());
%
% THE TWO DIRECTORIES BEHAVE DIFFERENTLY ON A RERUN, which is worth knowing
% before you run this twice:
%
%   figs/fast_4   overwritten. save_figure_stable deletes <out_dir>/<tag>*
%                 before saving, so each entry holds only the current set.
%   data/fast_4   must be ABSENT OR EMPTY. run_all_paper_analyses errors
%                 otherwise, because a reused run directory accumulates sweep
%                 folders rather than replacing them. Delete or move it between
%                 runs. Figures are cheap to regenerate; a run is not, and is
%                 not worth silently mixing with another.
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
    opts.run_dir     (1,:) char = 'data/fast_4_test'
    opts.fig_root    (1,:) char = 'figs/fast_4_test'
    % Off, like paper_config's own default: figures are built and saved but do
    % not pop up and steal focus. Especially wanted here -- a smoke test is the
    % run you are most likely to start and then carry on working through.
    opts.visible_figures (1,1) logical = false
end

cfg = paper_config( ...
    'preset_name', opts.preset_name, ...
    'run_mode',    opts.run_mode, ...
    'run_dir',     opts.run_dir, ...
    'fig_root',    opts.fig_root, ...
    'visible_figures', opts.visible_figures);
end

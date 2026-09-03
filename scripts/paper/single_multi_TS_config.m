function cfg = single_multi_TS_config(opts)
% SINGLE_MULTI_TS_CONFIG Single- vs multi-timescale adaptation: 3 conditions.
%
%   run_dir = run_all_paper_analyses(single_multi_TS_config());
%   results = make_all_paper_figures(single_multi_TS_config());
%
% THE COMPARISON THIS IS FOR (TR, 2026-09-03, after discussing the manuscript
% with his PI): no adaptation vs adaptation on ONE timescale vs adaptation on
% MANY, with SFA and STD always present together.
%
%   no_adaptation   -                      -
%   sfa1_std1       1 tau_a (0.25 s)       1 STD pair (2 / 0.25 s)
%   sfa3_std2       3 tau_a (0.25-10 s)    2 STD pairs
%
% The 7-condition set separates the mechanisms -- sfa_only, std_only and their
% one-timescale variants -- which is a more complicated story than this paper
% wants and is headed for a later one. Here the axis is TIMESCALE COUNT and only
% that.
%
% THIS MAY BECOME THE PAPER'S CONFIG. It is not yet: paper_config still names
% the 7-condition preset, and nothing here changes that. Same pipeline, same
% figure registry, three knobs turned:
%
%   preset_name  celltype_pairs_Sc0p2_noise0p025_dualStd_3cond -- the SAME
%                network as the 4- and 7-condition presets, differing only in
%                which regimes are swept, so runs stay comparable. Note this
%                set names the full-adaptation regime sfa3_std2 where the others
%                call it sfa_and_std; identical physics, different label, so
%                each set can title its own comparison.
%   run_mode     'fast' FOR NOW. This is the knob to change when the comparison
%                is worth real compute -- 'medium' for readable numbers,
%                'production' for the paper. Editing this line is the whole
%                change; nothing else needs to move.
%   run_dir      data/single_multi_TS  } both FIXED, so this cannot touch the
%   fig_root     figs/single_multi_TS  } paper's output whatever run mode it is
%                                        pointed at.
%
% COST: three conditions rather than seven is 57% less compute at every sweep
% stage, so a 'medium' run here is cheaper than a 'fast' seven-condition one at
% the grid stages.
%
% run_dir is one field used by BOTH entry points -- the analyses write there,
% the figures read there -- so the two halves need no argument passed between
% them. See paper_config for the convention.
%
% RERUNNING: delete data/single_multi_TS first. run_all_paper_analyses refuses a
% run directory that is not absent or empty, because a reused one accumulates
% sweep folders rather than replacing them. figs/single_multi_TS needs no such
% care; it is overwritten in place.
%
% THE FIGURE REGISTRY IS NOT DUPLICATED. This delegates to paper_config, so the
% figure list, the per-figure preset overrides and panelA_gammas stay defined in
% one place.
%
% THREE CONDITIONS, AND WHAT THAT DOES TO THE FIGURES. Condition-keyed plotters
% get three columns instead of seven. Anything that hardcodes a seven-condition
% layout will show up as a squashed or empty panel -- which is a reason to run
% this, not a reason not to.
%
% See also: paper_config, fast4_config, single_multi_TS_run,
%           srnn_condition_titles, srnn_param_preset

arguments
    opts.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond'
    opts.run_mode    (1,:) char = 'fast'
    opts.run_dir     (1,:) char = 'data/single_multi_TS'
    opts.fig_root    (1,:) char = 'figs/single_multi_TS'
    % Off, like paper_config's own default: figures are built and saved but do
    % not pop up and steal focus.
    opts.visible_figures (1,1) logical = false
end

cfg = paper_config( ...
    'preset_name', opts.preset_name, ...
    'run_mode',    opts.run_mode, ...
    'run_dir',     opts.run_dir, ...
    'fig_root',    opts.fig_root, ...
    'visible_figures', opts.visible_figures);
end

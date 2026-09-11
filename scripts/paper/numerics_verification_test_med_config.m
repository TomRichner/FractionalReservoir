function cfg = numerics_verification_test_med_config()
% NUMERICS_VERIFICATION_TEST_MED_CONFIG The numerics-verification stage alone, at 'medium'.
%
%   cfg = numerics_verification_test_med_config();
%   numerics_verification_test_med_run          % runs the stage, then the two figures
%
% A mini config for exercising ONE stage and its two figures without the
% sweep pipeline. It exists because run_all_paper_analyses runs every stage
% unconditionally (its stage list is a literal, not config-driven), so the
% only way to iterate on one stage is to call it directly -- which the
% matching _run script does -- and then hand make_all_paper_figures a registry
% holding just the figures that read it.
%
% EVERY SETTING IS STATED HERE, like every *_config.m since 2026-09-09; this
% does not call paper_config. The preset is the same network as
% single_multi_TS_independent_med_config, so the numbers here are the ones
% that run's supplemental figures will show at 'medium'.
%
% The medium sibling of numerics_verification_test_config: 10 s of LLE
% accumulation on the full network and n = 40 for the QR comparison, in ITS
% OWN roots so the fast run is not overwritten.
%
% Both roots are private to this config, so it cannot touch any other run.
% Unlike the full pipeline the stage does not require an empty run_dir, so
% this can be rerun in place.
%
% See also: numerics_verification_test_med_run, run_numerics_verification,
%           fig_numerics_verification, single_multi_TS_independent_med_config

cfg = struct();

cfg.preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
cfg.run_mode    = 'medium';

cfg.run_dir  = 'data/numerics_verification_test_med';   % the stage writes <run_dir>/numerics_verification
cfg.fig_root = 'figs/numerics_verification_test_med';   % overwritten in place

cfg.visible_figures = false;

% Only the two entries that read this stage. Names are the output folders,
% and the two fig_tags are not prefixes of each other (save_figure_stable
% deletes <tag>* before saving).
F = {};
F = add(F, 'fig_numerics_solver',     @fig_numerics_verification, false, {'variant', 'solver'});
F = add(F, 'fig_numerics_lya_method', @fig_numerics_verification, false, {'variant', 'lya_method'});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end

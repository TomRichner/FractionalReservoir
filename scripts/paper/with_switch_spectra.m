function cfg = with_switch_spectra(cfg, K, n_samples, first_state_s, opts)
% WITH_SWITCH_SPECTRA Figure 3 with its heat map on the switching run, one top-K pass throughout (FD-034).
%
%   cfg = WITH_SWITCH_SPECTRA(cfg)                 K = 100, 40 states
%   cfg = WITH_SWITCH_SPECTRA(cfg, K, n_samples)
%   cfg = WITH_SWITCH_SPECTRA(cfg, K, n_samples, first_state_s)
%   cfg = WITH_SWITCH_SPECTRA(..., 'as_default', true)   the paper's Figure 3 (below)
%
% THE mu7chain BUNDLES' DEFAULT since FD-046 (Tom, 2026-11-20: "we will use
% the switching version"; docs/Plan_Figure3_update.md, step 6):
% mu7chain_config ends with with_switch_spectra(cfg, 100, 40, 1,
% 'as_default', true). Without 'as_default' it is the opt-in it was
% (FD-034), adding versions beside the registered entries. In the registered Figure 3 the expansion
% row comes from the local-Lyapunov stage (a switching-input run) and the
% eigenvalue density, with its optional Ritz overlay, from the eig-heatmap
% stage (a separate baseline run with a top-K pass of its own). This turns
% on the version in which both are one trajectory and one pass:
%
%   cfg = with_switch_spectra(mu7chain_config('essential'));
%   run_all_paper_analyses(cfg);      % the local_lyapunov stage saves the spectra
%   make_all_paper_figures(cfg);      % the two entries below are drawn too
%
% What it does:
%   * cfg.local_lyapunov_K = K: the number of directions the stage tracks,
%     and so the number the expansion row sums, the size of H and the
%     number of exponent ticks. NOTE that this changes the stage's saved
%     results, and with them the registered fig_main3_grouped and
%     fig_local_lyapunov, whenever K is not the bundle's own value.
%   * cfg.local_lyapunov_spectra_n_samples = n_samples:
%     run_all_paper_analyses passes it to run_local_lyapunov as
%     'spectra_n_samples' (the Jacobian's spectra and H = Q'JQ of the
%     stage's own top-K pass at that many states of its run; see
%     run_local_lyapunov).
%   * cfg.local_lyapunov_spectra_start_s = first_state_s (default 1 s;
%     run_local_lyapunov's 'spectra_start_s'), only when the config's
%     accumulation window opens before that time. The states are evenly
%     spaced from there to the end of the run, so 40 of them are t = 1, 2,
%     ..., 40 s whether the window is [1, 40] s (the first version, where
%     the field is not set and nothing changes) or [0, 40] s (the timing of
%     docs/Plan_Figure3_update.md; Tom has not decided on moving the
%     states, 2026-11-20).
%   * registry entries, each after the entry it is a version of, in a
%     directory of its own, none in_paper:
%       fig_eig_heatmap_switch     the density of the switching run with the
%                                  Ritz overlay and its inset
%       fig_main3_grouped_switch   the composite: input, expansion, heat map
%                                  (no lambda_1 row), same overlay and inset
%       fig_main3_grouped_switch_hist  that composite with a fourth row, (D):
%                                  the eigenvalues projected on the real axis
%                                  against a uniform disk (Tom, 2026-11-20)
%       fig_eig_real_histogram_switch  the stand-alone histogram of the
%                                  switching run (only when the config
%                                  registers fig_eig_real_histogram)
%     The entries they are versions of are untouched.
%   * WITH 'as_default', true, the paper's Figure 3 is the switching one:
%     the entry fig_main3_grouped (its name, its in_paper) draws the
%     four-row composite (input, expansion, heat map with the Ritz overlay
%     and inset, and row (D), the real-axis histogram), and the baseline
%     composite it drew until then is kept, right after it, as
%     fig_main3_grouped_baseline, not in_paper (it still reads the
%     eig_heatmap stage). fig_main3_grouped_switch and _switch_hist are not
%     added (the first is not chosen, the second would be fig_main3_grouped
%     itself); fig_eig_heatmap_switch and fig_eig_real_histogram_switch are
%     added as above. The composite saves under the tag
%     Fig_Main3_Grouped_noC_switch_ritz_inset_realhist (fig_grouped_main
%     names it by its options) in <fig_root>/fig_main3_grouped/; the
%     baseline under Fig_Main3_Grouped in <fig_root>/fig_main3_grouped_baseline/.
%
% The config must register fig_eig_heatmap and fig_main3_grouped and name a
% local-Lyapunov stage with input steps (the mu7chain bundles do).
% mu7chain_run rebuilds the config per figure set, so a config changed by
% hand does not reach its figure passes; call this at the end of a bundle's
% config function to make it that bundle's default.
%
% See also: with_ritz_overlay, run_local_lyapunov, fig_eig_heatmap,
%           fig_main3_stability, run_all_paper_analyses
arguments
    cfg (1,1) struct
    K (1,1) double {mustBeInteger, mustBePositive} = 100
    n_samples (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(n_samples, 2)} = 40
    first_state_s (1,1) double = 1
    opts.as_default (1,1) logical = false
end
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
need = {'fig_eig_heatmap', 'fig_main3_grouped'};
missing = setdiff(need, names);
if ~isempty(missing)
    error('with_switch_spectra:MissingEntry', ...
        'The config''s registry lacks %s, which the added entries are versions of.', strjoin(missing, ', '));
end
added = {'fig_eig_heatmap_switch', 'fig_main3_grouped_switch'};
if any(ismember([added, {'fig_main3_grouped_baseline'}], names))
    error('with_switch_spectra:AlreadyAdded', 'This config already carries the switching-run entries.');
end
cfg.local_lyapunov_K = K;
cfg.local_lyapunov_spectra_n_samples = n_samples;
if isfield(cfg, 'local_lyapunov_accumulation_start_s') && cfg.local_lyapunov_accumulation_start_s < first_state_s
    % Only when the window opens before it: otherwise the stage's default
    % (the first state at the accumulation start) is that state already.
    cfg.local_lyapunov_spectra_start_s = first_state_s;
end

heat = cfg.figures{strcmp(names, 'fig_eig_heatmap')};
main3 = cfg.figures{strcmp(names, 'fig_main3_grouped')};
overlay = {'spectra_source', 'local_lyapunov', 'ritz_overlay', true, 'ritz_inset', true};
added_heat = struct('name', added{1}, 'fn', heat.fn, 'in_paper', false, 'args', {[heat.args, overlay]});
added_main3 = struct('name', added{2}, 'fn', main3.fn, 'in_paper', false, ...
    'args', {[main3.args, overlay, {'show_lambda_row', false}]});
added_main3_hist = struct('name', [added{2} '_hist'], 'fn', main3.fn, 'in_paper', false, ...
    'args', {[main3.args, overlay, {'show_lambda_row', false, 'real_histogram_row', true}]});
has_hist = ismember('fig_eig_real_histogram', names);
if has_hist
    hist0 = cfg.figures{strcmp(names, 'fig_eig_real_histogram')};
    added_hist = struct('name', 'fig_eig_real_histogram_switch', 'fn', hist0.fn, 'in_paper', false, ...
        'args', {[hist0.args, {'spectra_source', 'local_lyapunov'}]});
end

if opts.as_default
    % The paper's name and flag on the four-row switching composite; the
    % baseline composite kept beside it under a name of its own.
    paper_main3 = added_main3_hist;
    paper_main3.name = main3.name;
    paper_main3.in_paper = main3.in_paper;
    baseline = main3;
    baseline.name = 'fig_main3_grouped_baseline';
    baseline.in_paper = false;
    main3_entries = {paper_main3, baseline};
else
    main3_entries = {main3, added_main3, added_main3_hist};
end
F = {};
for k = 1:numel(cfg.figures)
    if strcmp(names{k}, 'fig_main3_grouped')
        F = [F, main3_entries]; %#ok<AGROW>
        continue;
    end
    F{end + 1} = cfg.figures{k}; %#ok<AGROW>
    if strcmp(names{k}, 'fig_eig_heatmap');   F{end + 1} = added_heat;  end %#ok<AGROW>
    if has_hist && strcmp(names{k}, 'fig_eig_real_histogram'); F{end + 1} = added_hist; end %#ok<AGROW>
end
cfg.figures = F;
end

function cfg = with_switch_spectra(cfg, K, n_samples)
% WITH_SWITCH_SPECTRA Figure 3 with its heat map on the switching run, one top-K pass throughout (FD-034).
%
%   cfg = WITH_SWITCH_SPECTRA(cfg)                 K = 100, 40 states
%   cfg = WITH_SWITCH_SPECTRA(cfg, K, n_samples)
%
% OPT-IN, like with_ritz_overlay: no bundle sets it, so no registered run,
% stage output or figure changes. In the registered Figure 3 the expansion
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
%   * two registry entries, each after the entry it is a version of, in a
%     directory of its own, neither in_paper:
%       fig_eig_heatmap_switch     the density of the switching run with the
%                                  Ritz overlay and its inset
%       fig_main3_grouped_switch   the composite: input, expansion, heat map
%                                  (no lambda_1 row), same overlay and inset
%     The entries they are versions of are untouched.
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
end
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
need = {'fig_eig_heatmap', 'fig_main3_grouped'};
missing = setdiff(need, names);
if ~isempty(missing)
    error('with_switch_spectra:MissingEntry', ...
        'The config''s registry lacks %s, which the added entries are versions of.', strjoin(missing, ', '));
end
added = {'fig_eig_heatmap_switch', 'fig_main3_grouped_switch'};
if any(ismember(added, names))
    error('with_switch_spectra:AlreadyAdded', 'This config already carries the switching-run entries.');
end
cfg.local_lyapunov_K = K;
cfg.local_lyapunov_spectra_n_samples = n_samples;

heat = cfg.figures{strcmp(names, 'fig_eig_heatmap')};
main3 = cfg.figures{strcmp(names, 'fig_main3_grouped')};
overlay = {'spectra_source', 'local_lyapunov', 'ritz_overlay', true, 'ritz_inset', true};
added_heat = struct('name', added{1}, 'fn', heat.fn, 'in_paper', false, 'args', {[heat.args, overlay]});
added_main3 = struct('name', added{2}, 'fn', main3.fn, 'in_paper', false, ...
    'args', {[main3.args, overlay, {'show_lambda_row', false}]});

F = {};
for k = 1:numel(cfg.figures)
    F{end + 1} = cfg.figures{k}; %#ok<AGROW>
    if strcmp(names{k}, 'fig_eig_heatmap');   F{end + 1} = added_heat;  end %#ok<AGROW>
    if strcmp(names{k}, 'fig_main3_grouped'); F{end + 1} = added_main3; end %#ok<AGROW>
end
cfg.figures = F;
end

function cfg = with_ritz_overlay(cfg, K)
% WITH_RITZ_OVERLAY Add the Ritz pass and its figures to a bundle config (FD-033).
%
%   cfg = WITH_RITZ_OVERLAY(cfg)        K = 50
%   cfg = WITH_RITZ_OVERLAY(cfg, K)
%
% OPT-IN. The Ritz overlay of Figure 3 (the eigenvalues of the Jacobian
% projected on the leading Lyapunov directions, drawn on the eigenvalue
% density) is not part of any bundle's registered default: no config sets it,
% so no run pays for the pass and no registered figure changes. This turns
% it on for one config:
%
%   cfg = with_ritz_overlay(mu7chain_config('essential'));
%   run_all_paper_analyses(cfg);      % the eig_heatmap stage makes the pass
%   make_all_paper_figures(cfg);      % the three entries below are drawn too
%
% What it does:
%   * cfg.eig_heatmap_ritz_K = K: run_all_paper_analyses passes it to
%     run_eig_heatmap as 'ritz_K' (a top-K pass per example, condition and
%     state scaling at the sampled states; see run_eig_heatmap). K = 50 is
%     the number of directions Figure 3's expansion panel tracks.
%   * three registry entries, each after the entry it is a version of, in a
%     directory of its own (the tags share a prefix and a save clears by
%     prefix):
%       fig_eig_heatmap_ritz         the density with the points
%       fig_eig_heatmap_ritz_inset   the same with the inset near zero
%       fig_main3_grouped_ritz       the composite Figure 3 with both
%     None is marked in_paper. The entries they are versions of are
%     untouched, so the registered figures are what they were.
%
% The config must register fig_eig_heatmap and fig_main3_grouped (the
% mu7chain bundles do). To make the overlay a bundle's default, call this at
% the end of its config function; mu7chain_run rebuilds the config per
% figure set, so a config changed by hand does not reach its figure passes.
%
% See also: run_eig_heatmap, fig_eig_heatmap, fig_main3_stability,
%           run_all_paper_analyses, make_all_paper_figures
arguments
    cfg (1,1) struct
    K (1,1) double {mustBeInteger, mustBePositive} = 50
end
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
need = {'fig_eig_heatmap', 'fig_main3_grouped'};
missing = setdiff(need, names);
if ~isempty(missing)
    error('with_ritz_overlay:MissingEntry', ...
        'The config''s registry lacks %s, which the overlay''s entries are versions of.', strjoin(missing, ', '));
end
if any(startsWith(names, 'fig_eig_heatmap_ritz')) || any(strcmp(names, 'fig_main3_grouped_ritz'))
    error('with_ritz_overlay:AlreadyAdded', 'This config already carries the Ritz overlay''s entries.');
end
cfg.eig_heatmap_ritz_K = K;

heat = cfg.figures{strcmp(names, 'fig_eig_heatmap')};
main3 = cfg.figures{strcmp(names, 'fig_main3_grouped')};
overlay = {'ritz_overlay', true};
added_heat = { ...
    struct('name', 'fig_eig_heatmap_ritz', 'fn', heat.fn, 'in_paper', false, 'args', {[heat.args, overlay]}), ...
    struct('name', 'fig_eig_heatmap_ritz_inset', 'fn', heat.fn, 'in_paper', false, ...
        'args', {[heat.args, overlay, {'ritz_inset', true}]})};
added_main3 = {struct('name', 'fig_main3_grouped_ritz', 'fn', main3.fn, 'in_paper', false, ...
    'args', {[main3.args, overlay, {'ritz_inset', true}]})};

F = {};
for k = 1:numel(cfg.figures)
    F{end + 1} = cfg.figures{k}; %#ok<AGROW>
    if strcmp(names{k}, 'fig_eig_heatmap');   F = [F, added_heat];  end %#ok<AGROW>
    if strcmp(names{k}, 'fig_main3_grouped'); F = [F, added_main3]; end %#ok<AGROW>
end
cfg.figures = F;
end

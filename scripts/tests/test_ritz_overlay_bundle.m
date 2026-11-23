% test_ritz_overlay_bundle.m - with_ritz_overlay: the opt-in that adds the
% Ritz pass and its figures to a bundle config (FD-033). Configs only; no
% simulation and no figure.
%
% Checks:
%   1. No bundle carries the overlay by default: mu7chain_config, in
%      'essential' and 'fast', has no cfg.eig_heatmap_ritz_K and none of the
%      opt-in's entries (fig_eig_heatmap_ritz*, fig_main3_grouped_baseline_ritz;
%      the bundles' own fig_main3_grouped_ritz, the switching composite with
%      the overlay of with_switch_spectra, is another entry, 2026-11-21).
%   2. with_ritz_overlay(cfg) sets cfg.eig_heatmap_ritz_K = 50 (or the K
%      given) and adds exactly three entries, each right after the entry it
%      is a version of, all 'appendix' (none in_paper): fig_eig_heatmap_ritz and
%      fig_eig_heatmap_ritz_inset (fig_eig_heatmap's function and args plus
%      the overlay's) and fig_main3_grouped_baseline_ritz (the baseline composite's:
%      fig_main3_grouped_baseline in the mu7chain bundles since FD-046, whose
%      fig_main3_grouped is the switching composite; fig_main3_grouped in a
%      config without that entry). Every other entry, and every other field
%      of the config, is unchanged.
%   3. It errors :AlreadyAdded when applied twice and :MissingEntry on a
%      config without the entries it builds on.
%   (That the functions accept the added args is exercised by
%   test_eig_heatmap_stage, which draws with them.)
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: with_ritz_overlay, mu7chain_config, run_eig_heatmap, fig_eig_heatmap

fprintf('=== Testing with_ritz_overlay ===\n\n');
all_passed = true;

ok = true;
for mode = {'essential', 'fast'}
    c = mu7chain_config(mode{1});
    n = cellfun(@(f) f.name, c.figures, 'UniformOutput', false);
    ok = ok && ~isfield(c, 'eig_heatmap_ritz_K') && ~any(startsWith(n, 'fig_eig_heatmap_ritz')) && ~any(strcmp(n, 'fig_main3_grouped_baseline_ritz'));
end
all_passed = check('mu7chain_config carries no Ritz pass and none of the opt-in''s entries by default', ok) && all_passed;

base = mu7chain_config('essential');
cfg  = with_ritz_overlay(base);
names0 = cellfun(@(f) f.name, base.figures, 'UniformOutput', false);
names  = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
added  = {'fig_eig_heatmap_ritz', 'fig_eig_heatmap_ritz_inset', 'fig_main3_grouped_baseline_ritz'};
all_passed = check('K defaults to 50, and a given K is kept', ...
    cfg.eig_heatmap_ritz_K == 50 && getfield(with_ritz_overlay(base, 20), 'eig_heatmap_ritz_K') == 20) && all_passed; %#ok<GFLD>
ih = find(strcmp(names, 'fig_eig_heatmap')); im = find(strcmp(names, 'fig_main3_grouped_baseline'));
all_passed = check('three entries added, each right after the entry it is a version of', ...
    numel(names) == numel(names0) + 3 && isequal(setdiff(names, names0), sort(added)) && ...
    isequal(names(ih + 1:ih + 2), added(1:2)) && strcmp(names{im + 1}, added{3})) && all_passed;
kept = cfg.figures(~ismember(names, added));
rest = rmfield(cfg, {'figures', 'eig_heatmap_ritz_K'});
all_passed = check('every other entry and every other field of the config is unchanged', ...
    isequal(kept, base.figures) && isequal(rest, rmfield(base, 'figures'))) && all_passed;
e1 = cfg.figures{ih + 1}; e2 = cfg.figures{ih + 2}; e3 = cfg.figures{im + 1};
h0 = base.figures{strcmp(names0, 'fig_eig_heatmap')}; m0 = base.figures{strcmp(names0, 'fig_main3_grouped_baseline')};
all_passed = check('the added entries are the originals'' functions and args plus the overlay''s, none in_paper', ...
    isequal(e1.fn, h0.fn) && isequal(e1.args, [h0.args, {'ritz_overlay', true}]) && ...
    isequal(e2.args, [h0.args, {'ritz_overlay', true, 'ritz_inset', true}]) && ...
    isequal(e3.fn, m0.fn) && isequal(e3.args, [m0.args, {'ritz_overlay', true, 'ritz_inset', true}]) && ...
    ~e1.in_paper && ~e2.in_paper && ~e3.in_paper) && all_passed;

all_passed = check('applied twice errors :AlreadyAdded', ...
    throws_id(@() with_ritz_overlay(cfg), 'with_ritz_overlay:AlreadyAdded')) && all_passed;
% Without the baseline entry, the composite version follows fig_main3_grouped.
nob = base; nob.figures = base.figures(~strcmp(names0, 'fig_main3_grouped_baseline'));
c2 = with_ritz_overlay(nob); n2 = cellfun(@(f) f.name, c2.figures, 'UniformOutput', false);
m2 = nob.figures{strcmp(cellfun(@(f) f.name, nob.figures, 'UniformOutput', false), 'fig_main3_grouped')};
all_passed = check('without fig_main3_grouped_baseline the composite version follows fig_main3_grouped', ...
    strcmp(n2{find(strcmp(n2, 'fig_main3_grouped')) + 1}, 'fig_main3_grouped_baseline_ritz') && ...
    isequal(c2.figures{strcmp(n2, 'fig_main3_grouped_baseline_ritz')}.args, [m2.args, {'ritz_overlay', true, 'ritz_inset', true}])) && all_passed;
bare = base; bare.figures = base.figures(~strcmp(names0, 'fig_eig_heatmap'));
all_passed = check('a config without fig_eig_heatmap errors :MissingEntry', ...
    throws_id(@() with_ritz_overlay(bare), 'with_ritz_overlay:MissingEntry')) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL with_ritz_overlay TESTS PASSED ===\n');
else
    fprintf('=== SOME with_ritz_overlay TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end

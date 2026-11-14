% test_noise_levels_timeseries.m - fig_noise_levels_timeseries on a tiny network.
%
% Two conditions (one without depression, one with the chain), the
% mu7chain preset at n = 60, 1 s settling and 1 s shown. Checks:
%   1. the figure draws and writes its PNG under the tag given;
%   2. the series differ between two different levels, in every condition;
%   3. two equal levels give bit-identical series (same network, same
%      seeds, noise seed included), within one call and across two calls;
%   4. the panel of a condition without depression holds no trace.
% A few seconds; no pool needed.

p = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
common = {'preset_name', p, 'conditions', {'no_adaptation', 'sfa3_stdc2'}, ...
    'model_args', {'n', 60, 'indegree', 20}, 'T_settle', 1, 'T_show', 1, ...
    'n_show_E', 4, 'n_show_I', 4, 'visible', false, 'verbose', 'near-none'};
out_dir = fullfile(tempdir, 'test_noise_levels_timeseries');
if isfolder(out_dir), rmdir(out_dir, 's'); end
c1 = onCleanup(@() rmdir(out_dir, 's'));

A = fig_noise_levels_timeseries(common{:}, 'levels', [0.02 0.04], 'out_dir', out_dir, ...
    'level_label', '\sigma_u', 'tag', 'Fig_test_noise_levels');
c2 = onCleanup(@() close(A.figs));
assert(isfile(fullfile(out_dir, 'Fig_test_noise_levels.png')), 'the PNG was not written');
assert(any(strcmp(A.files, 'Fig_test_noise_levels.png')), 'out.files does not list the PNG');
assert(isequal(size(A.series), [2 2]), 'series is not levels x conditions');

fields = {'x', 'r', 'synaptic_output', 'depression'};
for i = 1:2
    for f = fields
        a = A.series{1, i}.(f{1}); b = A.series{2, i}.(f{1});
        if isempty(a), continue; end
        assert(~isequal(a, b), 'condition %d: %s does not change with the level', i, f{1});
    end
end
assert(isempty(A.series{1, 1}.depression), 'no_adaptation should have no depression factor');
assert(~isempty(A.series{1, 2}.depression), 'sfa3_stdc2 should have a depression factor');

B = fig_noise_levels_timeseries(common{:}, 'levels', [0.02 0.02], 'save', false);
close(B.figs);
for i = 1:2
    assert(isequal(B.series{1, i}, B.series{2, i}), 'condition %d: equal levels differ within a call', i);
    assert(isequal(A.series{1, i}, B.series{1, i}), 'condition %d: the same level differs across calls', i);
end
fprintf('PASS test_noise_levels_timeseries\n');

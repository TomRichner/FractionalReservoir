% TEST_MC_GAIN_BIAS_FIGURES Block C's figures and statistics on a synthetic result (FD-010, Plan C section 5).
%
% Quick, no simulation. On make_synthetic_mc_gain_sweep's table (the shape
% Plan C section 2 fixes for run_mc_gain_sweep), checks that
%   * read_mc_gain_sweep reads a table, a struct and a .mat file alike,
%     accepts the stage's likely alternative column names, fills a missing
%     optional column with NaN and lists it, and refuses a missing required one;
%   * fig_mc_gain_bias_maps, fig_mc_gain_curves (one result and a list of
%     three as rows) and fig_mc_robustness build, invisible and unsaved, for
%     the paper's three conditions and for four, with the expected panels;
%   * a result is found by run_dir in <run_dir>/mc_gain_sweep;
%   * mc_gain_bias_statistics recovers what the generator put in: each
%     condition's peak and its location, the failed job, the order of the
%     shares of the map at half the best of any condition (chained MTS >
%     1TS > none), chained MTS higher at most grid points, a good Gaussian fit
%     for the adapting conditions, and NaN paired columns when the reference
%     condition is absent; and writes its .csv and .md.
%
% Asserts; silent until it throws, then prints one PASS line.
% See also: make_synthetic_mc_gain_sweep, read_mc_gain_sweep, mc_gain_bias_statistics

four  = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
three = four(1:3);
[T, settings, truth] = make_synthetic_mc_gain_sweep(four);
D = read_mc_gain_sweep(struct('T', T, 'settings', settings));
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));

%% The reader
assert(isequal(D.conditions, four) && numel(D.gains) == 9 && numel(D.biases) == 6 && numel(D.seeds) == 3);
assert(nnz(D.T.failed) == truth.n_failed && all(isnan(D.T.lle_input(D.T.failed))));
assert(D.preset_gain == 1 && D.preset_bias == 0.2 && isempty(D.missing));
% A file, with the stage's possible column names
T2 = renamevars(T, {'gain', 'bias', 'mc_total', 'lle_input'}, {'level_of_chaos', 'S_c', 'MC', 'lambda1_input'});
T2 = removevars(T2, 'lle_noinput');
results = T2;
f = fullfile(tmp, 'mc_gain_sweep_results.mat');
save(f, 'results', 'settings');
D2 = read_mc_gain_sweep(f);
assert(isequaln(D2.T.mc_total, D.T.mc_total) && isequaln(D2.T.lle_input, D.T.lle_input));
assert(isequal(D2.missing, {'lle_noinput'}) && all(isnan(D2.T.lle_noinput)));
threw = false;
try
    read_mc_gain_sweep(removevars(T, 'seed'));
catch err
    threw = strcmp(err.identifier, 'read_mc_gain_sweep:MissingColumn');
end
assert(threw, 'a missing required column must be refused');

%% The figures, three and four conditions
for sel = {three, four}
    nc = numel(sel{1});
    o = fig_mc_gain_bias_maps('data', D, 'conditions', sel{1}, 'save', false, 'visible', false);
    n_ax = numel(findall(o.figs, 'Type', 'axes'));
    assert(n_ax == 3 * nc, 'maps: %d axes for %d conditions', n_ax, nc);
    assert(numel(findall(o.figs, 'Type', 'colorbar')) == 3);
    assert(isempty(o.files)); close(o.figs);
    o = fig_mc_gain_bias_maps('data', D, 'conditions', sel{1}, 'show_iqr', false, 'save', false, 'visible', false);
    assert(numel(findall(o.figs, 'Type', 'axes')) == 2 * nc); close(o.figs);

    o = fig_mc_gain_curves('data', D, 'conditions', sel{1}, 'save', false, 'visible', false);
    assert(numel(findall(o.figs, 'Type', 'axes')) == 3); close(o.figs);

    o = fig_mc_robustness('data', D, 'conditions', sel{1}, 'save', false, 'visible', false);
    ax = findall(o.figs, 'Type', 'axes');
    assert(numel(ax) == 2);
    assert(numel(findall(o.figs, 'Type', 'stair')) == nc); close(o.figs);
end
% Three results as rows, one of them read from a file, one swept along gain only
D_gain_only = read_mc_gain_sweep(struct('T', T(T.bias == 0.2, :), 'settings', settings));
o = fig_mc_gain_curves('data', {D, D2, D_gain_only}, 'conditions', three, ...
    'row_labels', {'1 : 1', '3 : 2', '2 : 3'}, 'bias', 0.25, 'save', false, 'visible', false);
assert(numel(findall(o.figs, 'Type', 'axes')) == 9);
ttl = arrayfun(@(a) strjoin(cellstr(a.Title.String), ' '), findall(o.figs, 'Type', 'axes'), ...
    'UniformOutput', false);
assert(nnz(contains(ttl, 'nearest to 0.25')) == 3, 'each row title must state the nearest bias level');
close(o.figs);
% Located by run_dir
run_dir = fullfile(tmp, 'run'); mkdir(fullfile(run_dir, 'mc_gain_sweep'));
copyfile(f, fullfile(run_dir, 'mc_gain_sweep', 'mc_gain_sweep_results.mat'));
o = fig_mc_robustness('run_dir', run_dir, 'conditions', three, 'save', false, 'visible', false);
assert(contains(o.source, 'mc_gain_sweep_results.mat')); close(o.figs);

%% The statistics recover the generator's pattern
[S, files] = mc_gain_bias_statistics('data', D, 'out_dir', tmp, 'verbose', 'near-none');
assert(numel(files) == 2 && all(cellfun(@isfile, files)));
assert(isequal(S.condition', four) && height(S) == 4);
for c = 1:4
    p = truth.peak.(S.condition{c});
    assert(abs(S.peak_gain(c) - p(1)) < 1e-9 && abs(S.peak_bias(c) - p(2)) < 1e-9, ...
        '%s: peak at (%g, %g), generator put it at (%g, %g)', S.condition{c}, ...
        S.peak_gain(c), S.peak_bias(c), p(1), p(2));
end
row = @(n) find(strcmp(S.condition, n));
assert(S.n_failed(row('no_adaptation')) == 1 && sum(S.n_failed) == 1);
sa = S.share_half_any(cellfun(row, truth.share_any_order));
assert(all(diff(sa) < 0), 'share at half the best of any: %s', mat2str(sa', 3));
assert(S.peak_mc(row('no_adaptation')) == max(S.peak_mc), 'the narrow ridge holds the highest cell');
assert(all(S.share_ref_higher([row('no_adaptation'), row('sfa1_std1'), row('sfa3_std2')]) > 0.5));
assert(isnan(S.share_ref_higher(row('sfa3_stdc2'))));
assert(all(S.fit_r2([row('sfa1_std1'), row('sfa3_stdc2')]) > 0.95));
assert(S.fit_fwhm_log10(row('sfa3_stdc2')) > 2 * S.fit_fwhm_log10(row('no_adaptation')), ...
    'the chain''s width along gain must exceed the ridge''s');
% Absent reference: paired columns NaN, the rest computed
S3 = mc_gain_bias_statistics('data', D, 'conditions', {'no_adaptation', 'sfa1_std1'}, ...
    'reference', 'sfa3_stdc2_missing', 'verbose', 'near-none');
assert(height(S3) == 2 && all(isnan(S3.share_ref_higher)) && all(isfinite(S3.peak_mc)));
% A result holding three conditions only
D3 = read_mc_gain_sweep(struct('T', T(~strcmp(T.condition, 'sfa3_std2'), :)));
S4 = mc_gain_bias_statistics('data', D3, 'verbose', 'near-none');
assert(isequal(S4.condition', three));

disp('PASS: test_mc_gain_bias_figures -- reader, three figures at 3 and 4 conditions, statistics recover the synthetic pattern.');

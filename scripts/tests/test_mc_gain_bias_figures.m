% TEST_MC_GAIN_BIAS_FIGURES Block C's figures and statistics on a synthetic result (FD-010, Plan C section 5).
%
% Quick, no simulation. On make_synthetic_mc_gain_sweep's table (the shape
% Plan C section 2 fixes for run_mc_gain_sweep), checks that
%   * read_mc_gain_sweep reads a table, a struct and a .mat file alike,
%     returns its own output unchanged, fills a missing optional column
%     with NaN and lists it, and refuses a missing required one;
%   * fig_mc_gain_bias_maps, fig_mc_gain_curves (one result and a list of
%     three as rows) and fig_mc_robustness build, invisible and unsaved, for
%     the paper's three conditions and for four, with the expected panels;
%   * the maps' linear lambda_1 scale spans +-lle_clip;
%   * with one seed the maps drop the blank spread row and say "one seed";
%     the lambda_1 = 0 boundary is drawn;
%   * a result is found by run_dir in <run_dir>/mc_gain_sweep;
%   * fig_mc_vs_lambda1 draws four panels and one colorbar, writes its PNG,
%     and its Pearson r, Spearman rho and counts for one condition, in
%     out.stats and in the panel's text, match a direct computation;
%   * the general fig_mc_sweep_scatter behind it: rate_mean weights the
%     rates by the named preset's fractions (or a stage override of f_E)
%     and is refused when they are unknown, frac_on_slope is 1 - pinned -
%     silent; gain against rate_mean colored by frac_on_slope draws gain on
%     a log axis, an overlay on its one panel, no sign split, and a Pearson
%     r that matches a direct computation;
%   * read_mc_gain_sweep_jobs reads a folder of finished job files (half of
%     them, plus a *.tmp save in progress it must skip) into the same table,
%     without touching the folder;
%   * the tuning width at half the best of any condition, on a 4-point
%     example worked by hand (2 octaves closed; 3 octaves open);
%   * make_block_C_figures writes a set's figures and tables for a result
%     that lacks one of the set's conditions, with a note;
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
assert(isequaln(read_mc_gain_sweep(D), D), 'the reader must return its own output unchanged');
assert(isequaln(read_mc_gain_sweep(T).T, D.T), 'a bare table reads into the same rows');
% A file, as the stage saves it (T and settings), here without one optional column
T_full = T;
T = removevars(T_full, 'lle_noinput');
f = fullfile(tmp, 'mc_gain_sweep_results.mat');
save(f, 'T', 'settings');
T = T_full;
D2 = read_mc_gain_sweep(f);
assert(isequaln(D2.T.mc_total, D.T.mc_total) && isequaln(D2.T.lle_input, D.T.lle_input));
assert(isequal(D2.missing, {'lle_noinput'}) && all(isnan(D2.T.lle_noinput)));
assert(strcmp(D2.source, f) && isequal(D2.conditions, four));
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
    assert(numel(findall(o.figs, 'Type', 'axes')) == 2 * nc);
    % The lambda_1 = 0 boundary is drawn as cell edges where the generator
    % puts a sign change.
    assert(~isempty(findall(o.figs, 'Tag', 'lle_zero')), 'maps: no lambda_1 = 0 boundary drawn');
    close(o.figs);
    if nc == 3
        % Linear color scales: lambda_1's colorbar spans +-lle_clip
        o = fig_mc_gain_bias_maps('data', D, 'conditions', sel{1}, 'mc_scale', 'linear', ...
            'lle_scale', 'linear', 'lle_clip', 5, 'save', false, 'visible', false);
        cbs = findall(o.figs, 'Type', 'colorbar');
        lcb = cbs(contains(arrayfun(@(b) b.Label.String, cbs, 'UniformOutput', false), 'lambda'));
        assert(isscalar(lcb) && isequal(lcb.Limits, [-5 5]) && contains(lcb.Label.String, 'linear'), ...
            'maps: linear lambda_1 scale must span +-lle_clip');
        close(o.figs);
    end

    o = fig_mc_gain_curves('data', D, 'conditions', sel{1}, 'save', false, 'visible', false);
    assert(numel(findall(o.figs, 'Type', 'axes')) == 3); close(o.figs);

    o = fig_mc_robustness('data', D, 'conditions', sel{1}, 'save', false, 'visible', false);
    ax = findall(o.figs, 'Type', 'axes');
    assert(numel(ax) == 2);
    assert(numel(findall(o.figs, 'Type', 'stair')) == nc); close(o.figs);
end
% One seed: the spread row would be blank, so it is dropped and the caption says so
D_one = read_mc_gain_sweep(struct('T', T(T.seed == D.seeds(1), :), 'settings', settings));
o = fig_mc_gain_bias_maps('data', D_one, 'conditions', three, 'save', false, 'visible', false);
assert(numel(findall(o.figs, 'Type', 'axes')) == 2 * 3 && numel(findall(o.figs, 'Type', 'colorbar')) == 2);
tl = findall(o.figs, 'Type', 'tiledlayout');
cap = arrayfun(@(t) strjoin(cellstr(t.Title.String), ' '), tl, 'UniformOutput', false);
assert(any(contains(cap, 'One seed')), 'one-seed caption missing');
close(o.figs);
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

%% Memory capacity against lambda_1: draws, writes its file, its statistics
% match a direct computation (on the unclipped values of the finite runs)
o = fig_mc_vs_lambda1('data', D, 'conditions', four, 'out_dir', tmp, 'visible', false);
assert(numel(findall(o.figs, 'Type', 'axes')) == 4 && isscalar(findall(o.figs, 'Type', 'colorbar')));
assert(any(endsWith(o.files, 'Fig_MC_vs_lambda1.png')) && isfile(fullfile(tmp, 'Fig_MC_vs_lambda1.png')));
c = 'no_adaptation';
k = strcmp(D.T.condition, c) & ~D.T.failed & isfinite(D.T.lle_input) & isfinite(D.T.mc_total);
x = D.T.lle_input(k); y = D.T.mc_total(k);
pearson_direct = @(a, b) sum((a - mean(a)) .* (b - mean(b))) / ...
    sqrt(sum((a - mean(a)).^2) * sum((b - mean(b)).^2));
avg_rank = @(v) sum(v' < v, 2) + (sum(v' == v, 2) + 1) / 2;   % ties share their mean rank
r_direct = pearson_direct(x, y);
rho_direct = pearson_direct(avg_rank(x), avg_rank(y));
S = o.stats(strcmp(o.stats.condition, c), :);
assert(abs(S.pearson_r - r_direct) < 1e-12 && abs(S.spearman_rho - rho_direct) < 1e-12);
assert(S.n == nnz(k) && S.n_stable == nnz(x < 0) && S.n_unstable == nnz(x > 0) && ...
    S.n_failed == nnz(strcmp(D.T.condition, c)) - nnz(k) && S.n_clipped == nnz(abs(x) > 5));
ax = findall(o.figs, 'Type', 'axes');
panel = ax(strcmp(arrayfun(@(a) char(a.Title.String), ax, 'UniformOutput', false), titles_by_name({c})));
txt = findall(panel, 'Tag', 'mc_sweep_scatter_stats');
assert(isscalar(txt) && contains(txt.String{1}, sprintf('Pearson r = %.2f, n = %d', r_direct, nnz(k))) && ...
    contains(txt.String{2}, sprintf('%.2f', rho_direct)), 'the panel''s statistics must match');
close(o.figs);

%% The general scatter: derived quantities, a non-default x and y, an overlay
% rate_mean weights rate_E and rate_I by the fractions of the preset the
% result names (here the main sweep's, f = [0.5 0.5]); without a preset or
% an f override the fractions are unknown and the call refuses.
threw = false;
try
    mc_sweep_quantity('rate_mean', D);
catch err
    threw = strcmp(err.identifier, 'mc_sweep_population_fractions:NoFractions');
end
assert(threw, 'rate_mean without known fractions must be refused');
s_named = settings;
s_named.preset_name = mu7chain_config('essential').mc_preset;
Dn = read_mc_gain_sweep(struct('T', T, 'settings', s_named));
fr = srnn_param_preset(s_named.preset_name).f;
Qr = mc_sweep_quantity('rate_mean', Dn);
assert(max(abs(Qr.values - (fr(1) * Dn.T.rate_E + fr(2) * Dn.T.rate_I)), [], 'omitnan') < 1e-15);
s_named.model_args = {'f_E', 0.8};                % a stage override wins over the preset
Q8 = mc_sweep_quantity('rate_mean', read_mc_gain_sweep(struct('T', T, 'settings', s_named)));
assert(max(abs(Q8.values - (0.8 * T.rate_E + 0.2 * T.rate_I)), [], 'omitnan') < 1e-12);
Qs = mc_sweep_quantity('frac_on_slope', Dn);
assert(isequaln(Qs.values, 1 - Dn.T.frac_pinned - Dn.T.frac_silent));

curve = struct('x', [0.1 1], 'y', [0.5 0.5], 'label', 'ref', 'conditions', {{'no_adaptation'}});
o = fig_mc_sweep_scatter('data', Dn, 'conditions', three, 'x', 'gain', 'y', 'rate_mean', ...
    'color_by', 'frac_on_slope', 'overlay', curve, 'save', false, 'visible', false);
ax = findall(o.figs, 'Type', 'axes');
assert(numel(ax) == 3 && all(strcmp({ax.XScale}, 'log')), 'gain is drawn on a log axis');
assert(isscalar(findall(o.figs, 'Tag', 'mc_sweep_scatter_overlay')), 'the overlay belongs to one panel');
assert(~ismember('r_stable', o.stats.Properties.VariableNames), 'no sign split unless x is lambda_1');
k = strcmp(Dn.T.condition, 'sfa1_std1') & ~Dn.T.failed & isfinite(Qr.values);
S = o.stats(strcmp(o.stats.condition, 'sfa1_std1'), :);
assert(S.n == nnz(k) && abs(S.pearson_r - pearson_direct(Dn.T.gain(k), Qr.values(k))) < 1e-12);
close(o.figs);

%% The statistics recover the generator's pattern
[S, files, SB, SG] = mc_gain_bias_statistics('data', D, 'out_dir', tmp, 'verbose', 'near-none');
assert(numel(files) == 4 && all(cellfun(@isfile, files)));
assert(height(SB) == 4 * numel(D.biases) && height(SG) == 4 * numel(D.gains));
assert(all(strcmp(S.width_source([row_of(S, 'sfa1_std1'), row_of(S, 'sfa3_stdc2')]), 'fit')));
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

%% The tuning width at half the best of any condition, by hand
% A: MC 2 8 8 2 at gains 1 2 4 8; B: 10 everywhere. Half the best of any is
% 5; A crosses it half way between 1 and 2 and between 4 and 8 (linear in
% MC, log2 gain): from log2 gain 0.5 to 2.5, so 2 octaves, closed; B
% covers the whole row, 3 octaves, open.
Th = table(repmat([1; 2; 4; 8], 2, 1), zeros(8, 1), ones(8, 1), ...
    [repmat({'A'}, 4, 1); repmat({'B'}, 4, 1)], [2; 8; 8; 2; 10; 10; 10; 10], ...
    'VariableNames', {'gain', 'bias', 'seed', 'condition', 'mc_total'});
Sh = mc_gain_bias_statistics('data', Th, 'reference', 'B', 'verbose', 'near-none');
assert(abs(Sh.gain_oct_fitbias(1) - 2) < 1e-12 && ~Sh.gain_oct_fitbias_open(1));
assert(abs(Sh.gain_oct_best(2) - 3) < 1e-12 && Sh.gain_oct_best_open(2));
assert(Sh.best_any(1) == 10 && Sh.share_half_any(1) == 0.5);

%% A run in progress: its job files read back into the same table
% The signature carries what the reader takes the delays and the preset's
% point from: the protocol, the preset and its overrides (here a bias of
% 0.25 where the preset has 0.2).
job_dir = fullfile(tmp, 'partial', 'mc_gain_sweep', 'jobs'); mkdir(job_dir);
sig = struct('preset_name', mu7chain_config('essential').mc_preset, 'protocol', mc_protocol(), ...
    'model_args', {{'S_c', 0.25}});
[grp, gg] = findgroups(T(:, {'gain', 'bias', 'seed'}));
keep = 1:2:height(gg);                    % half the jobs "finished"
for k = keep
    R = T(grp == k, :);
    nr = height(R);
    job = struct('gain', R.gain(1), 'bias', R.bias(1), 'seed', R.seed(1), 'seed_stim', R.seed(1) + 1, ...
        'conditions', {R.condition'}, 'signature', sig, 'mc_total', R.mc_total', 'horizon', R.horizon', ...
        'R2_delay', R.R2_delay, 'lle_input', R.lle_input', 'lle_noinput', R.lle_noinput', ...
        'rate_E', R.rate_E', 'rate_I', R.rate_I', 'frac_pinned', R.frac_pinned', ...
        'frac_silent', R.frac_silent', 'seconds', ones(1, nr), 'failed', R.failed', ...
        'message', {repmat({''}, 1, nr)});
    save(fullfile(job_dir, sprintf('job_%03d.mat', k)), 'job');
end
fclose(fopen(fullfile(job_dir, 'job_999.mat.tmp'), 'w'));   % a save in progress: not read
before = dir(job_dir);
[Dj, info] = read_mc_gain_sweep_jobs(fullfile(tmp, 'partial'));
assert(info.n_read == numel(keep) && info.n_signatures == 1);
after = dir(job_dir);
assert(isequal({after.name}, {before.name}) && isequal([after.datenum], [before.datenum]), ...
    'the reader must not touch the job folder');
want = sortrows(D.T(ismember(grp, keep), :), {'gain', 'bias', 'seed', 'condition'});
got  = sortrows(Dj.T, {'gain', 'bias', 'seed', 'condition'});
for v = {'mc_total', 'lle_input', 'frac_silent', 'R2_delay', 'failed'}
    assert(isequaln(got.(v{1}), want.(v{1})), 'partial reader: %s differs', v{1});
end
assert(numel(Dj.delay_s) == numel(settings.delay_s) && max(abs(Dj.delay_s - settings.delay_s)) < 1e-12);
assert(Dj.preset_gain == 1 && Dj.preset_bias == 0.25, 'the preset''s point, with the signature''s overrides');
assert(isequal(Dj.conditions, four) && Dj.settings.partial);

%% One call for a run's whole set, on a result missing a condition
Dm = read_mc_gain_sweep(struct('T', T(~strcmp(T.condition, 'sfa1_std1') & T.gain <= 1, :), 'settings', settings));
ob = make_block_C_figures(Dm, fullfile(tmp, 'blockC'), 'sets', {'three', three}, 'verbose', 'near-none');
assert(isequal(ob.three.conditions, {'no_adaptation', 'sfa3_stdc2'}) && isscalar(ob.notes));
assert(all(cellfun(@isfile, ob.three.files)) && numel(ob.three.files) >= 7);
assert(any(endsWith(ob.three.files, 'MC_gain_bias_statistics.md')));

disp('PASS: test_mc_gain_bias_figures -- reader, partial reader, three figures at 3 and 4 conditions and one seed, statistics and tuning widths, the one-call set.');

function r = row_of(S, name)
r = find(strcmp(S.condition, name));
end

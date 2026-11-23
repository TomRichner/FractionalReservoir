%% test_outlier_figures.m
% The four figures on the outlier_spectra stage's output (FD-045):
% fig_outlier_plane, fig_outlier_dose, fig_outlier_slow_variables,
% fig_outlier_dwell, their one reader outlier_fig_data, and the censored
% dwell estimate censored_survival.
%
% Checks:
%   1. censored_survival against answers worked by hand: no censoring is
%      the empirical survival; a textbook case with two cut durations; a
%      longest duration that was cut leaves the estimate open (NaN beyond
%      it, quartiles not reached); everything cut; a cut and a complete
%      duration of the same length; nothing at all; and a made-up SERIES
%      through ospec_excursions, cut at both ends, whose median dwell is
%      0.5 s where counting the cut excursions as complete gives 0.3 s.
%   2. The stage ONCE, tiny: 60 neurons, 2 s of settling, a 3 s window,
%      3 calibration states, 2 seeds, three networks (reference and
%      E-to-E x 1.5 with the full-Jacobian passes; E-to-E x 1.25 without),
%      the three conditions.
%   3. Each figure draws from it and writes its PNG and its _table.md;
%      takes its options (seed, networks, conditions, outlier_factor, a
%      given window with a break, inset off; spread; x_thr and a given
%      window; band); reads the index file given as data_file; leaves the
%      numbers it reports equal to what the run files hold (the marked
%      eigenvalues per state counted again in units of W against each
%      state's own radius; the dose values against ospec_measures; the
%      dwell counts against ospec_excursions).
%   4. Each errors on an EMPTY run directory and on no run directory
%      (resolve_data_file: no fallback), on a network the index lacks, on
%      a seed it lacks; the slow-variables figure on a network without
%      dense spectra.
%   5. THE WHOLE-WINDOW CASE of the dwell figure, on a copy of the stage's
%      folder in which one condition's frozen-gain series is beyond zero
%      throughout in every seed: the curve is a level line at 1 that ends
%      at the window's length (no drop to zero there), and the table names
%      the condition.
%
% Asserts; prints one PASS line. About a minute. Assumes setup_paths has
% run and a pool allocation is set (pool_set(0) will do).
%
% NOT TESTED HERE: anything at 500 neurons (limits, marker sizes and the
% inset were looked at on the tiny run only), the figures through
% make_all_paper_figures (wired by FD-046), and how the pictures look.
%
% See also: outlier_fig_data, censored_survival, test_outlier_spectra_stage

%% 1. The censored dwell estimate
K = censored_survival([1 2 2 4]);
assert(isequal(K.t', [1 2 4]) && max(abs(K.S' - [0.75 0.25 0])) < 1e-12 && ~K.open && K.n == 4 && K.n_censored == 0);
assert(isequal(K.at([0.5 1 2 3 5]), [1 0.75 0.25 0.25 0]) && K.quantile(0.5) == 2 && K.quantile(0.25) == 1 && K.quantile(0.75) == 2);

K = censored_survival([1 2 3 4 5], logical([0 1 0 1 0]));
assert(isequal(K.t', [1 3 5]) && max(abs(K.S' - [4/5, 4/5 * 2/3, 0])) < 1e-12 && ~K.open && K.n_censored == 2, 'the textbook case');
assert(K.quantile(0.5) == 5 && abs(K.at(3.5) - 8/15) < 1e-12);

K = censored_survival([1 3 5], logical([0 0 1]));       % the longest was cut: open
assert(max(abs(K.S' - [2/3 1/3])) < 1e-12 && K.open && K.t_max == 5);
assert(abs(K.at(5) - 1/3) < 1e-12 && isnan(K.at(6)) && K.quantile(0.5) == 3 && isnan(K.quantile(0.75)), 'open beyond the longest');

K = censored_survival(20, true);                        % beyond zero for the whole window
assert(isempty(K.t) && K.open && isequaln(K.at([0.01 10 20 21]), [1 1 1 NaN]) && isnan(K.quantile(0.5)), 'everything cut');

K = censored_survival([2 2], logical([0 1]));           % a cut one of the same length counts as still going
assert(isequal(K.t, 2) && K.S == 0.5 && K.open);

K = censored_survival([], []);
assert(K.n == 0 && isnan(K.t_max) && all(isnan(K.at([0.1 1]))) && isnan(K.quantile(0.5)), 'no excursions');
assert(throws(@() censored_survival([1 2], true), 'censored_survival:Size'));

% A series: beyond zero for 3 samples at the start (cut), 2, 2, 5, and 4 at the end (cut); step 0.1 s.
x = [1 1 1 -1 1 1 -1 -1 1 1 -1 1 1 1 1 1 -1 1 1 1 1];
E = ospec_excursions(5 + 0.1 * (0:numel(x) - 1), x);
assert(max(abs(E.duration - [0.3 0.2 0.2 0.5 0.4])) < 1e-12 && isequal(E.censored, logical([1 0 0 0 1])));
K = censored_survival(E.duration, E.censored);
assert(max(abs(K.t' - [0.2 0.5])) < 1e-12 && max(abs(K.S' - [0.6 0])) < 1e-12 && ~K.open, 'the series: S');
assert(abs(K.quantile(0.5) - 0.5) < 1e-12 && abs(median(E.duration) - 0.3) < 1e-12, ...
    'the censored median is 0.5 s; counting the cut excursions as complete gives 0.3 s');

%% 2. The stage, once, tiny
tmp = tempname;
cleanup = onCleanup(@() rmdir_if(tmp));
cfg = tiny_cfg();
t0 = tic;
run_outlier_spectra_stage(cfg, tmp, outlier_spectra_stage_plan(cfg));
t_stage = toc(t0);
runs = ospec_read_index(tmp);
assert(height(runs) == 18 && all(runs.ok));
fig_dir = fullfile(tmp, 'figs');
common = {'run_dir', tmp, 'visible', false, 'verbose', 'near-none'};
open_before = findall(groot, 'Type', 'figure');

%% 3a. fig_outlier_plane
o = draw(@fig_outlier_plane, common, fullfile(fig_dir, 'plane'), 'Fig_Outlier_Plane');
assert(isequal(size(o.numbers), [2 3]));
% The marked eigenvalues per state, counted again from the run file in units of W.
s = o.numbers{2, 3};                                    % E-to-E x 1.5, multiple timescales
assert(isequal(s.seeds, [7 8; 1 2]), 'the seeds in the order of the index');
for k = 1:2
    R = ospec_load_run(tmp, struct('network', 'EEx1p5', 'condition', 'sfa3_stdc2', 'seed', s.seeds(k, :)));
    n_marked = mean(sum(abs(R.cal.ev_M) > 1.03 .* R.cal.R_eff, 1));
    assert(abs(s.n_marked(k) - n_marked) < 1e-12, 'marked against each state''s OWN radius, seed %d', k);
    assert(abs(s.rad(k) - median(R.cal.R_eff) / R.tau_d) < 1e-9 && s.kappa(k) == R.kappa);
    is_E = R.cal.cell_type == 1;
    assert(abs(s.E_slope(k) - mean(R.cal.slope(is_E, :) > 0, 'all')) < 1e-12 && ...
        abs(s.E_slope(k) + s.E_floor(k) + s.E_ceil(k) - 1) < 1e-12 && abs(s.I_slope(k) + s.I_floor(k) + s.I_ceil(k) - 1) < 1e-12);
end
insets = findall(o.figs, 'Tag', 'fig_outlier_plane_inset');
assert(numel(insets) == 6, 'one inset per panel');
dots = findall(insets(1), 'Type', 'scatter');
assert(numel(dots.XData) == 60 * 3, 'one dot per neuron per state in the inset');
% Polish of 2026-11-22: no eigenvalue dot lies on a tick label (a tick whose
% label would be covered in any panel is left out of every panel).
assert(isempty(intersect(o.window.re_ticks, o.window.hidden_re_ticks)) && ...
    isempty(intersect(o.window.im_ticks, o.window.hidden_im_ticks)), 'a hidden tick is not among the ticks drawn');
txt = fileread(fullfile(fig_dir, 'plane', 'Fig_Outlier_Plane_table.md'));
assert(contains(txt, 'TICKS LEFT OUT') ~= (isempty(o.window.hidden_re_ticks) && isempty(o.window.hidden_im_ticks)), ...
    'the sidecar says which ticks were left out, or that none was');
n_covered = 0;
panels = setdiff(findall(o.figs, 'Type', 'axes'), insets);
for a = panels'
    sc = findall(a, 'Type', 'scatter');
    x = [sc.XData]; y = [sc.YData];
    for t = findall(a, 'Type', 'text', 'Tag', 'complex_plane_axes')'
        if iscell(t.String) || isempty(regexp(t.String, ['^[\d.' char(8722) ']+$'], 'once')); continue; end   % tick labels only
        e = t.Extent; pad = 0.15 * e(4);                   % the glyphs are a little inside the text box
        n_covered = n_covered + nnz(x > e(1) & x < e(1) + e(3) & y > e(2) + pad & y < e(2) + e(4) - pad);
    end
end
assert(n_covered == 0, '%d eigenvalue dot(s) lie on a tick label', n_covered);
close(o.figs);

o2 = draw(@fig_outlier_plane, [common, {'seed', [1 2], 'outlier_factor', 1.3, 'inset', false, 'networks', {'EEx1p5'}, ...
    'conditions', {'sfa3_stdc2', 'no_adaptation'}}], fullfile(fig_dir, 'plane_opts'), 'Fig_Outlier_Plane');
assert(isequal(size(o2.numbers), [1 2]) && isempty(findall(o2.figs, 'Tag', 'fig_outlier_plane_inset')));
assert(all(o2.numbers{1, 1}.n_marked <= s.n_marked), 'a larger factor marks no more');
assert(contains(fileread(fullfile(fig_dir, 'plane_opts', 'Fig_Outlier_Plane_table.md')), 'network seed [1 2]'));
close(o2.figs);
o3 = draw(@fig_outlier_plane, [common, {'re_lim', [-200 60], 'im_lim', [-150 150], 'far_lim', [100 140], 'tick_step', 50}], ...
    fullfile(fig_dir, 'plane_window'), 'Fig_Outlier_Plane');
assert(isequal(o3.window.re_lim, [-200 60]) && isequal(o3.window.far_lim, [100 140]) && isequal(o3.window.y_lim, [-150 150]), ...
    'a given window with a break is taken');
assert(contains(fileread(fullfile(fig_dir, 'plane_window', 'Fig_Outlier_Plane_table.md')), 'THE REAL AXIS IS BROKEN'));
close(o3.figs);
o4 = fig_outlier_plane('data_file', fullfile(tmp, 'outlier_spectra', 'ospec_index.mat'), 'save', false, 'visible', false);
assert(isempty(o4.files) && isequaln(o4.numbers, o.numbers), 'the index file given as data_file is the same run');
close(o4.figs);
assert(throws(@() fig_outlier_plane(common{:}, 'save', false, 'networks', {'reference', 'EEx2p0'}), 'fig_outlier_plane:NoNetwork'));
assert(throws(@() fig_outlier_plane(common{:}, 'save', false, 'seed', 99), 'fig_outlier_plane:NoSeed'));
assert(throws(@() fig_outlier_plane(common{:}, 'save', false, 're_lim', [-200 5], 'im_lim', [-150 150]), 'fig_outlier_plane:MarkOutside'), ...
    'a diamond outside a given window is an error, not a silent clip');

%% 3b. fig_outlier_dose
o = draw(@fig_outlier_dose, common, fullfile(fig_dir, 'dose'), 'Fig_Outlier_Dose');
assert(isequal(sort(o.networks), {'EEx1p25', 'EEx1p5', 'reference'}) && issorted(o.x), 'the levels, placed by their structural ratio');
assert(isequal(size(o.values.ratio_theory_med), [3 3 2]));
w = find(strcmp(o.networks, 'EEx1p25'));
R = ospec_load_run(tmp, struct('network', 'EEx1p25', 'condition', 'sfa1_std1', 'seed', [1 2]));
m = ospec_measures(R);
assert(o.values.ratio_theory_med(w, 2, 2) == m.row.ratio_theory_med && o.values.frac_on_E(w, 2, 2) == m.row.frac_on_E && ...
    o.values.gain_ratio_EI(w, 2, 2) == m.row.gain_ratio_EI && o.values.structural_ratio(w, 2, 2) == m.row.structural_ratio, ...
    'the dose values are ospec_measures''');
assert(abs(o.x(w) - median(squeeze(o.values.structural_ratio(w, 1, :)))) < 1e-12, 'each level at the median across seeds');
n_thin = numel(findall(o.figs, 'Type', 'line', 'LineWidth', 0.75));
assert(n_thin == 3 * 3 * 2, 'one thin line per panel, condition and seed (found %d)', n_thin);
ax = findall(o.figs, 'Type', 'axes');
yl = arrayfun(@(a) a.YLim(2), ax);
assert(max(yl) >= 1.3 * max(o.diamond) - 1e-9, 'headroom above the top diamond');
close(o.figs);
o2 = draw(@fig_outlier_dose, [common, {'spread', 'band', 'conditions', {'sfa3_stdc2'}}], fullfile(fig_dir, 'dose_band'), 'Fig_Outlier_Dose');
assert(isempty(findall(o2.figs, 'Type', 'line', 'LineWidth', 0.75)) && numel(findall(o2.figs, 'Type', 'patch')) == 3 && ...
    size(o2.values.ratio_theory_med, 2) == 1);
close(o2.figs);
o2 = fig_outlier_dose(common{:}, 'save', false, 'spread', 'none');
assert(isempty(findall(o2.figs, 'Type', 'patch')) && isempty(findall(o2.figs, 'Type', 'line', 'LineWidth', 0.75)));
close(o2.figs);

%% 3c. fig_outlier_slow_variables
o = draw(@fig_outlier_slow_variables, common, fullfile(fig_dir, 'slow'), 'Fig_Outlier_Slow_Variables');
s = o.numbers{2, 3};
R1 = ospec_load_run(tmp, struct('network', 'EEx1p5', 'condition', 'sfa3_stdc2', 'seed', [7 8]));
R2 = ospec_load_run(tmp, struct('network', 'EEx1p5', 'condition', 'sfa3_stdc2', 'seed', [1 2]));
assert(s.n_seeds == 2 && isequal(s.cont_ov, [median([R1.dense.cont_ov]), median([R2.dense.cont_ov])]) && ...
    isequal(s.cont_ov_all, [R1.dense.cont_ov, R2.dense.cont_ov]), 'the continuation''s overlap, across seeds');
s0 = o.numbers{1, 1};                                   % no adaptation: the full Jacobian is the frozen one
assert(max(abs(s0.cont_ov - 1)) < 1e-6 && max(abs(s0.cont_re - s0.frozen_re)) < 1e-3, 'without slow variables the two coincide');
assert(o.window.re_lim(1) < 0 && o.window.re_lim(2) > 0 && all(o.window.re_ticks > 0) && isempty(o.window.far_lim), ...
    'through the origin, no break, real ticks right of zero only');
% The common canvas: canvas_w wide unless a composite passes fig_w; a tick
% whose label a mark would lie on is left out, never the last of its axis.
sty = manuscript_style(); o.figs.Units = 'inches';
assert(abs(o.figs.Position(3) - sty.canvas_w) < 1e-6, 'the stand-alone figure is canvas_w wide (%.2f in)', o.figs.Position(3));
assert(isempty(intersect(o.window.re_ticks, o.window.hidden_re_ticks)) && ~isempty(o.window.re_ticks) && ...
    isempty(intersect(o.window.im_ticks, o.window.hidden_im_ticks)) && ~isempty(o.window.im_ticks), ...
    'ticks under a mark are left out; one stays on each axis');
close(o.figs);
o2 = fig_outlier_slow_variables(common{:}, 'save', false, 'visible', false, 'fig_w', 6);
o2.figs.Units = 'inches';
assert(abs(o2.figs.Position(3) - 6) < 1e-6, 'fig_w sets the width');
close(o2.figs);
o2 = draw(@fig_outlier_slow_variables, [common, {'seed', 1, 'x_thr', 0.3, 're_lim', [-14 10], 'im_lim', [-9 9], ...
    'networks', {'reference'}, 'conditions', {'sfa1_std1'}}], fullfile(fig_dir, 'slow_opts'), 'Fig_Outlier_Slow_Variables');
assert(isequal(o2.window.re_lim, [-14 10]) && isequal(o2.window.y_lim, [-9 9]) && isequal(size(o2.numbers), [1 1]));
assert(contains(fileread(fullfile(fig_dir, 'slow_opts', 'Fig_Outlier_Slow_Variables_table.md')), 'at least 0.3 of its squared norm'));
close(o2.figs);
assert(throws(@() fig_outlier_slow_variables(common{:}, 'save', false, 'networks', {'reference', 'EEx1p25'}), ...
    'fig_outlier_slow_variables:NeedsDense'));

%% 3d. fig_outlier_dwell
o = draw(@fig_outlier_dwell, common, fullfile(fig_dir, 'dwell'), 'Fig_Outlier_Dwell');
assert(isequal(size(o.curves), [2 3 2]));
q = o.curves{1, 3, 1};                                  % reference, multiple timescales, frozen gains
Ra = ospec_load_run(tmp, struct('network', 'reference', 'condition', 'sfa3_stdc2', 'seed', [7 8]));
Ea = ospec_excursions(Ra.t, Ra.rate_top);
Rb = ospec_load_run(tmp, struct('network', 'reference', 'condition', 'sfa3_stdc2', 'seed', [1 2]));
Eb = ospec_excursions(Rb.t, Rb.rate_top);
assert(q.n_seeds == 2 && q.n_exc == numel(Ea.duration) + numel(Eb.duration) && ...
    q.n_cut == nnz(Ea.censored) + nnz(Eb.censored) && isequal(q.fraction, [Ea.fraction, Eb.fraction]) && ...
    abs(q.dt - 0.01) < 1e-12, 'the dwell counts are ospec_excursions''');
qf = o.curves{1, 3, 2};
Ef = ospec_excursions(Rb.full.t, real(Rb.full.net));
assert(abs(qf.dt - 0.1) < 1e-12 && qf.fraction(2) == Ef.fraction, 'the full-Jacobian series at its own step');
for i = 1:numel(o.curves)
    q = o.curves{i};
    if isempty(q) || isempty(q.x); continue; end
    assert(q.x(1) == q.dt && q.med(1) == 1 && all(q.med >= 0 & q.med <= 1) && all(q.lo <= q.med + 1e-12) && ...
        all(q.hi >= q.med - 1e-12) && issorted(q.x), 'curve %d: a survival curve from one sample on', i);
end
close(o.figs);
o2 = draw(@fig_outlier_dwell, [common, {'band', 'range', 'networks', {'EEx1p5'}, 'conditions', {'sfa1_std1'}}], ...
    fullfile(fig_dir, 'dwell_opts'), 'Fig_Outlier_Dwell');
assert(isequal(size(o2.curves), [1 1 2]));
close(o2.figs);
o2 = fig_outlier_dwell(common{:}, 'save', false, 'band', 'none');
assert(isempty(findall(o2.figs, 'Type', 'patch')));
close(o2.figs);
% A network without the full-Jacobian series draws its frozen-gain curves alone.
o2 = draw(@fig_outlier_dwell, [common, {'networks', {'EEx1p25'}}], fullfile(fig_dir, 'dwell_nofull'), 'Fig_Outlier_Dwell');
assert(~isempty(o2.curves{1, 1, 1}) && isempty(o2.curves{1, 1, 2}) && ...
    contains(fileread(fullfile(fig_dir, 'dwell_nofull', 'Fig_Outlier_Dwell_table.md')), 'No full-Jacobian series'));
close(o2.figs);

%% 4. No run, no figure
empty_dir = tempname; mkdir(empty_dir);
cleanup_empty = onCleanup(@() rmdir_if(empty_dir));
figs = {@fig_outlier_plane, @fig_outlier_dose, @fig_outlier_slow_variables, @fig_outlier_dwell};
for k = 1:numel(figs)
    assert(throws(@() figs{k}('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:NotFound'), ...
        'figure %d on an empty run directory', k);
    assert(throws(@() figs{k}('save', false, 'visible', false), 'resolve_data_file:NoRunDir'), 'figure %d with no run directory', k);
end

%% 5. Beyond zero for the whole window in every seed
tmp2 = tempname;
cleanup2 = onCleanup(@() rmdir_if(tmp2));
copyfile(fullfile(tmp, 'outlier_spectra'), fullfile(tmp2, 'outlier_spectra'));
for k = find(strcmp(runs.network, 'reference') & strcmp(runs.condition, 'sfa1_std1'))'
    f = fullfile(tmp2, 'outlier_spectra', runs.file{k});
    L = load(f, 'R'); R = L.R;
    R.rate_top = abs(R.rate_top) + 0.5;                 % made up: beyond zero at every sample
    save(f, 'R');
end
o = draw(@fig_outlier_dwell, {'run_dir', tmp2, 'visible', false, 'verbose', 'near-none'}, fullfile(tmp2, 'figs'), 'Fig_Outlier_Dwell');
q = o.curves{1, 2, 1};
assert(q.n_exc == 2 && q.n_cut_both == 2 && all(q.fraction == 1), 'one excursion per seed, cut at both ends');
assert(abs(q.x(end) - q.window_s) < 1e-9 && numel(q.x) == 2 && isequal(q.med, [1 1]) && all(isnan(q.pooled_q)), ...
    'a level line at 1 up to the window''s length: no step down, no quartile');
txt = fileread(fullfile(tmp2, 'figs', 'Fig_Outlier_Dwell_table.md'));
assert(contains(txt, 'BEYOND ZERO FOR THE WHOLE WINDOW IN EVERY SEED') && contains(txt, 'reference, sfa1_std1, frozen gains'), ...
    'the table names the condition');
tri = findall(o.figs, 'Type', 'line', 'Marker', '>');
assert(any(arrayfun(@(h) isscalar(h.YData) && h.YData == 1 && abs(h.XData - q.window_s) < 1e-9, tri)), ...
    'the triangle at the end of the level line');
close(o.figs);

open_after = findall(groot, 'Type', 'figure');
assert(numel(open_after) == numel(open_before), 'the test closed every figure it drew, and an erroring figure left none open');

fprintf('PASS: test_outlier_figures -- censored_survival by hand, the four figures on one tiny stage run (%.0f s), their options, their numbers against the run files, no run no figure, the whole-window dwell case.\n', t_stage);

%% ------------------------------------------------------------------------
function cfg = tiny_cfg()
% The bundle's essential config with the stage shrunk to seconds.
cfg = mu7chain_config('essential');
cfg.verbose = 'near-none';
cfg.outlier_spectra_n_override = 60;
cfg.outlier_spectra_networks = {'reference', 'EEx1p25', 'EEx1p5'};
cfg.outlier_spectra_full_networks = {'reference', 'EEx1p5'};
cfg.outlier_spectra_conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
cfg.outlier_spectra_n_seeds = 2;
cfg.outlier_spectra_T_settle = 2;
cfg.outlier_spectra_window_s = 3;
cfg.outlier_spectra_full_window_s = 3;
cfg.outlier_spectra_lya_window_s = 1.5;
cfg.outlier_spectra_n_cal = 3;
cfg.outlier_spectra_n_null = 4;
cfg.outlier_spectra_n_shuffle = 3;
cfg.outlier_spectra_chunk_s = 0.7;
end

function o = draw(fn, args, out_dir, tag)
% Draw and save one figure; assert its PNG and its sidecar were written.
o = fn(args{:}, 'out_dir', out_dir);
assert(isscalar(o.figs) && isgraphics(o.figs, 'figure') && strcmp(o.figs.Visible, 'off'), '%s: the figure', tag);
assert(any(strcmp(o.files, [tag '.png'])) && any(strcmp(o.files, [tag '_table.md'])), '%s: the files', tag);
png = dir(fullfile(out_dir, [tag '.png']));
txt = fileread(fullfile(out_dir, [tag '_table.md']));
assert(~isempty(png) && png.bytes > 5000, '%s: the PNG', tag);
assert(startsWith(txt, ['# ' tag]) && contains(txt, '|---|') && numel(txt) > 800 && ~contains(txt, 'NaN'), '%s: the sidecar', tag);
assert(endsWith(o.source, 'ospec_index.mat'), '%s: the source', tag);
end

function tf = throws(fn, id)
try
    fn();
    tf = false;
catch err
    tf = strcmp(err.identifier, id);
    if ~tf; fprintf(2, 'expected %s, got %s: %s\n', id, err.identifier, err.message); end
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end

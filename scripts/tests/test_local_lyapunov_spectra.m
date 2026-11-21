% test_local_lyapunov_spectra.m - run_local_lyapunov's 'spectra_n_samples'
% (FD-034): the Jacobian's spectra and the projected Jacobian of the stage's
% own top-K pass at states of its switching run, the figures drawn from them,
% and the opt-in config switch. A 40-neuron network, 10 s, two input steps.
%
% Checks:
%   1. Off by default: no local_lyapunov_spectra.mat, and settings has no
%      spectra_n_samples.
%   2. On (10 states): results (top-K and Benettin, input, indices) are
%      identical to the run without it, timings apart; the stage's file
%      differs only by settings.spectra_n_samples.
%   3. The spectra file: 10 states at 1, 2, ..., 10 s, five per input step;
%      per condition N eigenvalues per state, H of size K x K x states (raw
%      norm only), every state valid, the pass's exponents equal to the
%      stage's bit for bit, real parts of the Ritz values summing to
%      trace(H), and the diagonal of H close to the pass's interval rates
%      (correlation above 0.9; NaN only after the last state).
%   4. fig_eig_heatmap with 'spectra_source', 'local_lyapunov' and the
%      overlay: saved as fig_eig_heatmap_switch_ritz_inset with a table that
%      carries the shared-pass section; fig_grouped_main group 3 saves
%      Fig_Main3_Grouped_switch_ritz_inset with one density image and one
%      inset per panel. On a run without the spectra both error through
%      resolve_data_file, and the default source is untouched by the option.
%      The fourth row ('real_histogram_row'): one histogram and one dotted
%      reference per column, under the plot box of the heat map above it and
%      on its x limits, no text on it, a taller canvas of the same width, a
%      _table.md with the row's caption sentence; fig_eig_real_histogram
%      with 'spectra_source', 'local_lyapunov'; the gains the stage saves
%      (disk_by_cond). With the lambda_1 row it errors :HistogramRow.
%      'spectra_start_s': a later first state gives the same H, exponents
%      and gains at the shared states, bit for bit.
%   5. with_switch_spectra: the mu7chain bundles carry it by default since
%      FD-046 ('essential' and 'fast' are exactly the config without it, put
%      back by undo_switch below, with with_switch_spectra(cfg, 150, 40, 1,
%      'as_default', true) applied: fig_main3_grouped the four-row
%      switching composite, the baseline composite beside it). As the
%      opt-in (on the config without it): it sets K = 150 (100 until 2026-11-21) and 40 states (or
%      those given) and adds exactly four entries, each after the entry it
%      is a version of, all 'appendix' (none in_paper), every other entry unchanged; it
%      names the first state only when the accumulation window opens
%      before it; twice errors :AlreadyAdded, a config without the entries
%      :MissingEntry.
%
% ~1.5 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: run_local_lyapunov, fig_eig_heatmap, fig_main3_stability, with_switch_spectra

fprintf('=== Testing run_local_lyapunov spectra_n_samples (FD-034) ===\n\n');
all_passed = true;

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_steps5s_noise0p025_chainStd_4cond_mu7chain';
C = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
K = 12; n_s = 10;
base = fullfile(tempdir, 'local_lyapunov_spectra_test');
if exist(base, 'dir'); rmdir(base, 's'); end
off_dir = fullfile(base, 'off'); on_dir = fullfile(base, 'on');
args = {'preset_name', P, 'run_mode', 'essential', 'K', K, 'T', 10, 'step_duration_s', 5, ...
    'accumulation_start_s', 1, 'warmup_s', 3, 'simulation_start_s', -3, 'n_override', 40, 'conditions', C, ...
    'use_parallel', false}; %#ok<NASGU>   serial: the test needs no pool
evalc('f_off = run_local_lyapunov(args{:}, ''out_dir'', fullfile(off_dir, ''local_lyapunov''));');
evalc('f_on = run_local_lyapunov(args{:}, ''out_dir'', fullfile(on_dir, ''local_lyapunov''), ''spectra_n_samples'', n_s);');
A = load(f_off); B = load(f_on);
spectra_file = fullfile(on_dir, 'local_lyapunov', 'local_lyapunov_spectra.mat');

all_passed = check('off by default: no spectra file, no settings.spectra_n_samples', ...
    ~isfile(fullfile(off_dir, 'local_lyapunov', 'local_lyapunov_spectra.mat')) && ~isfield(A.settings, 'spectra_n_samples')) && all_passed;

strip = @(r) rmfield(r, 'seconds');
strip_topk = @(r) arrayfun(@(x) setfield(setfield(x, 'topk', rmfield(x.topk, 'seconds')), 'ben', drop_seconds(x.ben)), r); %#ok<SFLD>
all_passed = check('on: the stage''s results are identical to the run without it (timings apart)', ...
    isequaln(strip_topk(strip(A.results)), strip_topk(strip(B.results))) && isequal(A.cond_names, B.cond_names) && ...
    isequaln(rmfield(A.settings, 'minutes'), rmfield(B.settings, {'minutes', 'spectra_n_samples'})) && ...
    B.settings.spectra_n_samples == n_s) && all_passed;

ok = isfile(spectra_file);
if ok
    D = load(spectra_file);
    ok = isequal(D.cond_names, C) && max(abs(D.state_times(:)' - (1:10))) < 1e-9 && ...
        isequal(D.state_step(:)', [1 1 1 1 1 2 2 2 2 2]) && D.lle_window == 9;
    worst_trace = 0; worst_corr = 1;
    for i = 1:numel(C)
        R = D.ritz_by_cond{i}; tk = B.results(i).topk;
        dg = cell2mat(arrayfun(@(s) diag(R.H(:, :, s))', (1:n_s)', 'UniformOutput', false));
        z = cell2mat(arrayfun(@(s) eig(R.H(:, :, s)), 1:n_s, 'UniformOutput', false));
        worst_trace = max(worst_trace, max(abs(sum(real(z), 1) - sum(dg, 2)')));
        avg = mean(cat(3, R.local_rates_before, R.local_rates_after), 3, 'omitnan');
        cc = corrcoef(dg(:), avg(:)); worst_corr = min(worst_corr, cc(1, 2));
        ok = ok && numel(D.evals_by_cond{i}) == R.N * n_s && isequal(size(R.H), [K K n_s]) && all(R.valid) && ...
            isequal(R.scalings, {'raw'}) && isequal(R.LE_spectrum, tk.LE_spectrum) && isequal(R.sort_idx(:), tk.sort_idx(:)) && ...
            all(isfinite(R.local_rates_before(:))) && all(all(isfinite(R.local_rates_after(1:end - 1, :)))) && ...
            all(isnan(R.local_rates_after(end, :))) && D.lle_by_cond(i) == tk.LLE && ...
            isequal(sort(R.local_rates_before(3, :), 'descend'), sort(tk.local_LE_spectrum_t(find(tk.t_lya < 3 - 1e-6, 1, 'last'), :), 'descend'));
    end
    ok = ok && worst_trace < 1e-9 && worst_corr > 0.9;
    fprintf('    sum Re(Ritz) - trace(H): %.2g; smallest correlation of diag(H) with the interval rates: %.3f\n', worst_trace, worst_corr);
end
all_passed = check('the spectra file: states, steps, H of the stage''s own pass, trace identity, diag(H) near the interval rates', ok) && all_passed;

evalc(['h = fig_eig_heatmap(''run_dir'', on_dir, ''out_dir'', fullfile(on_dir, ''heat''), ''visible'', false, ' ...
    '''spectra_source'', ''local_lyapunov'', ''ritz_overlay'', true, ''ritz_inset'', true);']);
tbl_file = fullfile(on_dir, 'heat', 'fig_eig_heatmap_switch_ritz_inset_table.md');
ok = any(strcmp(h.files, 'fig_eig_heatmap_switch_ritz_inset.png')) && isfile(tbl_file) && ...
    n_tagged(h.figs, 'ritz_density') == numel(C) && n_tagged(h.figs, 'ritz_inset') == numel(C);
if ok
    tbl = fileread(tbl_file);
    ok = contains(tbl, '## The shared pass') && contains(tbl, 'By input step') && contains(tbl, 'switching-input run');
end
close(h.figs);
evalc(['g = fig_grouped_main(3, ''run_dir'', on_dir, ''out_dir'', fullfile(on_dir, ''g3''), ''show_lambda_row'', false, ' ...
    '''spectra_source'', ''local_lyapunov'', ''ritz_overlay'', true, ''ritz_inset'', true);']);
ok = ok && any(strcmp(g.files, 'Fig_Main3_Grouped_switch_ritz_inset.png')) && ...
    n_tagged(g.figs, 'ritz_density') == numel(C) && n_tagged(g.figs, 'ritz_inset') == numel(C);
g_pos = g.figs.Position;
close(g.figs);
all_passed = check('figures from the switching run: tags with _switch, the overlay and inset per panel, the shared-pass table', ok) && all_passed;

% The fourth row (the real-axis histogram) and the stand-alone histogram.
evalc(['g4 = fig_grouped_main(3, ''run_dir'', on_dir, ''out_dir'', fullfile(on_dir, ''g4''), ''show_lambda_row'', false, ' ...
    '''spectra_source'', ''local_lyapunov'', ''ritz_overlay'', true, ''ritz_inset'', true, ''real_histogram_row'', true);']);
tag4 = 'Fig_Main3_Grouped_switch_ritz_inset_realhist';   % no _noC since the merge with FD-031: no lambda_1 row is the default
ok = any(strcmp(g4.files, [tag4 '.png'])) && isfile(fullfile(on_dir, 'g4', [tag4 '_table.md'])) && ...
    n_tagged(g4.figs, 'eig_real_histogram') == numel(C) && n_tagged(g4.figs, 'uniform_disk_reference') == numel(C) && ...
    g4.figs.Position(3) == g_pos(3) && g4.figs.Position(4) > g_pos(4) && numel(findall(g4.figs, 'Tag', 'row_letter')) == 4;
for c = 1:numel(C)
    heat = findall(g4.figs, 'Tag', sprintf('stability_r3_c%d', c)); hx = findall(g4.figs, 'Tag', sprintf('stability_r4_c%d', c));
    ok = ok && isscalar(heat) && isscalar(hx);
    if ok
        pb = plot_box_position(heat);
        ok = isequal(hx.XLim, heat.XLim) && max(abs(hx.Position([1 3]) - pb([1 3]))) < 1e-6 && ...
            hx.Position(2) + hx.Position(4) < pb(2) && isempty(findobj(hx, 'Type', 'text'));
    end
end
if ok
    tbl = fileread(fullfile(on_dir, 'g4', [tag4 '_table.md']));
    ok = contains(tbl, 'Row (D)') && contains(tbl, 'radius drawn: `gain`') && contains(tbl, 'switching-input run');
end
close(g4.figs);
evalc(['rh = fig_eig_real_histogram(''run_dir'', on_dir, ''out_dir'', fullfile(on_dir, ''rh''), ''visible'', false, ' ...
    '''spectra_source'', ''local_lyapunov'');']);
ok = ok && any(strcmp(rh.files, 'fig_eig_real_histogram_switch.png')) && strcmp(rh.reference(1).mode, 'gain') && ...
    isfile(fullfile(on_dir, 'rh', 'fig_eig_real_histogram_switch_table.md'));
close(rh.figs);
G = D.disk_by_cond{end};
ok = ok && numel(D.disk_by_cond) == numel(C) && G.n == 40 && numel(G.rms_gain) == n_s && ...
    all(G.rms_gain > 0 & G.rms_gain < 1) && max(abs(G.t(:)' - (1:10))) < 1e-9;
all_passed = check('the histogram row: (D) under each heat map''s plot box on its x limits, no text on it, table written; the stand-alone histogram; the saved gains', ok) && all_passed;
all_passed = check('the histogram row needs the layout without the lambda_1 row', ...
    throws_prefix(@() fig_main3_stability(on_dir, {}, 'real_histogram_row', true, 'show_lambda_row', true), ...
    'fig_main3_stability:HistogramRow')) && all_passed;   % named: no lambda_1 row is the default since the merge with FD-031

% A first state after the accumulation start: the same pass, the same H.
late_dir = fullfile(base, 'late');
evalc(['run_local_lyapunov(args{:}, ''out_dir'', fullfile(late_dir, ''local_lyapunov''), ''spectra_n_samples'', n_s - 1, ' ...
    '''spectra_start_s'', 2);']);
L = load(fullfile(late_dir, 'local_lyapunov', 'local_lyapunov_spectra.mat'));
ok = max(abs(L.state_times(:)' - (2:10))) < 1e-9 && L.settings.spectra_start_s == 2 && ~isfield(D.settings, 'spectra_start_s');
for i = 1:numel(C)
    R = D.ritz_by_cond{i}; RL = L.ritz_by_cond{i};
    ok = ok && isequal(RL.H, R.H(:, :, 2:end)) && isequal(RL.LE_spectrum, R.LE_spectrum) && isequal(RL.interval, R.interval) && ...
        RL.warmup_s == R.warmup_s && isequal(L.evals_by_cond{i}, D.evals_by_cond{i}(R.N + 1:end)) && ...
        isequal(L.disk_by_cond{i}.rms_gain, D.disk_by_cond{i}.rms_gain(2:end));
end
all_passed = check('spectra_start_s: later first state, the pass still over the accumulation window (H, exponents and gains bit for bit)', ok) && all_passed;
all_passed = check('a run without the spectra errors through resolve_data_file', ...
    throws_prefix(@() fig_eig_heatmap('run_dir', off_dir, 'save', false, 'visible', false, 'spectra_source', 'local_lyapunov'), ...
    'resolve_data_file:')) && all_passed;
rmdir(base, 's');

%% The config switch: the bundles' default (FD-046), and the opt-in
ok = true;
for mode = {'essential', 'fast'}
    c = mu7chain_config(mode{1});
    ok = ok && isequal(with_switch_spectra(undo_switch(c), 150, 40, 1, 'as_default', true), c) && ...
        c.local_lyapunov_K == 150 && c.local_lyapunov_spectra_n_samples == 40 && c.local_lyapunov_spectra_start_s == 1;
end
all_passed = check('mu7chain_config carries the switching Figure 3 by default (K = 150, 40 states from 1 s)', ok) && all_passed;
base_cfg = undo_switch(mu7chain_config('essential'));
cfg = with_switch_spectra(base_cfg);
names0 = cellfun(@(f) f.name, base_cfg.figures, 'UniformOutput', false);
names  = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
added  = {'fig_eig_heatmap_switch', 'fig_main3_grouped_switch', 'fig_main3_grouped_switch_hist', 'fig_eig_real_histogram_switch'};
c2 = with_switch_spectra(base_cfg, 30, 8);
ih = find(strcmp(names, 'fig_eig_heatmap')); im = find(strcmp(names, 'fig_main3_grouped'));
extra = {'spectra_source', 'local_lyapunov', 'ritz_overlay', true, 'ritz_inset', true};
h0 = base_cfg.figures{strcmp(names0, 'fig_eig_heatmap')}; m0 = base_cfg.figures{strcmp(names0, 'fig_main3_grouped')};
all_passed = check('with_switch_spectra: K = 150 and 40 states by default, given values kept', ...
    cfg.local_lyapunov_K == 150 && cfg.local_lyapunov_spectra_n_samples == 40 && ...
    c2.local_lyapunov_K == 30 && c2.local_lyapunov_spectra_n_samples == 8) && all_passed;
ir = find(strcmp(names, 'fig_eig_real_histogram')); r0 = base_cfg.figures{strcmp(names0, 'fig_eig_real_histogram')};
all_passed = check('four entries added after the entries they are versions of, all appendix, the rest unchanged', ...
    numel(names) == numel(names0) + 4 && strcmp(names{ih + 1}, added{1}) && strcmp(names{im + 1}, added{2}) && ...
    strcmp(names{im + 2}, added{3}) && strcmp(names{ir + 1}, added{4}) && ...
    isequal(cfg.figures(~ismember(names, added)), base_cfg.figures) && ...
    isequal(cfg.figures{ih + 1}.args, [h0.args, extra]) && isequal(cfg.figures{ih + 1}.fn, h0.fn) && ...
    isequal(cfg.figures{im + 1}.args, [m0.args, extra, {'show_lambda_row', false}]) && isequal(cfg.figures{im + 1}.fn, m0.fn) && ...
    isequal(cfg.figures{im + 2}.args, [m0.args, extra, {'show_lambda_row', false, 'real_histogram_row', true}]) && ...
    isequal(cfg.figures{im + 2}.fn, m0.fn) && ...
    isequal(cfg.figures{ir + 1}.args, [r0.args, {'spectra_source', 'local_lyapunov'}]) && isequal(cfg.figures{ir + 1}.fn, r0.fn) && ...
    ~any(cellfun(@(f) f.in_paper, cfg.figures(ismember(names, added)))) && ...
    isequal(rmfield(drop_start(cfg), {'figures', 'local_lyapunov_K', 'local_lyapunov_spectra_n_samples'}), ...
            rmfield(base_cfg, {'figures', 'local_lyapunov_K'}))) && all_passed;
all_passed = check('as_default: the paper''s name on the four-row composite, the baseline beside it', ...
    isequal(with_switch_spectra(base_cfg, 150, 40, 1, 'as_default', true), mu7chain_config('essential'))) && all_passed;
% The first state: left to the stage when the window opens at or after it, named when it opens before.
% (Whatever the bundle's own accumulation start is.)
early = base_cfg; early.local_lyapunov_accumulation_start_s = 0;
late = base_cfg; late.local_lyapunov_accumulation_start_s = 1;
c3 = with_switch_spectra(early); c4 = with_switch_spectra(early, 150, 40, 2); c5 = with_switch_spectra(late);
all_passed = check('first state: not set for a window from 1 s, 1 s (or the value given) for a window from 0 s', ...
    ~isfield(c5, 'local_lyapunov_spectra_start_s') && ...
    c3.local_lyapunov_spectra_start_s == 1 && c4.local_lyapunov_spectra_start_s == 2) && all_passed;
bare = base_cfg; bare.figures = base_cfg.figures(~strcmp(names0, 'fig_eig_heatmap'));
all_passed = check('applied twice errors :AlreadyAdded; without fig_eig_heatmap :MissingEntry', ...
    throws_prefix(@() with_switch_spectra(cfg), 'with_switch_spectra:AlreadyAdded') && ...
    throws_prefix(@() with_switch_spectra(bare), 'with_switch_spectra:MissingEntry')) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL local_lyapunov spectra TESTS PASSED ===\n');
else
    fprintf('=== SOME local_lyapunov spectra TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function b = drop_seconds(b)
if isfield(b, 'seconds'); b = rmfield(b, 'seconds'); end
end

function c = undo_switch(c)
% A mu7chain config as it was before with_switch_spectra became its default:
% without the switching entries, fig_main3_grouped the baseline composite,
% and no spectra settings (K stays the bundle's, 150 since 2026-11-21).
names = cellfun(@(f) f.name, c.figures, 'UniformOutput', false);
base = c.figures{strcmp(names, 'fig_main3_grouped_baseline')};
base.name = 'fig_main3_grouped';
base.in_paper = true;
base.group = 'main';
% as_default names the baseline's file (fig_tag, last in its args); before it, the entry had none.
assert(isequal(base.args(end - 1:end), {'fig_tag', 'Fig_Main3_Grouped_baseline'}));
base.args = base.args(1:end - 2);
c.figures{strcmp(names, 'fig_main3_grouped')} = base;
c.figures = c.figures(~ismember(names, {'fig_main3_grouped_baseline', 'fig_eig_heatmap_switch', ...
    'fig_eig_real_histogram_switch'}));
c = rmfield(c, {'local_lyapunov_spectra_n_samples', 'local_lyapunov_spectra_start_s'});
end

function c = drop_start(c)
if isfield(c, 'local_lyapunov_spectra_start_s'); c = rmfield(c, 'local_lyapunov_spectra_start_s'); end
end

function n = n_tagged(fig, tag)
n = numel(findall(fig, 'Tag', tag));
end

function ok = throws_prefix(fn, prefix)
ok = false;
try
    fn();
catch ME
    ok = startsWith(ME.identifier, prefix);
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

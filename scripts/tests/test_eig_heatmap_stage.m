% test_eig_heatmap_stage.m - run_eig_heatmap with the E:I imbalance examples,
% fig_eig_heatmap and fig_eig_heatmap_imbalance, end to end on a small network.
%
% Checks:
%   1. A bad run mode errors :badMode naming every mode (also in test_run_modes).
%   2. 'fast' with n_override 60 (so the examples run at 60 too), 5 states
%      per example, into a temp run directory: the .mat keeps every legacy
%      variable (evals_by_cond, lle_by_cond, lya_by_cond, num/spec abscissa,
%      J_times_by_cond, lle_window, lya_T_interval, settings, cond_names,
%      condition_titles) AND carries `examples`: 3 entries (inhibition-
%      dominant, reference, excitation-dominant) x 3 conditions with finite
%      lambda_1, mean rate in [0, 1] and B_E in [0, 1]; the reference entry's
%      eigenvalues equal the legacy variable; the two imbalanced examples
%      carry mu_EE_relative at 0.5x and 1.5x the preset's value and B_E
%      ordered inhibition < reference < excitation.
%   3. fig_eig_heatmap and fig_eig_heatmap_imbalance write files and the
%      latter its table; the imbalance figure errors :NoExamples on a .mat
%      without examples, and both error via resolve_data_file on a run
%      directory without the stage.
%
%   4. THE RITZ PASS (FD-033). Off by default: the .mat holds exactly the
%      variables it always did and `examples` has no ritz field. With
%      'ritz_K', 10 (reference example only, same seeds and states): every
%      legacy variable is identical to the run without it; ritz_by_cond
%      holds, per condition, H of size K x K x states x 2 scalings with
%      K = min(10, state dimension), every state valid and at the times of
%      J_times_by_cond, values = eig(H) with real parts summing to trace(H),
%      block fractions summing to 1, and the scales by kind of variable
%      (x alone without adaptation, where the scaled pass must then give the
%      raw Ritz values); and for the diagnostics JQ_gram (residuals in
%      [0, 1]) and the full Jacobian's spectral abscissa, count of
%      eigenvalues right of zero and numerical abscissa per state, the last
%      bounding the projected one. A scaling name it does not know errors
%      :badRitzScaling.
%   5. fig_eig_heatmap with 'ritz_overlay', in each ritz_style: the default
%      'density' (an image tagged ritz_density per panel), 'marginal',
%      'fine' and 'scatter' (points tagged ritz_values), with the exponent
%      ticks; no text on the panels and one title line unless the
%      annotation and the tick label are asked for; saved under the tag
%      fig_eig_heatmap_ritz (plus the style, scaling and inset) with its
%      table, which carries the diagnostics; fig_ritz_by_state draws; without the option the figure has no such object and the tag
%      is unchanged; on a run directory without the pass it draws as before
%      and errors :NoRitz when the overlay is asked for.
%   7. fig_eig_real_histogram (Re of the eigenvalues against the projection
%      of a uniform disk): draws from the run with the saved gains (radius
%      from the rms gain), falls back to the unit-gain disk on a run without
%      them, and its reference sums to n per state over a disk inside the
%      window.
%   6. The composite Figure 3 (fig_main3_stability, fig_grouped_main group
%      3) on a small local_lyapunov run beside the Ritz run: with the
%      overlay every heat-map panel carries the points and its inset and the
%      file is Fig_Main3_Grouped_ritz_inset; without it there are none and
%      the tag is Fig_Main3_Grouped.
%
% ~3 min. Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: run_eig_heatmap, fig_eig_heatmap, fig_eig_heatmap_imbalance

fprintf('=== Testing run_eig_heatmap + figures (imbalance examples) ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';

all_passed = check('bad run mode errors :badMode', ...
    throws_id(@() run_eig_heatmap('run_mode', 'nope'), 'run_eig_heatmap:badMode')) && all_passed;

run_dir = fullfile(tempdir, 'eig_heatmap_stage_test');
if exist(run_dir, 'dir'); rmdir(run_dir, 's'); end
mkdir(run_dir);
t0 = tic;
out = evalc(['mat_file = run_eig_heatmap(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(run_dir, ''eig_heatmap''), ''n_override'', 60, ''n_samples'', 5, ' ...
    '''n_samples_examples'', 5, ''use_parallel'', false);']);
fprintf('  (stage ran in %.0f s)\n', toc(t0));
all_passed = check('stage returns the .mat path and it exists', ischar(mat_file) && isfile(mat_file)) && all_passed;
D = load(mat_file);
legacy = {'evals_by_cond', 'lle_by_cond', 'lya_by_cond', 'num_abscissa_by_cond', ...
    'spec_abscissa_by_cond', 'J_times_by_cond', 'lle_window', 'lya_T_interval', ...
    'settings', 'cond_names', 'condition_titles'};
all_passed = check('legacy variables all present', all(isfield(D, legacy))) && all_passed;
all_passed = check('no Ritz pass by default: the variables are the legacy ones and examples, no ritz field', ...
    isequal(sort(fieldnames(D)), sort([legacy, {'examples'}]')) && ~isfield(D.examples, 'ritz_by_cond') && ...
    ~isfield(D.settings, 'ritz_K')) && all_passed;
all_passed = check('examples present: 3 entries, labels in order', isfield(D, 'examples') && ...
    numel(D.examples) == 3 && isequal({D.examples.label}, {'inhibition-dominant', 'reference', 'excitation-dominant'})) && all_passed;
E = D.examples;
n_cond = numel(D.cond_names);
fin = true; rng_ok = true; sz = true;
for e = 1:numel(E)
    fin = fin && all(isfinite(E(e).lle_by_cond)) && all(isfinite(E(e).mean_rate_by_cond)) && all(isfinite(E(e).B_E_by_cond));
    rng_ok = rng_ok && all(E(e).mean_rate_by_cond >= 0 & E(e).mean_rate_by_cond <= 1) && ...
        all(E(e).B_E_by_cond >= 0 & E(e).B_E_by_cond <= 1);
    sz = sz && numel(E(e).evals_by_cond) == n_cond && numel(E(e).lle_by_cond) == n_cond && ...
        all(cellfun(@(v) ~isempty(v), E(e).evals_by_cond)) && E(e).n == 60;
    fprintf('  %-20s mu_EE %6.3f  lambda_1 %s  <r> %s  B_E %s\n', E(e).label, E(e).mu_tilde_relative(1, 1), ...
        mat2str(E(e).lle_by_cond, 3), mat2str(E(e).mean_rate_by_cond, 3), mat2str(E(e).B_E_by_cond, 3));
end
all_passed = check('every example x condition: finite lambda_1, mean rate, B_E', fin) && all_passed;
all_passed = check('mean rate and B_E in [0, 1]; n = 60; eigenvalues pooled', rng_ok && sz) && all_passed;
i_ref = find(strcmp({E.label}, 'reference'));
all_passed = check('reference example equals the legacy variables', ...
    isequal(E(i_ref).evals_by_cond, D.evals_by_cond) && isequal(E(i_ref).lle_by_cond, D.lle_by_cond)) && all_passed;
d = srnn_param_preset(P);
mu_ref = d.mu_tilde_relative(1, 1);
i_inh = find(strcmp({E.label}, 'inhibition-dominant')); i_exc = find(strcmp({E.label}, 'excitation-dominant'));
all_passed = check('imbalanced examples carry mu_EE at 0.5x and 1.5x the preset', ...
    abs(E(i_inh).overrides.mu_EE_relative - 0.5 * mu_ref) < 1e-12 && abs(E(i_exc).overrides.mu_EE_relative - 1.5 * mu_ref) < 1e-12 && ...
    abs(E(i_inh).mu_tilde_relative(1, 1) - 0.5 * mu_ref) < 1e-9 && abs(E(i_exc).mu_tilde_relative(1, 1) - 1.5 * mu_ref) < 1e-9 && ...
    abs(E(i_ref).mu_tilde_relative(1, 1) - mu_ref) < 1e-9) && all_passed;
all_passed = check('B_E ordered inhibition < reference < excitation (every condition)', ...
    all(E(i_inh).B_E_by_cond < E(i_ref).B_E_by_cond) && all(E(i_ref).B_E_by_cond < E(i_exc).B_E_by_cond)) && all_passed;
all_passed = check('settings record the example labels and the override', ...
    isequal(D.settings.examples, {E.label}) && D.settings.n_samples_examples == 5 && D.settings.n_examples_override == 60) && all_passed;

%% The Ritz pass (FD-033)
ritz_dir = fullfile(tempdir, 'eig_heatmap_stage_test_ritz');
if exist(ritz_dir, 'dir'); rmdir(ritz_dir, 's'); end
mkdir(ritz_dir);
ref_only = struct('label', 'reference', 'overrides', struct()); %#ok<NASGU>  used inside evalc
t0 = tic;
evalc(['mat_ritz = run_eig_heatmap(''preset_name'', P, ''run_mode'', ''fast'', ' ...
    '''out_dir'', fullfile(ritz_dir, ''eig_heatmap''), ''n_override'', 60, ''n_samples'', 5, ' ...
    '''examples'', ref_only, ''use_parallel'', false, ''ritz_K'', 10);']);
fprintf('  (stage with the Ritz pass, reference only, ran in %.0f s)\n', toc(t0));
DR = load(mat_ritz);
same = true;
for k = 1:numel(legacy)
    if ~ismember(legacy{k}, {'settings'}); same = same && isequaln(DR.(legacy{k}), D.(legacy{k})); end
end
all_passed = check('with the Ritz pass every legacy variable is identical to the run without it', same) && all_passed;
all_passed = check('settings record ritz_K and the scalings; the reference example carries the same structs', ...
    DR.settings.ritz_K == 10 && isequal(DR.settings.ritz_scalings, {'raw', 'type_sd'}) && ...
    isequaln(DR.examples(1).ritz_by_cond, DR.ritz_by_cond)) && all_passed;
shape = true; ident = true; kinds_ok = true; worst_trace = 0; diag_ok = true;
for i = 1:n_cond
    R = DR.ritz_by_cond{i};
    n_s = numel(DR.J_times_by_cond{i});
    shape = shape && R.K == min(10, R.N) && isequal(size(R.H), [R.K, R.K, n_s, 2]) && ...
        isequal(size(R.values), [R.K, n_s, 2]) && all(R.valid) && isequal(R.t(:), DR.J_times_by_cond{i}(:)) && ...
        all(isfinite(R.H), 'all') && all(isfinite(R.LE_spectrum), 'all');
    for sc = 1:2
        for k = 1:n_s
            worst_trace = max(worst_trace, abs(sum(real(R.values(:, k, sc))) - trace(R.H(:, :, k, sc))));
        end
    end
    ident = ident && max(abs(sum(R.block_fractions, 2) - 1), [], 'all') < 1e-12;
    for sc = 1:2
        for k = 1:n_s
            Dg = ritz_diagnostics(R.H(:, :, k, sc), R.JQ_gram(:, :, k, sc));
            diag_ok = diag_ok && Dg.residual >= 0 && Dg.residual <= 1 && all(Dg.pair_residual >= 0 & Dg.pair_residual <= 1 + 1e-9) && ...
                Dg.max_re <= Dg.num_abscissa + 1e-10 && Dg.num_abscissa <= R.full_num_abscissa(k, sc) + 1e-6;
        end
    end
    diag_ok = diag_ok && isequal(size(R.full_max_re), [n_s, 1]) && all(R.full_n_positive >= 0) && ...
        max(R.full_max_re) == max(real(DR.evals_by_cond{i}));
    if strcmp(DR.cond_names{i}, 'no_adaptation')
        d_raw = sort(real(R.values(:, :, 1))); d_sc = sort(real(R.values(:, :, 2)));
        kinds_ok = kinds_ok && isequal(R.scales.kinds, {'x'}) && max(abs(d_raw - d_sc), [], 'all') < 1e-8;
    else
        kinds_ok = kinds_ok && numel(R.scales.kinds) > 1 && all(R.scales.sd > 0);
    end
    fprintf('  %-14s N = %3d, K = %2d: lambda_1 %+.3f (raw) %+.3f (type_sd); x share of the subspace %.2f / %.2f; kinds %s\n', ...
        DR.cond_names{i}, R.N, R.K, R.LE_spectrum(1, 1), R.LE_spectrum(1, 2), ...
        mean(R.block_fractions(:, 1, 1)), mean(R.block_fractions(:, 1, 2)), strjoin(R.scales.kinds, ' '));
end
all_passed = check('per condition: H is K x K x states x 2, all states valid, at the sampled times', shape) && all_passed;
all_passed = check(sprintf('real parts of the Ritz values sum to trace(H) (worst %.1e)', worst_trace), worst_trace < 1e-9) && all_passed;
all_passed = check('block fractions sum to 1', ident) && all_passed;
all_passed = check(['diagnostics: residuals in [0, 1]; max Re Ritz <= num. abscissa of H <= that of the full J in the same norm; ' ...
    'the full spectral abscissa per state is that of the pooled spectra'], diag_ok) && all_passed;
all_passed = check('scales by kind: x alone without adaptation (scaled = raw there), several kinds with it', kinds_ok) && all_passed;
all_passed = check('an unknown scaling errors :badRitzScaling', ...
    throws_id(@() run_eig_heatmap('preset_name', P, 'run_mode', 'fast', 'ritz_K', 5, 'ritz_scalings', {'zscore'}), ...
    'run_eig_heatmap:badRitzScaling')) && all_passed;

%% Figures
evalc('res1 = fig_eig_heatmap(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig1''), ''visible'', false);');
all_passed = check(sprintf('fig_eig_heatmap writes files (%d)', numel(res1.files)), numel(res1.files) >= 1) && all_passed;
close(res1.figs);
evalc('res2 = fig_eig_heatmap_imbalance(''run_dir'', run_dir, ''out_dir'', fullfile(run_dir, ''fig2''), ''visible'', false);');
all_passed = check(sprintf('fig_eig_heatmap_imbalance writes files (%d) and the table', numel(res2.files)), ...
    numel(res2.files) >= 1 && isfile(fullfile(run_dir, 'fig2', 'Fig_Eig_Heatmap_Imbalance_table.md'))) && all_passed;
close(res2.figs);
old_dir = fullfile(tempdir, 'eig_heatmap_stage_test_old'); if ~exist(old_dir, 'dir'); mkdir(fullfile(old_dir, 'eig_heatmap')); end
Dold = rmfield(D, 'examples'); save(fullfile(old_dir, 'eig_heatmap', 'eig_heatmap_data.mat'), '-struct', 'Dold');
all_passed = check('imbalance figure errors :NoExamples on a pre-examples .mat', ...
    throws_id(@() fig_eig_heatmap_imbalance('run_dir', old_dir, 'save', false, 'visible', false), 'fig_eig_heatmap_imbalance:NoExamples')) && all_passed;

%% The Ritz overlay (FD-033)
n_tagged = @(fh, tag) numel(findall(fh, 'Tag', tag));
evalc('res0 = fig_eig_heatmap(''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''fig0''), ''visible'', false);');
all_passed = check('without the option: no Ritz object, tag unchanged, no table', ...
    n_tagged(res0.figs, 'ritz_values') == 0 && n_tagged(res0.figs, 'ritz_density') == 0 && n_tagged(res0.figs, 'ritz_inset') == 0 && ...
    any(strcmp(res0.files, 'fig_eig_heatmap.png')) && isempty(dir(fullfile(ritz_dir, 'fig0', '*_table.md')))) && all_passed;
close(res0.figs);
evalc('res3 = fig_eig_heatmap(''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''fig3''), ''visible'', false, ''ritz_overlay'', true, ''ritz_style'', ''scatter'', ''ritz_annotate'', true, ''ritz_tick_label'', true);');
sc3 = findall(res3.figs, 'Tag', 'ritz_values');
all_passed = check('overlay, scatter style, text asked for: one scatter of K x states points, exponent ticks, their label and the annotation per panel, no inset', ...
    numel(sc3) == n_cond && n_tagged(res3.figs, 'ritz_exponent_ticks') == n_cond && n_tagged(res3.figs, 'ritz_inset') == 0 && ...
    n_tagged(res3.figs, 'ritz_tick_label') == n_cond && n_tagged(res3.figs, 'ritz_annotation') == n_cond && ...
    sum(arrayfun(@(h) numel(h.XData), sc3)) == sum(cellfun(@(R) R.K * nnz(R.valid), DR.ritz_by_cond))) && all_passed;
all_passed = check('overlay, scatter style: saved as fig_eig_heatmap_ritz_scatter with its table', ...
    any(strcmp(res3.files, 'fig_eig_heatmap_ritz_scatter.png')) && isfile(fullfile(ritz_dir, 'fig3', 'fig_eig_heatmap_ritz_scatter_table.md'))) && all_passed;
close(res3.figs);
ok5 = true;
for sty = {'density', 'marginal'}
    evalc(['res5 = fig_eig_heatmap(''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''fig5'', sty{1}), ''visible'', false, ' ...
        '''ritz_overlay'', true, ''ritz_style'', sty{1});']);
    tl5 = findall(res5.figs, 'Type', 'tiledlayout');
    ok5 = ok5 && n_tagged(res5.figs, ['ritz_' sty{1}]) == n_cond && n_tagged(res5.figs, 'ritz_annotation') == 0 && ...
        n_tagged(res5.figs, 'ritz_tick_label') == 0 && isscalar(cellstr(tl5.Title.String));
    if strcmp(sty{1}, 'density')
        tbl = fileread(fullfile(ritz_dir, 'fig5', 'density', 'fig_eig_heatmap_ritz_table.md'));
        ok5 = ok5 && contains(tbl, 'Caption paragraph') && any(strcmp(res5.files, 'fig_eig_heatmap_ritz.png')) && n_tagged(res5.figs, 'ritz_values') == 0 && ...
            contains(tbl, 'What the points are: diagnostics') && contains(tbl, 'Re > 0 per state');
    else
        ok5 = ok5 && any(strcmp(res5.files, 'fig_eig_heatmap_ritz_marginal.png'));
    end
    close(res5.figs);
end
all_passed = check('overlay, density (the default) and marginal styles: a density image or a histogram per panel, no text on the panels, one title line; the table has the diagnostics and a caption', ok5) && all_passed;
evalc('bs = fig_ritz_by_state(''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''bs''), ''visible'', false);');
all_passed = check('fig_ritz_by_state draws three panels per condition and saves', ...
    numel(findall(bs.figs, 'Type', 'axes')) == 3 * n_cond && any(strcmp(bs.files, 'fig_ritz_by_state.png'))) && all_passed;
close(bs.figs);
evalc(['res4 = fig_eig_heatmap(''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''fig4''), ''visible'', false, ' ...
    '''ritz_overlay'', true, ''ritz_inset'', true, ''ritz_scaling'', ''type_sd'', ''ritz_k'', 4, ''ritz_exponent_ticks'', false, ''ritz_style'', ''fine'');']);
sc4 = findall(res4.figs, 'Tag', 'ritz_values');
all_passed = check('overlay with inset, fine style, scaled, k = 4: 4 x states points per panel, an inset per panel, no ticks, own tag', ...
    sum(arrayfun(@(h) numel(h.XData), sc4)) == sum(cellfun(@(R) 4 * nnz(R.valid), DR.ritz_by_cond)) && ...
    n_tagged(res4.figs, 'ritz_inset') == n_cond && n_tagged(res4.figs, 'ritz_exponent_ticks') == 0 && ...
    any(strcmp(res4.files, 'fig_eig_heatmap_ritz_fine_typesd_inset.png')) && ...
    isfile(fullfile(ritz_dir, 'fig4', 'fig_eig_heatmap_ritz_fine_typesd_inset_table.md'))) && all_passed;
close(res4.figs);
all_passed = check('a run directory without the Ritz pass: the overlay errors :NoRitz', ...
    throws_id(@() fig_eig_heatmap('run_dir', run_dir, 'save', false, 'visible', false, 'ritz_overlay', true), ...
    'fig_eig_heatmap:NoRitz')) && all_passed;

%% The real-axis histogram against a uniform disk
evalc('rh = fig_eig_real_histogram(''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''rh''), ''visible'', false);');
rh_children = get(findall(rh.figs, 'Type', 'axes'), 'Children'); rh_children = vertcat(rh_children{:});   % no text on the panels
all_passed = check('fig_eig_real_histogram: a histogram and one dotted reference per panel, rms-gain radius from the saved gains, no text on the panels, table written', ...
    n_tagged(rh.figs, 'eig_real_histogram') == n_cond && n_tagged(rh.figs, 'uniform_disk_reference') == n_cond && ...
    all(strcmp({rh.reference.mode}, 'gain')) && all([rh.reference.rho] > 0) && all([rh.reference.rho] < [rh.reference.rho_unit]) && ...
    isempty(findobj(rh_children, 'flat', 'Type', 'text')) && ...
    isfile(fullfile(ritz_dir, 'rh', 'fig_eig_real_histogram_table.md'))) && all_passed;
close(rh.figs);
evalc('rh0 = fig_eig_real_histogram(''run_dir'', run_dir, ''save'', false, ''visible'', false, ''count_scale'', ''clipped'');');
all_passed = check('a run without saved gains: no error, the unit-gain disk is drawn (R / tau_d)', ...
    all(strcmp({rh0.reference.mode}, 'unit')) && all([rh0.reference.rho] == [rh0.reference.rho_unit]) && ...
    all(isnan([rh0.reference.rho_gain]))) && all_passed;
close(rh0.figs);
evalc('rh3 = fig_eig_real_histogram(''run_dir'', ritz_dir, ''save'', false, ''visible'', false, ''disk_radius'', 3);');
r3 = rh3.reference(1);
all_passed = check('a disk inside the window: the reference sums to n per state over its bins', ...
    rh3.edges(1) < r3.c - 3 && rh3.edges(end) > r3.c + 3 && abs(sum(r3.expected) - r3.n) < 1e-9 * r3.n && ...
    all(r3.expected(rh3.edges(2:end) < r3.c - 3) == 0)) && all_passed;
close(rh3.figs);
all_passed = check('with the gains saved, the stage records n, tau_d, R and one rms gain per state', ...
    isfield(DR, 'disk_by_cond') && ~isfield(D, 'disk_by_cond') && numel(DR.disk_by_cond{1}.rms_gain) == numel(DR.J_times_by_cond{1}) && ...
    DR.disk_by_cond{1}.n == 60 && all(DR.disk_by_cond{end}.rms_gain > 0 & DR.disk_by_cond{end}.rms_gain < 1)) && all_passed;

%% The composite Figure 3 with the overlay
evalc(['run_local_lyapunov(''preset_name'', P, ''run_mode'', ''fast'', ''out_dir'', fullfile(ritz_dir, ''local_lyapunov''), ' ...
    '''K'', 5, ''T'', 8, ''n_override'', 60);']);
evalc('m3 = fig_main3_stability(ritz_dir, {}, ''ritz_overlay'', true, ''ritz_inset'', true);');
panels = findall(m3.figs, 'Type', 'axes', '-regexp', 'Tag', '^stability_r1_c\d+$');
in_panel = arrayfun(@(ax) numel(findall(ax, 'Tag', 'ritz_density')), panels);
all_passed = check('fig_main3_stability: every heat-map panel carries the points; one inset per panel', ...
    numel(panels) == n_cond && all(in_panel == 1) && n_tagged(m3.figs, 'ritz_inset') == n_cond) && all_passed;
close(m3.figs);
evalc('m0 = fig_main3_stability(ritz_dir, {});');
all_passed = check('fig_main3_stability without the option: no Ritz object', ...
    n_tagged(m0.figs, 'ritz_density') == 0 && n_tagged(m0.figs, 'ritz_inset') == 0) && all_passed;
close(m0.figs);
evalc(['g3 = fig_grouped_main(3, ''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''g3''), ' ...
    '''ritz_overlay'', true, ''ritz_inset'', true);']);
evalc('g0 = fig_grouped_main(3, ''run_dir'', ritz_dir, ''out_dir'', fullfile(ritz_dir, ''g0''));');
all_passed = check('fig_grouped_main group 3: Fig_Main3_Grouped_ritz_inset with the points; Fig_Main3_Grouped without', ...
    any(strcmp(g3.files, 'Fig_Main3_Grouped_ritz_inset.png')) && n_tagged(g3.figs, 'ritz_density') == n_cond && ...
    n_tagged(g3.figs, 'ritz_inset') == n_cond && any(strcmp(g0.files, 'Fig_Main3_Grouped.png')) && ...
    n_tagged(g0.figs, 'ritz_density') == 0) && all_passed;
close([g3.figs, g0.figs]);

empty_dir = fullfile(tempdir, 'eig_heatmap_stage_test_empty'); if ~exist(empty_dir, 'dir'); mkdir(empty_dir); end
all_passed = check('both figures error on a run directory without the stage', ...
    throws_prefix(@() fig_eig_heatmap('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:') && ...
    throws_prefix(@() fig_eig_heatmap_imbalance('run_dir', empty_dir, 'save', false, 'visible', false), 'resolve_data_file:')) && all_passed;
rmdir(empty_dir, 's'); rmdir(old_dir, 's'); rmdir(run_dir, 's'); rmdir(ritz_dir, 's');

fprintf('\n');
if all_passed
    fprintf('=== ALL eig_heatmap stage TESTS PASSED ===\n');
else
    fprintf('=== SOME eig_heatmap stage TESTS FAILED ===\n');
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

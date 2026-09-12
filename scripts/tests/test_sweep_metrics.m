% test_sweep_metrics.m - The metric registry, the per-job Lyapunov fields the
% sweeps now store, and the plotters that read them.
%
% Checks:
%   1. sweep_metrics is self-consistent: unique keys and fields, every field
%      is a stored result field (lya_summary_fields or mean_rate /
%      mean_synaptic_output), lookup by key and by field, unknown -> error,
%      the four sheet measures are lle / r / hks / dky.
%   2. analysis_run_config carries the top-K settings for every analysis and
%      mode (lya_method 'topk', K 15, auto, cap 30, lya_dt 0.05), and a
%      config can override lya_method to 'benettin' after the fact.
%   3. A tiny Pairs sweep (2 gain levels x 2 reps, n = 30, two conditions)
%      through the real driver: every result carries every registry field
%      and every lya_summary field as a scalar; LLE, h_KS_bits, D_KY,
%      K_used, n_positive finite on successful jobs; K_used in {15, 30};
%      D_KY_resolved 0/1; the leading local series present;
%      collect_level_values works for h_KS_bits and D_KY.
%   4. The same sweep with lya_method = 'benettin' stores NaN for the
%      spectrum-derived fields and finite LLE and transient scalars.
%   5. psa.plot, plot_sensitivity and plot_unit_histograms run for every
%      in_sheets measure, with figure Names '<field> Distribution',
%      '<field> Sensitivity - <param>', '<field> Unit Histogram'.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: sweep_metrics, ParamSpaceAnalysis2, SRNNCellTypePairs.lya_summary

fprintf('=== Testing sweep_metrics and the per-job Lyapunov fields ===\n\n');
all_passed = true;

%% 1. Registry
M = sweep_metrics();
keys_ = {M.key}; fields_ = {M.field};
all_passed = check('unique keys and fields', numel(unique(keys_)) == numel(keys_) && numel(unique(fields_)) == numel(fields_)) && all_passed;
known = [SRNNCellTypePairs.lya_summary_fields(), {'mean_rate', 'mean_synaptic_output'}];
all_passed = check('every registry field is a stored result field', all(ismember(fields_, known))) && all_passed;
all_passed = check('lookup by key and by field agree', isequal(sweep_metrics('hks'), sweep_metrics('h_KS_bits'))) && all_passed;
all_passed = check('unknown metric errors', throws_id(@() sweep_metrics('nope'), 'sweep_metrics:UnknownMetric')) && all_passed;
sheets = {M([M.in_sheets]).key};
all_passed = check(sprintf('sheet measures are lle, r, hks, dky (%s)', strjoin(sheets, ', ')), ...
    isequal(sort(sheets), sort({'lle', 'r', 'hks', 'dky'}))) && all_passed;

%% 2. analysis_run_config
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
preset = srnn_param_preset(P);
ok = true;
for an = {'sensitivity', 'tau_sensitivity', 'param_space'}
    for md = run_mode_names()
        c = analysis_run_config(an{1}, md{1}, preset);
        ok = ok && strcmp(c.model.lya_method, 'topk') && c.model.lya_K == 15 && ...
            c.model.lya_K_auto && c.model.lya_K_max == 30 && c.model.lya_dt == 0.05;
    end
end
all_passed = check('every analysis x mode carries topk / K 15 / auto / cap 30 / lya_dt 0.05', ok) && all_passed;
c = analysis_run_config('sensitivity', 'fast', preset); c.model.lya_method = 'benettin';
all_passed = check('a config can override lya_method afterwards', strcmp(c.model.lya_method, 'benettin')) && all_passed;

%% 3. A tiny sweep through the driver
[d, cls, conds] = srnn_param_preset(P);
cfg = analysis_run_config('sensitivity', 'fast', d);
d.n = 30; d.indegree = 6; d.F_tracks_network = true; d.sigma_u_noise = 0;
cfg.model.ode_solver = 'rk4'; cfg.model.T_range = [0 8]; cfg.model.lya_T_interval = [3 8];
tmp = fullfile(tempdir, 'psa_sweep_metrics_test');
if exist(tmp, 'dir'); rmdir(tmp, 's'); end
psa = make_psa(cls, merge_struct(d, cfg.model), conds([1 3]), tmp, 'topk');
evalc('psa.run();');
res = [psa.results.no_adaptation; psa.results.sfa3_std2];
res = res(~cellfun(@isempty, res));
succ = cellfun(@(r) r.success, res);
all_passed = check(sprintf('%d of %d jobs succeeded', nnz(succ), numel(res)), all(succ)) && all_passed;
need = [SRNNCellTypePairs.lya_summary_fields(), fields_];
has_all = all(cellfun(@(r) all(cellfun(@(f) isfield(r, f) && isscalar(r.(f)), need)), res));
all_passed = check('every result carries every registry and lya_summary field as a scalar', has_all) && all_passed;
fin = all(cellfun(@(r) all(isfinite([r.LLE r.h_KS_bits r.n_positive r.K_used r.mean_rate])), res));
all_passed = check('LLE, h_KS_bits, n_positive, K_used, mean_rate finite', fin) && all_passed;
ku = cellfun(@(r) r.K_used, res);
all_passed = check(sprintf('K_used in {15, 30} (%s)', mat2str(unique(ku)')), all(ismember(ku, [15 30]))) && all_passed;
resolved = cellfun(@(r) r.D_KY_resolved, res);
all_passed = check('D_KY_resolved is 0/1 and D_KY finite where resolved', all(ismember(resolved, [0 1])) && ...
    all(cellfun(@(r) ~r.D_KY_resolved || isfinite(r.D_KY), res))) && all_passed;
all_passed = check('leading local series stored per job', all(cellfun(@(r) ~isempty(r.local_rate_lead) && numel(r.local_rate_lead) == numel(r.t_lya_lead), res))) && all_passed;
all_passed = check('block fractions sum to 1 per job', all(cellfun(@(r) abs(r.lead_frac_x + r.lead_frac_sfa + r.lead_frac_std + r.lead_frac_stf - 1) < 1e-9, res))) && all_passed;
v_h = ParamSpaceAnalysis2.collect_level_values(psa, 'level_of_chaos', 1, 'no_adaptation', 'h_KS_bits');
v_d = ParamSpaceAnalysis2.collect_level_values(psa, 'level_of_chaos', 2, 'sfa3_std2', 'D_KY');
all_passed = check('collect_level_values works for h_KS_bits and D_KY', numel(v_h) == 2 && numel(v_d) <= 2) && all_passed;
fprintf('  no_adaptation: LLE %s, h_KS %s bit/s, D_KY %s, K_used %s\n', ...
    mat2str(cellfun(@(r) r.LLE, psa.results.no_adaptation)', 3), mat2str(cellfun(@(r) r.h_KS_bits, psa.results.no_adaptation)', 3), ...
    mat2str(cellfun(@(r) r.D_KY, psa.results.no_adaptation)', 3), mat2str(cellfun(@(r) r.K_used, psa.results.no_adaptation)'));

%% 4. Benettin comparison run on the same networks
tmpb = fullfile(tempdir, 'psa_sweep_metrics_test_benettin');
if exist(tmpb, 'dir'); rmdir(tmpb, 's'); end
psb = make_psa(cls, merge_struct(d, cfg.model), conds([1 3]), tmpb, 'benettin');
evalc('psb.run();');
rb = psb.results.no_adaptation; rb = rb(~cellfun(@isempty, rb));
all_passed = check('benettin: LLE and transient scalars finite, spectrum fields NaN', ...
    all(cellfun(@(r) isfinite(r.LLE) && isfinite(r.frac_local_positive) && isnan(r.h_KS_bits) && isnan(r.D_KY) && isnan(r.K_used), rb))) && all_passed;
la = cellfun(@(r) r.LLE, psa.results.no_adaptation); lb = cellfun(@(r) r.LLE, psb.results.no_adaptation);
fprintf('  paired lambda_1, topk vs benettin, no_adaptation: %s vs %s\n', mat2str(la', 3), mat2str(lb', 3));
all_passed = check('topk and benettin lambda_1 agree within 0.3 on every job (same networks)', max(abs(la - lb)) < 0.3) && all_passed;

%% 5. Plotters for every sheet measure
names_ok = true;
for spec = M([M.in_sheets])
    evalc('psa.plot(''metric'', spec.field);');
    f = findobj(0, 'Type', 'figure', 'Name', sprintf('%s Distribution', spec.field));
    names_ok = names_ok && ~isempty(f); close(f);
    evalc('psa.plot_sensitivity(''metric'', spec.field);');
    f = findobj(0, 'Type', 'figure', '-regexp', 'Name', sprintf('^%s Sensitivity - ', spec.field));
    names_ok = names_ok && ~isempty(f); close(f);
end
all_passed = check('plot and plot_sensitivity run for every sheet measure with the Name convention', names_ok) && all_passed;
evalc('psa.plot_unit_histograms(''Metrics'', {M([M.in_sheets]).key}, ''color_by'', ''level_of_chaos'');');
uh_ok = true;
for spec = M([M.in_sheets])
    f = findobj(0, 'Type', 'figure', 'Name', sprintf('%s Unit Histogram', spec.field));
    uh_ok = uh_ok && ~isempty(f); close(f);
end
close(findobj(0, 'Type', 'figure', 'Name', 'f Value Colorbar'));
all_passed = check('plot_unit_histograms makes one figure per sheet measure', uh_ok) && all_passed;
[~, fh] = load_and_make_unit_histograms(psa.output_dir, 'Metrics', {'lle', 'hks', 'dky'}, 'ColorBy', 'level_of_chaos');
all_passed = check('load_and_make_unit_histograms accepts registry keys', numel(fh) == 3) && all_passed;
close(fh); close(findobj(0, 'Type', 'figure', 'Name', 'f Value Colorbar'));

rmdir(tmp, 's'); rmdir(tmpb, 's');
fprintf('\n');
if all_passed
    fprintf('=== ALL sweep_metrics TESTS PASSED ===\n');
else
    fprintf('=== SOME sweep_metrics TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function psa = make_psa(cls, defaults, conds, out_dir, method)
defaults.lya_method = method;
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 10, 'note', 'sweep_metrics_test', ...
    'randomize_order', false, 'verbose', false, 'use_parallel', false);
psa.model_class = cls;
psa.integer_params = {'n', 'indegree'};
psa.model_defaults = defaults;
psa.output_dir = out_dir;
psa.set_conditions(conds);
psa.add_grid_parameter('level_of_chaos', [1.0, 2.0]);
psa.add_grid_parameter('reps', 1:2);
end

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

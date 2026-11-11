% test_psa_pairs_end_to_end.m - a tiny ParamSpaceAnalysis2 sweep of
% SRNNCellTypePairs, configured the way the paper's sweeps are, checked
% from the context struct to the files on disk.
%
% Written 2026-11-10 (docs/ISSUES.md, FD-001). The sweep stages
% (run_sensitivity_analysis and friends) are sized by run mode and take
% minutes even at 'fast', so this builds the same thing they build, from
% the same resolve_run_context output, and shrinks it: the manuscript's
% preset and its three conditions, noise on, the integrator and the top-K
% Lyapunov settings analysis_run_config chose, but n = 30, a 4-s run, and
% a 2 x 2 grid (level_of_chaos x reps). 12 jobs, serial.
%
% Checks:
%   1. resolve_run_context: Pairs, three conditions, f_param 'f_E', the
%      stochastic integrator for a noisy preset, model_defaults = preset
%      overlaid by the run mode's model settings.
%   2. The sweep: every job of every condition succeeds and carries the
%      expected fields (identity, mean_rate, mean_synaptic_output, every
%      lya_summary scalar), with finite LLE and a rate in (0, 1).
%   3. Pairing: one grid point uses one network seed in all three
%      conditions, and different grid points use different seeds; the
%      condition's own parameters reach the model (effective_param).
%   4. A sweep result equals a direct run: rebuild_model + build + run
%      reproduces the stored LLE and mean rate exactly, for a job with one
%      depression component per route and a job with two.
%   5. On disk: psa_object.mat, param_space_summary.mat and one result
%      file per condition (in a subfolder named after the condition);
%      from_dir gives back the same results;
%      resolved_defaults records the class, the noise and the integrator.
%   6. Reproducibility: the same sweep run again into another directory
%      gives identical LLE, mean rate and D_KY in every job, and
%      same_config accepts the pair.
%
% test_psa_model_class and test_sweep_metrics already run small Pairs
% sweeps (with noise off); this adds the noisy configuration, the pairing
% across conditions, the sweep-equals-direct-run check and the rerun.
%
% Prints PASS/FAIL per check and a final banner. About 6 s. Assumes
% setup_paths has run.
%
% See also: ParamSpaceAnalysis2, resolve_run_context, test_sweep_metrics

fprintf('=== Testing a tiny end-to-end Pairs sweep ===\n\n');
all_passed = true;

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium';
cond_names = {'no_adaptation', 'sfa1_std1', 'sfa3_std2'};

%% 1. The context the stages receive
ctx = resolve_run_context('param_space', 'preset_name', P, 'run_mode', 'fast', 'verbose', 'near-none');
all_passed = check('ctx: SRNNCellTypePairs, f_param f_E, integer_params {n, indegree}', ...
    strcmp(ctx.model_class, 'SRNNCellTypePairs') && strcmp(ctx.f_param, 'f_E') && ...
    isequal(ctx.integer_params, {'n', 'indegree'})) && all_passed;
all_passed = check('ctx: the preset''s three conditions, in order', ...
    isequal(cellfun(@(c) c.name, ctx.conditions, 'UniformOutput', false), cond_names)) && all_passed;
all_passed = check('ctx: a noisy preset gets the stochastic integrator (sra1) and top-K', ...
    ctx.cfg.is_stochastic && strcmp(ctx.model_defaults.ode_solver, 'sra1') && ...
    strcmp(ctx.model_defaults.lya_method, 'topk') && ctx.model_defaults.sigma_u_noise == 0.025) && all_passed;
all_passed = check('ctx: model_defaults is the preset overlaid by the run mode''s model settings', ...
    isequaln(ctx.model_defaults, merge_struct(ctx.preset_defaults, ctx.cfg.model))) && all_passed;

%% 2. The sweep
root = fullfile(tempdir, 'psa_pairs_end_to_end_test');
if exist(root, 'dir'); rmdir(root, 's'); end
psa = make_psa(ctx, fullfile(root, 'a'));
evalc('psa.run();');

n_points = 4;                                   % 2 gains x 2 reps
res = struct();
ok_count = true; ok_success = true;
for c = 1:numel(cond_names)
    rc = psa.results.(cond_names{c});
    ok_count = ok_count && numel(rc) == n_points && ~any(cellfun(@isempty, rc));
    res.(cond_names{c}) = rc;
    ok_success = ok_success && all(cellfun(@(r) isfield(r, 'success') && r.success, rc));
end
all_passed = check(sprintf('%d jobs returned per condition, none empty', n_points), ok_count) && all_passed;
all_passed = check('every job succeeded', ok_success) && all_passed;
if ~ok_success
    for c = 1:numel(cond_names)
        for k = 1:numel(res.(cond_names{c}))
            r = res.(cond_names{c}){k};
            if ~r.success; fprintf('      %s job %d: %s\n', cond_names{c}, k, r.error_message); end
        end
    end
end

need = [{'success', 'config', 'config_idx', 'condition_name', 'network_seed', 'run_duration', ...
    'mean_rate', 'mean_synaptic_output', 'local_rate_lead', 't_lya_lead'}, ...
    SRNNCellTypePairs.lya_summary_fields()];
ok_fields = true; ok_names = true; ok_values = true;
for c = 1:numel(cond_names)
    for k = 1:n_points
        r = res.(cond_names{c}){k};
        missing = need(~isfield(r, need));
        if ~isempty(missing)
            ok_fields = false;
            fprintf('      %s job %d lacks: %s\n', cond_names{c}, k, strjoin(missing, ', '));
            continue
        end
        ok_names = ok_names && strcmp(r.condition_name, cond_names{c}) && r.config_idx == k && ...
            all(isfield(r.config, {'level_of_chaos', 'reps'}));
        ok_values = ok_values && isfinite(r.LLE) && r.mean_rate > 0 && r.mean_rate < 1 && ...
            isfinite(r.mean_synaptic_output) && r.mean_synaptic_output > 0 && ...
            r.mean_synaptic_output <= r.mean_rate + 1e-12 && r.K_used >= 15;
    end
end
all_passed = check(sprintf('every result carries the %d expected fields', numel(need)), ok_fields) && all_passed;
all_passed = check('condition_name, config_idx and the grid values are recorded per job', ok_names) && all_passed;
all_passed = check('LLE finite, 0 < synaptic output <= rate < 1, K_used >= 15', ok_values) && all_passed;
theta_ratio = cellfun(@(r) r.mean_synaptic_output / r.mean_rate, res.sfa3_std2);
all_passed = check(sprintf('no depression: synaptic output == rate; two components: clearly below it (ratio %.2f to %.2f)', ...
    min(theta_ratio), max(theta_ratio)), ...
    all(cellfun(@(r) abs(r.mean_synaptic_output - r.mean_rate) < 1e-12, res.no_adaptation)) && ...
    all(theta_ratio < 0.9)) && all_passed;

%% 3. Pairing across conditions
seeds = zeros(n_points, numel(cond_names));
for c = 1:numel(cond_names)
    seeds(:, c) = cellfun(@(r) r.network_seed, res.(cond_names{c}));
end
all_passed = check('one grid point, one network seed in all three conditions', ...
    all(seeds == seeds(:, 1), 'all')) && all_passed;
all_passed = check('different grid points (and reps) use different seeds', ...
    numel(unique(seeds(:, 1))) == n_points) && all_passed;
tau1 = psa.effective_param(res.sfa1_std1{1}, 'tau_a');
tau3 = psa.effective_param(res.sfa3_std2{1}, 'tau_a');
gains = cellfun(@(r) psa.effective_param(r, 'level_of_chaos'), res.sfa3_std2);
all_passed = check('effective_param reads the condition''s SFA ladder (1 rung vs 3) and the swept gain', ...
    numel(tau1{1}) == 1 && numel(tau3{1}) == 3 && isequal(unique(gains)', [0.5 1.5])) && all_passed;
all_passed = check('effective_param reads a run-wide default (n = 30)', ...
    psa.effective_param([], 'n') == 30) && all_passed;

%% 4. A sweep result equals a direct run
for cn = {'sfa1_std1', 'sfa3_std2'}
    r = res.(cn{1}){end};
    mdl = psa.rebuild_model(r);
    mdl.verbose = 'near-none';
    evalc('mdl.build(); mdl.run();');
    s = mdl.lya_summary();
    % mean_rate is the mean over every stored sample of every neuron, the
    % types pooled in cell_type_names order.
    rates = [mdl.plot_data.r.E(:); mdl.plot_data.r.I(:)];
    all_passed = check(sprintf('%s: rebuild_model + run reproduces the stored LLE (%+.4f) and mean rate exactly', ...
        cn{1}, r.LLE), ...
        s.LLE == r.LLE && mean(rates) == r.mean_rate) && all_passed;
    all_passed = check(sprintf('%s: the rebuilt model has %d depression component(s) on all four routes', ...
        cn{1}, mdl.n_b_pairs(1)), ...
        all(mdl.n_b_pairs(:) == str2double(cn{1}(end)))) && all_passed;
end

%% 5. On disk
% The two run-level files at the top, each condition's results in a
% subfolder named after the condition.
files = [{'psa_object.mat', 'param_space_summary.mat'}, ...
    cellfun(@(c) fullfile(c, sprintf('param_space_results_%s.mat', c)), cond_names, 'UniformOutput', false)];
on_disk = cellfun(@(f) exist(fullfile(psa.output_dir, f), 'file') == 2, files);
all_passed = check('the run directory holds psa_object.mat, param_space_summary.mat and <condition>/param_space_results_<condition>.mat', ...
    all(on_disk)) && all_passed;
all_passed = check('the run directory is a dated folder inside the output_dir that was set', ...
    startsWith(psa.output_dir, fullfile(root, 'a')) && ~strcmp(psa.output_dir, fullfile(root, 'a'))) && all_passed;
loaded = ParamSpaceAnalysis2.from_dir(psa.output_dir);
all_passed = check('from_dir returns the same results, conditions and grid', ...
    isequaln(loaded.results, psa.results) && isequaln(loaded.conditions, psa.conditions) && ...
    isequal(loaded.grid_params, psa.grid_params) && strcmp(loaded.model_class, 'SRNNCellTypePairs')) && all_passed;
rd = psa.resolved_defaults;
all_passed = check('resolved_defaults records the noise, the integrator and the Lyapunov method', ...
    rd.sigma_u_noise == 0.025 && strcmp(rd.ode_solver, 'sra1') && strcmp(rd.lya_method, 'topk') && ...
    rd.n == 30) && all_passed;
all_passed = check('resolved_defaults leaves out the grid axes and the condition fields', ...
    ~any(isfield(rd, {'level_of_chaos', 'reps', 'tau_a', 'synapse_config'}))) && all_passed;

%% 6. The same sweep again
psb = make_psa(ctx, fullfile(root, 'b'));
evalc('psb.run();');
same = true;
for c = 1:numel(cond_names)
    for f = {'LLE', 'mean_rate', 'mean_synaptic_output', 'D_KY', 'h_KS_bits', 'network_seed'}
        va = cellfun(@(r) r.(f{1}), psa.results.(cond_names{c}));
        vb = cellfun(@(r) r.(f{1}), psb.results.(cond_names{c}));
        if ~isequaln(va, vb)
            same = false;
            fprintf('      %s %s differs: %s vs %s\n', cond_names{c}, f{1}, mat2str(va', 6), mat2str(vb', 6));
        end
    end
end
all_passed = check('a second sweep gives identical LLE, rate, synaptic output, D_KY and h_KS in every job', same) && all_passed;
[tf, why] = psa.same_config(psb);
all_passed = check(sprintf('same_config accepts the two runs%s', why_text(why)), tf) && all_passed;

rmdir(root, 's');
fprintf('\n');
if all_passed
    fprintf('=== ALL end-to-end Pairs sweep TESTS PASSED ===\n');
else
    fprintf('=== SOME end-to-end Pairs sweep TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function psa = make_psa(ctx, out_dir)
% What run_param_space_analysis / run_sensitivity_analysis set up from ctx,
% on a small network and a short run.
d = ctx.model_defaults;
d.n = 30; d.indegree = 6; d.F_tracks_network = true;
d.T_range = [0 4]; d.lya_T_interval = [2 4]; d.lya_warmup = 1;
psa = ParamSpaceAnalysis2('n_levels', 2, 'batch_size', 6, 'note', 'end_to_end_test', ...
    'randomize_order', false, 'verbose', 'near-none', 'use_parallel', false);
psa.model_class    = ctx.model_class;
psa.integer_params = ctx.integer_params;
psa.model_defaults = d;
psa.output_dir     = out_dir;
psa.set_conditions(ctx.conditions);
psa.add_grid_parameter('level_of_chaos', [0.5, 1.5]);
psa.add_grid_parameter('reps', 1:2);
end

function s = why_text(why)
if isempty(why); s = ''; else; s = sprintf(' (%s)', why); end
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

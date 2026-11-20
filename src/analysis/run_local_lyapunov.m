function mat_file = run_local_lyapunov(cfg)
% RUN_LOCAL_LYAPUNOV Local Lyapunov exponents and local KS entropy under a stimulus staircase.
%
%   mat_file = RUN_LOCAL_LYAPUNOV('preset_name', p, 'run_mode', m, 'out_dir', d)
%
% One network (rng_seeds [1 2]) of the preset under each of its adaptation
% conditions, run TWICE on the same Brownian path: top-K QR (exactly K, no
% retry) and Benettin (K = 1, the manuscript's estimator). Legacy defaults
% use T=60 s and accumulation[30,60] after15-s alignment. Explicit duration,
% accumulation and stimulus-cadence options support revised bundles (mu7revised:
% K50,T40,8 five-second input intervals). K defaults by mode: fast30,
% medium100, production100; explicit K overrides these.
%
% What it is for (TR, 2026-09-14): how often the local rates are positive, and
% whether the local KS entropy rate sum_k max(local_k, 0) spikes at every
% change of the input and returns to zero -- the dynamical signature of a
% stable adapting network encoding a new input. Grown out of
% scripts/explorations/explore_local_lyapunov_exponents.m.
%
% Output: <out_dir>/local_lyapunov_data.mat with results(i) per condition
% (name, title, topk = lya_results, ben = lya_results, u_ex, t_ex,
% type_indices, cell_type_names, seconds) and settings. A standalone run lands
% in <root>/data/local_lyapunov.
%
% THE SPECTRA ALONG THE SAME RUN (FD-034; cfg.spectra_n_samples > 0, off by
% default, and then nothing below is computed and the saved file is what it
% was). Figure 3 draws the local rates of this stage beside an eigenvalue
% density with a Ritz overlay; from the eig_heatmap stage those belong to a
% separate baseline run with a top-K pass of its own. With the option, both
% are taken on THIS switching run, and the Ritz values come from the pass
% whose local rates the figure plots:
%   * spectra_n_samples states evenly spaced over the accumulation window
%     [accumulation start, T], both ends included (40 over [1, 40] s: one
%     per second, 5 per 5-s input step);
%   * the eigenvalues of the full Jacobian at those states
%     (jacobian_spectra_at_states);
%   * H = Q' J Q at those states (ritz_pass_at_states, raw norm) from a
%     top-K pass with the model's own arguments (K, seed, warm-up,
%     accumulation window, re-orthonormalization interval). The model's
%     pass cannot return H itself, so the pass is made a second time on the
%     stored trajectory; lyapunov_topk's outputs do not depend on whether H
%     is asked for, and the stage errors :PassMismatch unless the exponents
%     of the two passes are equal bit for bit. The diagonal of H is then the
%     instantaneous value of the local rates saved in results(i).topk, in
%     Gram-Schmidt order (sorted column j is raw column sort_idx(j)).
% Saved beside the stage's file, as <out_dir>/local_lyapunov_spectra.mat, in
% the layout fig_eig_heatmap reads (evals_by_cond, ritz_by_cond, cond_names,
% condition_titles, lle_by_cond, lle_window and the per-state abscissae), so
% that figure and the composite Figure 3 draw it with
% 'spectra_source', 'local_lyapunov'. Added for this run: state_times and
% state_step (the input step each state belongs to: a state at a switch
% instant has not yet felt the new input and is counted with the step that
% ends there), and in each ritz_by_cond{i} the pass's sort_idx and its local
% rates over the QR intervals ending and starting at each state
% (local_rates_before, local_rates_after; states x K, Gram-Schmidt order;
% NaN where the run has no such interval). The trajectory is held in full
% while a condition is processed (store_full_state: steps x state dimension
% doubles), one condition at a time. local_lyapunov_data.mat is unchanged
% but for settings.spectra_n_samples.
%
% See also: fig_local_lyapunov, lyapunov_topk, SRNNCellTypePairs.lya_summary,
%           fig_eig_heatmap, fig_main3_stability, with_switch_spectra

arguments
    cfg.preset_name (1,:) char = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_steps5s_noise0p025_dualStd_3cond_mu7'
    cfg.run_mode    (1,:) char = 'production'
    cfg.out_dir     (1,:) char = ''
    cfg.verbose                 = 'minimal'
    cfg.K           (1,1) double = 0     % 0 -> per run_mode
    cfg.T           (1,1) double = 60
    cfg.step_duration_s (1,1) double = 0 %0 preserves preset input config
    cfg.accumulation_start_s (1,1) double = NaN % NaN preserves T/2
    cfg.warmup_s (1,1) double = NaN % NaN preserves T/4
    cfg.simulation_start_s (1,1) double = 0
    cfg.n_override  (1,1) double = 0     % tests only
    cfg.conditions               = {}    % the preset's conditions to run, by name; {} = all (FD-030)
    cfg.spectra_n_samples (1,1) double {mustBeInteger, mustBeNonnegative} = 0   % 0 = off; else states at which to save the spectra and H (FD-034)
    cfg.use_parallel (1,1) logical = true    % the spectra's eigenvalue problems in a parfor (spectra_n_samples > 0 only)
end

setup_paths();
switch cfg.run_mode
    case {'essential', 'fast'}, K = 30;      % 'essential' takes 'fast's sizes
    % 'fastmedium' takes medium's K (Tom, 2026-11-14). The mu7chain bundles
    % name K = 50 themselves (cfg.local_lyapunov_K), which ran 8.0 min at
    % 'fast' (2026-11-13), so there the mode changes nothing.
    case {'fastmedium', 'medium', 'medium2'}, K = 100;
    case 'production',          K = 100;
    otherwise
        error('run_local_lyapunov:badMode', 'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.K > 0; K = cfg.K; end
T = cfg.T;
acc_start=cfg.accumulation_start_s; if isnan(acc_start), acc_start=T/2; end
warmup=cfg.warmup_s; if isnan(warmup), warmup=T/4; end
assert(isfinite(acc_start) && acc_start>=0 && acc_start<T && ...
    isfinite(warmup) && warmup>=0 && isfinite(cfg.simulation_start_s) && ...
    cfg.simulation_start_s<=acc_start-warmup, ...
    'run_local_lyapunov:Window','Simulation must include the full alignment before accumulation.');
if isempty(cfg.out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'local_lyapunov');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir); mkdir(out_dir); end

[preset_defaults, ~, conditions] = preset_with_conditions(cfg.preset_name, cfg.conditions);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
title_map  = srnn_condition_titles();
titles     = cellfun(@(n) title_map(n), cond_names, 'UniformOutput', false);
seeds  = [1 2];
common = {'rng_seeds', seeds, 'fs', 400, 'T_range', [cfg.simulation_start_s T], 'lya_T_interval', [acc_start T], ...
          'lya_warmup', warmup, 'verbose', cfg.verbose};
if cfg.step_duration_s>0
    n_steps=round(T/cfg.step_duration_s);
    assert(n_steps>=1 && abs(n_steps*cfg.step_duration_s-T)<1e-9, ...
        'run_local_lyapunov:Cadence','Duration must contain an integer number of stimulus steps.');
    input=preset_defaults.input_config;
    input.n_steps=n_steps; input.no_stim_pattern=false(1,n_steps);
    common=[common,{'input_config',input}];
end
if cfg.n_override > 0
    common = [common, {'n', cfg.n_override, 'indegree', max(2, round(0.2 * cfg.n_override)), ...
        'F_tracks_network', true}];
end
vprintf(cfg.verbose, 'minimal', '[local_lyapunov] preset=%s run_mode=%s: %d conditions, T = %g s, K = %d\n', ...
    cfg.preset_name, cfg.run_mode, numel(cond_names), T, K);

with_spectra = cfg.spectra_n_samples > 0;
topk_extra = {};
if with_spectra
    assert(cfg.spectra_n_samples >= 2, 'run_local_lyapunov:SpectraSamples', ...
        'spectra_n_samples must be 0 (off) or at least 2.');
    topk_extra = {'store_full_state', true};     % the trajectory is needed after the run
    spectra = struct('evals', {cell(1, numel(cond_names))}, 'ritz', {cell(1, numel(cond_names))}, ...
        'num_abscissa', {cell(1, numel(cond_names))}, 'spec_abscissa', {cell(1, numel(cond_names))}, ...
        'times', {cell(1, numel(cond_names))}, 'seconds', nan(2, numel(cond_names)));
end

t_stage = tic;
results = struct('name', cond_names, 'title', titles, 'topk', [], 'ben', [], 'u_ex', [], 't_ex', [], ...
    'type_indices', [], 'cell_type_names', [], 'seconds', []);
for i = 1:numel(cond_names)
    t0 = tic;
    rng(0, 'twister');
    m = build_from_preset(cfg.preset_name, cond_names{i}, common{:}, ...
        'lya_method', 'topk', 'lya_K', K, 'lya_K_auto', false, 'lya_dt', 0.05, topk_extra{:});
    m.run();
    results(i).topk = m.lya_results;
    results(i).u_ex = m.u_ex; results(i).t_ex = m.t_ex;
    results(i).type_indices = m.type_indices; results(i).cell_type_names = m.cell_type_names;
    if with_spectra
        [spectra.evals{i}, spectra.num_abscissa{i}, spectra.spec_abscissa{i}, spectra.times{i}, ...
            spectra.ritz{i}, spectra.seconds(:, i)] = spectra_along_run(m, cfg, acc_start, T, warmup);
        vprintf(cfg.verbose, 'minimal', '  %-14s spectra at %d states: eigenvalues %.0f s, H from the top-%d pass %.0f s\n', ...
            cond_names{i}, numel(spectra.times{i}), spectra.seconds(1, i), spectra.ritz{i}.K, spectra.seconds(2, i));
    end
    if with_spectra; clear m; end                % a stored trajectory is not held across the Benettin run
    rng(0, 'twister');
    m = build_from_preset(cfg.preset_name, cond_names{i}, common{:}, 'lya_method', 'benettin');
    m.run();
    results(i).ben = m.lya_results;
    results(i).seconds = toc(t0);
    vprintf(cfg.verbose, 'minimal', '  %-14s top-%d lambda_1 %+.4f | Benettin %+.4f | h_KS %.2f bit/s | %.0f s\n', ...
        cond_names{i}, K, results(i).topk.LE_spectrum(1), results(i).ben.LLE, ...
        results(i).topk.h_KS_bits, results(i).seconds);
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, 'seeds', seeds, ...
    'T', T, 'K', K, 'lya_T_interval', [acc_start T], 'lya_warmup', warmup, 'fs', 400, ...
    'simulation_start_s', cfg.simulation_start_s, 'step_duration_s',cfg.step_duration_s, 'n_override', cfg.n_override, 'minutes', toc(t_stage) / 60);
condition_titles = titles; %#ok<NASGU>
mat_file = fullfile(out_dir, 'local_lyapunov_data.mat');
if with_spectra
    % Only when asked for, so a run without it saves what it always did.
    settings.spectra_n_samples = cfg.spectra_n_samples;
    save_spectra(out_dir, spectra, results, cond_names, settings, cfg);
end
save(mat_file, 'results', 'cond_names', 'condition_titles', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '[local_lyapunov] %.1f min -> %s\n', settings.minutes, mat_file);
end

%% ------------------------------------------------------------------------
function [evals, omega, alpha, t_used, R, seconds] = spectra_along_run(m, cfg, acc_start, T, warmup)
% The full-Jacobian spectra and H = Q'JQ at states evenly spaced over the
% accumulation window of the run just made (FD-034; see the header). m holds
% its trajectory (store_full_state) and its own top-K result.
t_req = linspace(acc_start, T, cfg.spectra_n_samples);
half_step = 0.5 / m.fs;
idx = arrayfun(@(tt) find(m.t_out >= tt - half_step, 1, 'first'), t_req);   % the stored row at each time
assert(numel(unique(idx)) == numel(idx), 'run_local_lyapunov:SpectraSamples', ...
    'spectra_n_samples = %d asks for more states than the run stores in [%g, %g] s.', cfg.spectra_n_samples, acc_start, T);
t0 = tic;
[evals, omega, alpha, t_used] = jacobian_spectra_at_states(m, idx, cfg.use_parallel);
seconds = [toc(t0); NaN];
t0 = tic;
tk = m.lya_results;
R = ritz_pass_at_states(m, idx, tk.K, {'raw'}, evals, struct('warmup_s', warmup, 'lya_dt', tk.lya_dt));
seconds(2) = toc(t0);
% The same pass as the model's: same function and arguments, so the same bits.
if ~isequal(R.LE_spectrum(:, 1), tk.LE_spectrum)
    error('run_local_lyapunov:PassMismatch', ...
        ['The top-%d pass that returns H does not reproduce the model''s own pass ' ...
         '(largest difference of the exponents %.3g 1/s); H would not belong to the saved local rates.'], ...
        tk.K, max(abs(R.LE_spectrum(:, 1) - tk.LE_spectrum)));
end
% The pass's local rates over the QR intervals ending and starting at each
% state, in H's (Gram-Schmidt) column order: t_lya holds segment starts.
raw = zeros(1, tk.K); raw(tk.sort_idx) = 1:tk.K;          % raw column j is sorted column raw(j)
R.sort_idx = tk.sort_idx(:)';
R.local_rates_before = nan(numel(idx), tk.K);
R.local_rates_after  = nan(numel(idx), tk.K);
for s = 1:numel(idx)
    % Interval k covers [t_lya(k), t_lya(k) + lya_dt]. "Before": the last one
    % that starts before the state (it ends at the state or contains it);
    % "after": the first one that ends after it. The same interval for a
    % state strictly inside one.
    before = find(tk.t_lya < t_used(s) - half_step, 1, 'last');
    after  = find(tk.t_lya + tk.lya_dt > t_used(s) + half_step, 1, 'first');
    if ~isempty(before); R.local_rates_before(s, :) = tk.local_LE_spectrum_t(before, raw); end
    if ~isempty(after);  R.local_rates_after(s, :)  = tk.local_LE_spectrum_t(after, raw);  end
end
end

function save_spectra(out_dir, spectra, results, cond_names, settings, cfg)
% local_lyapunov_spectra.mat, in the layout fig_eig_heatmap reads.
st = manuscript_style();
condition_titles = cond_names;
for i = 1:numel(cond_names)
    if st.condition_title.isKey(cond_names{i}); condition_titles{i} = st.condition_title(cond_names{i}); end
end
evals_by_cond         = spectra.evals;
ritz_by_cond          = spectra.ritz;
num_abscissa_by_cond  = spectra.num_abscissa;
spec_abscissa_by_cond = spectra.spec_abscissa;
J_times_by_cond       = spectra.times;
lle_by_cond    = arrayfun(@(r) r.topk.LLE, results);
lya_T_interval = settings.lya_T_interval;
lle_window     = diff(lya_T_interval);
state_times    = spectra.times{1};
% The step whose input shaped the state: (0, d] is step 1, (d, 2d] step 2, ...
if settings.step_duration_s > 0
    state_step = max(1, ceil(state_times / settings.step_duration_s - 1e-9));
else
    state_step = nan(size(state_times));
end
settings.trajectory = 'the switching-input run of the local_lyapunov stage';
settings.spectra_seconds = spectra.seconds;
settings.ritz_K = settings.K;
settings.ritz_scalings = {'raw'};
settings.n_samples = cfg.spectra_n_samples;
save(fullfile(out_dir, 'local_lyapunov_spectra.mat'), 'evals_by_cond', 'ritz_by_cond', 'cond_names', ...
    'condition_titles', 'lle_by_cond', 'lle_window', 'lya_T_interval', 'num_abscissa_by_cond', ...
    'spec_abscissa_by_cond', 'J_times_by_cond', 'state_times', 'state_step', 'settings', '-v7.3');
end

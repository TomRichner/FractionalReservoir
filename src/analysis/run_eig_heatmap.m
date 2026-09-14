function mat_file = run_eig_heatmap(cfg)
% RUN_EIG_HEATMAP Sample Jacobian eigenvalues through a run, per adaptation regime.
%
%   mat_file = RUN_EIG_HEATMAP()
%   mat_file = RUN_EIG_HEATMAP('run_mode', 'fast')
%
% The COMPUTE half of the eigenvalue-occupancy figure. Runs the four adaptation
% regimes on a shared network, samples the instantaneous Jacobian at fixed
% intervals through each run, pools the eigenvalues, and saves them to
% eig_heatmap_data.mat. Plotting is fig_eig_heatmap, so the look can be iterated
% without re-simulating.
%
% Because the network is nonlinear, the eigenvalues of the instantaneous
% Jacobian move around the complex plane as the state evolves. Pooling them over
% many sampled times shows how much time they spend in each region, and in
% particular to the RIGHT of the imaginary axis (Re > 0, locally unstable).
%
% All four conditions share rng_seeds, so W is identical and they are directly
% comparable.
%
% PORTED to the paper's preset. Two decisions worth stating:
%
%   * level_of_chaos COMES FROM THE PRESET (1.0), not from the 3.0 the original
%     script set. That 3.0 existed because the original was DETERMINISTIC and
%     needed high gain to make the eigenvalues wander at all; the paper's preset
%     is stochastic (sigma_u_noise = 0.025), and the noise is what moves the
%     state -- and therefore the Jacobian -- around. Using the preset's own gain
%     keeps this panel showing the same network as every other figure.
%   * THE INTEGRATOR IS NAMED EXPLICITLY. sigma_u_noise > 0 requires a
%     stochastic scheme, and this function does not go through
%     analysis_run_config, which is what selects one for the sweeps. Left to the
%     class default it would fail outright with 'requires a stochastic
%     integrator'. build_from_preset picks sra1.
%
% Eigenvalues are computed through the CLASS's own static compute_Jacobian_fast,
% resolved by name, so this works for either model class.
% SRNNModel2.eigenvalue_time_series is not used: SRNNCellTypePairs has no such
% method, and its b-states are per route rather than per population, so the two
% classes do not share a state layout.
%
% See also: fig_eig_heatmap, build_from_preset, srnn_param_preset

arguments
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.run_mode    (1,:) char    = 'production'
    cfg.out_dir     (1,:) char    = ''
    cfg.verbose                    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (or a logical); see verbose_level
    cfg.n_samples   (1,1) double  = 0      % 0 -> per run_mode
    cfg.use_parallel (1,1) logical = true
end

setup_paths();
% A standalone run writes into data/, NOT next to this file. The old default
% dropped eig_heatmap_data.mat into the figure folder, where fig_eig_heatmap
% then read it forever -- so a pipeline run could write fresh data into the run
% directory and the figure would go on plotting the standalone copy. That is
% what happened between Aug 22 and Aug 26.
%
% The pipeline still passes out_dir explicitly (run_all_paper_analyses), landing
% the .mat in <run_dir>/eig_heatmap/.
if isempty(cfg.out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'eig_heatmap');
else
    out_dir = cfg.out_dir;
end
if ~isfolder(out_dir); mkdir(out_dir); end

% Cost/fidelity. T_range buys a longer trajectory to sample; n_samples buys
% resolution of the occupancy density. n is NOT reduced in fast mode: the
% eigenvalue cloud's shape depends on network size, so shrinking it would change
% the thing being measured rather than just measuring it less well.
% medium2 runs at medium effort: it differs from medium only in sweep dimensions
% (levels, reps, fs for the stochastic integrator), and this stage has no sweep.
% Same collapse as run_memory_capacity and run_dc_lle_analysis.
switch cfg.run_mode
    case 'fast',                   T_range = [0 40];  n_samples = 40;   fs = 200;
    case {'medium', 'medium2'},    T_range = [0 100]; n_samples = 150;  fs = 400;
    case 'production',             T_range = [0 200]; n_samples = 300;  fs = 400;
    otherwise
        % Every name in run_mode_names() must be handled above; test_run_modes
        % asserts it. Reaching here means a mode was added to the sweeps and this
        % stage was not taught about it -- which is how a 'fast2' run lost it.
        error('run_eig_heatmap:badMode', ...
            'Unknown run_mode ''%s'' (expected %s).', ...
            cfg.run_mode, strjoin(run_mode_names(), ', '));
end
if cfg.n_samples > 0; n_samples = cfg.n_samples; end

lle_window     = min(30, diff(T_range) / 2);
lya_T_interval = [T_range(2) - lle_window, T_range(2)];

% The conditions come FROM THE PRESET, in its order. Hardcoding the four
% original names here meant a 7-regime preset would have been sampled for only
% four of them, silently, and the figure would have looked complete.
[~, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
titles     = cellfun(@(n) pretty(n), cond_names, 'UniformOutput', false);

vprintf(cfg.verbose, 'minimal', '[eig_heatmap] preset=%s run_mode=%s T=%g s n_samples=%d\n', ...
    cfg.preset_name, cfg.run_mode, T_range(2), n_samples);

n_cond        = numel(cond_names);
evals_by_cond = cell(1, n_cond);
lle_by_cond   = nan(1, n_cond);
lya_by_cond   = cell(1, n_cond);
num_abscissa_by_cond  = cell(1, n_cond);   % max eig((J + J')/2) per sampled state
spec_abscissa_by_cond = cell(1, n_cond);   % max real eig(J) per sampled state
J_times_by_cond       = cell(1, n_cond);

for i = 1:n_cond
    vprintf(cfg.verbose, 'verbose', '\n=== %d/%d %s ===\n', i, n_cond, titles{i});
    model = build_from_preset(cfg.preset_name, cond_names{i}, 'verbose', cfg.verbose, ...
        'T_range',          T_range, ...
        'fs',               fs, ...
        'rng_seeds',        [1 2], ...      % same W across conditions
        'lya_method',       'topk', ...     % the sweeps' estimator (K 15, retry to 30)
        'lya_K',            15, 'lya_K_auto', true, 'lya_K_max', 30, 'lya_dt', 0.05, ...
        'lya_T_interval',   lya_T_interval, ...
        'store_full_state', true);          % required to read S_out below
    model.run();

    lle_by_cond(i) = model.lya_results.LLE;
    lya_by_cond{i} = model.lya_summary();

    % Sample after the LLE warmup window opens, so the transient is excluded.
    t_start = lya_T_interval(1);
    t_end   = T_range(2);
    J_times = linspace(t_start, t_end, n_samples);
    [evals_by_cond{i}, num_abscissa_by_cond{i}, spec_abscissa_by_cond{i}, J_times_by_cond{i}] = ...
        sample_eigenvalues(model, J_times, cfg.use_parallel);
    vprintf(cfg.verbose, 'minimal', '  %-14s LLE = %+.4f | %d eigenvalues pooled | numerical abscissa median %+.3f (spectral %+.3f)\n', ...
        lle_by_cond(i), numel(evals_by_cond{i}), ...
        median(num_abscissa_by_cond{i}), median(spec_abscissa_by_cond{i}));
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'T_range', T_range, 'fs', fs, 'n_samples', n_samples, ...
    'lle_window', lle_window, 'lya_T_interval', lya_T_interval, ...
    'model_class', class(model), 'level_of_chaos', model.level_of_chaos, ...
    'n', model.n, 'ode_solver', model.ode_solver, 'lya_method', model.lya_method);

condition_titles = titles;                      %#ok<NASGU>  saved name
mat_file = fullfile(out_dir, 'eig_heatmap_data.mat');
save(mat_file, 'evals_by_cond', 'condition_titles', 'cond_names', ...
    'lle_by_cond', 'lya_by_cond', 'num_abscissa_by_cond', 'spec_abscissa_by_cond', ...
    'J_times_by_cond', 'lle_window', 'lya_T_interval', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '[eig_heatmap] saved %s\n', mat_file);
end

%% ------------------------------------------------------------------------
function [ev_all, omega, alpha, t_used] = sample_eigenvalues(model, J_times_sec, use_parallel)
% Pool the Jacobian eigenvalues at the requested times, and per state the
% NUMERICAL ABSCISSA omega = max eig((J_xx + J_xx')/2) beside the SPECTRAL
% ABSCISSA alpha = max real eig(J_xx) of the DENDRITIC BLOCK
% J_xx = (-I + W diag(theta'(x)))/tau_d -- the rate-network Jacobian of
% Hennequin et al., with the adaptation and depression states held fixed.
% omega bounds the instantaneous growth rate of a perturbation of x
% (d/dt ||dx|| <= omega ||dx||), so omega > alpha is the margin by which the
% non-normal recurrent coupling can amplify transiently even when every
% eigenvalue decays (Trefethen & Embree). omega >= alpha always.
%
% WHY THE x BLOCK ONLY (measured 2026-09-12): on the FULL Jacobian the
% numerical abscissa is 50-180 /s against a spectral abscissa of 0-10, and it
% ranks the regimes by how their state coordinates are scaled rather than by
% dynamics -- the symmetric part mixes rows in units of 1/tau_d = 10 (x) with
% rows in units of 1/tau_a and 1/tau_rel (a, b), and the (x,a) vs (a,x)
% blocks differ by a factor W/(c tau_d) against 1/tau_a. The numerical
% abscissa is not invariant to a diagonal rescaling of the state, so it is
% only meaningful on a block with one unit. The full-J eigenvalues are still
% pooled for the heatmap as before.
%
% Resolves compute_Jacobian_fast BY CLASS NAME, so one implementation serves
% both model classes. They do not share a state layout -- SRNNCellTypePairs
% carries b-states per ROUTE where SRNNModel2 carries them per population -- so
% indexing S_out by hand would be class-specific and fragile; the x rows are
% the LAST n of the state on both classes.
cls    = class(model);
params = model.get_params();
S_out  = model.S_out;
t_out  = model.t_out;
n      = model.n;

idx = arrayfun(@(tt) find(t_out >= tt, 1, 'first'), J_times_sec, ...
    'UniformOutput', false);
idx = unique([idx{~cellfun(@isempty, idx)}]);
t_used = t_out(idx);
t_used = t_used(:);

n_idx = numel(idx);
ev = cell(1, n_idx);
omega = zeros(n_idx, 1);
alpha = zeros(n_idx, 1);
if use_parallel
    parfor k = 1:n_idx
        J = full(feval([cls '.compute_Jacobian_fast'], S_out(idx(k), :)', params));
        ev{k} = eig(J);
        Jxx = J(end - n + 1:end, end - n + 1:end);
        alpha(k) = max(real(eig(Jxx)));
        omega(k) = max(eig((Jxx + Jxx') / 2));
    end
else
    for k = 1:n_idx
        J = full(feval([cls '.compute_Jacobian_fast'], S_out(idx(k), :)', params));
        ev{k} = eig(J);
        Jxx = J(end - n + 1:end, end - n + 1:end);
        alpha(k) = max(real(eig(Jxx)));
        omega(k) = max(eig((Jxx + Jxx') / 2));
    end
end
ev_all = vertcat(ev{:});
end

function s = pretty(name)
% Display name for a condition. Delegates to manuscript_style so this figure and
% every other one title the same regime identically, rather than keeping a
% private list that silently falls back to the raw snake_case name whenever a
% regime is added.
st = manuscript_style();
if st.condition_title.isKey(name)
    s = st.condition_title(name);
else
    s = strrep(name, '_', ' ');
end
end

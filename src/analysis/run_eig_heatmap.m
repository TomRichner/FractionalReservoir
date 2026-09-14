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
% EXAMPLES (2026-09-14). Besides the reference network the stage runs a
% prespecified set of E:I-imbalanced examples -- mu_EE_relative at 0.5x and
% 1.5x the preset's value by default -- with the SAME seeds, and stores per
% (example, condition) the pooled eigenvalues, lambda_1, the mean rate and the
% realised E:I weight balance, so fig_eig_heatmap_imbalance is a defined
% comparison annotated with matched numbers rather than a chosen seed. The
% legacy top-level variables describe the reference example exactly as before;
% `examples` (struct array) carries all of them. In fast mode the non-reference
% examples run at n = 250 (a dense eig at N = 4000 costs ~30 s a state).
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
    cfg.examples                   = []       % struct array (label, overrides); [] -> mu_EE x 0.5 / 1 / 1.5
    cfg.n_override  (1,1) double  = 0        % 0 -> preset n (fast mode runs the NON-reference examples at 250)
    cfg.n_samples_examples (1,1) double = 0  % states per non-reference example; 0 -> per run_mode (40 / 60 / 150)
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

n_cond = numel(cond_names);

%% The examples: the reference network plus prespecified E:I imbalances
% (TR / Codex audit sec. 7, 2026-09-14). One interpretable block-mean
% coordinate, mu_EE_relative, at 0.5x / 1x / 1.5x the preset's value, so the
% Jacobian-occupancy picture is a DEFINED comparison rather than a chosen
% seed: every example x condition shares rng_seeds [1 2], and each panel is
% annotated with its matched finite-time lambda_1, mean rate and realised E:I
% weight balance B_E (fig_eig_heatmap_imbalance). Pass cfg.examples to draw
% a different set; the reference example is always run and is what the
% legacy variables below (fig_eig_heatmap, fig_transient_amplification)
% describe.
examples = cfg.examples;
if isempty(examples)
    mu_EE_ref = preset_mu_EE(cfg.preset_name);
    examples = struct('label', {'inhibition-dominant', 'reference', 'excitation-dominant'}, ...
        'overrides', {struct('mu_EE_relative', 0.5 * mu_EE_ref), struct(), ...
                      struct('mu_EE_relative', 1.5 * mu_EE_ref)});
end
if ~any(strcmp({examples.label}, 'reference'))
    examples(end + 1) = struct('label', 'reference', 'overrides', struct());
end
n_ex = numel(examples);
n_examples_override = cfg.n_override;
if strcmp(cfg.run_mode, 'fast') && cfg.n_override == 0
    n_examples_override = 250;      % a dense eig at N = 4000 is ~30 s; fast keeps the examples affordable
end
n_samples_ex = cfg.n_samples_examples;
if n_samples_ex == 0
    switch cfg.run_mode
        case 'fast',                n_samples_ex = 40;
        case {'medium', 'medium2'}, n_samples_ex = 60;
        case 'production',          n_samples_ex = 150;
    end
end

ex_out = repmat(struct('label', '', 'overrides', struct(), 'mu_tilde_relative', [], 'n', NaN, ...
    'evals_by_cond', {cell(1, n_cond)}, 'lle_by_cond', nan(1, n_cond), ...
    'mean_rate_by_cond', nan(1, n_cond), 'B_E_by_cond', nan(1, n_cond), ...
    'lya_by_cond', {cell(1, n_cond)}, 'num_abscissa_by_cond', {cell(1, n_cond)}, ...
    'spec_abscissa_by_cond', {cell(1, n_cond)}, 'J_times_by_cond', {cell(1, n_cond)}), 1, n_ex);
model = [];
for e = 1:n_ex
    is_ref = strcmp(examples(e).label, 'reference');
    if is_ref
        n_ov = cfg.n_override; n_smp = n_samples;      % the legacy loop, unchanged
    else
        n_ov = n_examples_override; n_smp = n_samples_ex;
    end
    ex_out(e).label     = examples(e).label;
    ex_out(e).overrides = examples(e).overrides;
    for i = 1:n_cond
        vprintf(cfg.verbose, 'verbose', '\n=== %s | %d/%d %s ===\n', examples(e).label, i, n_cond, titles{i});
        r = sample_condition(cfg, cond_names{i}, examples(e).overrides, n_ov, T_range, fs, ...
            lya_T_interval, n_smp);
        ex_out(e).evals_by_cond{i}         = r.evals;
        ex_out(e).lle_by_cond(i)           = r.lle;
        ex_out(e).mean_rate_by_cond(i)     = r.mean_rate;
        ex_out(e).B_E_by_cond(i)           = r.B_E;
        ex_out(e).lya_by_cond{i}           = r.lya;
        ex_out(e).num_abscissa_by_cond{i}  = r.num_abscissa;
        ex_out(e).spec_abscissa_by_cond{i} = r.spec_abscissa;
        ex_out(e).J_times_by_cond{i}       = r.J_times;
        ex_out(e).mu_tilde_relative        = r.mu_tilde_relative;
        ex_out(e).n                        = r.n;
        model = r.model;
        vprintf(cfg.verbose, 'minimal', '  %-20s %-14s LLE = %+.4f | <r> %.3f | B_E %.2f | %d eigenvalues | num. abscissa median %+.3f (spectral %+.3f)\n', ...
            examples(e).label, cond_names{i}, r.lle, r.mean_rate, r.B_E, numel(r.evals), ...
            median(r.num_abscissa), median(r.spec_abscissa));
    end
end

% Legacy top-level variables = the reference example, exactly as before.
i_ref = find(strcmp({ex_out.label}, 'reference'), 1);
evals_by_cond         = ex_out(i_ref).evals_by_cond;
lle_by_cond           = ex_out(i_ref).lle_by_cond;
lya_by_cond           = ex_out(i_ref).lya_by_cond;
num_abscissa_by_cond  = ex_out(i_ref).num_abscissa_by_cond;
spec_abscissa_by_cond = ex_out(i_ref).spec_abscissa_by_cond;
J_times_by_cond       = ex_out(i_ref).J_times_by_cond;

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'T_range', T_range, 'fs', fs, 'n_samples', n_samples, ...
    'lle_window', lle_window, 'lya_T_interval', lya_T_interval, ...
    'model_class', class(model), 'level_of_chaos', model.level_of_chaos, ...
    'n', ex_out(i_ref).n, 'ode_solver', model.ode_solver, 'lya_method', model.lya_method, ...
    'examples', {{ex_out.label}}, 'example_overrides', {{ex_out.overrides}}, ...
    'n_examples_override', n_examples_override, 'n_samples_examples', n_samples_ex);

condition_titles = titles;                      %#ok<NASGU>  saved name
examples = ex_out;                              %#ok<NASGU>  saved name
mat_file = fullfile(out_dir, 'eig_heatmap_data.mat');
save(mat_file, 'evals_by_cond', 'condition_titles', 'cond_names', ...
    'lle_by_cond', 'lya_by_cond', 'num_abscissa_by_cond', 'spec_abscissa_by_cond', ...
    'J_times_by_cond', 'lle_window', 'lya_T_interval', 'settings', 'examples', '-v7.3');
vprintf(cfg.verbose, 'minimal', '[eig_heatmap] saved %s\n', mat_file);
end

%% ------------------------------------------------------------------------
function r = sample_condition(cfg, cond_name, overrides, n_override, T_range, fs, lya_T_interval, n_samples)
% One (example, condition): build on the shared seeds with the example's
% overrides, run with the sweeps' top-K estimator, sample the Jacobian after
% the LLE window opens, and read the mean rate and the realised E:I weight
% balance off the same model.
args = {'T_range', T_range, 'fs', fs, 'rng_seeds', [1 2], ...   % same W across conditions AND examples
    'lya_method', 'topk', 'lya_K', 15, 'lya_K_auto', true, 'lya_K_max', 30, 'lya_dt', 0.05, ...
    'lya_T_interval', lya_T_interval, 'store_full_state', true, 'verbose', cfg.verbose};
if n_override > 0
    args = [args, {'n', n_override, 'indegree', max(2, round(0.2 * n_override)), 'F_tracks_network', true}];
end
args = [args, struct2namevalue(overrides)];
model = build_from_preset(cfg.preset_name, cond_name, args{:});
model.run();

J_times = linspace(lya_T_interval(1), T_range(2), n_samples);
[evals, omega, alpha, t_used] = sample_eigenvalues(model, J_times, cfg.use_parallel);

% Mean rate over the sampled window (the plot_data grid may be decimated).
pd = model.plot_data;
keep = pd.t >= lya_T_interval(1);
names = fieldnames(pd.r);
acc = [];
for q = 1:numel(names)
    v = pd.r.(names{q});
    acc = [acc; reshape(v(:, keep), [], 1)]; %#ok<AGROW>
end
mean_rate = mean(acc(~isnan(acc)));

r = struct('evals', evals, 'lle', model.lya_results.LLE, 'lya', model.lya_summary(), ...
    'num_abscissa', omega, 'spec_abscissa', alpha, 'J_times', t_used, ...
    'mean_rate', mean_rate, 'B_E', ei_weight_balance(model), ...
    'mu_tilde_relative', model.mu_tilde_relative, 'n', model.n, 'model', model);
end

function v = ei_weight_balance(m)
% The realised E:I weight balance of the drawn W: the excitatory share of the
% total summed weight, |sum W_E| / (|sum W_E| + |sum W_I|). Same definition
% as ei_weight_fraction in fig_EI_weights_param_space (kept in step by hand).
ti = m.type_indices;
S_E = full(sum(sum(m.W(:, ti{1}))));
S_I = full(sum(sum(m.W(:, ti{2}))));
denom = abs(S_E) + abs(S_I);
if denom == 0
    v = 0.5;
else
    v = abs(S_E) / denom;
end
end

function mu = preset_mu_EE(preset_name)
% The preset's mu_EE_relative, read through the class's own alias (index
% (post, pre) = (1, 1); a 1 x C row is broadcast down the columns), by
% reading the preset struct: (post, pre) = (1, 1) is E->E on both layouts.
d = srnn_param_preset(preset_name);
mu = d.mu_tilde_relative;
mu = mu(1, 1);                      % (post, pre) = (E, E); a 1 x C row broadcasts, so (1) is E too
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

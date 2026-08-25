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
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond'
    cfg.run_mode    (1,:) char    = 'production'
    cfg.out_dir     (1,:) char    = ''
    cfg.n_samples   (1,1) double  = 0      % 0 -> per run_mode
    cfg.use_parallel (1,1) logical = true
end

setup_paths();
this_dir = fileparts(mfilename('fullpath'));
if isempty(cfg.out_dir); out_dir = this_dir; else; out_dir = cfg.out_dir; end
if ~isfolder(out_dir); mkdir(out_dir); end

% Cost/fidelity. T_range buys a longer trajectory to sample; n_samples buys
% resolution of the occupancy density. n is NOT reduced in fast mode: the
% eigenvalue cloud's shape depends on network size, so shrinking it would change
% the thing being measured rather than just measuring it less well.
switch cfg.run_mode
    case 'fast',       T_range = [0 40];  n_samples = 40;   fs = 200;
    case 'medium',     T_range = [0 100]; n_samples = 150;  fs = 400;
    case 'production', T_range = [0 200]; n_samples = 300;  fs = 400;
    otherwise
        error('run_eig_heatmap:badMode', 'Unknown run_mode ''%s''.', cfg.run_mode);
end
if cfg.n_samples > 0; n_samples = cfg.n_samples; end

lle_window     = min(30, diff(T_range) / 2);
lya_T_interval = [T_range(2) - lle_window, T_range(2)];

cond_names  = {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'};
[~, ~, conditions] = srnn_param_preset(cfg.preset_name);
titles = cellfun(@(n) pretty(n), cond_names, 'UniformOutput', false);

fprintf('[eig_heatmap] preset=%s run_mode=%s T=%g s n_samples=%d\n', ...
    cfg.preset_name, cfg.run_mode, T_range(2), n_samples);

n_cond        = numel(cond_names);
evals_by_cond = cell(1, n_cond);
lle_by_cond   = nan(1, n_cond);

for i = 1:n_cond
    fprintf('\n=== %d/%d %s ===\n', i, n_cond, titles{i});
    model = build_from_preset(cfg.preset_name, cond_names{i}, ...
        'T_range',          T_range, ...
        'fs',               fs, ...
        'rng_seeds',        [1 2], ...      % same W across conditions
        'lya_method',       'benettin', ...
        'lya_T_interval',   lya_T_interval, ...
        'store_full_state', true);          % required to read S_out below
    model.run();

    lle_by_cond(i) = model.lya_results.LLE;

    % Sample after the LLE warmup window opens, so the transient is excluded.
    t_start = lya_T_interval(1);
    t_end   = T_range(2);
    J_times = linspace(t_start, t_end, n_samples);
    evals_by_cond{i} = sample_eigenvalues(model, J_times, cfg.use_parallel);
    fprintf('  LLE = %+.4f | %d eigenvalues pooled\n', ...
        lle_by_cond(i), numel(evals_by_cond{i}));
end

settings = struct('preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'T_range', T_range, 'fs', fs, 'n_samples', n_samples, ...
    'lle_window', lle_window, 'lya_T_interval', lya_T_interval, ...
    'model_class', class(model), 'level_of_chaos', model.level_of_chaos, ...
    'n', model.n, 'ode_solver', model.ode_solver);

condition_titles = titles;                      %#ok<NASGU>  saved name
mat_file = fullfile(out_dir, 'eig_heatmap_data.mat');
save(mat_file, 'evals_by_cond', 'condition_titles', 'cond_names', ...
    'lle_by_cond', 'lle_window', 'lya_T_interval', 'settings', '-v7.3');
fprintf('\nSaved: %s\n', mat_file);
end

%% ------------------------------------------------------------------------
function ev_all = sample_eigenvalues(model, J_times_sec, use_parallel)
% Pool the Jacobian eigenvalues at the requested times.
%
% Resolves compute_Jacobian_fast BY CLASS NAME, so one implementation serves
% both model classes. They do not share a state layout -- SRNNCellTypePairs
% carries b-states per ROUTE where SRNNModel2 carries them per population -- so
% indexing S_out by hand would be class-specific and fragile.
cls    = class(model);
params = model.get_params();
S_out  = model.S_out;
t_out  = model.t_out;

idx = arrayfun(@(tt) find(t_out >= tt, 1, 'first'), J_times_sec, ...
    'UniformOutput', false);
idx = unique([idx{~cellfun(@isempty, idx)}]);

n_idx = numel(idx);
ev = cell(1, n_idx);
if use_parallel
    parfor k = 1:n_idx
        J = feval([cls '.compute_Jacobian_fast'], S_out(idx(k), :)', params);
        ev{k} = eig(full(J));
    end
else
    for k = 1:n_idx
        J = feval([cls '.compute_Jacobian_fast'], S_out(idx(k), :)', params);
        ev{k} = eig(full(J));
    end
end
ev_all = vertcat(ev{:});
end

function s = pretty(name)
switch name
    case 'no_adaptation', s = 'No adaptation';
    case 'sfa_only',      s = 'SFA only';
    case 'std_only',      s = 'STD only';
    case 'sfa_and_std',   s = 'SFA + STD';
    otherwise,            s = name;
end
end

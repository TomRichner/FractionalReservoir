function out_dir = run_dc_lle_analysis(opts)
% RUN_DC_LLE_ANALYSIS Multi-seed DC-staircase Lyapunov analysis, all conditions.
%
%   out_dir = RUN_DC_LLE_ANALYSIS()
%   out_dir = RUN_DC_LLE_ANALYSIS('run_mode', 'fast')
%   out_dir = RUN_DC_LLE_ANALYSIS('run_mode', 'fast', 'output_dir', d)
%
% Question: does tonic DC drive change the largest local Lyapunov exponent (LLE)
% of the network, does the effect hold across networks, and -- since 2026-09-02 --
% does the ADAPTATION REGIME change the answer?
%
% For each network seed and each adaptation condition we run ONE simulation in
% which a uniform DC is stepped through a staircase of levels, compute the local
% LLE with Benettin's method, and bin it by DC level, skipping psd_settle
% seconds after each step so SFA/STD have settled. Aggregating the per-level mean
% across seeds gives mean +/- std (seed variability) as a function of DC.
%
% PAIRED ACROSS CONDITIONS. Within one seed every condition is built from the
% same rng_seeds, so it sees the same network and the same stimulus and only the
% adaptation differs -- the same fairness design run_memory_capacity uses. Seeds
% are the parfor axis; conditions run serially inside a seed, which is what keeps
% them paired.
%
% TWO KNOBS, KEPT APART, as everywhere else:
%   preset_name  WHICH NETWORK. The paper's operating point by default, so this
%                measures the network the rest of the paper is about.
%   run_mode     HOW MUCH COMPUTE: seeds and the integrator. See
%                dc_lle_run_config below.
%
% WHY THE PROTOCOL IS NOT IN THE PRESET. input_config, T_range and fs are set
% here, and that is not the analysis fighting the preset -- a preset may not
% carry any of them:
%
%   * fs and T_range are the run_mode timing knobs, which CLAUDE.md's preset
%     contract explicitly excludes (they are analysis_run_config's cfg.model).
%   * T_range is DERIVED, [0, nL*hold_dur], so it depends on how many DC levels
%     the run mode chose. No fixed value is right anywhere.
%   * input_config carries a FUNCTION HANDLE, @dc_staircase_stimulus. The
%     bursting_pairs preset already ruled on this in its own comment: a handle
%     in a preset is what the "nonlinearity is named data" rule exists to avoid,
%     and it would break the property that presets are comparable data frozen
%     into resolved_defaults.
%
% This is the same split mc_esn's comment described for the MC protocol
% settings: the preset says what the NETWORK is, the analysis says what the
% EXPERIMENT is.
%
% NOTE input_config is REPLACED wholesale, not merged (see pairs_input_config's
% header), so every field the class needs is stated below even when the value
% matches the default.
%
% ONE NOISE SOURCE. input_config.noise_intensity is 0: the preset's additive
% Wiener noise on x is the noise, and it is the one the Lyapunov machinery is
% built around -- it keeps the diffusion constant and cancels in Benettin's
% trajectory difference, so the exponent stays measurable. Stimulus noise on top
% would be a second, redundant source the exponent is not designed for.
%
% NOT A ParamSpaceAnalysis2 SWEEP, deliberately: DC is a sub-field of
% input_config, which the grid framework cannot vary, and the staircase fixes
% T_range by construction. Hence a purpose-built parfor with its own run-mode
% table.
%
% PORTED to SRNNCellTypePairs on 2026-09-02, from a hardcoded SRNNModel2 network
% (n = 50, indegree = 4, mu -4/3, its own c). See the commit for what changed;
% no earlier output is comparable.
%
% See also: dc_staircase_stimulus, fig_dc_lle, replot_dc_lle, build_from_preset,
%           run_memory_capacity, confplot

arguments
    opts.preset_name    (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    opts.run_mode       (1,:) char    = 'production'
    opts.output_dir     (1,:) char    = ''
    opts.save_figs      (1,1) logical = false
    opts.use_parallel   (1,1) logical = true    % false for serial debugging
    % Integrator override. Empty means "decide from the preset" -- deterministic
    % at sigma_u_noise = 0, stochastic above it. Same rule and same reason as
    % run_memory_capacity.
    opts.ode_solver     (1,:) char    = ''
end

setup_paths();

%% Resolve the two knobs
[preset, model_class, conditions] = srnn_param_preset(opts.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_dc_lle_analysis:BadModelClass', ...
        ['Preset ''%s'' is written for %s. This analysis was ported to ' ...
         'SRNNCellTypePairs on 2026-09-02 and reads its conditions through ' ...
         'that class''s vocabulary.'], opts.preset_name, model_class);
end
cfg = dc_lle_run_config(opts.run_mode, preset, opts.ode_solver);

condition_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
n_cond = numel(condition_names);
nL     = numel(cfg.dc_levels);

fprintf('[dc_lle] preset=%s  run_mode=%s\n', opts.preset_name, opts.run_mode);
fprintf('[dc_lle] seeds=%d  conditions=%d  DC levels=%d  solver=%s  fs=%d\n', ...
    cfg.n_seeds, n_cond, nL, cfg.ode_solver, cfg.fs);
fprintf('[dc_lle] conditions: %s\n', strjoin(condition_names, ', '));

%% The staircase stimulus -- stated in full, since assigning input_config
%  REPLACES the struct the class built rather than merging into it.
T_range = [0, nL * cfg.hold_dur];
input_config = struct();
input_config.dc_levels       = cfg.dc_levels;
input_config.hold_dur        = cfg.hold_dur;
input_config.ramp_dur        = cfg.ramp_dur;
input_config.noise_intensity = cfg.noise_intensity;   % 0; see the header
input_config.intrinsic_drive = [];                    % required by the class
input_config.generator       = @dc_staircase_stimulus;

%% Output directory
project_root = fileparts(which('setup_paths'));
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
if isempty(opts.output_dir)
    base_dir = fullfile(project_root, 'data', 'dc_lle');
else
    base_dir = fullfile(opts.output_dir, 'dc_lle');
end
out_dir = fullfile(base_dir, sprintf('dc_lle_nSeeds_%d_%s', cfg.n_seeds, dt_str));
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
fprintf('[dc_lle] out_dir = %s\n', out_dir);

%% Run: parfor over SEEDS, conditions serial inside so they stay paired
n_seeds     = cfg.n_seeds;
seed_offset = 0;                      % deterministic: same networks every run
seeds       = seed_offset + (1:n_seeds);
seed_results = cell(n_seeds, 1);

% parfor-transparent locals
preset_name = opts.preset_name;
extra_args  = {'input_config', input_config, 'T_range', T_range, ...
               'fs', cfg.fs, 'ode_solver', cfg.ode_solver, ...
               'lya_method', 'benettin', 'store_full_state', false};

run_start = tic;
if opts.use_parallel
    parfor s = 1:n_seeds
        seed_results{s} = run_one_seed(s, seed_offset, preset_name, ...
            condition_names, extra_args, cfg.hold_dur, cfg.psd_settle, nL);
    end
else
    for s = 1:n_seeds
        seed_results{s} = run_one_seed(s, seed_offset, preset_name, ...
            condition_names, extra_args, cfg.hold_dur, cfg.psd_settle, nL);
    end
end
fprintf('[dc_lle] %d seeds x %d conditions finished in %.1f min\n', ...
    n_seeds, n_cond, toc(run_start)/60);

%% Aggregate: [n_seeds x nL x n_cond]
per_seed_level_mean = nan(n_seeds, nL, n_cond);
per_seed_level_std  = nan(n_seeds, nL, n_cond);
LLE_by_seed         = nan(n_seeds, n_cond);
for s = 1:n_seeds
    per_seed_level_mean(s, :, :) = seed_results{s}.lvl_mean;
    per_seed_level_std(s, :, :)  = seed_results{s}.lvl_std;
    LLE_by_seed(s, :)            = seed_results{s}.LLE;
end
level_mean = squeeze(mean(per_seed_level_mean, 1, 'omitnan'));   % nL x n_cond
level_std  = squeeze(std(per_seed_level_mean, 0, 1, 'omitnan')); % nL x n_cond

t_lya = seed_results{1}.t_lya;

%% Save
config = struct('preset_name', opts.preset_name, 'run_mode', opts.run_mode, ...
    'model_class', model_class, ...
    'hold_dur', cfg.hold_dur, 'ramp_dur', cfg.ramp_dur, ...
    'noise_intensity', cfg.noise_intensity, 'psd_settle', cfg.psd_settle, ...
    'fs', cfg.fs, 'T_range', T_range, 'ode_solver', cfg.ode_solver, ...
    'sigma_u_noise', preset_field(preset, 'sigma_u_noise', 0), ...
    'n', preset_field(preset, 'n', []), 'f', preset_field(preset, 'f', []));

dc_lle_results = struct( ...
    'dc_levels', cfg.dc_levels, 'conditions', {condition_names}, ...
    'hold_dur', cfg.hold_dur, 'psd_settle', cfg.psd_settle, ...
    'n_seeds', n_seeds, 'seeds', seeds, 'seed_offset', seed_offset, ...
    't_lya', t_lya, ...
    'per_seed_level_mean', per_seed_level_mean, ...
    'per_seed_level_std', per_seed_level_std, ...
    'level_mean', level_mean, 'level_std', level_std, ...
    'LLE_by_seed', LLE_by_seed, 'config', config);

save(fullfile(out_dir, 'dc_lle_results.mat'), 'dc_lle_results', '-v7.3');
fprintf('[dc_lle] saved dc_lle_results.mat\n');
copyfile([mfilename('fullpath') '.m'], out_dir);

%% Plot
% Not assigned: save_some_figs_to_folder_2 collects OPEN figures, so the handle
% is not needed here and a variable kept only to be unused is worse than none.
plot_dc_lle_summary(dc_lle_results);
if opts.save_figs
    fig_dir = fullfile(out_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'dc_lle', [], {'fig', 'png'});
    fprintf('[dc_lle] figures saved to %s\n', fig_dir);
end

%% Summary
fprintf('\n=== DC LLE Analysis Summary ===\n');
fprintf('Output: %s\n', out_dir);
fprintf('%-10s', 'DC');
fprintf('%16s', condition_names{:}); fprintf('\n');
for k = 1:nL
    fprintf('%-10.4g', cfg.dc_levels(k));
    for c = 1:n_cond
        fprintf('%16s', sprintf('%+.3f+/-%.3f', level_mean(k, c), level_std(k, c)));
    end
    fprintf('\n');
end
fprintf('Done.\n');
end

%% ======================================================================
%  ONE SEED: every condition on the same network
%  ======================================================================
function res = run_one_seed(s, seed_offset, preset_name, condition_names, ...
        extra_args, hold_dur, psd_settle, nL)
% All conditions share rng_seeds, so they share W and the stimulus. That is what
% makes the per-condition comparison paired rather than confounded by the draw.
rng_seed = s + seed_offset;
n_cond = numel(condition_names);

lvl_mean = nan(nL, n_cond);
lvl_std  = nan(nL, n_cond);
LLE      = nan(1, n_cond);
t_lya    = [];

for c = 1:n_cond
    model = build_from_preset(preset_name, condition_names{c}, ...
        extra_args{:}, 'rng_seeds', [rng_seed rng_seed]);
    model.run();

    ll = model.lya_results.local_lya;
    t  = model.lya_results.t_lya;
    if isempty(t_lya); t_lya = t; end
    LLE(c) = model.lya_results.LLE;

    for k = 1:nL
        lo  = (k-1)*hold_dur + psd_settle;
        hi  = k*hold_dur;
        sel = t > lo & t <= hi;
        if any(sel)
            lvl_mean(k, c) = mean(ll(sel));
            lvl_std(k, c)  = std(ll(sel));
        end
    end
end

res = struct('seed', s, 't_lya', t_lya, 'lvl_mean', lvl_mean, ...
    'lvl_std', lvl_std, 'LLE', LLE);
end

%% ======================================================================
%  RUN-MODE TABLE
%  ======================================================================
function cfg = dc_lle_run_config(run_mode, preset_defaults, ode_solver_override)
% Cost and protocol for one DC-LLE run. The analogue of mc_run_config, and
% separate from analysis_run_config for the reason the header gives: the
% staircase fixes T_range by construction, so the fs/T_range/LLE tuning the
% sweeps share does not apply.
%
% SEEDS ARE THE COST AXIS, and they were cut hard when this moved onto the
% paper's network: n = 500 with indegree = 100 against the old n = 50 with
% indegree = 4, and 7 conditions per seed instead of 1. 25/38/50 seeds would
% have made even 'fast' an hour.
switch run_mode
    case {'fast', 'fast2'}
        cfg = pack(5);
    case {'medium', 'medium2'}
        cfg = pack(15);
    case 'production'
        cfg = pack(30);
    otherwise
        error('run_dc_lle_analysis:badMode', ...
            ['Unknown run_mode ''%s'' (expected ''fast'', ''fast2'', ' ...
             '''medium'', ''medium2'' or ''production'').'], run_mode);
end

cfg.ode_solver = select_solver(cfg, preset_defaults, ode_solver_override);
end

function cfg = pack(n_seeds)
cfg = struct();
cfg.n_seeds = n_seeds;

% Protocol constants -- identical in every mode, so they are not knobs.
% The levels were tuned on the pre-2026-09-02 network (S_c = 0.5, mu 3/-4) and
% have NOT been retuned for the paper's operating point (S_c = 0.20, mu
% 5.5/-5.5). Whether they still span anything interesting is the first thing to
% look at in a 'fast' run.
cfg.dc_levels       = [0 0.0125 0.025 0.05 0.075 0.1 0.15 0.2 0.3];
cfg.hold_dur        = 30;      % s each DC level is held
cfg.ramp_dur        = 10;      % s ramping 0 -> dc_levels(1)
cfg.psd_settle      = 15;      % s skipped after each step before the LLE window
cfg.fs              = 400;
cfg.noise_intensity = 0;       % ONE noise source; the preset's Wiener noise

% The two integrators this analysis would use. select_solver picks between them
% from the preset, or takes an explicit override.
cfg.det_solver = 'rk4';
cfg.sde_solver = 'sra1';
end

function solver = select_solver(cfg, preset_defaults, override)
% Mirrors run_memory_capacity's rule, and analysis_run_config's before that.
%
% The old table named 'ode45' for production, which the paper preset would
% REJECT outright: sigma_u_noise = 0.025 requires a stochastic scheme. sra1
% costs two drift evaluations against rk4's four, so this is also cheaper.
if ~isempty(override)
    solver = override;
else
    is_stochastic = isstruct(preset_defaults) && ...
        isfield(preset_defaults, 'sigma_u_noise') && ...
        any(preset_defaults.sigma_u_noise(:) > 0);
    if is_stochastic
        solver = cfg.sde_solver;
    else
        solver = cfg.det_solver;
    end
end

sigma = preset_field(preset_defaults, 'sigma_u_noise', 0);
check_noise_settings(sigma, solver, 'run_dc_lle_analysis');
end

%% ======================================================================
function v = preset_field(s, name, default_val)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default_val;
end
end

function fig = plot_dc_lle_summary(r)
% Mean local LLE vs DC level, one band per adaptation condition.
titles = srnn_condition_titles();
cond   = r.conditions;
colors = lines(numel(cond));

fig = figure('Name', 'DC LLE: local Lyapunov vs DC level (across seeds)');
hold on;
h = gobjects(1, numel(cond));
for c = 1:numel(cond)
    shade = 1 - 0.35 * (1 - colors(c, :));      % lighter band, same hue
    confplot(r.dc_levels, r.level_mean(:, c), r.level_std(:, c), ...
        r.level_std(:, c), [colors(c, :); shade]);
    h(c) = plot(nan, nan, '-', 'Color', colors(c, :), 'LineWidth', 2);
end
yline(0, '--k', 'edge of chaos', 'HandleVisibility', 'off');
hold off;

labels = cond;
for c = 1:numel(cond)
    if isKey(titles, cond{c}); labels{c} = titles(cond{c}); end
end
legend(h, labels, 'Location', 'best', 'Box', 'off');
xlabel('DC level (input units)');
ylabel('mean local \lambda_1  (\pm std across seeds)');
title(sprintf('Local Lyapunov exponent vs DC stim level  (n_{seeds} = %d)', ...
    r.n_seeds));
grid off;
end

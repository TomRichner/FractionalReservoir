function plan = burst_mode_stage_plan(cfg)
%BURST_MODE_STAGE_PLAN The runs of the opt-in stage burst_mode, validated before any compute.
%
%   plan = BURST_MODE_STAGE_PLAN(cfg)
%
% The pre-flight half of the stage 'burst_mode' (Plan 2, package E; the
% data of Figure 8 A to C and of S10 B to D). Config in, plan out; nothing
% is computed and no file is touched. run_burst_mode_stage carries it out;
% bm_run_network is one network; bm_load_runs reads what it wrote.
%
% THE QUESTION (Report 5.1, the pilot of 2026-11-21): is a population burst
% of the bursting network the moment the population mode of its
% linearization crosses the stability boundary, and does sustained input
% keep that mode stable? Per network, through a baseline period and a
% sustained step of input:
%   - the frozen-gain effective connectivity M = W diag(d) at every
%     dt_sample: its rightmost eigenvalue as a rate (lambda - 1)/tau_d
%     (1/s; the boundary is Re = 0), that mode's coherence across the
%     population, the bulk radius predicted from the gains (R_eff), the
%     mean gains and rates;
%   - the rightmost NETWORK mode of the full Jacobian at every dt_full, by
%     a dense eigendecomposition (full_jacobian_network_mode). The pilot
%     showed the frozen-gain boundary is not the criterion once the
%     adaptive variables respond: in three of eight further seeds the
%     frozen-gain mode sat at or beyond the boundary during the step with
%     no bursting, while the full Jacobian had it damped;
%   - the frozen-adaptation arm: from the last sample before the step, the
%     adaptive variables held at their values there (bm_frozen_rhs), the
%     same step applied;
%   - bursts, crossings before bursts and their lead, share of time beyond
%     the boundary before and during the step (gains frozen, full
%     Jacobian, frozen arm), population synchrony before and during.
%
% THE NETWORK, STIMULUS AND CONDITION are those of the bursting figures
% (fig_stim_engages_adaptation): cfg.bursting_preset (bursting_pairs:
% n = 50, 35 E / 15 I) under its most-adapted condition (build_from_preset
% with an empty condition: sfa3_std1), a DC staircase [0 dc_level] held
% hold_dur s per level from t = 0 to every neuron with that figure's white
% noise probe (dc_staircase_stimulus), u_ex_scale 1, at the bundle's fs and
% integrator. The step is at t = hold_dur. The pilot used the same network
% and step but noise-free, rk4 at 400 Hz and 10 s of settling before t = 0
% (bursting_condition_model); its numbers are therefore not this stage's.
%
% THE SIZES, by run mode: n_seeds 3 at 'essential', 9 from 'fast' up
% (the pilot's nine networks). Everything else is the same in every mode.
%
% THE KNOBS (cfg.burst_mode_<name>; default)
%   seeds        m x 2 rng_seeds: [42 42] (the bursting figures' seed; the
%                stage's first seed is the one Figure 8 draws), then
%                [1 1] ... [8 8] (the pilot's further seeds)
%   n_seeds      the table's
%   dc_level     0.4      the bursting figure's second level
%   hold_dur     50 s     per level (the figure's); the step is at hold_dur
%   noise_intensity  0.001  the figure's white-noise probe on u
%   transient_s  10 s     skipped at the start of each level: the windows
%                         are before = [transient_s, hold_dur] and during =
%                         [hold_dur + transient_s, 2 hold_dur]
%   dt_sample    0.01 s   frozen-gain series, rates, adaptation traces
%   dt_full      0.05 s   full-Jacobian series (a multiple of dt_sample);
%                         the keep_ev seeds also get it every dt_sample
%                         from 1 s before to 0.5 s after each burst onset
%                         before the step (Figure 8 B)
%   x_thr        0.1      dendritic fraction of a network mode
%   burst        struct('thresh', 0.2, 'min_dur', 0.1, 'merge', 0.2): the
%                bursting pilots' rule on the mean E rate
%   tol_pre, tol_post   0.6 s, 0.1 s: a crossing coincides with a burst when
%                it overlaps [onset - tol_pre, offset + tol_post] (pilot)
%   frozen_pre_s 2 s      the frozen arm's series starts this long before
%                         the step (it equals the intact run there)
%   keep_ev      1        the first keep_ev seeds keep every eigenvalue of
%                         M at every sample (Figure 8 C draws the first)
%   fs           the bundle's cfg.fs (stage_fs), else 500 Hz
% THE INTEGRATOR is the bundle's cfg.ode_solver (fixed step:
% ospec_bundle_solver), else build_from_preset's rule (rk4: the preset is
% noise-free). The frozen arm integrates on with the model's own
% integrate(), which needs a fixed-step scheme.
%
% COST (measured on the workstation, 2026-11-21): see run_burst_mode_stage.
%
% THE PLAN: plan.folder 'burst_mode'; plan.runs (seed, file, keep_ev);
% every knob as resolved; plan.preset_name, plan.condition, plan.run_mode,
% plan.ode_solver, plan.verbose.
%
% See also: run_burst_mode_stage, bm_run_network, bm_load_runs,
%           fig_burst_mode, fig_supp10_grouped, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

SIZES = {'essential', 3; 'fast', 9; 'fastmedium', 9; 'medium', 9; 'medium2', 9; 'production', 9};
KNOBS = {'seeds', 'n_seeds', 'dc_level', 'hold_dur', 'noise_intensity', 'transient_s', 'dt_sample', ...
    'dt_full', 'x_thr', 'burst', 'tol_pre', 'tol_post', 'frozen_pre_s', 'keep_ev', 'fs'};
PREFIX = 'burst_mode_';

assert(isequal(sort(SIZES(:, 1)), sort(run_mode_names())'), 'burst_mode_stage_plan:TableIncomplete', ...
    'The size table must name every run mode of run_mode_names.');
if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, SIZES(:, 1))
    bad_config('The burst_mode stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(SIZES(:, 1)', ', '));
end
names = fieldnames(cfg);
names = names(startsWith(names, PREFIX));
unknown = setdiff(names, strcat(PREFIX, KNOBS));
if ~isempty(unknown)
    bad_config('Unknown setting(s) %s. The knobs of the burst_mode stage are %s%s.', ...
        strjoin(unknown', ', '), PREFIX, strjoin(KNOBS, [', ' PREFIX]));
end
knob = @(name, default) override(cfg, [PREFIX name], default);

%% The network: the bundle's bursting preset under its most-adapted condition
preset_name = field_or_empty(cfg, 'bursting_preset');
if isempty(preset_name); preset_name = 'bursting_pairs'; end
[~, model_class, conditions] = srnn_param_preset(preset_name);   % errors on an unknown name
if ~strcmp(model_class, 'SRNNCellTypePairs')
    bad_config('The burst_mode stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        preset_name, model_class);
end
condition = full_adaptation_condition(conditions);

%% The knobs
row = SIZES(strcmp(SIZES(:, 1), cfg.run_mode), :);
P = struct();
P.n_seeds         = knob('n_seeds', row{2});
P.dc_level        = knob('dc_level', 0.4);
P.hold_dur        = knob('hold_dur', 50);
P.noise_intensity = knob('noise_intensity', 0.001);
P.transient_s     = knob('transient_s', 10);
P.dt_sample       = knob('dt_sample', 0.01);
P.dt_full         = knob('dt_full', 0.05);
P.x_thr           = knob('x_thr', 0.1);
P.burst           = knob('burst', struct('thresh', 0.2, 'min_dur', 0.1, 'merge', 0.2));
P.tol_pre         = knob('tol_pre', 0.6);
P.tol_post        = knob('tol_post', 0.1);
P.frozen_pre_s    = knob('frozen_pre_s', 2);
P.keep_ev         = knob('keep_ev', 1);
P.fs              = knob('fs', stage_fs(field_or_empty(cfg, 'fs'), 500, 'burst_mode_stage_plan'));

positive = {'dc_level', 'hold_dur', 'dt_sample', 'dt_full', 'tol_pre', 'fs'};
for k = 1:numel(positive)
    v = P.(positive{k});
    if ~(is_real_scalar(v) && v > 0)
        bad_config('cfg.%s%s must be a positive scalar.', PREFIX, positive{k});
    end
end
nonneg = {'noise_intensity', 'transient_s', 'tol_post', 'frozen_pre_s'};
for k = 1:numel(nonneg)
    v = P.(nonneg{k});
    if ~(is_real_scalar(v) && v >= 0)
        bad_config('cfg.%s%s must be a scalar of at least 0.', PREFIX, nonneg{k});
    end
end
if ~(is_whole(P.n_seeds) && P.n_seeds >= 1)
    bad_config('cfg.%sn_seeds must be a whole number of at least 1.', PREFIX);
end
if ~(is_whole(P.keep_ev) && P.keep_ev >= 0)
    bad_config('cfg.%skeep_ev must be a whole number of at least 0.', PREFIX);
end
if ~(is_real_scalar(P.x_thr) && P.x_thr > 0 && P.x_thr < 1)
    bad_config('cfg.%sx_thr must be a fraction between 0 and 1.', PREFIX);
end
if ~(isstruct(P.burst) && all(isfield(P.burst, {'thresh', 'min_dur', 'merge'})))
    bad_config('cfg.%sburst must be a struct with thresh, min_dur and merge.', PREFIX);
end
if P.transient_s >= P.hold_dur
    bad_config('cfg.%stransient_s (%g s) must be shorter than hold_dur (%g s).', PREFIX, P.transient_s, P.hold_dur);
end
if P.frozen_pre_s >= P.hold_dur
    bad_config('cfg.%sfrozen_pre_s must be shorter than hold_dur.', PREFIX);
end
step_s = P.dt_sample * P.fs; step_f = P.dt_full / P.dt_sample;
if abs(step_s - round(step_s)) > 1e-9 || round(step_s) < 1
    bad_config('cfg.%sdt_sample (%g s) must be a whole number of integration steps (1/fs = %g s).', ...
        PREFIX, P.dt_sample, 1 / P.fs);
end
if abs(step_f - round(step_f)) > 1e-9 || round(step_f) < 1
    bad_config('cfg.%sdt_full (%g s) must be a whole multiple of dt_sample (%g s).', PREFIX, P.dt_full, P.dt_sample);
end

%% The seeds
seeds = knob('seeds', [42 42; (1:8)' * [1 1]]);
if ~(isnumeric(seeds) && ismatrix(seeds) && size(seeds, 2) == 2 && all(seeds(:) == round(seeds(:))) ...
        && all(seeds(:) >= 0) && size(unique(seeds, 'rows'), 1) == size(seeds, 1))
    bad_config('cfg.%sseeds must be an m x 2 matrix of distinct pairs of whole numbers (rng_seeds).', PREFIX);
end
if size(seeds, 1) < P.n_seeds
    bad_config('cfg.%sn_seeds is %d but cfg.%sseeds has %d row(s).', PREFIX, P.n_seeds, PREFIX, size(seeds, 1));
end
P.seeds = seeds(1:P.n_seeds, :);

runs = struct('seed', {}, 'file', {}, 'keep_ev', {});
for s = 1:P.n_seeds
    runs(end + 1) = struct('seed', P.seeds(s, :), ...
        'file', sprintf('burst_mode_s%d-%d.mat', P.seeds(s, 1), P.seeds(s, 2)), 'keep_ev', s <= P.keep_ev); %#ok<AGROW>
end

[ode_solver, problem] = ospec_bundle_solver(cfg, 'burst_mode_stage_plan');
if ~isempty(problem); bad_config('%s', problem); end
verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end

plan = P;
plan.folder      = 'burst_mode';
plan.preset_name = preset_name;
plan.condition   = condition;
plan.run_mode    = cfg.run_mode;
plan.ode_solver  = ode_solver;
plan.verbose     = verbose;
plan.runs        = runs;
end

%% ------------------------------------------------------------------------
function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function tf = is_real_scalar(v)
tf = isscalar(v) && isnumeric(v) && isreal(v) && isfinite(v);
end

function tf = is_whole(v)
tf = is_real_scalar(v) && v == round(v);
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('burst_mode_stage_plan:BadConfig', varargin{:});
end

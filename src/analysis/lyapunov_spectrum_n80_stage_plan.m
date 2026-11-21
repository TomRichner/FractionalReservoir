function plan = lyapunov_spectrum_n80_stage_plan(cfg)
% LYAPUNOV_SPECTRUM_N80_STAGE_PLAN The run a bundle's lyapunov_spectrum_n80 stage makes, validated before any compute.
%
%   plan = LYAPUNOV_SPECTRUM_N80_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'lyapunov_spectrum_n80'
% (paper_opt_in_stages; FD-032, Tom, 2026-11-19: "it should get made when we
% do a run. It should use the 80 neuron network and use topk = 240"). The
% top 240 Lyapunov exponents per condition on the 80-neuron version of the
% bundle's network, by run_lyapunov_spectrum, for the ranked spectrum and
% its histogram (fig_lyapunov_ranked_spectrum). Config in, plan out; nothing
% is computed and no file is touched. run_lyapunov_spectrum_n80_stage
% carries the plan out.
%
% WHY 80 NEURONS. At 500 neurons the leading exponents of a stable adapting
% network are a band of slow single-neuron adaptation variables that top-K
% never gets past; at 80 neurons 240 exponents reach well into the
% spectrum (scripts/explorations/small_network/).
%
% THE CONFIG
%   cfg.lyapunov_spectrum_n80_preset   the 80-neuron network (in the
%                     mu7chain bundle, ..._mu7chain_n80); a noisy
%                     SRNNCellTypePairs preset. Required.
%   cfg.run_mode      sets the number of seeds (the table below); nothing
%                     else: the network, K and window are the same in every
%                     mode
%   cfg.run_conditions   the conditions simulated (bundle_run_conditions);
%                     empty = every condition the preset states
%   cfg.lyapunov_spectrum_n80   optional struct of overrides, for tests:
%     n_seeds, K, T, lya_T_interval, lya_warmup. A field the stage does not
%     know is refused, so a typo cannot be ignored.
%
% THE SETTINGS, in every mode (the lead's decisions under FD-032):
%   T_range          [0 140] s, from the network's initial state
%   lya_T_interval   [80 140] s: exponents accumulated over 60 s
%   lya_warmup       20 s: the tangent basis aligns over [60 80] s, so the
%                    network has settled for 60 s before iteration starts.
%                    In the exploration a 10 s window from a cold start put
%                    the top band of exponents about 0.04 1/s too low; this
%                    window agreed with its own half-window estimate to
%                    0.007 1/s.
%   K                240, capped at each condition's state dimension N
%                    (no adaptation has N = 80 and gets all 80)
%   K_max            = K: NO RETRY DOUBLING. The figure compares ranks across
%                    seeds and conditions at one K. The doubling exists to
%                    resolve D_KY when the top K exponents still sum
%                    positive; where N <= 240 the whole spectrum is taken
%                    and there is nothing to double into, and where N > 240
%                    the adapting conditions are stable (lambda_1 < 0, D_KY
%                    resolved at 0). A doubling to 480 would also double the
%                    cost of exactly the runs that are already the slowest.
%   fs               the bundle's cfg.fs (stage_fs; 500 Hz in the mu7chain
%                    bundles), else 500 Hz
%   integrator       the bundle's cfg.ode_solver (stage_solver; 'sra1' in
%                    the mu7chain bundles), else the one the preset selects
%                    (sra1: the preset is stochastic; build_from_preset).
%                    Both since the merge with the solver plumbing
%                    (2026-11-21); in the mu7chain bundles the values are
%                    the ones the stage always used.
%   lya_dt           0.05 s (run_lyapunov_spectrum's)
%   noise            the preset's only; no noise-off variant
% EVERY RUN USES THE ONE T_range. Runs of different length on one seed are
% different noise realizations (wiener_increments draws the second Gaussian
% of every step after the first over the whole run), so a run is comparable
% to another only at the same T.
%
% THE SEEDS, by run mode:
%   essential, fast                     1   (network seed pair [1 2])
%   fastmedium, medium(2), production   3   ([1 2], [2 3], [3 4])
%
% THE PLAN
%   plan.folder   where it writes, relative to the run directory
%                 ('lyapunov_spectrum_n80')
%   plan.args     name-value arguments for run_lyapunov_spectrum, all but
%                 out_dir
%   plan.n_seeds, plan.K, plan.T, plan.lya_T_interval, plan.lya_warmup,
%   plan.fs, plan.ode_solver ('' = the preset's), plan.conditions,
%   plan.preset_name
%                 what the run will do, for the stage's one-line report
%
% See also: run_lyapunov_spectrum_n80_stage, run_lyapunov_spectrum,
%           fig_lyapunov_ranked_spectrum, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

%            mode           seeds
SEEDS = {   'essential',    1; ...
            'fast',         1; ...
            'fastmedium',   3; ...
            'medium',       3; ...
            'medium2',      3; ...
            'production',   3};
K_DEFAULT = 240;
T_DEFAULT = 140;
WINDOW_DEFAULT = [80 140];
WARMUP_DEFAULT = 20;
FS = 500;
ALLOWED = {'n_seeds', 'K', 'T', 'lya_T_interval', 'lya_warmup'};

assert(isequal(sort(SEEDS(:, 1)), sort(run_mode_names())'), 'lyapunov_spectrum_n80_stage_plan:TableIncomplete', ...
    'The seed table must name every run mode of run_mode_names.');
if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, SEEDS(:, 1))
    bad_config('The lyapunov_spectrum_n80 stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(SEEDS(:, 1)', ', '));
end
preset_name = field_or_empty(cfg, 'lyapunov_spectrum_n80_preset');
if isempty(preset_name)
    bad_config(['The lyapunov_spectrum_n80 stage runs on cfg.lyapunov_spectrum_n80_preset, ' ...
        'which this bundle does not set.']);
end

opts = struct();
if isfield(cfg, 'lyapunov_spectrum_n80') && ~isempty(cfg.lyapunov_spectrum_n80)
    opts = cfg.lyapunov_spectrum_n80;
end
if ~(isstruct(opts) && isscalar(opts))
    bad_config('cfg.lyapunov_spectrum_n80 must be a scalar struct.');
end
unknown = setdiff(fieldnames(opts), ALLOWED);
if ~isempty(unknown)
    bad_config('cfg.lyapunov_spectrum_n80 has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(ALLOWED, ', '));
end

%% The size: the mode's seeds, the fixed settings, then the overrides
row = SEEDS(strcmp(SEEDS(:, 1), cfg.run_mode), :);
n_seeds    = override(opts, 'n_seeds', row{2});
K          = override(opts, 'K', K_DEFAULT);
T          = override(opts, 'T', T_DEFAULT);
window     = override(opts, 'lya_T_interval', WINDOW_DEFAULT);
warmup     = override(opts, 'lya_warmup', WARMUP_DEFAULT);
if ~(is_whole(n_seeds) && n_seeds >= 1) || ~(is_whole(K) && K >= 1)
    bad_config('cfg.lyapunov_spectrum_n80.n_seeds and K must be whole numbers of at least 1.');
end
if ~(isscalar(T) && isnumeric(T) && T > 0)
    bad_config('cfg.lyapunov_spectrum_n80.T must be a positive scalar (seconds).');
end
if ~(isnumeric(window) && numel(window) == 2 && window(1) >= 0 && window(1) < window(2) && window(2) <= T)
    bad_config('The accumulation window must be [t1, t2] with 0 <= t1 < t2 <= T = %g s.', T);
end
if ~(isscalar(warmup) && isnumeric(warmup) && warmup >= 0 && warmup <= window(1))
    bad_config('The warm-up must be a scalar in [0, t1] = [0, %g] s.', window(1));
end

%% The network: a noisy SRNNCellTypePairs preset, and the bundle's conditions
[preset, model_class, conditions] = preset_with_conditions(preset_name, bundle_run_conditions(cfg));
if ~strcmp(model_class, 'SRNNCellTypePairs')
    bad_config('The lyapunov_spectrum_n80 stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        preset_name, model_class);
end
if ~(isfield(preset, 'sigma_u_noise') && any(preset.sigma_u_noise(:) > 0))
    bad_config(['The lyapunov_spectrum_n80 stage computes the spectrum with the preset''s noise on; ' ...
        'preset ''%s'' has none.'], preset_name);
end
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
plan = struct();
plan.folder         = 'lyapunov_spectrum_n80';
plan.preset_name    = preset_name;
plan.n_seeds        = n_seeds;
plan.K              = K;
plan.T              = T;
plan.lya_T_interval = window;
plan.lya_warmup     = warmup;
% The bundle's rate and integrator (cfg.fs, cfg.ode_solver; Tom, 2026-11-17:
% every paper simulation at SRA1 and 500 Hz), else the stage's own.
fs = stage_fs(field_or_empty(cfg, 'fs'), FS, 'lyapunov_spectrum_n80_stage_plan');
solver = stage_solver(field_or_empty(cfg, 'ode_solver'), '', 'lyapunov_spectrum_n80_stage_plan');
plan.fs             = fs;
plan.ode_solver     = solver;
plan.conditions     = cond_names;
plan.args = {'preset_name', preset_name, 'run_mode', cfg.run_mode, 'verbose', verbose, ...
    'conditions', cond_names, 'n_seeds', n_seeds, 'K', K, 'K_max', K, ...
    'T', T, 'lya_T_interval', window, 'lya_warmup', warmup, 'fs', fs, ...
    'noise_off_too', false};
if ~isempty(solver); plan.args = [plan.args, {'ode_solver', solver}]; end
end

%% ------------------------------------------------------------------------
function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function tf = is_whole(v)
tf = isscalar(v) && isnumeric(v) && isfinite(v) && v == round(v);
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('lyapunov_spectrum_n80_stage_plan:BadConfig', varargin{:});
end

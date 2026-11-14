function plan = noise_switch_stage_plan(cfg)
% NOISE_SWITCH_STAGE_PLAN The run a bundle's noise_switch stage makes, validated before any compute.
%
%   plan = NOISE_SWITCH_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'noise_switch'
% (paper_opt_in_stages): copies of one network that share the noise for a
% while and then run on noise seeds of their own (run_noise_switch_rollouts,
% Tom, 2026-11-14). Config in, plan out; nothing is computed and no file is
% touched. run_noise_switch_stage carries the plan out.
%
% THE CONFIG
%   cfg.preset_name   the network: the bundle's NOISY preset (the stage has
%                     nothing to share or stop sharing at sigma_u_noise = 0)
%   cfg.run_mode      sets the size (the table below) and the integrator and
%                     rate (the sweeps', analysis_run_config)
%   cfg.noise_switch  optional struct; every field is optional, and one that
%                     is absent or empty leaves the default:
%     n_copies        copies of the network
%     t_rollout       seconds each copy runs on its own seed after the switch
%     t_shared        seconds of identical noise before the switch
%     conditions      names, in order; default every condition the preset
%                     states (the figures draw the set's own selection)
%     model_args      name-value overrides of the preset (for tests on a
%                     small network)
%   A field the stage does not know is refused, so a typo cannot be ignored.
%
% THE SIZES, one row per run mode (the stage's own table):
%   mode          copies   shared   rollout
%   essential        4      10 s      5 s
%   fast             8      10 s     10 s
%   medium(2)       16      10 s     20 s
%   production      32      10 s     20 s
% The first run (8 copies, 5 s) showed several blocks still rising at 5 s,
% hence the longer rollouts above 'essential'.
%
% WHICH CONDITIONS. Every one the preset states, as every analysis of the
% bundle runs them all: one run directory serves both figure sets, each of
% which draws its own selection (cfg.conditions) from it.
%
% THE PLAN
%   plan.folder   where it writes, relative to the run directory ('noise_switch')
%   plan.args     name-value arguments for run_noise_switch_rollouts, all but
%                 out_dir
%   plan.n_copies, plan.t_shared, plan.t_rollout, plan.conditions
%                 what the run will do, for the stage's one-line report
%
% See also: run_noise_switch_stage, run_noise_switch_rollouts,
%           fig_noise_switch_rollouts, paper_opt_in_stages

arguments
    cfg (1,1) struct
end

%            mode          copies  t_shared  t_rollout
SIZES = {   'essential',    4,      10,       5; ...
            'fast',         8,      10,       10; ...
            'medium',       16,     10,       20; ...
            'medium2',      16,     10,       20; ...
            'production',   32,     10,       20};
ALLOWED = {'n_copies', 't_rollout', 't_shared', 'conditions', 'model_args'};

if ~isfield(cfg, 'run_mode') || ~ismember(cfg.run_mode, SIZES(:, 1))
    bad_config('The noise_switch stage has no size for run mode ''%s''; it has %s.', ...
        char(string(field_or_empty(cfg, 'run_mode'))), strjoin(SIZES(:, 1)', ', '));
end
assert(isequal(sort(SIZES(:, 1)), sort(run_mode_names())'), 'noise_switch_stage_plan:TableIncomplete', ...
    'The size table must name every run mode of run_mode_names.');
if ~isfield(cfg, 'preset_name') || isempty(cfg.preset_name)
    bad_config('The noise_switch stage runs on cfg.preset_name, which this bundle does not set.');
end

opts = struct();
if isfield(cfg, 'noise_switch') && ~isempty(cfg.noise_switch); opts = cfg.noise_switch; end
if ~(isstruct(opts) && isscalar(opts))
    bad_config('cfg.noise_switch must be a scalar struct.');
end
unknown = setdiff(fieldnames(opts), ALLOWED);
if ~isempty(unknown)
    bad_config('cfg.noise_switch has unknown field(s) %s. Allowed: %s.', ...
        strjoin(unknown', ', '), strjoin(ALLOWED, ', '));
end

%% The size: the mode's row, then the bundle's overrides
row = SIZES(strcmp(SIZES(:, 1), cfg.run_mode), :);
n_copies  = override(opts, 'n_copies',  row{2});
t_shared  = override(opts, 't_shared',  row{3});
t_rollout = override(opts, 't_rollout', row{4});
if ~(isscalar(n_copies) && isnumeric(n_copies) && n_copies >= 2 && n_copies == round(n_copies))
    bad_config('cfg.noise_switch.n_copies must be a whole number of at least 2.');
end
if ~(isscalar(t_shared) && isnumeric(t_shared) && t_shared > 0) || ...
        ~(isscalar(t_rollout) && isnumeric(t_rollout) && t_rollout > 0)
    bad_config('cfg.noise_switch.t_shared and t_rollout must be positive scalars (seconds).');
end

%% The network: a noisy SRNNCellTypePairs preset, and conditions it states
[preset, model_class, conditions] = srnn_param_preset(cfg.preset_name);   % errors on an unknown name
if ~strcmp(model_class, 'SRNNCellTypePairs')
    bad_config('The noise_switch stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        cfg.preset_name, model_class);
end
model_args = override(opts, 'model_args', {});
if ~iscell(model_args)
    bad_config('cfg.noise_switch.model_args must be a cell of name-value pairs.');
end
sigma = noise_level(preset, model_args);
if ~(sigma > 0)
    bad_config(['The noise_switch stage needs noise; preset ''%s'' has sigma_u_noise = %g ' ...
        '(the stage runs on cfg.preset_name, the noisy network, not the memory-capacity twin).'], ...
        cfg.preset_name, sigma);
end
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
cond_names = override(opts, 'conditions', names);
cond_names = cellstr(cond_names);
missing = setdiff(cond_names, names);
if ~isempty(missing)
    bad_config('cfg.noise_switch.conditions names %s, which preset ''%s'' does not state.', ...
        strjoin(missing, ', '), cfg.preset_name);
end

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
plan = struct();
plan.folder     = 'noise_switch';
plan.n_copies   = n_copies;
plan.t_shared   = t_shared;
plan.t_rollout  = t_rollout;
plan.conditions = cond_names;
plan.args = {'preset_name', cfg.preset_name, 'run_mode', cfg.run_mode, ...
    'conditions', cond_names, 'n_copies', n_copies, 't_shared', t_shared, ...
    't_rollout', t_rollout, 'model_args', model_args, 'verbose', verbose};
end

%% ------------------------------------------------------------------------
function sigma = noise_level(preset, model_args)
% sigma_u_noise as the run will see it: a model_args override wins over the
% preset, which wins over the class default (no noise).
sigma = 0;
if isfield(preset, 'sigma_u_noise'); sigma = preset.sigma_u_noise; end
k = find(strcmp(model_args(1:2:end), 'sigma_u_noise'), 1, 'last');
if ~isempty(k); sigma = model_args{2 * k}; end
end

function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('noise_switch_stage_plan:BadConfig', varargin{:});
end

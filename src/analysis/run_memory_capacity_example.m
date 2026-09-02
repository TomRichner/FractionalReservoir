function mat_file = run_memory_capacity_example(opts)
% RUN_MEMORY_CAPACITY_EXAMPLE Compute half of the example memory-capacity figure.
%
%   mat_file = RUN_MEMORY_CAPACITY_EXAMPLE()
%   mat_file = RUN_MEMORY_CAPACITY_EXAMPLE('run_mode', 'fast')
%
% Runs the memory-capacity protocol ONCE per adaptation condition (a single
% trial, not the paired ensemble that run_memory_capacity does) and saves the
% per-delay R^2 plus the reconstruction traces, so Fig_memory_capacity_example
% can render without re-simulating.
%
% It is the compute half of ONE figure rather than a general analysis: it keeps
% the per-delay predictions the reconstruction panels need, which the ensemble
% run (run_memory_capacity) deliberately discards.
%
% It used to live inside that figure's own folder for exactly that reason, and
% defaulted to writing its .mat beside its own .m -- which is how
% fig_memory_capacity_example came to plot four-day-old data while every sibling
% figure used the run it was handed. Being figure-specific is not a reason to sit
% outside the analysis layer; it writes data, so it lives with the analyses and
% writes into data/.
%
% Shares the network with run_memory_capacity via the 'mc_pairs_dualStd' preset,
% so the example and the ensemble figure describe the same reservoir.
%
% WAS A SCRIPT (compute_memory_capacity_example.m) with a hardcoded 60-line
% config block that duplicated looped_memory_capacity's settings by hand -- and
% had already drifted from them (level_of_chaos 2.0 vs 2.5 at one point, a
% different tau_a_E, and its own T_wash).
%
% See also: run_memory_capacity, Fig_memory_capacity_example, SRNN_ESN_reservoir

arguments
    opts.preset_name (1,:) char    = 'mc_pairs_dualStd'
    opts.run_mode    (1,:) char    = 'production'
    opts.output_dir  (1,:) char    = ''      % '' -> data/mc_example/
    opts.verbose     (1,1) logical = true
    % Empty -> deterministic at sigma = 0, stochastic above it. Name one to
    % force it; 'sra1' is legal at sigma = 0. Same rule as run_memory_capacity.
    opts.ode_solver  (1,:) char    = ''
end

setup_paths();

[preset, model_class, conditions] = srnn_param_preset(opts.preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('run_memory_capacity_example:BadModelClass', ...
        'Preset ''%s'' is written for %s; SRNN_ESN_reservoir needs SRNNCellTypePairs.', ...
        opts.preset_name, model_class);
end

% Protocol settings. Deliberately a SEPARATE, smaller table from
% run_memory_capacity's: this is one trial for a qualitative figure, so it does
% not need that function's trial count, bootstrap or permutation machinery.
switch opts.run_mode
    case 'fast',       T_train_sec = 60;  T_test_sec = 30;  d_max_sec = 5;
    case 'medium',     T_train_sec = 300; T_test_sec = 90;  d_max_sec = 10;
    case 'production', T_train_sec = 600; T_test_sec = 150; d_max_sec = 15;
    otherwise
        error('run_memory_capacity_example:badMode', ...
            'Unknown run_mode ''%s''.', opts.run_mode);
end
fs         = 200;
T_wash_sec = 10;
T_hold     = 0.3;
input_type = 'sample_hold';
readout_signal = 'synaptic';

% Integrator: the same rule run_memory_capacity applies, kept in step with it by
% deriving from the preset rather than hardcoding 'rk4' as this file used to.
if ~isempty(opts.ode_solver)
    solver = opts.ode_solver;
elseif isfield(preset, 'sigma_u_noise') && any(preset.sigma_u_noise(:) > 0)
    solver = 'sra1';
else
    solver = 'rk4';
end
sigma_probe = 0;
if isfield(preset, 'sigma_u_noise'); sigma_probe = preset.sigma_u_noise; end
check_noise_settings(sigma_probe, solver, 'run_memory_capacity_example');

fprintf('[mc_example] preset=%s run_mode=%s T_train=%gs d_max=%gs solver=%s\n', ...
    opts.preset_name, opts.run_mode, T_train_sec, d_max_sec, solver);

% Nothing stripped from the preset: on SRNNCellTypePairs tau_a is the settable
% property and n_a is Dependent on it, so the conditions carry tau_a directly.
% (On SRNNModel2 it was the other way round and tau_a_E had to be removed.)
base_args = [ ...
    namevalue(preset), ...
    {'fs',         fs, ...
     'ode_solver', solver, ...
     'rng_seeds',  [3, 4], ...     % fixed: one representative network + stimulus
     'input_type', input_type, ...
     'T_hold',     T_hold, ...
     'T_wash',     round(T_wash_sec  * fs), ...
     'T_train',    round(T_train_sec * fs), ...
     'T_test',     round(T_test_sec  * fs), ...
     'd_max',      round(d_max_sec   * fs)}];

% From the preset, via srnn_adaptation_conditions -- the same four the ensemble
% run uses, under the project's snake_case names. srnn_condition_titles supplies
% display text at plot time.
condition_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
condition_args  = cellfun(@(c) namevalue(rmfield(c, 'name')), conditions, ...
    'UniformOutput', false);
n_cond = numel(condition_names);

%% Build every condition on one network
esn = cell(1, n_cond);
for i = 1:n_cond
    esn{i} = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
    esn{i}.build();
end
% tau_a and synapse_config are what the conditions vary, so they are expected to
% differ; the network, input weights and stimulus must be shared.
SRNN_ESN_reservoir.verify_shared_build(esn, ...
    {'tau_a', 'synapse_config'}, ...
    {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});

%% Run
MC = zeros(1, n_cond);
R2 = cell(1, n_cond);
results = cell(1, n_cond);
for i = 1:n_cond
    fprintf('\n===== CONDITION %d/%d: %s =====\n', i, n_cond, condition_names{i});
    [MC(i), R2{i}, results{i}] = esn{i}.run_memory_capacity( ...
        'readout_signal', readout_signal, 'verbose', opts.verbose);
end

fprintf('\nMemory Capacity:\n');
for i = 1:n_cond
    fprintf('  %-10s MC = %.2f\n', condition_names{i}, MC(i));
end

%% Save only what the figure needs
% The full mc_results also carries the complete state time series (u_ex, x, r,
% b, br -- each [n x n_samples], ~0.6 GB per condition). The figure uses none of
% it, so stripping keeps the .mat at a few MB instead of ~9 GB.
delay_s = results{1}.delay_s;
plot_fields = {'predictions', 't_pred', 'R2_d', 'd', 'delay_s', 'hold_len', 'MC'};
results = cellfun(@(r) keep_fields(r, plot_fields), results, 'UniformOutput', false);

settings = struct('preset_name', opts.preset_name, 'run_mode', opts.run_mode, ...
    'fs', fs, 'T_hold', T_hold, 'input_type', input_type, ...
    'readout_signal', readout_signal, 'T_train_sec', T_train_sec, ...
    'T_test_sec', T_test_sec, 'T_wash_sec', T_wash_sec, 'd_max_sec', d_max_sec);

% A standalone run writes into data/, NOT next to this file. The old default
% dropped mc_example_data.mat into the figure folder, where
% fig_memory_capacity_example then read it forever, ignoring whatever the
% pipeline had written into the run directory.
if isempty(opts.output_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'data', 'mc_example');
else
    out_dir = opts.output_dir;
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
mat_file = fullfile(out_dir, 'mc_example_data.mat');
save(mat_file, 'results', 'MC', 'R2', 'delay_s', 'condition_names', ...
    'base_args', 'condition_args', 'settings');
fprintf('\nSaved: %s\n', mat_file);
end

%% ------------------------------------------------------------------------
function nv = namevalue(s)
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f); nv{2*i-1} = f{i}; nv{2*i} = s.(f{i}); end
end

function out = keep_fields(s, fields)
% A struct with only the named fields of s that actually exist.
out = struct();
for i = 1:numel(fields)
    if isfield(s, fields{i}); out.(fields{i}) = s.(fields{i}); end
end
end

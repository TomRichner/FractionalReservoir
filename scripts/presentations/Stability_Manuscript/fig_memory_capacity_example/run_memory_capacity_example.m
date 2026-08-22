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
% Lives in the figure's own folder, not in scripts/memory_capacity/, because it
% is the compute half of ONE figure rather than a general analysis: it keeps the
% per-delay predictions that the reconstruction panels need, which the ensemble
% run deliberately discards.
%
% Shares the network with run_memory_capacity via the 'mc_esn' preset, so the
% example and the ensemble figure describe the same reservoir.
%
% WAS A SCRIPT (compute_memory_capacity_example.m) with a hardcoded 60-line
% config block that duplicated looped_memory_capacity's settings by hand -- and
% had already drifted from them (level_of_chaos 2.0 vs 2.5 at one point, a
% different tau_a_E, and its own T_wash).
%
% See also: run_memory_capacity, Fig_memory_capacity_example, SRNN_ESN_reservoir

arguments
    opts.preset_name (1,:) char    = 'mc_esn'
    opts.run_mode    (1,:) char    = 'production'
    opts.output_dir  (1,:) char    = ''      % '' -> next to this file
    opts.verbose     (1,1) logical = true
end

setup_paths();
this_dir = fileparts(mfilename('fullpath'));

[preset, model_class] = srnn_param_preset(opts.preset_name);
if ~strcmp(model_class, 'SRNNModel2')
    error('run_memory_capacity_example:BadModelClass', ...
        'Preset ''%s'' is written for %s; SRNN_ESN_reservoir needs SRNNModel2.', ...
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

fprintf('[mc_example] preset=%s run_mode=%s T_train=%gs d_max=%gs\n', ...
    opts.preset_name, opts.run_mode, T_train_sec, d_max_sec);

base_args = [ ...
    namevalue(rmfield_if(preset, {'tau_a_E'})), ...
    {'fs',         fs, ...
     'ode_solver', 'rk4', ...
     'rng_seeds',  [3, 4], ...     % fixed: one representative network + stimulus
     'input_type', input_type, ...
     'T_hold',     T_hold, ...
     'T_wash',     round(T_wash_sec  * fs), ...
     'T_train',    round(T_train_sec * fs), ...
     'T_test',     round(T_test_sec  * fs), ...
     'd_max',      round(d_max_sec   * fs)}];

condition_names = {'Baseline', 'SFA', 'STD', 'SFA+STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...
    {'n_a_E', 3, 'n_b_E', 0}, ...
    {'n_a_E', 0, 'n_b_E', 1}, ...
    {'n_a_E', 3, 'n_b_E', 1} };
n_cond = numel(condition_names);

%% Build every condition on one network
esn = cell(1, n_cond);
for i = 1:n_cond
    esn{i} = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
    esn{i}.build();
end
% tau_a_E is auto-filled per n_a_E, so it legitimately differs across conditions.
SRNN_ESN_reservoir.verify_shared_build(esn, ...
    {'n_a_E', 'n_b_E', 'tau_a_E'}, ...
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

if isempty(opts.output_dir)
    out_dir = this_dir;
else
    out_dir = opts.output_dir;
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
end
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

function s = rmfield_if(s, names)
for i = 1:numel(names)
    if isfield(s, names{i}); s = rmfield(s, names{i}); end
end
end

function out = keep_fields(s, fields)
% A struct with only the named fields of s that actually exist.
out = struct();
for i = 1:numel(fields)
    if isfield(s, fields{i}); out.(fields{i}) = s.(fields{i}); end
end
end

function idx = stim_target_indices(params, rng_seed, config)
% STIM_TARGET_INDICES The neurons a pulse stimulus drives: a seeded random subset.
%
%   idx = STIM_TARGET_INDICES(params, rng_seed, config)
%
% Inputs:
%   params   - model params (get_params): .n, and for target_types also
%              .cell_type_names and .type_indices
%   rng_seed - seed of the draw (the model passes rng_seeds(2))
%   config   - struct with optional fields
%              .neuron_frac   fraction of the target pool driven (default 0.2)
%              .target_types  cell of cell-type names forming the pool
%                             (default {}: every neuron)
%
% Output:
%   idx - sorted row vector of driven neuron indices
%
% One function for both pulse generators (generate_am_pulse_train and
% generate_paired_pulse_input), so the sine-modulated trains and the
% paired-pulse protocol drive THE SAME neurons for the same seed: the
% "stimulated site" is a fixed set, as an electrode is.
%
% Why target_types exists: in a network whose E and I populations are
% statistically identical except for the sign of their output (the sfaEI
% presets), driving E and I alike makes the evoked recurrent input to the
% rest of the network cancel to first order. Driving E only is closer to
% what a cortical pulse does at a distance (long-range axons are
% excitatory) and gives the undriven neurons a real evoked response.
%
% The global RNG state is saved and restored, so the draw leaves W, the
% noise and the initial state untouched.
%
% See also: generate_am_pulse_train, generate_paired_pulse_input

neuron_frac = 0.2;
if isfield(config, 'neuron_frac'), neuron_frac = config.neuron_frac; end
target_types = {};
if isfield(config, 'target_types'), target_types = config.target_types; end

if isempty(target_types)
    pool = 1:params.n;
else
    pool = [];
    for k = 1:numel(target_types)
        q = find(strcmp(params.cell_type_names, target_types{k}), 1);
        if isempty(q)
            error('stim_target_indices:NoSuchType', ...
                'No cell type ''%s'' (have %s).', target_types{k}, ...
                strjoin(params.cell_type_names, ', '));
        end
        pool = [pool, reshape(params.type_indices{q}, 1, [])]; %#ok<AGROW>
    end
end

state = rng;
rng(rng_seed);
pick = randperm(numel(pool), round(numel(pool) * neuron_frac));
rng(state);
idx = sort(pool(pick));
end

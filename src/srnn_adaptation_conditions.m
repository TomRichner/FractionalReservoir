function conds = srnn_adaptation_conditions(model_class, synapse_config)
% SRNN_ADAPTATION_CONDITIONS The four adaptation regimes, per model class.
%
% Returns a cell array of condition structs for ParamSpaceAnalysis2.set_conditions.
% The four NAMES are the same whichever class is being swept -- no_adaptation,
% sfa_only, std_only, sfa_and_std -- so every downstream consumer keyed on the
% name keeps working, including the condition_titles maps inside
% ParamSpaceAnalysis2's plotters. What differs is how each regime is spelled:
%
%   SRNNModel2         n_a_E / n_b_E counts. n_b_E = 1 depresses EVERY outgoing
%                      excitatory synapse; there is no way to say "E->E only".
%   SRNNCellTypePairs  n_a as a 1 x C row plus a whole synapse_config struct,
%                      which names the routes individually. That is the point of
%                      running the sweep on this class: here STD is on E->E AND
%                      I->I, and nothing else.
%
% A condition may carry any model property, not just these, and PSA excludes
% whatever the conditions actually set from resolved_defaults and from the legal
% grid axes -- so n_a and synapse_config become condition-owned automatically,
% the same way n_a_E and n_b_E are.
%
% The optional SYNAPSE_CONFIG replaces the routes the STD conditions use
% (SRNNCellTypePairs only). It exists because the depression timescales are
% physics that belongs with a parameter preset, yet synapse_config can only
% reach the model through a condition -- so a preset hands its routes in here
% rather than trying to put them in model_defaults, where they could never take
% effect. Omit it for the default E->E and I->I routes.
%
% See also: ParamSpaceAnalysis2/set_conditions, srnn_param_preset

arguments
    model_class (1,:) char
    synapse_config struct = default_std_routes()
end

switch model_class
    case 'SRNNModel2'
        if ~isequal(synapse_config, default_std_routes())
            error('srnn_adaptation_conditions:NoSynapseConfig', ...
                ['SRNNModel2 has no synapse_config -- its depression is the ' ...
                'n_b_E count plus tau_b_E_rec/tau_b_E_rel in model_defaults. ' ...
                'Named routes require SRNNCellTypePairs.']);
        end
        conds = { ...
            struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
            struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
            struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', 1), ...
            struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1) ...
            };

    case 'SRNNCellTypePairs'
        sc = synapse_config;

        % An empty struct is the documented "no routes" value: compile_synapse_config
        % treats an absent or empty mechanism as absent, so no b states are created.
        no_synapses = struct();

        % SFA on E only, three timescales; tau_a{1} is then filled by
        % complete_type_defaults as logspace(0.25, 10, 3) unless swept.
        conds = { ...
            struct('name', 'no_adaptation', 'n_a', [0 0], 'synapse_config', no_synapses), ...
            struct('name', 'sfa_only',      'n_a', [3 0], 'synapse_config', no_synapses), ...
            struct('name', 'std_only',      'n_a', [0 0], 'synapse_config', sc), ...
            struct('name', 'sfa_and_std',   'n_a', [3 0], 'synapse_config', sc) ...
            };

    otherwise
        error('srnn_adaptation_conditions:UnknownModelClass', ...
            ['No adaptation conditions defined for model class ''%s''. ' ...
            'Valid: SRNNModel2, SRNNCellTypePairs.'], model_class);
end
end

function sc = default_std_routes()
% STD on the two recurrent within-type routes. The I->I recovery is slower and
% its release constant larger than E->E's, so the two populations depress on
% genuinely different timescales.
sc = struct();
sc.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.I.I.std = struct('tau_rec', 4, 'tau_rel', 1);
end

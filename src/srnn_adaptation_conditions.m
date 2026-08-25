function conds = srnn_adaptation_conditions(model_class, synapse_config, sfa_timescales, n_cell_types)
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

% The optional SFA_TIMESCALES are the adaptation timescales the SFA regimes
% switch on, as a plain vector in seconds. It defaults to the standard ladder
% srnn_sfa_timescales(3), which is the paper's operating point. It exists for the
% single-neuron METHODS figures, which show one exaggerated timescale so the rate
% decay is legible on a single trace.
%
% IT IS THE TIMESCALES, NOT A COUNT, and that is the point. This argument used to
% be n_a_sfa, an integer, with tau_a left to the preset or to the model's
% auto-fill -- and the two then had to agree in length or validate() rejected the
% pair. Two consequences, both of which this spelling removes:
%
%   * The auto-fill was logspace(log10(0.25), log10(10), n_a), and MATLAB's
%     logspace(a,b,1) returns 10^b. Asking for ONE timescale silently gave the
%     SLOW 10 s end. Nothing reported it; the model just ran a different
%     experiment.
%   * A count cannot say WHICH timescales, so single-timescale regimes were not
%     expressible at all.
%
% The conditions therefore carry tau_a EXPLICITLY, exactly as they already carry
% synapse_config rather than a count of depression timescales. That symmetry is
% deliberate: on the synapse side n_b has always been numel(tau_rec), derived
% from the timescales the condition states.
%
% MIGRATION HAZARD, and there is no safe automatic guard for it. The third
% argument used to be a COUNT, so an old positional call
%
%   srnn_adaptation_conditions(cls, routes, 3)     % once: THREE timescales
%                                                  % now:  ONE 3-second timescale
%
% still passes validation and means something different. A heuristic ("a scalar
% integer >= 2 is probably a count") cannot be added, because a single 3 s
% timescale is exactly what the single_neuron_stf preset legitimately asks for.
% Every in-repo caller was updated with the rename; anything outside the repo
% must be checked by hand. Write srnn_sfa_timescales(K) when you mean K of them.

% The optional N_CELL_TYPES is how long the n_a rows are (SRNNCellTypePairs
% only). It defaults to 2, which every sweep preset uses. It exists because
% SRNNCellTypePairs can build a ONE-cell-type network -- the Sompolinsky panel
% and the single-neuron methods figures are all C = 1 -- and validate() requires
% numel(n_a) == n_cellTypes exactly, so a hardcoded [0 0] would reject them.
%
% SFA always lands on the FIRST type and the rest get zero, which is what the
% two-element form already meant.

arguments
    model_class (1,:) char
    synapse_config = []          % [] -> default_std_routes(); else a struct
    sfa_timescales (1,:) double {mustBePositive} = srnn_sfa_timescales(3)
    n_cell_types (1,1) double {mustBeInteger, mustBePositive} = 2
end

% [] means "the default routes". Accepting it here is what lets a caller pass
% n_a_sfa without also having to restate the routes -- and so keeps the default
% defined in exactly one place, namely default_std_routes() below.
if isempty(synapse_config)
    synapse_config = default_std_routes();
elseif ~isstruct(synapse_config)
    error('srnn_adaptation_conditions:BadSynapseConfig', ...
        'synapse_config must be a struct or [] (got %s).', class(synapse_config));
end

switch model_class
    case 'SRNNModel2'
        % SRNNModel2 is intrinsically two populations -- its whole parameter
        % vocabulary is _E / _I -- so accepting any other count would let the
        % argument lie about the model being configured.
        if n_cell_types ~= 2
            error('srnn_adaptation_conditions:NotTwoTypes', ...
                ['SRNNModel2 is a two-population class; n_cell_types = %d is ' ...
                'meaningless for it. Use SRNNCellTypePairs for any other ' ...
                'number of cell types.'], n_cell_types);
        end
        if ~isequal(synapse_config, default_std_routes())
            error('srnn_adaptation_conditions:NoSynapseConfig', ...
                ['SRNNModel2 has no synapse_config -- its depression is the ' ...
                'n_b_E count plus tau_b_E_rec/tau_b_E_rel in model_defaults. ' ...
                'Named routes require SRNNCellTypePairs.']);
        end
        % SRNNModel2 still speaks COUNTS -- n_a_E with tau_a_E auto-filled -- and
        % is deliberately left that way: it is slated for deletion, and porting
        % it to explicit timescales would be work on code that is going away.
        % The count is therefore derived from the timescales here, so both
        % classes are driven by the same argument.
        n_a_sfa = numel(sfa_timescales);
        conds = { ...
            struct('name', 'no_adaptation', 'n_a_E', 0,       'n_b_E', 0), ...
            struct('name', 'sfa_only',      'n_a_E', n_a_sfa, 'n_b_E', 0), ...
            struct('name', 'std_only',      'n_a_E', 0,       'n_b_E', 1), ...
            struct('name', 'sfa_and_std',   'n_a_E', n_a_sfa, 'n_b_E', 1) ...
            };

    case 'SRNNCellTypePairs'
        sc = synapse_config;

        % An empty struct is the documented "no routes" value: compile_synapse_config
        % treats an absent or empty mechanism as absent, so no b states are created.
        no_synapses = struct();

        % SFA on the FIRST cell type only; every other type gets no adaptation.
        % tau_a is stated EXPLICITLY rather than left to the model's auto-fill,
        % so the condition records which timescales it ran, not merely how many.
        %
        % The double brace is load-bearing: struct() given a bare cell value
        % builds a struct ARRAY with one element per cell entry. Wrapping makes
        % the cell the value of a scalar struct's field, which is what a
        % condition is.
        % n_a is NOT set: it is Dependent on tau_a and read-only, so setting it
        % would throw. The timescales are the whole statement.
        sfa_taus  = [{sfa_timescales}, repmat({zeros(1,0)}, 1, n_cell_types - 1)];
        no_taus   = repmat({zeros(1,0)}, 1, n_cell_types);
        conds = { ...
            struct('name','no_adaptation', 'tau_a',{no_taus},  'synapse_config',no_synapses), ...
            struct('name','sfa_only',      'tau_a',{sfa_taus}, 'synapse_config',no_synapses), ...
            struct('name','std_only',      'tau_a',{no_taus},  'synapse_config',sc), ...
            struct('name','sfa_and_std',   'tau_a',{sfa_taus}, 'synapse_config',sc) ...
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

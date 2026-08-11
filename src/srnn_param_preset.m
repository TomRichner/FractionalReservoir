function d = srnn_param_preset(name)
% SRNN_PARAM_PRESET Named sets of SRNNModel2 parameter overrides.
%
% Returns a struct suitable for assigning to ParamSpaceAnalysis2.model_defaults
% (or splatting into an SRNNModel2 constructor). The NAME is only a lookup key;
% it is the returned struct that reaches the model.
%
% Usage:
%   psa.model_defaults = srnn_param_preset('overconnected');
%   psa.model_defaults.fs = 200;                 % layer a tweak on top
%
% What belongs in a preset: the physics -- which network is being simulated.
% Three things deliberately do NOT:
%
%   * n_levels / n_reps -- not SRNNModel2 properties at all; they size the sweep
%     and would be rejected by validate_model_defaults. See analysis_run_config.
%   * fs / T_range / ode_solver / lya_T_interval -- these are cost/fidelity
%     knobs owned by run_mode, again via analysis_run_config. A preset that set
%     them would silently redefine what 'fast' and 'production' mean.
%   * n_a_E / n_a_I / n_b_E / n_b_I -- owned by the adaptation conditions.
%
% Swept parameters (n, f, level_of_chaos) DO belong here. A grid axis always
% overrides the preset for the sweep that varies it, while sweeps that hold that
% axis fixed use the preset's value. Since run_sensitivity_analysis builds one
% PSA per swept parameter, a preset carrying all three makes each 1-D sweep hold
% the other two at this operating point rather than at class defaults.
%
% See also: analysis_run_config, merge_struct, ParamSpaceAnalysis2, SRNNModel2

arguments
    name (1,:) char
end

switch name
    case 'default'
        % Everything at SRNNModel2's class defaults.
        d = struct();

    case 'overconnected'
        % The parameter set from scripts/tests/test_SRNN2_defaults_overconnected.m:
        % E:I of 2:1 with per-synapse inhibition doubled to compensate, a slower
        % dendritic time constant, stronger adaptation, and the piecewise sigmoid
        % centred at 0.
        S_a = 0.9;
        S_c = 0.0;
        d = struct( ...
            'n',                   300, ...
            'f',                   2/3, ...
            'indegree',            100, ...
            'level_of_chaos',      1, ...
            'S_a',                 S_a, ...
            'S_c',                 S_c, ...
            'tau_d',               1, ...
            'tau_b_E_rec',         1, ...
            'c_E',                 1/3, ...
            'mu_E_tilde_relative',  3, ...   % class default
            'mu_I_tilde_relative', -6);      % class default -4; doubled, half as many I neurons
        % Activation handles capture S_a/S_c as literals, so they must be built
        % from the same values set above.
        d.activation_function = @(x) SRNNModel2.piecewiseSigmoid(x, S_a, S_c);
        d.activation_function_derivative = @(x) SRNNModel2.piecewiseSigmoidDerivative(x, S_a, S_c);

    otherwise
        error('srnn_param_preset:UnknownPreset', ...
            'Unknown preset ''%s''. Valid presets: %s.', ...
            name, strjoin(srnn_param_preset_names(), ', '));
end
end

function names = srnn_param_preset_names()
% The valid preset names, kept next to the switch above so they stay in sync.
names = {'default', 'overconnected'};
end

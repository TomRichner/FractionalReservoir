function [d, model_class, conditions] = srnn_param_preset(name)
% SRNN_PARAM_PRESET Named sets of model parameter overrides.
%
% Returns a struct suitable for assigning to ParamSpaceAnalysis2.model_defaults
% (or splatting into a model constructor), plus the MODEL CLASS those overrides
% are written for and the ADAPTATION CONDITIONS to sweep them under.
% The NAME is only a lookup key; it is the returned struct that reaches the model.
%
% Usage:
%   [psa.model_defaults, psa.model_class] = srnn_param_preset('celltype_pairs');
%   psa.model_defaults.fs = 200;                 % layer a tweak on top
%
% The second output exists because the two model classes have disjoint parameter
% vocabularies -- 'overconnected' is meaningless to SRNNCellTypePairs and
% 'celltype_pairs' is meaningless to SRNNModel2 -- so a preset that did not carry
% its class would only fail later, inside validate_model_defaults. It defaults to
% 'SRNNModel2', so existing single-output callers are unaffected.
%
% The third is srnn_adaptation_conditions(model_class) unless the preset needs
% different depression routes. Those timescales are physics and would belong in
% the struct above, except that synapse_config can only reach the model through
% a condition -- so the preset passes them to srnn_adaptation_conditions instead.
%
% What belongs in a preset: the physics -- which network is being simulated.
% Three things deliberately do NOT:
%
%   * n_levels / n_reps -- not SRNNModel2 properties at all; they size the sweep
%     and would be rejected by validate_model_defaults. See analysis_run_config.
%   * fs / T_range / ode_solver / lya_T_interval -- these are cost/fidelity
%     knobs owned by run_mode, again via analysis_run_config. A preset that set
%     them would silently redefine what 'fast' and 'production' mean.
%   * n_a_E / n_a_I / n_b_E / n_b_I -- owned by the adaptation conditions, as
%     are their SRNNCellTypePairs counterparts n_a and synapse_config. See
%     srnn_adaptation_conditions.
%
% The nonlinearity is named data (`activation` plus S_a/S_c), not a function
% handle, so a preset cannot end up with a handle whose captured constants
% disagree with its own S_a/S_c fields. Use `activation_custom` only for a
% nonlinearity outside the three named ones.
%
% Swept parameters (n, f, level_of_chaos) DO belong here. A grid axis always
% overrides the preset for the sweep that varies it, while sweeps that hold that
% axis fixed use the preset's value. Since run_sensitivity_analysis builds one
% PSA per swept parameter, a preset carrying all three makes each 1-D sweep hold
% the other two at this operating point rather than at class defaults.
%
% See also: analysis_run_config, merge_struct, srnn_adaptation_conditions,
%           ParamSpaceAnalysis2, SRNNModel2, SRNNCellTypePairs

arguments
    name (1,:) char
end

model_class = 'SRNNModel2';
std_routes = [];        % [] = whatever srnn_adaptation_conditions defaults to

switch name
    case 'default'
        % Everything at SRNNModel2's class defaults.
        d = struct();

    case 'overconnected'
        % The parameter set from scripts/tests/test_SRNN2_defaults_overconnected.m:
        % E:I of 2:1 with per-synapse inhibition doubled to compensate, a slower
        % dendritic time constant, stronger adaptation, and the piecewise sigmoid
        % centred at 0.
        d = struct( ...
            'n',                   300, ...
            'f',                   2/3, ...
            'indegree',            100, ...
            'level_of_chaos',      1, ...
            'activation',          'piecewise', ...
            'S_a',                 0.9, ...
            'S_c',                 0.0, ...
            'tau_d',               1, ...
            'tau_b_E_rec',         1, ...
            'c_E',                 1/3, ...
            'mu_E_tilde_relative',  3, ...   % class default
            'mu_I_tilde_relative', -6);      % class default -4; doubled, half as many I neurons

    case 'celltype_pairs'
        % The operating point from scripts/tests/test_SRNNCellTypePairs_defaults.m:
        % a two-type E/I network run through SRNNCellTypePairs so that depression
        % can be put on ONE ROUTE at a time. The question it exists to ask is what
        % changes when I->I synapses depress as well as E->E, which SRNNModel2
        % cannot express (its n_b_E = 1 depresses every outgoing E synapse).
        %
        % Note the asymmetric mu block: E receives more excitation than I does
        % (4 vs 3) and both are inhibited equally (-3), which is what pushes this
        % network away from the symmetric class default.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    300, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [4 -3; 3 -3], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1 1; 1 1], ...
            'level_of_chaos',       1.3, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.1));
        % Per-neuron setpoint heterogeneity (mu_S_c / sigma_S_c) is deliberately
        % left off, matching the commented-out lines in the source script: with
        % both empty, S_c_vec stays empty and every path takes the scalar branch.

    case 'celltype_pairs_S_c_by_type'
        % scripts/tests/example_SRNNCellTypePairs_S_c_by_type.m: 'celltype_pairs'
        % with a separate nonlinearity setpoint per cell type, plus a stronger,
        % more asymmetric connectivity.
        %
        % The setpoint is where the interest is. phi is centred on S_c, so
        % raising I's to 0.25 slides its curve RIGHT and makes the I population
        % LESS excitable -- and the example's own comparison shows the network
        % does not simply turn inhibition down: the I rate barely moves while the
        % E rate more than doubles, so the whole E/I operating point shifts.
        % sigma_S_c = 0 means no cell-to-cell spread, only the two type means.
        %
        % One deliberate divergence from that script: level_of_chaos is 1.4 here,
        % not the script's 1.5. The sweeps vary it anyway; this is only the value
        % the analyses that hold it fixed operate at.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    300, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [4.5 -3; 4 -3], ...  % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1 1; 1 1], ...
            'level_of_chaos',       1.4, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % unused: mu_S_c below overrides it
            'mu_S_c',               [0.0 0.25], ...
            'sigma_S_c',            [0.0 0.0], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.1));
        % This preset's I->I release constant is 0.5, not the default 1.
        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.I.I.std = struct('tau_rec', 4, 'tau_rel', 0.5);

    case 'celltype_pairs_S_c_by_type_n500'
        % 'celltype_pairs_S_c_by_type' at n = 500 instead of 300. Derived from it
        % rather than copied, so the two cannot drift apart -- the whole point of
        % this preset is that n is the ONLY difference.
        %
        % indegree stays at 100, so alpha drops from 1/3 to 0.2. That does not
        % move the bulk radius: F_tracks_network is true by default, which makes
        % R exactly independent of n (the n*alpha cancels in get.R). What changes
        % is the network, not its scale -- sparser connectivity and a larger
        % state vector.
        %
        % Note this is the fixed-n value only. The n sensitivity sweep and the
        % param_space grid both vary n over [100, 1000] and override it; it is
        % the OTHER sweeps (f_E, level_of_chaos, mu_EE_relative, tau_a_E) that
        % actually run at 500.
        [d, model_class, conditions] = srnn_param_preset('celltype_pairs_S_c_by_type');
        d.n = 500;
        return;     % conditions already resolved by the recursive call

    case 'celltype_pairs_S_c_by_type_n500_fixedF'
        % ..._n500 with the weight scale F PINNED instead of tracking n.
        %
        % Why this exists. F = 1/sqrt(n*alpha*(2-alpha)) and n*alpha = indegree,
        % so F = 1/sqrt(indegree*(2 - indegree/n)): with F_tracks_network = true,
        % n enters only through alpha, and F drifts ~4% from n = 300 to n = 500.
        % The bigger problem is the bulk radius. R is independent of n only when
        % the relative multipliers are unity -- the n*alpha cancellation leaves
        % f*[sigma_rel^2 + (1-alpha)*mu_rel^2]/(2-alpha), which is alpha-free
        % only if mu_rel^2 = sigma_rel^2 = 1 (then 1 + (1-alpha) = 2-alpha).
        % This preset's mu_rel runs 3 to 4.5, so R moves 1.40 -> 3.73 across the
        % n sweep and "vary n" is really "vary n AND criticality".
        %
        % F_ref_n is 500, not the class default 300, so that AT the preset's own
        % operating point F is exactly what the tracking version computed. The
        % fixed-n sweeps (f_E, level_of_chaos, the mu blocks, tau_a_E) are then
        % numerically identical to ..._n500 and the change is isolated to the n
        % sweep and the param_space grid -- which is the comparison worth having.
        %
        % Pinning F does NOT freeze the network: build() still passes the real
        % alpha to the generator, so connectivity tracks the grid point. Only the
        % weight SCALE is held.
        [d, model_class, conditions] = srnn_param_preset('celltype_pairs_S_c_by_type_n500');
        d.F_tracks_network = false;
        d.F_ref_n          = 500;
        d.F_ref_indegree   = 100;
        return;     % conditions already resolved by the recursive call

    case 'celltype_pairs_all_std_n500'
        % Depression on ALL FOUR routes, SFA on E only, and a setpoint split of
        % E at 0.15 / I at 0.25. Everything else is ..._n500_fixedF: n = 500,
        % pinned F, level_of_chaos 1.4, mu_tilde_relative [4.5 -3; 4 -3].
        %
        % Where the earlier presets depressed only the two within-type routes,
        % here every synapse in the network depresses -- n_b becomes [1 1; 1 1]
        % rather than [1 0; 0 1], so the state vector gains the two cross-route
        % b variables. This is the case SRNNModel2 cannot express at all: its
        % n_b_E/n_b_I are per-PREsynaptic-population counts, so it can say "all
        % outgoing E synapses depress" but never distinguish E->E from E->I.
        %
        % The two STD time constants are set independently, and they are not the
        % same knob:
        %   tau_rec  recovery, b -> 1 at rate (1-b)/tau_rec. SMALLER = recovers
        %            sooner. E->E is 1 s against 3 s elsewhere, so excitatory
        %            recurrence bounces back while the rest stays depressed.
        %   tau_rel  release, the depression term b*r/tau_rel. SMALLER = depresses
        %            HARDER for a given rate. I->I is 0.15 against 0.25
        %            elsewhere, so inhibition-onto-inhibition is the quickest to
        %            give way -- which the mu_II sweep flagged as the strongest
        %            single driver of chaos once STD is present.
        [d, model_class] = srnn_param_preset('celltype_pairs_S_c_by_type_n500_fixedF');
        d.mu_S_c    = [0.15 0.25];   % E less excitable than before, I unchanged
        d.sigma_S_c = [0.0  0.0];    % no cell-to-cell spread, only the type means

        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', 1, 'tau_rel', 0.25);
        std_routes.E.I.std = struct('tau_rec', 3, 'tau_rel', 0.25);
        std_routes.I.E.std = struct('tau_rec', 3, 'tau_rel', 0.25);
        std_routes.I.I.std = struct('tau_rec', 3, 'tau_rel', 0.15);
        % NOTE the nesting is synapse_config.<PRE>.<POST>, so .E.I is E -> I.
        % That is the opposite order from the mu_EI naming, which is (post, pre).

    case 'celltype_pairs_uniform_std_n500'
        % The homogeneous control for celltype_pairs_all_std_n500: same four
        % depressing routes, but nothing distinguishes any of them from any
        % other, and nothing distinguishes E from I except its sign.
        %
        % Three heterogeneities are removed at once:
        %   * STD is identical on every route, tau_rec = 2 and tau_rel = 0.25,
        %     rather than E->E recovering faster and I->I depressing harder.
        %   * Both types share S_c = 0, rather than 0.15 / 0.25.
        %   * The weights are symmetric, mu_tilde_relative = [4 -4; 4 -4],
        %     rather than the asymmetric [4.5 -3; 4 -3]. Every population sends
        %     the same magnitude and receives the same magnitude; only the sign
        %     differs. Note this leaves the RMT unity condition unmet (mu_rel is
        %     4, not 1), so R still varies across the n sweep.
        %
        % The setpoint is cleared rather than set to [0 0]. Both give every
        % neuron a centre of zero, but leaving mu_S_c/sigma_S_c EMPTY keeps
        % S_c_vec empty, which holds the whole model on the scalar-S_c code
        % path -- the branch that is bit-identical to the pre-heterogeneity
        % code, and cheaper. Setting [0 0] would draw a vector of zeros and take
        % the per-neuron path to compute the same thing.
        [d, model_class] = srnn_param_preset('celltype_pairs_S_c_by_type_n500_fixedF');
        d.mu_tilde_relative = [4 -4; 4 -4];   % multiples of F, (post <- pre)
        d.S_c       = 0.0;
        d.mu_S_c    = [];
        d.sigma_S_c = [];

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_uniform_std_n500_mu5p5'
        % celltype_pairs_uniform_std_n500 with the weight means raised to a
        % magnitude of 5.5 and level_of_chaos returned to 1, so the scale lives
        % entirely in the weights and nothing multiplies W afterwards.
        %
        % NOT the same network as its parent at level_of_chaos = 1.4. That
        % equivalence needs BOTH tildes scaled -- level_of_chaos multiplies the
        % assembled W, and a weight is drawn as mu + sigma*randn, so
        % 1.4*(mu + sigma*randn) = 1.4*mu + (1.4*sigma)*randn. Here sigma stays
        % at 1, and mu goes to 5.5 rather than the 5.6 that 4 * 1.4 would give.
        % The mean-to-spread ratio therefore changes: this network has a
        % relatively stronger deterministic block structure and a relatively
        % weaker random bulk than the parent does.
        [d, model_class, conditions] = srnn_param_preset('celltype_pairs_uniform_std_n500');
        d.mu_tilde_relative = [5.5 -5.5; 5.5 -5.5];   % multiples of F, (post <- pre)
        d.level_of_chaos    = 1.0;
        return;     % conditions already resolved by the recursive call

    otherwise
        error('srnn_param_preset:UnknownPreset', ...
            'Unknown preset ''%s''. Valid presets: %s.', ...
            name, strjoin(srnn_param_preset_names(), ', '));
end

if isempty(std_routes)
    conditions = srnn_adaptation_conditions(model_class);
else
    conditions = srnn_adaptation_conditions(model_class, std_routes);
end
end

function names = srnn_param_preset_names()
% The valid preset names, kept next to the switch above so they stay in sync.
names = {'default', 'overconnected', 'celltype_pairs', ...
    'celltype_pairs_S_c_by_type', 'celltype_pairs_S_c_by_type_n500', ...
    'celltype_pairs_S_c_by_type_n500_fixedF', 'celltype_pairs_all_std_n500', ...
    'celltype_pairs_uniform_std_n500', 'celltype_pairs_uniform_std_n500_mu5p5'};
end

function ic = pairs_input_config(intrinsic_drive)
% Assigning input_config REPLACES the struct SRNNCellTypePairs.set_defaults
% built, so a preset that wants to change one field has to restate all of them.
% This mirrors that default exactly apart from intrinsic_drive.
%
% intrinsic_drive is a SCALAR here, not 0.1*ones(n,1) as in the source script:
% n is a sweep axis, so no fixed-length vector is right at more than one grid
% point. generate_external_input accepts "scalar or n-by-1".
%
% step_density is left empty; complete_type_defaults fills one 0.2 entry per
% cell type, which is what the class default would have done anyway.
ic = struct( ...
    'n_steps',          3, ...
    'step_density',     struct(), ...
    'amp',              0.5, ...
    'no_stim_pattern',  logical([1 0 1]), ...
    'intrinsic_drive',  intrinsic_drive, ...
    'positive_only',    false);
end

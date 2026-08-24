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
% its class would only fail later, inside validate_model_defaults.
%
% The third is srnn_adaptation_conditions(model_class) unless the preset needs
% different depression routes. Those timescales are physics and would belong in
% the struct above, except that synapse_config can only reach the model through
% a condition -- so the preset passes them to srnn_adaptation_conditions instead.
%
% EVERY CASE IS SELF-CONTAINED. A preset states its own model_class, its own
% complete override struct, and its own std_routes -- nothing is inherited from
% another case at run time. This file used to build ten of its fourteen presets
% by CALLING ITSELF, so reading off `level_of_chaos` for the preset behind a
% finished run meant walking up to seven recursive hops, and editing a mid-chain
% preset silently moved every descendant. Worse, whether a descendant inherited
% its conditions depended on whether its recursive call asked for two outputs or
% three -- load-bearing, and invisible at a glance.
%
% Lineage is preserved as PROSE instead: each derived preset opens with a
% "Copied from X. Changed: ..." block naming its parent and listing exactly what
% differs. That is the part worth keeping; the coupling is not. When you change
% a preset, the descendants listed under "Derived presets" in its comment are
% the ones to consider changing too -- deliberately, one at a time.
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

% No default model_class. Every case names its own, so a new SRNNCellTypePairs
% preset cannot silently inherit 'SRNNModel2' and fail much later inside
% validate_model_defaults with a wall of "not a property" messages.
std_routes = [];        % [] = whatever srnn_adaptation_conditions defaults to
n_a_sfa    = 3;         % SFA timescales the two SFA conditions switch on.
                        % Override ONLY when the preset carries its own tau_a of
                        % a different length -- the count and the timescales must
                        % agree or validate() rejects the pair.

switch name
    case 'default'
        % Everything at SRNNModel2's class defaults.
        model_class = 'SRNNModel2';
        d = struct();

    case 'overconnected'
        % The parameter set from scripts/tests/test_SRNN2_defaults_overconnected.m:
        % E:I of 2:1 with per-synapse inhibition doubled to compensate, a slower
        % dendritic time constant, stronger adaptation, and the piecewise sigmoid
        % centred at 0.
        model_class = 'SRNNModel2';
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
        % std_routes stays [] -- the default E->E and I->I routes.

    case 'celltype_pairs_S_c_by_type'
        % Copied from celltype_pairs. Changed:
        %   mu_tilde_relative  [4 -3; 3 -3] -> [4.5 -3; 4 -3]
        %   level_of_chaos     1.3 -> 1.4
        %   mu_S_c             (absent)     -> [0.0 0.25]   NEW
        %   sigma_S_c          (absent)     -> [0.0 0.0]    NEW
        %   std_routes         (default)    -> I->I tau_rel 1 -> 0.5
        % Derived presets: celltype_pairs_S_c_by_type_n500.
        %
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
        % Copied from celltype_pairs_S_c_by_type. Changed:
        %   n   300 -> 500
        % Derived presets: celltype_pairs_S_c_by_type_n500_fixedF.
        %
        % n is the ONLY difference from the parent -- keep it that way.
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
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
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
        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.I.I.std = struct('tau_rec', 4, 'tau_rel', 0.5);

    case 'celltype_pairs_S_c_by_type_n500_fixedF'
        % Copied from celltype_pairs_S_c_by_type_n500. Changed:
        %   F_tracks_network  (default true) -> false   NEW
        %   F_ref_n           (absent)       -> 500     NEW
        %   F_ref_indegree    (absent)       -> 100     NEW
        % Derived presets: celltype_pairs_all_std_n500,
        %                  celltype_pairs_uniform_std_n500.
        %
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
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
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
            'input_config',         pairs_input_config(0.1), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100);
        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.I.I.std = struct('tau_rec', 4, 'tau_rel', 0.5);

    case 'celltype_pairs_all_std_n500'
        % Copied from celltype_pairs_S_c_by_type_n500_fixedF. Changed:
        %   mu_S_c      [0.0 0.25] -> [0.15 0.25]
        %   std_routes  E->E, I->I -> all four routes, own timescales
        % Derived presets: none.
        %
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
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
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
            'mu_S_c',               [0.15 0.25], ...  % E less excitable than the parent, I unchanged
            'sigma_S_c',            [0.0  0.0], ...   % no cell-to-cell spread, only the type means
            'c',                    [0.5/3, 0], ...   % SFA scaling, E only
            'input_config',         pairs_input_config(0.1), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100);

        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', 1, 'tau_rel', 0.25);
        std_routes.E.I.std = struct('tau_rec', 3, 'tau_rel', 0.25);
        std_routes.I.E.std = struct('tau_rec', 3, 'tau_rel', 0.25);
        std_routes.I.I.std = struct('tau_rec', 3, 'tau_rel', 0.15);
        % NOTE the nesting is synapse_config.<PRE>.<POST>, so .E.I is E -> I.
        % That is the opposite order from the mu_EI naming, which is (post, pre).

    case 'celltype_pairs_uniform_std_n500'
        % Copied from celltype_pairs_S_c_by_type_n500_fixedF. Changed:
        %   mu_tilde_relative  [4.5 -3; 4 -3] -> [4 -4; 4 -4]
        %   mu_S_c             [0.0 0.25]     -> []   (cleared)
        %   sigma_S_c          [0.0 0.0]      -> []   (cleared)
        %   std_routes         E->E, I->I     -> all four routes, all 2 / 0.25
        %   (S_c restated as 0.0, unchanged in value but now the operative one)
        % Derived presets: celltype_pairs_uniform_std_n500_mu5p5.
        %
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
        % The setpoint is CLEARED rather than set to [0 0]. Both give every
        % neuron a centre of zero, but leaving mu_S_c/sigma_S_c EMPTY keeps
        % S_c_vec empty, which holds the whole model on the scalar-S_c code
        % path -- the branch that is bit-identical to the pre-heterogeneity
        % code, and cheaper. Setting [0 0] would draw a vector of zeros and take
        % the per-neuron path to compute the same thing. The empty fields are
        % still WRITTEN OUT below rather than omitted, because omitting them
        % would let a stale mu_S_c survive from elsewhere in model_defaults.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [4 -4; 4 -4], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1 1; 1 1], ...
            'level_of_chaos',       1.4, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.1), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_uniform_std_n500_mu5p5'
        % Copied from celltype_pairs_uniform_std_n500. Changed:
        %   mu_tilde_relative  [4 -4; 4 -4] -> [5.5 -5.5; 5.5 -5.5]
        %   level_of_chaos     1.4 -> 1.0
        % Derived presets: celltype_pairs_uniform_std_n500_mu5p5_nodrive.
        %
        % The weight means are raised to a magnitude of 5.5 and level_of_chaos
        % returned to 1, so the scale lives entirely in the weights and nothing
        % multiplies W afterwards.
        %
        % NOT the same network as its parent at level_of_chaos = 1.4. That
        % equivalence needs BOTH tildes scaled -- level_of_chaos multiplies the
        % assembled W, and a weight is drawn as mu + sigma*randn, so
        % 1.4*(mu + sigma*randn) = 1.4*mu + (1.4*sigma)*randn. Here sigma stays
        % at 1, and mu goes to 5.5 rather than the 5.6 that 4 * 1.4 would give.
        % The mean-to-spread ratio therefore changes: this network has a
        % relatively stronger deterministic block structure and a relatively
        % weaker random bulk than the parent does.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1 1; 1 1], ...
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.1), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_uniform_std_n500_mu5p5_nodrive'
        % Copied from celltype_pairs_uniform_std_n500_mu5p5. Changed:
        %   input_config  pairs_input_config(0.1) -> pairs_input_config(0.0)
        % Derived presets: celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5.
        %
        % ..._mu5p5 with the tonic drive removed: intrinsic_drive 0.1 -> 0.
        %
        % intrinsic_drive is the constant term added to every neuron's input for
        % the whole simulation, distinct from the stepped stimulus that
        % input_config's n_steps / amp / step_density build. Taking it to zero
        % leaves the network driven ONLY by that stimulus, so between steps it
        % runs autonomously. Whether it still does anything there is the point of
        % the comparison: at S_c = 0 the piecewise nonlinearity sits at
        % phi(0) = 0.5, so an undriven neuron is at mid-range rather than silent,
        % and the network need not fall quiet just because the drive is gone.
        %
        % input_config is stated whole -- assigning the property replaces the
        % struct, so there is no such thing as changing one field of it.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1 1; 1 1], ...
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5'
        % Copied from celltype_pairs_uniform_std_n500_mu5p5_nodrive. Changed:
        %   sigma_tilde_relative  [1 1; 1 1] -> [1.5 1.5; 1.5 1.5]
        % Derived presets: ..._sig1p5_noise0p02, ..._sig1p5_noise0p01,
        %                  celltype_pairs_Sc0p2_noise0p025.
        %
        % The weight SPREAD is raised, mu left at 5.5.
        %
        % This is the opposite lever from the mu change. mu sets the block means
        % -- the deterministic E/I structure, which is what produces the outlier
        % eigenvalues -- while sigma sets the scatter around them, which is what
        % sets the bulk radius R. Raising sigma alone therefore inflates the
        % random bulk against a fixed deterministic skeleton, the reverse of the
        % mu5p5 step, and pushes R up substantially.
        %
        % Note the property is sigma_tilde_relative, a multiplier of
        % F = 1/sqrt(n*alpha*(2-alpha)). The absolute sigma_tilde is Dependent
        % and read-only; assigning it raises SRNNModel:RenamedProperty.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1.5 1.5; 1.5 1.5], ...     % multiples of F
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p02'
        % Copied from celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5.
        % Changed:
        %   sigma_u_noise  (absent) -> 0.02   NEW
        % Derived presets: none.
        %
        % ..._sig1p5 run as an SDE: additive Wiener noise on the dendritic
        % state, dx = (...)/tau_d dt + sigma_u_noise/tau_d dW.
        %
        % WHAT 0.02 MEANS HERE. sigma_u_noise is INPUT-REFERRED -- the units of
        % u -- which usually makes it readable against intrinsic_drive. That
        % framing does not apply to this preset: the drive is 0. The meaningful
        % comparison is against the nonlinearity instead. With S_c = 0 and
        % S_a = 0.8 the piecewise sigmoid spans +/-0.6, and
        %
        %   x_noise_std = sigma_u_noise / sqrt(2*tau_d) = 0.02/sqrt(0.2) = 0.0447
        %
        % so the noise moves a neuron over about 7.5% of its input range -- a
        % perturbation of the operating point rather than a redefinition of it.
        %
        % The INTEGRATOR IS NOT NAMED HERE, deliberately. ode_solver is a
        % run_mode knob (analysis_run_config), it is on the list of names a
        % preset may not carry, and merge_struct gives cfg.model precedence so
        % naming it here would be inert anyway. What happens instead is that
        % analysis_run_config sees sigma_u_noise > 0 and selects that mode's
        % STOCHASTIC integrator ('sra1') in place of its deterministic one.
        % Carrying sigma_u_noise is therefore the whole of what marks this
        % preset stochastic.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1.5 1.5; 1.5 1.5], ...     % multiples of F
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.02);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01'
        % Copied from celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5.
        % Changed:
        %   sigma_u_noise  (absent) -> 0.01   NEW
        % Derived presets: none. Sibling of ..._noise0p02, which uses 0.02.
        %
        % HALF the noise of the noise0p02 sibling, so the two bracket the
        % amplitude rather than testing a single value. Also the preset the
        % longer/finer 'medium2' runs use.
        %
        % x_noise_std = 0.01/sqrt(2*tau_d) = 0.0224, about 3.7% of the piecewise
        % sigmoid's 0.6 half-width at S_c = 0, against 7.5% at 0.02 -- a light
        % perturbation of the operating point. At 0.02 the same network crossed
        % from chaotic to contracting over a fast2 sweep, so half of that is
        % deliberately on the gentler side of the transition. Do not expect the
        % effect to halve with sigma, though: noise enters the exponent
        % quadratically more often than linearly, which is the reason to keep
        % both rather than interpolate between them.
        %
        % Like its sibling it names no integrator: analysis_run_config sees
        % sigma_u_noise > 0 and picks the mode's stochastic scheme ('sra1').
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1.5 1.5; 1.5 1.5], ...     % multiples of F
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.0, ...     % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.01);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_Sc0p2_noise0p025'
        % Copied from celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5.
        % Changed:
        %   S_c            0.0     -> 0.20
        %   sigma_u_noise  (absent) -> 0.025   NEW
        % Derived presets: celltype_pairs_Sc0p2_noise0p025_dualStd.
        %
        % The hand-tuned operating point from
        % scripts/examples/example_SRNNCellTypePairs_from_preset.m, promoted out
        % of that script's pre-build overrides so the sweeps can use it. This is
        % the preset behind the production run in
        % data/param_space/run_all_aug_14_26_17_25.
        %
        % The example also assigned level_of_chaos = 1.0 and reasserted
        % tau_rec = 2 / tau_rel = 0.25 on every route; those match the values
        % above already, so they were no-ops there. plot_deci was set there too
        % and is deliberately NOT carried: it decimates plot_data only, changes
        % no physics, and is one of the fields same_config ignores when deciding
        % whether two runs may be pooled.
        %
        % S_c = 0.2 RAISES THE THRESHOLD with the drive still at zero. phi is
        % centred on S_c, so an unstimulated neuron sits at phi(0) rather than
        % phi(S_c) = 0.5: with S_a = 0.8 the linear branch runs over
        % [S_c - 0.4, S_c + 0.4], x = 0 falls inside it, and
        % phi(0) = (0 - S_c) + 0.5 = 0.3. Still active, so the network does not
        % go quiet the way it would at a setpoint above 0.5 -- but each neuron
        % now starts further from saturation than the S_c = 0 parent.
        %
        % sigma_u_noise = 0.025 is above BOTH sibling noise presets:
        % x_noise_std = 0.025/sqrt(2*tau_d) = 0.0559, about 9.3% of the
        % piecewise half-width, against 7.5% at 0.02 and 3.7% at 0.01. The three
        % together span a factor of 2.5 in amplitude.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1.5 1.5; 1.5 1.5], ...     % multiples of F
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.20, ...    % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.025);

        std_routes = struct();
        uniform_std = struct('tau_rec', 2, 'tau_rel', 0.25);
        std_routes.E.E.std = uniform_std;
        std_routes.E.I.std = uniform_std;
        std_routes.I.E.std = uniform_std;
        std_routes.I.I.std = uniform_std;

    case 'celltype_pairs_Sc0p2_noise0p025_dualStd'
        % Copied from celltype_pairs_Sc0p2_noise0p025. Changed:
        %   STD, every route   tau_rec  2     -> [2 4]        (n_b 1 -> 2)
        %                      tau_rel  0.25  -> [0.25 0.5]
        % Derived presets: none.
        %
        % A SECOND, SLOWER DEPRESSION TIMESCALE on each of the four routes. The
        % original pair (2, 0.25) is kept as the first element, so this preset
        % adds a slow component to the depression already in use rather than
        % replacing it. compile_synapse_config reads n_b off numel(tau_rec), so
        % a 2-element row IS the request for two timescales; the b-states are
        % integrated as columns and enter the synapse as prod(b, 2).
        %
        % DEEPER DEPRESSION IS THE POINT, and the second timescale is how it is
        % obtained. Both pairs share the same ratio tau_rec/tau_rel = 8, so both
        % b-states relax toward the same steady state 1/(1 + 8r) and differ only
        % in how fast they get there; the synapse multiplies them, so the
        % steady-state depression is the SQUARE of the parent's:
        %
        %   parent  b     = 1/(1 + 8r)        = 0.29  at r = 0.3
        %   here    b1*b2 = 1/(1 + 8r)^2      = 0.086 at r = 0.3
        %
        % i.e. about 3.4x deeper as well as slower. The ratios are deliberately
        % left equal rather than retuned to hold the depth fixed -- retuning
        % them would move the steady state the parent was tuned at, and the
        % extra depth is wanted here.
        %
        % Cost: the std_only and sfa_and_std conditions carry 2000 b-states
        % rather than 1000, so N_sys_eqs goes 2250 -> 3250 at n = 500.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    500, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.5 0.5], ...
            'mu_tilde_relative',    [5.5 -5.5; 5.5 -5.5], ...   % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1.5 1.5; 1.5 1.5], ...     % multiples of F
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.20, ...    % operative: mu_S_c is empty
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'c',                    [0.5/3, 0], ...     % SFA scaling, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.025);

        std_routes = struct();
        dual_std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
        std_routes.E.E.std = dual_std;
        std_routes.E.I.std = dual_std;
        std_routes.I.E.std = dual_std;
        std_routes.I.I.std = dual_std;

    % ====================================================================
    %  FIGURE PRESETS. Each exists because one manuscript figure is
    %  DELIBERATELY a different network from the paper's operating point,
    %  and hardcoding that network inside the figure script is what let the
    %  paper drift into showing several unrelated models.
    % ====================================================================

    case 'bursting_pairs'
        % The hand-tuned bursting network from
        % scripts/presentations/Stability_Manuscript/fig_stim_engages_adaptation/
        % bursting_SRNN_example.m, ported from SRNNModel2 to SRNNCellTypePairs.
        % Derived presets: none.
        %
        % Numerically equivalent to that script by construction -- every value
        % below is copied from it. Only the CLASS changes, so the whole paper
        % speaks one model vocabulary.
        %
        % WHY THIS NETWORK IS DIFFERENT, and must stay different: it is small
        % (n = 50) and sparse (indegree = 10, so alpha = 0.2) precisely so that
        % individual neurons are visible as traces and the population is small
        % enough to burst coherently. The paper's operating point (n = 500,
        % indegree = 100) averages that away.
        %
        % TRANSLATION FROM THE SRNNModel2 SPELLING:
        %   f 0.7                      -> f = [0.7 0.3]
        %   mu_E_tilde_relative  3     -> mu_tilde_relative(:,1) = 3   (E is PRE)
        %   mu_I_tilde_relative -2     -> mu_tilde_relative(:,2) = -2  (I is PRE)
        %   sigma_{E,I}_tilde_relative -> sigma_tilde_relative = ones(2)
        %   c_E 0.5/3, c_I unused      -> c = [0.5/3, 0]   (SFA on E only)
        % The mu block is column-uniform because SRNNModel2 has no notion of a
        % postsynaptic-dependent mean: every neuron sees the same statistics
        % from a given presynaptic population. Writing it as a full 2x2 with
        % equal rows is the exact Pairs statement of that.
        %
        % The DC-staircase stimulus is NOT here. input_config carries a
        % generator FUNCTION HANDLE (@dc_staircase_stimulus) plus the sweep of
        % DC levels, which is stimulus protocol rather than network physics, and
        % a handle in a preset is exactly what the "nonlinearity is named data"
        % rule exists to avoid. The figure builds it.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    50, ...
            'indegree',             10, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.7 0.3], ...
            'mu_tilde_relative',    [3 -2; 3 -2], ...   % (post <- pre), column-uniform
            'sigma_tilde_relative', [1 1; 1 1], ...
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.9, ...
            'S_c',                  0.5, ...
            'tau_d',                0.1, ...
            'c',                    [0.5/3, 0]);        % SFA scaling, E only
        % STD on E->E only. SRNNModel2's n_b_E = 1 depresses every OUTGOING
        % excitatory synapse, i.e. E->E and E->I alike; the source script had
        % n_b_I = 0, so I synapses did not depress. E->E and E->I with the same
        % constants is the faithful translation.
        std_routes = struct();
        bursting_std = struct('tau_rec', 1, 'tau_rel', 0.25);
        std_routes.E.E.std = bursting_std;
        std_routes.E.I.std = bursting_std;

    case 'sompolinsky_pairs'
        % The Sompolinsky-Crisanti-Sommers (1988) reproduction behind Figure 1
        % panel A (traces + eigenspectra), on SRNNCellTypePairs.
        % Derived presets: none.
        %
        % A purely RANDOM, Dale-free network: zero-mean Gaussian weights, tanh,
        % no adaptation, no external input. The spectral radius is then
        % R = level_of_chaos exactly, so chaos onset sits at gain 1 -- which is
        % the entire point of the figure.
        %
        % ONE CELL TYPE, named 'all'. This network has a single undifferentiated
        % population, and now says so: the weights take both signs, so 'E'/'I'
        % would be a lie, and there is no second type to name. Nothing here needs
        % the f_E / mu_EE_relative aliases, which require two types by design.
        %
        % (Until 2026-08-23 this was TWO statistically identical types named
        % 'A'/'B', purely because SRNNCellTypePairs could not build a one-type
        % model -- build_network configured RMTBlocks piecemeal where set_types
        % is required to change D. That is fixed; see singleCellTypeRefactor.md.
        % W is BIT-IDENTICAL across that change, which was not expected: the two
        % types had identical zero-mean statistics, so the per-block scaling was
        % uniform and the underlying draw never depended on how it was
        % partitioned. Figure 1 panel A is unchanged, LLEs included, and
        % test_pairs_single_celltype.m carries the checksum that keeps it so.)
        %
        % level_of_chaos is the figure's variable (gammas 0.9 / 1.6 / 2.5); the
        % value here is only the operating point for anything holding it fixed.
        % A grid axis or an explicit override wins, as always.
        %
        % tau_d = 1, not the 0.1 used everywhere else: the reference result is
        % stated in units of the membrane time constant.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    200, ...
            'indegree',             200, ...    % fully connected -> alpha = 1
            'n_cellTypes',          1, ...
            'cell_type_names',      {{'all'}}, ...
            'f',                    1, ...
            'mu_tilde_relative',    0, ...      % zero mean -> both signs, Dale-free
            'sigma_tilde_relative', 1, ...
            'level_of_chaos',       1.6, ...
            'activation',           'tanh', ...       % uses neither S_a nor S_c
            'mu_S_c',               [], ...           % MUST stay empty: per-type setpoints
            'sigma_S_c',            [], ...           % error under 'tanh' (no centre to vary)
            'tau_d',                1.0, ...
            'c',                    0, ...            % no SFA
            'x0_std',               1.0);             % visible relaxation in the stable panel
        % No synapses depress or facilitate: this is the bare random network.
        std_routes = struct();

    case 'single_neuron_stf'
        % The SFA / STD / STF single-neuron methods figure
        % (fig_adaptation_methods figure B), rebuilt on SRNNCellTypePairs.
        % Derived presets: none.
        %
        % Restores the figure produced by the deleted
        % fig_adaptation_methods/test_single_neuron_stf.m (last version at commit
        % 60c2992), which called SRNNModelCellTypes.dynamics_fast_ct directly --
        % a class that no longer exists.
        %
        % THE FACILITATION PARAMETERS MAP EXACTLY, which is not obvious. The old
        % model carried a release probability p resting at p0 with gain p/p0:
        %
        %   dp/dt = (p0 - p)/tau_f + kappa*(1-p)*r ,   gain = p/p0
        %
        % Substituting the gain variable u = p/p0 (so p = u*p0) gives
        %
        %   du/dt = (1-u)/tau_f + kappa*(1/p0 - u)*r
        %
        % which IS this class's  dg/dt = (1-g)/tau_dec + (G-g)*r/tau_fac  with
        %
        %   tau_dec = tau_f    = 6      (relaxation back to rest)
        %   tau_fac = 1/kappa  = 2.5    (kappa was 0.4)
        %   G       = 1/p0     = 2.857  (p0 was 0.35)
        %
        % Both rest at gain 1 and rise toward the same ceiling.
        %
        % WHAT DOES NOT MAP: the old STD depleted in proportion to p,
        %   db/dt = (1-b)/tau_rec - (p*b*r)/tau_rel,
        % where this class's depletion is independent of g. That is the
        % Tsodyks-Markram coupling, and no value of tau_rel reproduces it -- at
        % rest it is a constant factor (p0 = 0.35) but it ACCELERATES as p rises.
        % tau_rel is carried LITERALLY as 0.3 by decision, which makes depression
        % about 2.9x stronger at rest than the archived PNG shows. The rebuilt
        % figure will therefore not match that image, and that is expected.
        % Whether the decoupling is right at all is parked in UserNotes.md.
        %
        % Exaggerated on purpose: c = 1.0 and a single tau_a = 3 s, so the rate
        % decay is legible on one neuron. This is a mechanism cartoon, not the
        % paper's operating point.
        %
        % ONE cell type, ONE neuron, zero weights: n = 1, indegree = 1,
        % mu = sigma = 0, so W = 0 (1x1) and there is no recurrence at all.
        % (Before 2026-08-23 this had to be n = 2 with two types, because the
        % class could not build a one-type model and enforces n >= n_cellTypes.
        % The second neuron was inert but had to be explained away.) indegree
        % must be >= 1, so the single neuron nominally connects to itself --
        % with mu = sigma = 0 that weight is zero, so it makes no difference.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    1, ...
            'indegree',             1, ...
            'n_cellTypes',          1, ...
            'cell_type_names',      {{'E'}}, ...
            'f',                    1, ...
            'mu_tilde_relative',    0, ...            % W == 0: no recurrence at all
            'sigma_tilde_relative', 0, ...
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  1.0, ...          % hard sigmoid (piecewise linear)
            'S_c',                  0.5, ...
            'tau_d',                0.1, ...
            'c',                    1.0, ...          % exaggerated SFA
            'tau_a',                {{3}}, ...        % single 3 s timescale
            'x0_std',               0);               % deterministic x(0) = 0
        % ONE SFA timescale, not the usual three. The conditions must agree with
        % the tau_a above or validate() rejects the pair ("tau_a{1} must contain
        % n_a(1) positive values") -- which is precisely what this argument is
        % for. It is the only preset in the file that needs it.
        n_a_sfa = 1;
        % Both mechanisms on the E->E route. The figure switches them on and off
        % per column; these are the timescales it switches ON.
        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.3);
        std_routes.E.E.stf = struct('tau_dec', 6, 'tau_fac', 2.5, 'G', 1/0.35);

    case 'single_neuron_dualStd'
        % The SFA / STD single-neuron methods figure (fig_adaptation_methods
        % figure A), on the PAPER'S timescales.
        % Copied from celltype_pairs_Sc0p2_noise0p025_dualStd. Changed:
        %   n_cellTypes  2            -> 1
        %   n            500          -> 1
        %   indegree     100          -> 1
        %   f            [0.5 0.5]    -> 1
        %   mu_tilde_relative    2x2  -> 0     (W == 0)
        %   sigma_tilde_relative 2x2  -> 0     (W == 0)
        %   c            [0.5/3, 0]   -> 0.5/3
        %   input_config, F_ref_*     -> dropped (no recurrence, no sweep)
        % Derived presets: none.
        %
        % WHY THIS EXISTS. The figure is a mechanism cartoon: one unconnected
        % neuron given a step, with SFA and STD switched on and off per column.
        % It must carry the paper's timescales -- so that the cartoon explains
        % the network figures rather than some other model -- WITHOUT the
        % paper's recurrence, which would drown the mechanisms in network
        % dynamics.
        %
        % It exists because until 2026-08-23 the figure did exactly that wrong
        % thing: it named the dualStd preset and overrode no network parameters,
        % so it built the whole n = 500, indegree = 100 recurrent network and
        % plotted neuron 1. Its own generated README printed "One unconnected
        % neuron" above "n 500". The result showed none of the mechanisms -- the
        % no-adaptation column was network chaos and the STD columns were
        % silent with b ~ 1. Naming the network as a preset is what stops that
        % recurring; see singleCellTypeRefactor.md section 3c.
        %
        % NOISE IS DELIBERATELY KEPT (sigma_u_noise = 0.025, exactly the
        % paper's). On one neuron there is no population averaging, so the
        % jitter is fully visible: x_noise_std = sigma_u/sqrt(2*tau_d) = 0.056
        % against a 0.5 step, about 11%. That is accepted as honest about the
        % model the paper characterises (TR's decision) rather than smoothed
        % away.
        %
        % n_a_sfa is NOT overridden: unlike single_neuron_stf this shows all
        % three of the paper's SFA timescales, auto-filled by
        % complete_type_defaults, at the paper's per-timescale c = 0.5/3.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    1, ...
            'indegree',             1, ...
            'n_cellTypes',          1, ...
            'cell_type_names',      {{'E'}}, ...
            'f',                    1, ...
            'mu_tilde_relative',    0, ...      % W == 0: no recurrence at all
            'sigma_tilde_relative', 0, ...
            'level_of_chaos',       1.0, ...
            'activation',           'piecewise', ...
            'S_a',                  0.8, ...
            'S_c',                  0.20, ...
            'mu_S_c',               [], ...
            'sigma_S_c',            [], ...
            'tau_d',                0.1, ...
            'c',                    0.5/3, ...  % SFA scaling, per timescale
            'x0_std',               0, ...      % deterministic x(0) = 0
            'sigma_u_noise',        0.025);
        % The paper's dual depression, on the one route this network has.
        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);

    case 'mc_esn'
        % The memory-capacity network, from scripts/memory_capacity/
        % run_memory_capacity.m (was looped_memory_capacity.m).
        % Derived presets: none.
        %
        % THE ONE SRNNModel2 PRESET LEFT IN THE PAPER, and not by choice:
        % SRNN_ESN_reservoir subclasses SRNNModel2, so the memory-capacity
        % protocol cannot run on SRNNCellTypePairs at all. Porting the ESN
        % readout onto the Pairs class is tracked follow-up work; until then the
        % MC figures show a different network from every other figure, and the
        % methods section has to say so.
        %
        % WHAT IS DELIBERATELY ABSENT. The MC protocol settings -- input_type,
        % T_hold, T_wash, T_train, T_test, d_max, u_f_cutoff, u_alpha -- are
        % NOT here, and neither are fs / ode_solver. They size the experiment
        % rather than describing the network, which makes them run_mode knobs in
        % this project's vocabulary; run_memory_capacity owns them. Keeping them
        % out is what lets 'fast' and 'production' MC runs share one network.
        %
        % f = 0.6, off perfect balance, so the no-adaptation condition is not
        % accidentally favoured. level_of_chaos = 2 sits above the edge of chaos
        % (the logistic's mean slope < 1 pushes that edge past 1).
        model_class = 'SRNNModel2';
        d = struct( ...
            'n',              300, ...
            'f',              0.6, ...
            'level_of_chaos', 2.0, ...
            'tau_d',          0.1, ...
            'activation',     'logistic', ...
            'S_c',            0.35, ...
            'S_a',            0.9, ...      % unused by the logistic; recorded for parity
            ... % n_a_I / n_b_I are NOT here. looped_memory_capacity set them to 0
            ... % explicitly, but they are condition-owned fields and a preset may
            ... % not carry any of n_a_E / n_a_I / n_b_E / n_b_I. Both are already
            ... % the SRNNModel2 default (0), so nothing changes -- I neurons still
            ... % get no SFA and no STD in every condition.
            'c_E',            0.5/3, ...
            'tau_a_E',        logspace(log10(0.1), log10(10), 3), ...
            'tau_b_E_rec',    1.0, ...
            'tau_b_E_rel',    0.25, ...
            'std_zero_floor', false);

    otherwise
        error('srnn_param_preset:UnknownPreset', ...
            'Unknown preset ''%s''. Valid presets: %s.', ...
            name, strjoin(srnn_param_preset_names(), ', '));
end

% std_routes stays [] for a preset that wants the default routes;
% srnn_adaptation_conditions reads [] as "use your own default", so the default
% lives in exactly one place rather than being restated here.
%
% The cell-type count comes from n_cellTypes, which every SRNNCellTypePairs
% preset states explicitly, and NOT from numel(d.f): the 'default' preset has no
% f field at all, so that would throw. SRNNModel2 presets have no n_cellTypes
% and take the 2 that class is fixed at.
if isfield(d, 'n_cellTypes')
    n_cell_types = d.n_cellTypes;
else
    n_cell_types = 2;
end
conditions = srnn_adaptation_conditions(model_class, std_routes, n_a_sfa, n_cell_types);
end

function names = srnn_param_preset_names()
% The valid preset names, kept next to the switch above so they stay in sync.
names = {'default', 'overconnected', 'celltype_pairs', ...
    'celltype_pairs_S_c_by_type', 'celltype_pairs_S_c_by_type_n500', ...
    'celltype_pairs_S_c_by_type_n500_fixedF', 'celltype_pairs_all_std_n500', ...
    'celltype_pairs_uniform_std_n500', 'celltype_pairs_uniform_std_n500_mu5p5', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p02', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01', ...
    'celltype_pairs_Sc0p2_noise0p025', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd', ...
    ... % figure presets -- networks that are deliberately not the paper's
    ... % operating point, named so the figures stop hardcoding them
    'bursting_pairs', 'sompolinsky_pairs', 'single_neuron_stf', ...
    'single_neuron_dualStd', 'mc_esn'};
end

function ic = pairs_input_config(intrinsic_drive)
% Assigning input_config REPLACES the struct SRNNCellTypePairs.set_defaults
% built, so a preset that wants to change one field has to restate all of them.
% This mirrors that default exactly apart from intrinsic_drive.
%
% This helper is field-level, not preset-level: it spells out one property, not
% one preset, so it does not reintroduce the inter-preset coupling this file was
% flattened to remove. Every case still states its own call explicitly.
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

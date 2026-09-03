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
std_routes     = [];    % [] = whatever srnn_adaptation_conditions defaults to
sfa_timescales = []; %#ok<NASGU> -- a SENTINEL, and checkcode is right that every
                        % case currently overwrites it. That is the point: it
                        % exists so a case that FORGETS reaches the isempty()
                        % guard after the switch and gets an error naming the
                        % preset, instead of MATLAB's "Unrecognized function or
                        % variable 'sfa_timescales'".
                        % EVERY case must set this; see that check.
                        % The adaptation timescales the SFA conditions switch on,
                        % in seconds, normally log_ladder(lo, hi, K).
                        %
                        % NO SHARED DEFAULT, deliberately. This used to be a call
                        % to srnn_sfa_timescales(3), so a preset that said nothing
                        % inherited the paper's ladder. Timescales are PHYSICS,
                        % and physics belongs in the preset that describes the
                        % network -- every other physics value here is restated
                        % per case, including the depression timescales, which
                        % are the direct analogue. Three different ladders
                        % already existed in the repo, so the shared default was
                        % never the single source it looked like.
                        %
                        % tau_a itself must NOT appear in the `d` struct below:
                        % it is condition-owned, exactly like synapse_config, and
                        % a preset carrying it would be silently overridden by
                        % whichever condition is applied (and warned about by
                        % validate_model_defaults).
regimes        = 'standard';
                        % 'standard' = the four original regimes.
                        % 'timescales' adds one-timescale variants, for a preset
                        % that exists to compare timescale STRUCTURE. Costs 75%
                        % more compute per sweep, so it is opt-in.

switch name
    case 'default'
        % Everything at SRNNModel2's class defaults.
        model_class = 'SRNNModel2';
        d = struct();
        sfa_timescales = log_ladder(0.25, 10, 3);

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
            'c_E',                 1.0, ...   % TOTAL budget (was 1/3, hand-divided by n_a_E=3)
            'mu_E_tilde_relative',  3, ...   % class default
            'mu_I_tilde_relative', -6);      % class default -4; doubled, half as many I neurons
        sfa_timescales = log_ladder(0.25, 10, 3);

    case retired_presets()
        % THE EXPLORATORY celltype_pairs_* FAMILY, RETIRED 2026-08-25.
        %
        % These 13 presets carried a HAND-DIVIDED c (typically 0.5/3) written on
        % the assumption that every condition runs n_a = 3. The model now divides
        % c by the number of timescales actually in use, so those values would be
        % divided twice and every one of these networks would silently run at one
        % third of its intended adaptation. Nothing would report it.
        %
        % They are retired rather than rescaled because they are historical
        % exploration -- the path that produced the paper's operating point, not
        % the operating point itself. Their parameter values remain in git
        % (last live at commit cbcf637) and in the run directories they produced,
        % both of which record what they were far better than a case body kept
        % alive for reading.
        %
        % The successors carry c as a TOTAL budget:
        %   celltype_pairs_Sc0p2_noise0p025_dualStd_4cond   the same 4 regimes
        %   celltype_pairs_Sc0p2_noise0p025_dualStd_7cond   plus 3 timescale ones
        error('srnn_param_preset:RetiredPreset', ...
            ['Preset ''%s'' was retired on 2026-08-25 by the c/K adaptation ' ...
             'change: its c is a per-timescale value and would now be divided ' ...
             'twice, giving one third the intended adaptation with no error. ' ...
             'Use celltype_pairs_Sc0p2_noise0p025_dualStd_4cond (same four ' ...
             'regimes) or celltype_pairs_Sc0p2_noise0p025_dualStd_7cond (adds ' ...
             'the single-timescale regimes). Values remain in git at cbcf637.'], ...
            name);

    case 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond'
        % THE PAPER'S NETWORK, four adaptation regimes.
        % Copied from the retired celltype_pairs_Sc0p2_noise0p025_dualStd.
        % Changed:  c  [0.5/3, 0] -> [0.5, 0]
        % Derived presets: celltype_pairs_Sc0p2_noise0p025_dualStd_7cond,
        %                  single_neuron_dualStd.
        %
        % NUMERICALLY IDENTICAL to its retired parent. c is now the TOTAL
        % adaptation budget and the model divides it by the number of timescales
        % in use, so 0.5 with n_a = 3 gives exactly the 0.5/3 per timescale the
        % parent hand-wrote. The rename is what makes that safe: a value meaning
        % something new gets a name meaning something new, rather than silently
        % changing under the old one.
        %
        % DUAL DEPRESSION on all four routes. Both pairs share
        % tau_rec/tau_rel = 8, so each b relaxes toward 1/(1 + 8r) and the
        % synapse sees their PRODUCT -- the square of a single-timescale
        % synapse's depression, 0.086 against 0.29 at r = 0.3. Deeper depression
        % is the point.
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
            'c',                    [0.5, 0], ...       % TOTAL SFA budget, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.025);

        % SFA on the standard ladder: 3 timescales spanning 0.25 s to 10 s.
        sfa_timescales = log_ladder(0.25, 10, 3);

        std_routes = struct();
        dual_std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
        std_routes.E.E.std = dual_std;
        std_routes.E.I.std = dual_std;
        std_routes.I.E.std = dual_std;
        std_routes.I.I.std = dual_std;

    case 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
        % THE PAPER'S NETWORK, seven adaptation regimes. THE MAIN PRESET.
        % Copied from celltype_pairs_Sc0p2_noise0p025_dualStd_4cond.
        % Changed:  regimes  'standard' -> 'timescales'
        % Derived presets: none.
        %
        % THE NETWORK IS IDENTICAL to the 4-condition preset -- same n, same
        % connectivity, same nonlinearity, same noise, same c, same depression
        % timescales. The only difference is which adaptation regimes get swept,
        % so results are directly comparable between the two and the four shared
        % regime names mean exactly the same thing in both.
        %
        % WHAT THE THREE EXTRA REGIMES BUY. sfa_only_oneTS, std_only_oneTS and
        % sfa3_std1 hold the mechanism fixed and vary only how many timescales
        % carry it, which isolates timescale STRUCTURE from mechanism STRENGTH.
        % That comparison is only meaningful because c is now the total
        % adaptation budget: the model divides by the number of timescales, so
        % one-timescale SFA adapts exactly as hard in total as three-timescale
        % SFA. Under the old per-timescale c it would have been three times
        % weaker and the comparison would have measured two things at once.
        %
        % Depression needs no such normalization -- it enters as a PRODUCT, each
        % b rests at 1, and the dual-timescale steady state is deliberately the
        % square of the single one -- so std_only_oneTS IS genuinely weaker
        % depression than std_only, and that is the intended contrast.
        %
        % COST: 7 regimes rather than 4 is 75% more compute at every sweep stage.
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
            'c',                    [0.5, 0], ...       % TOTAL SFA budget, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.025);

        % SFA on the standard ladder: 3 timescales spanning 0.25 s to 10 s.
        % sfa_only_oneTS takes the FIRST of these, so the fast end is what the
        % one-timescale SFA regime runs.
        sfa_timescales = log_ladder(0.25, 10, 3);

        % The one-timescale regimes take their route from these by keeping the
        % FIRST tau_rec/tau_rel pair, so they cannot describe an older network.
        std_routes = struct();
        dual_std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
        std_routes.E.E.std = dual_std;
        std_routes.E.I.std = dual_std;
        std_routes.I.E.std = dual_std;
        std_routes.I.I.std = dual_std;
        regimes = 'timescales';

    case 'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond'
        % THE PAPER'S NETWORK, three adaptation regimes: none / one timescale /
        % many. Copied from celltype_pairs_Sc0p2_noise0p025_dualStd_7cond.
        % Changed:  regimes  'timescales' -> 'single_multi'
        % Derived presets: none.
        %
        % THE NETWORK IS IDENTICAL to the 4- and 7-condition presets -- same n,
        % connectivity, nonlinearity, noise, c and depression timescales. Only
        % the regime SET differs, so a 3-condition run is directly comparable
        % with the exploratory 7-condition ones.
        %
        % ONE CAVEAT ON THAT COMPARISON, and it is a naming one, not a physics
        % one: this set calls the full-adaptation regime sfa3_std2 where the
        % other sets call it sfa_and_std. The two are byte-identical -- see
        % srnn_adaptation_conditions -- but a human lining up two run
        % directories has to know that. Only no_adaptation shares its name
        % across all three sets.
        %
        % WHY THREE, AND WHY THESE THREE (TR, 2026-09-03, after discussing the
        % manuscript with his PI). The 7-regime set separates SFA from STD --
        % sfa_only, std_only and their one-timescale variants -- which turned out
        % to be too complicated a story for this paper and is headed for a later
        % one. The interesting contrast is NO adaptation vs adaptation on ONE
        % timescale vs adaptation on MANY, with both mechanisms always present
        % together:
        %
        %   no_adaptation   -                      -
        %   sfa1_std1       1 tau_a (0.25 s)       1 STD pair (2 / 0.25 s)
        %   sfa3_std2       3 tau_a (0.25-10 s)    2 STD pairs
        %
        % So the axis is TIMESCALE COUNT and only that: every regime with
        % adaptation has SFA and STD in the same proportion. The names state the
        % counts because the counts are the content, and that is what lets the
        % figures read "Single-Timescale Adaptation" against
        % "Multiple-Timescale Adaptation" rather than "SFA + STD" twice.
        %
        % THIS COMPARISON IS ONLY CLEAN BECAUSE c IS THE TOTAL SFA BUDGET. The
        % model divides c by the number of timescales in use, so one-timescale
        % SFA adapts exactly as hard IN TOTAL as three-timescale SFA and the
        % contrast is structure, not strength. Depression is deliberately NOT
        % normalized that way -- it enters as a product, each b rests at 1, so
        % two timescales genuinely depress harder than one. The single/multi
        % contrast therefore holds SFA strength fixed while STD strength grows
        % with its timescale count, which is the intended asymmetry rather than
        % an oversight: see Equations_stability_paper.md.
        %
        % COST: 3 regimes against 7 is 57% less compute at every sweep stage.
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
            'c',                    [0.5, 0], ...       % TOTAL SFA budget, E only
            'input_config',         pairs_input_config(0.0), ...
            'F_tracks_network',     false, ...
            'F_ref_n',              500, ...
            'F_ref_indegree',       100, ...
            'sigma_u_noise',        0.025);

        % SFA on the standard ladder: 3 timescales spanning 0.25 s to 10 s.
        % sfa1_std1 takes the FIRST of these -- 0.25 s, the fast end -- which is
        % the whole content of the single-timescale regime.
        sfa_timescales = log_ladder(0.25, 10, 3);

        % sfa1_std1 takes its route from these by keeping the FIRST
        % tau_rec/tau_rel pair, so it cannot describe an older network.
        std_routes = struct();
        dual_std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
        std_routes.E.E.std = dual_std;
        std_routes.E.I.std = dual_std;
        std_routes.I.E.std = dual_std;
        std_routes.I.I.std = dual_std;
        regimes = 'single_multi';

    % ====================================================================
    %  FIGURE PRESETS. Each exists because one manuscript figure is
    %  DELIBERATELY a different network from the paper's operating point,
    %  and hardcoding that network inside the figure script is what let the
    %  paper drift into showing several unrelated models.
    % ====================================================================

    case 'bursting_pairs'
        % The hand-tuned bursting network from
        % From the bursting figure (fig_stim_engages_adaptation),
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
            'c',                    [0.5, 0]);          % TOTAL SFA budget, E only
        % The paper's standard ladder, as in the source script.
        sfa_timescales = log_ladder(0.25, 10, 3);

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
        % is required to change D. That is fixed.
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
        % Stated for completeness only: c = 0, so no condition here adapts
        % whatever timescales it is handed. The standard ladder keeps the value
        % honest rather than arbitrary.
        sfa_timescales = log_ladder(0.25, 10, 3);

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
            'x0_std',               0);               % deterministic x(0) = 0
        % ONE exaggerated 3 s timescale, so the rate decay is legible on a single
        % trace. This lives here rather than as a `tau_a` field above because
        % tau_a is condition-owned: the SFA-off columns need it EMPTY, and only a
        % condition can say that. Carrying it in the preset is what used to make
        % build_from_preset('single_neuron_stf','no_adaptation') fail outright
        % with "tau_a{1} must contain n_a(1) positive values".
        % ONE exaggerated 3 s timescale, not a ladder -- so log_ladder is not
        % what this wants. The figure shows a single unconnected neuron, and a
        % slow constant is what makes the rate decay legible on one trace.
        sfa_timescales = 3;
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
        % recurring.
        %
        % NOISE IS DELIBERATELY KEPT (sigma_u_noise = 0.025, exactly the
        % paper's). On one neuron there is no population averaging, so the
        % jitter is fully visible: x_noise_std = sigma_u/sqrt(2*tau_d) = 0.056
        % against a 0.5 step, about 11%. That is accepted as honest about the
        % model the paper characterises (TR's decision) rather than smoothed
        % away.
        %
        % sfa_timescales is NOT overridden: unlike single_neuron_stf this shows
        % all three of the paper's SFA timescales, at the paper's c.
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
            'c',                    0.5, ...    % TOTAL SFA budget
            'x0_std',               0, ...      % deterministic x(0) = 0
            'sigma_u_noise',        0.025);
        % All three of the paper's SFA timescales -- unlike single_neuron_stf,
        % which exaggerates a single one so the decay is legible on one trace.
        sfa_timescales = log_ladder(0.25, 10, 3);

        % The paper's dual depression, on the one route this network has.
        std_routes = struct();
        std_routes.E.E.std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);

    case 'mc_pairs_dualStd'
        % The memory-capacity network, on SRNNCellTypePairs.
        % Replaces 'mc_esn', deleted 2026-09-02 when SRNN_ESN_reservoir was
        % re-parented from SRNNModel2 onto SRNNCellTypePairs. Memory capacity
        % was the last part of the paper on a different model class; it no
        % longer is, and the methods section no longer has to say so.
        % Derived presets: none.
        %
        % PHYSICS CARRIED OVER FROM mc_esn unchanged: n = 300, level_of_chaos =
        % 2.0 (above the edge of chaos -- the logistic's mean slope < 1 pushes
        % that edge past 1), S_c = 0.35, a TOTAL SFA budget of 0.5 over three
        % timescales, and f off perfect balance so the no-adaptation condition
        % is not accidentally favoured. mc_esn's scalar f = 0.6 becomes the row
        % [0.6 0.4]; that is the same partition, since the two cell-type
        % partitioning algorithms provably agree at C = 2 (see commit 1ed6789).
        %
        % DUAL DEPRESSION ON ALL FOUR ROUTES, matching the paper's preset rather
        % than mc_esn's route set. This is not cosmetic. mc_esn ran STD via
        % n_b_E = 1, which on SRNNModel2 depresses every outgoing E synapse --
        % E->E and E->I, but nothing from I. Putting identical STD on all four
        % routes is what makes the 'synaptic' readout well defined here: a
        % presynaptic neuron's outgoing routes then share their (tau_rec,
        % tau_rel), so its b trajectories are bit-identical down every route and
        % it has ONE synaptic output. SRNN_ESN_reservoir.assert_route_redundancy
        % enforces exactly that and errors if a future edit breaks it.
        %
        % NOISE IS OFF, deliberately and for now. mc_run_config selects the
        % integrator the same way analysis_run_config does, so sigma_u_noise > 0
        % would switch MC to 'sra1' automatically and work -- but the MC protocol
        % has not been re-validated with noise, and TR wants the first ported run
        % comparable on one axis at a time. Set sigma_u_noise here to turn it on;
        % nothing else needs to change.
        model_class = 'SRNNCellTypePairs';
        d = struct( ...
            'n',                    300, ...
            'indegree',             100, ...
            'n_cellTypes',          2, ...
            'cell_type_names',      {{'E', 'I'}}, ...
            'f',                    [0.6 0.4], ...
            'mu_tilde_relative',    [3 -4; 3 -4], ...    % multiples of F, (post <- pre)
            'sigma_tilde_relative', [1 1; 1 1], ...      % multiples of F
            'level_of_chaos',       2.0, ...
            'tau_d',                0.1, ...
            'activation',           'logistic', ...
            'S_c',                  0.35, ...
            'c',                    [0.5, 0], ...        % TOTAL SFA budget, E only
            'std_zero_floor',       false, ...
            'sigma_u_noise',        0);

        % A DIFFERENT LADDER from the paper's 0.25-10: this one starts at 0.1 s,
        % carried over from mc_esn. Presets state their own timescales precisely
        % so a second ladder like this is visible here rather than hidden behind
        % a shared default that only ever described the first one.
        sfa_timescales = log_ladder(0.1, 10, 3);

        std_routes = struct();
        mc_dual_std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
        std_routes.E.E.std = mc_dual_std;
        std_routes.E.I.std = mc_dual_std;
        std_routes.I.E.std = mc_dual_std;
        std_routes.I.I.std = mc_dual_std;

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

% Every case states its own ladder, so a new preset that forgets fails HERE with
% its own name rather than silently adopting whatever the last shared default
% happened to be. That silent adoption is the failure this replaced: presets
% inherited log_ladder(0.25, 10, 3) without saying so, while two of them already
% used different ladders.
if isempty(sfa_timescales)
    error('srnn_param_preset:NoSfaTimescales', ...
        ['Preset ''%s'' does not set sfa_timescales. Every preset states its ' ...
         'own adaptation timescales -- normally log_ladder(lo, hi, K) -- ' ...
         'because they are physics and physics lives in the preset.'], name);
end

conditions = srnn_adaptation_conditions(model_class, ...
    'synapse_config', std_routes, 'sfa_timescales', sfa_timescales, ...
    'n_cell_types', n_cell_types, 'regimes', regimes);
end

function names = retired_presets()
% The exploratory celltype_pairs_* family, retired by the c/K change.
%
% Named here rather than deleted from the file entirely so that asking for one
% gets a message explaining what happened and what replaced it, instead of
% "Unknown preset" with a list that no longer mentions them.
names = {'celltype_pairs', ...
    'celltype_pairs_S_c_by_type', 'celltype_pairs_S_c_by_type_n500', ...
    'celltype_pairs_S_c_by_type_n500_fixedF', 'celltype_pairs_all_std_n500', ...
    'celltype_pairs_uniform_std_n500', 'celltype_pairs_uniform_std_n500_mu5p5', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p02', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01', ...
    'celltype_pairs_Sc0p2_noise0p025', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd'};
end

function names = srnn_param_preset_names()
% The valid preset names, kept next to the switch above so they stay in sync.
% Retired names are NOT here: they are a separate list with their own error.
names = {'default', 'overconnected', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond', ...
    ... % figure presets -- networks that are deliberately not the paper's
    ... % operating point, named so the figures stop hardcoding them
    'bursting_pairs', 'sompolinsky_pairs', 'single_neuron_stf', ...
    'single_neuron_dualStd', 'mc_pairs_dualStd'};
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

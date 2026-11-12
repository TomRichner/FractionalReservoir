function R = p3_run_one(network, condition, ablation, ref, opts)
% P3_RUN_ONE Build, run and measure one network x condition x ablation of pilot 3.
%
%   R = P3_RUN_ONE(network, condition)                   the intact run
%   R = P3_RUN_ONE(network, condition, ablation, ref)    an ablation of it
%
% network     'bursting' (bursting_condition_model) or 'run13'
%             (run13_condition_model). The builders are called, never
%             re-derived; the ablation reaches them through 'model_args'.
% condition   'original' | 'chain_mts' | 'product_mts'
% ablation    see p3_ablation_args (default 'intact')
% ref         the intact run's R.ref (same network and condition)
%
% Options: 'T_show' (100 s analysed), 'seed' (42; the bursting network's
% pair is [seed seed], run13's noise seed follows the builder's rule).
% Settling is 10 s in both networks.
%
% A CONSTANT INPUT (the mean-matched ablations) is given as the
% input_config's n x 1 intrinsic_drive. The model applies no input before
% t = 0 (generate_stimulus fills the settling interval with zeros), so
% these runs use T_range = [0, 10 + T_show] and are analysed over the last
% T_show seconds: the same number of steps, hence the same noise path, as
% the intact run's [-10, T_show].
%
% R.t                1 x nt, s, 0 .. T_show (100 Hz bursting, 200 Hz run13)
% R.rE, R.rI         neurons x nt rates;  R.x  n x nt (network order)
% R.groups           n x 1 module label (ones for the bursting network)
% R.sfa_E, R.sfa_I   1 x nt mean adaptation term (c/K) sum_k a_k
% R.dep              1 x nt mean depression factor over the depressed
%                    routes (std_route_factor); NaN if there is none
% R.meas             the measures (p3_measures)
% R.spec             the ablation (p3_ablation_args)
% R.ref              what an ablation of THIS run needs: tau_a,
%                    synapse_config, c, n_a, names, type_indices, n, and
%                    the time averages abar (n x 1 adaptation term),
%                    D_route / D_route_rw (C x C, pre by post) and
%                    D_neuron (C x C cell of n_pre x 1)
%
% See also: p3_ablation_args, p3_measures, run_p3_ablation

arguments
    network   (1,:) char {mustBeMember(network, {'bursting', 'run13'})}
    condition (1,:) char
    ablation  (1,:) char = 'intact'
    ref = []
    opts.T_show (1,1) double = 100
    opts.seed   (1,1) double = 42
end
T_settle = 10;
Tw = opts.T_show;
spec = p3_ablation_args(ablation, ref);
margs = spec.model_args;
started = tic;

switch network
    case 'bursting'
        ic = [];
        if ~isempty(spec.drive)
            ic = struct('n_steps', 1, 'step_density', struct(), 'amp', 0, ...
                'no_stim_pattern', true, 'intrinsic_drive', spec.drive, 'positive_only', false);
            margs = [margs, {'T_range', [0, T_settle + Tw]}];
        end
        assert(isempty(spec.col_scale), 'A per-neuron weight factor needs replace_connectivity (run13 only).');
        m = bursting_condition_model(condition, 'rng_seeds', [opts.seed opts.seed], ...
            'T_settle', T_settle, 'T_show', Tw, 'input_config', ic, 'model_args', margs);
        groups = ones(m.n, 1);
    case 'run13'
        ic = [];
        if ~isempty(spec.drive)
            preset_args = srnn_param_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium');
            ic = preset_args.input_config;
            ic.intrinsic_drive = spec.drive;
            ic.no_stim_pattern = true(1, 3);
            margs = [margs, {'T_range', [0, T_settle + Tw], 'T_plot', [T_settle, T_settle + Tw]}];
        end
        [m, groups] = run13_condition_model(condition, 'duration', Tw, 'seed', opts.seed, ...
            'input_config', ic, 'model_args', margs);
        if ~isempty(spec.col_scale)
            m.replace_connectivity(m.W * spdiags(spec.col_scale, 0, m.n, m.n));
        end
end
if isempty(spec.drive)
    assert(all(m.u_ex == 0, 'all'), 'The input must be zero.');
else
    assert(max(abs(m.u_ex - spec.drive), [], 'all') == 0, 'The constant input is not what was asked for.');
end
m.run();

%% Traces over the last T_show seconds, time re-zeroed
p = m.plot_data;
params = m.get_params();
names = m.cell_type_names;
C = numel(names);
t = p.t(:)';
keep = t >= t(end) - Tw - 1e-9;
t = t(keep) - (t(end) - Tw);

R = struct('network', network, 'condition', condition, 'ablation', ablation, 'spec', spec);
R.t = t;
R.rE = p.r.E(:, keep);
R.rI = p.r.I(:, keep);
R.x = [p.x.E(:, keep); p.x.I(:, keep)];
R.groups = groups;
R.E_idx = m.type_indices{1};
R.I_idx = m.type_indices{2};

% Adaptation term per neuron, (c/K) sum_k a_k, with the model's own c_eff.
abar = zeros(m.n, 1);
sfa = zeros(C, numel(t));
for q = 1:C
    a = p.a.(names{q});
    if isempty(a) || size(a, 2) == 0, continue; end
    term = params.c_eff(q) .* reshape(sum(a(:, :, keep), 2), size(a, 1), []);   % n_q x nt
    abar(m.type_indices{q}) = mean(term, 2);
    sfa(q, :) = mean(term, 1);
end
R.sfa_E = sfa(1, :);
R.sfa_I = sfa(2, :);

% Depression factor per route, from the model.
rates = {R.rE, R.rI};
D_route = nan(C); D_route_rw = nan(C); D_neuron = cell(C);
acc = zeros(0, numel(t));
for pre = 1:C
    for post = 1:C
        b = [];
        if isfield(p.b, names{pre}) && isfield(p.b.(names{pre}), names{post})
            b = p.b.(names{pre}).(names{post});
        end
        if isempty(b), continue; end
        D = SRNNCellTypePairs.std_route_factor(b(:, :, keep), params, pre, post);
        D = reshape(D, size(D, 1), []);                       % n_pre x nt
        r = rates{pre};
        D_neuron{pre, post} = mean(D, 2);
        D_route(pre, post) = mean(D, 'all');
        D_route_rw(pre, post) = sum(r .* D, 'all') / sum(r, 'all');
        acc = [acc; D]; %#ok<AGROW>
    end
end
if isempty(acc), R.dep = nan(1, numel(t)); else, R.dep = mean(acc, 1); end

R.ref = struct('names', {names}, 'n', m.n, 'type_indices', {m.type_indices}, ...
    'tau_a', {m.tau_a}, 'n_a', m.n_a, 'c', m.c, 'synapse_config', m.synapse_config, ...
    'abar', abar, 'D_route', D_route, 'D_route_rw', D_route_rw, 'D_neuron', {D_neuron});
R.meas = p3_measures(R);
R.run_seconds = toc(started);
end
